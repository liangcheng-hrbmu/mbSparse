import time
import os
import argparse
import sys
import numpy as np
import pickle as pkl
import networkx as nx
import scipy.sparse as sp
import resource
import datetime
import torch
from torch.utils.data import Dataset, DataLoader
from torch import nn, optim
from torch.nn import functional as F
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_samples, silhouette_score, adjusted_rand_score
from sklearn.cluster import KMeans, SpectralClustering, AffinityPropagation, AgglomerativeClustering, Birch, DBSCAN, FeatureAgglomeration, OPTICS, MeanShift
from .model import AE, VAE
from .util_function import *
from .graph_function import *
from .benchmark_util import *
from .gae_embedding import GAEembedding, measure_clustering_results, test_clustering_benchmark_results
import torch.multiprocessing as mp

def train(epoch, train_loader, model, optimizer, args, device, EMFlag=False):
    model.train()
    train_loss = 0
    for batch_idx, (data, dataindex) in enumerate(train_loader):
        if args.precisionModel == 'Double':
            data = data.type(torch.DoubleTensor)
        elif args.precisionModel == 'Float':
            data = data.type(torch.FloatTensor)
        data = data.to(device)
        regulationMatrixBatch = None

        optimizer.zero_grad()
        if args.model == 'VAE':
            recon_batch, mu, logvar, z = model(data)
            if EMFlag and (not args.EMreguTag):
                loss = loss_function_graph(recon_batch, data.view(-1, recon_batch.shape[1]), mu, logvar, gammaPara=args.gammaPara, regulationMatrix=regulationMatrixBatch,
                                           regularizer_type='noregu', reguPara=args.alphaRegularizePara, modelusage=args.model, reduction=args.reduction)
            else:
                loss = loss_function_graph(recon_batch, data.view(-1, recon_batch.shape[1]), mu, logvar, gammaPara=args.gammaPara, regulationMatrix=regulationMatrixBatch,
                                           regularizer_type=args.regulized_type, reguPara=args.alphaRegularizePara, modelusage=args.model, reduction=args.reduction)
        elif args.model == 'AE':
            recon_batch, z = model(data)
            mu_dummy = ''
            logvar_dummy = ''
            if EMFlag and (not args.EMreguTag):
                loss = loss_function_graph(recon_batch, data.view(-1, recon_batch.shape[1]), mu_dummy, logvar_dummy, gammaPara=args.gammaPara,
                                           regulationMatrix=regulationMatrixBatch, regularizer_type='noregu', reguPara=args.alphaRegularizePara, modelusage=args.model, reduction=args.reduction)
            else:
                loss = loss_function_graph(recon_batch, data.view(-1, recon_batch.shape[1]), mu_dummy, logvar_dummy, gammaPara=args.gammaPara, regulationMatrix=regulationMatrixBatch,
                                           regularizer_type=args.regulized_type, reguPara=args.alphaRegularizePara, modelusage=args.model, reduction=args.reduction)

        # L1 and L2 regularization
        # 0.0 for no regularization
        l1 = 0.0
        l2 = 0.0
        for p in model.parameters():
            l1 = l1 + p.abs().sum()
            l2 = l2 + p.pow(2).sum()
        loss = loss + args.L1Para * l1 + args.L2Para * l2

        loss.backward()
        train_loss += loss.item()
        optimizer.step()
        if batch_idx % args.log_interval == 0:
            print('Train Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(
                epoch, batch_idx * len(data), len(train_loader.dataset),
                100. * batch_idx / len(train_loader),
                loss.item() / len(data)))

        if batch_idx == 0:
            recon_batch_all = recon_batch
            data_all = data
            z_all = z
        else:
            recon_batch_all = torch.cat((recon_batch_all, recon_batch), 0)
            data_all = torch.cat((data_all, data), 0)
            z_all = torch.cat((z_all, z), 0)

    print('====> Epoch: {} Average loss: {:.4f}'.format(
          epoch, train_loss / len(train_loader.dataset)))

    return recon_batch_all, data_all, z_all


def generateAdj_Cluster(data, k=10, feature_train_epochs=10, feature_pretrain_epochs=10):
    parser = argparse.ArgumentParser(description='Main entrance of scGNN')
    parser.add_argument('--batch-size', type=int, default=64, metavar='N',
                        help='input batch size for training (default: 12800)')
    parser.add_argument('--knn-distance', type=str, default='euclidean',
                        help='KNN graph distance type: euclidean/cosine/correlation (default: euclidean)')
    parser.add_argument('--converge-graphratio', type=float, default=0.01,
                        help='converge condition: ratio of graph ratio change in EM iteration (default: 0.01), 0-1')
    parser.add_argument('--model', type=str, default='AE',
                        help='VAE/AE (default: AE)')
    parser.add_argument('--gammaPara', type=float, default=0.1,
                        help='regulized intensity (default: 0.1)')
    parser.add_argument('--alphaRegularizePara', type=float, default=0.9,
                        help='regulized parameter (default: 0.9)')    
    parser.add_argument('--useGAEembedding', action='store_true', default=True,
                        help='whether use GAE embedding for clustering(default: False)')
    parser.add_argument('--useBothembedding', action='store_true', default=False,
                        help='whether use both embedding and Graph embedding for clustering(default: False)')
    parser.add_argument('--alpha', type=float, default=0.5,
                        help='iteration alpha (default: 0.5) to control the converge rate, should be a number between 0~1')
    parser.add_argument('--GAEmodel', type=str,
                        default='gcn_vae', help="models used")
    parser.add_argument('--GAEepochs', type=int, default=200,
                        help='Number of epochs to train.')
    parser.add_argument('--GAEhidden1', type=int, default=32,
                        help='Number of units in hidden layer 1.')
    parser.add_argument('--GAEhidden2', type=int, default=16,
                        help='Number of units in hidden layer 2.')
    parser.add_argument('--GAElr', type=float, default=0.01,
                        help='Initial learning rate.')
    parser.add_argument('--GAEdropout', type=float, default=0.,
                        help='Dropout rate (1 - keep probability).')
    parser.add_argument('--GAElr_dw', type=float, default=0.001,
                        help='Initial learning rate for regularization.')
    
    parser.add_argument('--no-cuda', action='store_true', default=False,
                        help='Disable GPU training. If you only have CPU, add --no-cuda in the command line')
    parser.add_argument('--seed', type=int, default=1, metavar='S',
                        help='random seed (default: 1)')
    parser.add_argument('--regulized-type', type=str, default='noregu',
                        help='regulized type (default: LTMG) in EM, otherwise: noregu/LTMG/LTMG01')
    parser.add_argument('--reduction', type=str, default='sum',
                        help='reduction type: mean/sum, default(sum)')
    parser.add_argument('--prunetype', type=str, default='KNNgraphStatsSingleThread',
                        help='prune type, KNNgraphStats/KNNgraphML/KNNgraphStatsSingleThread (default: KNNgraphStatsSingleThread)')
    parser.add_argument('--precisionModel', type=str, default='Float',
                        help='Single Precision/Double precision: Float/Double (default:Float)')
    parser.add_argument('--coresUsage', type=str, default='1',
                        help='how many cores used: all/1/... (default:1)')
    parser.add_argument('--log-interval', type=int, default=100, metavar='N',
                        help='how many batches to wait before logging training status')
    parser.add_argument('--L1Para', type=float, default=1.0,
                        help='L1 regulized parameter (default: 0.001)')
    parser.add_argument('--L2Para', type=float, default=0.0,
                        help='L2 regulized parameter (default: 0.001)')
    parser.add_argument('--EMreguTag', action='store_true', default=False,
                        help='whether regu in EM process')

    args = parser.parse_args()
    args.cuda = not args.no_cuda and torch.cuda.is_available()
    

    checkargs(args)

    torch.manual_seed(args.seed)
    device = torch.device("cuda" if args.cuda else "cpu")
    # device = "cpu"
    print('Using device:'+str(device))

    if not args.coresUsage == 'all':
        torch.set_num_threads(int(args.coresUsage))

    kwargs = {'num_workers': 1, 'pin_memory': True} if args.cuda else {}
    start_time = time.time()

    print('---0:00:00---scRNA starts loading.')
    scData = scDataset(data)
    train_loader = DataLoader(
        scData, batch_size=args.batch_size, shuffle=False, **kwargs)
    print('---'+str(datetime.timedelta(seconds=int(time.time()-start_time))) +
          '---TrainLoader has been successfully prepared.')

    regulationMatrix = None
    if args.model == 'VAE':
        model = VAE(dim=scData.features.shape[1]).to(device)
    elif args.model == 'AE':
        model = AE(dim=scData.features.shape[1]).to(device)
    if args.precisionModel == 'Double':
        model = model.double()
    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    print('---'+str(datetime.timedelta(seconds=int(time.time()-start_time))) +
          '---Pytorch model ready.')
    start_time = time.time()

    print('Start training...')
    for epoch in range(0, feature_pretrain_epochs):
        recon, original, z = train(0, train_loader, model, optimizer, args, device, EMFlag=False)

    zOut = z.detach().cpu().numpy()
    print('zOut ready at ' + str(time.time()-start_time))
    ptstatus = model.state_dict()
    reconOri = recon.clone()
    reconOri = reconOri.detach().cpu().numpy()
    print('---'+str(datetime.timedelta(seconds=int(time.time()-start_time)))+'---Start Prune')
    adj, edgeList = generateAdj(zOut, graphType=args.prunetype, para=args.knn_distance+':'+str(
        k), adjTag=(args.useGAEembedding or args.useBothembedding))
    print('---'+str(datetime.timedelta(seconds=int(time.time() -
                                                   start_time)))+'---Prune Finished')

    G0 = nx.Graph()
    G0.add_weighted_edges_from(edgeList)
    nlG0 = nx.normalized_laplacian_matrix(G0)
    adjOld = nlG0

    print('---'+str(datetime.timedelta(seconds=int(time.time()-start_time))
                    )+"---EM process starts")

    for bigepoch in range(0, feature_train_epochs):
        print('---'+str(datetime.timedelta(seconds=int(time.time() -
                                                       start_time)))+'---Start %sth iteration.' % (bigepoch))
        recon, original, z = train(bigepoch, train_loader, model, optimizer, args, device, EMFlag=True)
        zOut = z.detach().cpu().numpy()

        print('---'+str(datetime.timedelta(seconds=int(time.time()-start_time)))+'---Start Prune')
        adj, edgeList = generateAdj(zOut, graphType=args.prunetype, para=args.knn_distance+':'+str(k), adjTag=(args.useGAEembedding or args.useBothembedding or (bigepoch == int(total_train_epoch)-1)))
        print('---'+str(datetime.timedelta(seconds=int(time.time() -
                                                       start_time)))+'---Prune Finished')
        Gc = nx.Graph()
        Gc.add_weighted_edges_from(edgeList)
        adjGc = nx.adjacency_matrix(Gc)
        adjNew = args.alpha*nlG0 + \
            (1-args.alpha) * adjGc/np.sum(adjGc, axis=0)
        print('---'+str(datetime.timedelta(seconds=int(time.time() -
                                                       start_time)))+'---New adj ready')
        graphChange = np.mean(abs(adjNew-adjOld))
        graphChangeThreshold = args.converge_graphratio * \
            np.mean(abs(nlG0))
        adjOld = adjNew
        if graphChange < graphChangeThreshold:
            print('Converge now!')
            break
    print('---'+str(datetime.timedelta(seconds=int(time.time()-start_time))
                    )+'---All iterations finished, start output results.')
    return adj.toarray()