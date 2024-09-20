import sys
import gc
import numpy as np
import torch
import torch.nn.functional as F
import torch.optim as optim
from tqdm import tqdm, trange
from torch_geometric.nn import GCNConv
from .cvae_models import VAE
from torch.utils.data import TensorDataset, DataLoader, RandomSampler
import copy
import random

exc_path = sys.path[0]

def loss_fn(recon_x, x, mean, log_var):
    BCE = torch.nn.functional.binary_cross_entropy(recon_x, x, reduction='sum')
    KLD = -0.5 * torch.sum(1 + log_var - mean.pow(2) - log_var.exp())
    return (BCE + KLD) / x.size(0)

def generated_features(cuda, device, adj, features, cave_train_epochs=3, batch_size=64, latent_size=10, cave_lr=0.01, conditional=True, seed=100):
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)
    x_list, c_list = [], [] 
    for i in trange(adj.shape[0]):
        x = features[adj[i].nonzero()]
        c = np.tile(features[i], (x.shape[0], 1))
        x_list.append(x)
        c_list.append(c)
    features_x = np.vstack(x_list)
    features_c = np.vstack(c_list)
    del x_list
    del c_list
    gc.collect()
    
    features_x = torch.tensor(features_x, dtype=torch.float32)
    features_c = torch.tensor(features_c, dtype=torch.float32)

    cvae_features = torch.tensor(features, dtype=torch.float32)

    cvae_dataset = TensorDataset(features_x, features_c)
    cvae_dataset_sampler = RandomSampler(cvae_dataset)
    cvae_dataset_dataloader = DataLoader(cvae_dataset, sampler=cvae_dataset_sampler, batch_size=batch_size)
    
    cvae = VAE(encoder_layer_sizes=[features.shape[1], 256], 
               latent_size=latent_size, 
               decoder_layer_sizes=[256, features.shape[1]],
               conditional=conditional, 
               conditional_size=features.shape[1])
    cvae_optimizer = optim.Adam(cvae.parameters(), lr=cave_lr)

    if cuda:
        cvae_features = cvae_features.to(device)
        cvae.to(device)

    for _ in trange(cave_train_epochs, desc='Run CVAE Train'):
        for _, (x, c) in enumerate(tqdm(cvae_dataset_dataloader)):
            cvae.train()
            x, c = x.to(device), c.to(device)
            if conditional:
                recon_x, mean, log_var, _ = cvae(x, c)
            else:
                recon_x, mean, log_var, _ = cvae(x)
            cvae_loss = loss_fn(recon_x, x, mean, log_var)
            cvae_optimizer.zero_grad()
            cvae_loss.backward()
            cvae_optimizer.step()
            z = torch.randn([cvae_features.size(0), latent_size]).to(device)
            augmented_features = cvae.inference(z, cvae_features)
    return augmented_features
