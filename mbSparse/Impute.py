import pandas as pd
import numpy as np
from .vae.vae_train import *
import torch
from .cave.cvae_train import generated_features

def normalize_matrix(X):
    scale = np.sum(X, axis=0) / 10 ** 6
    X = X / scale
    X = np.log10(X + 1)
    return X

def imputation(scfile, k, feature_train_epochs, feature_pretrain_epochs, cave_train_epochs, unnormalized=True):  
    
    np_scfile = np.array(scfile)
    if unnormalized == True: 
        np_scfile = normalize_matrix(np_scfile)
    imputation_file=np.array(np_scfile)
    
    adj = generateAdj_Cluster(np_scfile, k, feature_train_epochs, feature_pretrain_epochs)

    features = np.array(np_scfile)    
    features_min = features.min(axis=0)
    features_max = features.max(axis=0)
    features = (features - features_min) / (features_max - features_min)
    features = features.T
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    cuda = torch.cuda.is_available()
    augmented_features = generated_features(cuda, device, adj, features, cave_train_epochs)
    augmented_features = augmented_features.detach().cpu().numpy().T
    augmented_features = augmented_features * (features_max - features_min) + features_min
    
    if len(scfile.columns) > 500:
        for k in range(adj.shape[1]):
            scfile_col = np.where(adj[:,k]==1)[0]
            imputation_file_row = np.where(imputation_file[:,k]==0)[0]
            imputation_file[np.where(imputation_file[:,k]==0)[0],k] = np.mean(np.add(np_scfile[:,scfile_col][imputation_file_row,:],augmented_features[:,scfile_col][imputation_file_row,:]),axis=1)
    else:
        for k in range(adj.shape[1]):
            scfile_col = np.where(adj[:,k]==1)[0]
            imputation_file_row = np.where(imputation_file[:,k]==0)[0]
            imputation_file[np.where(imputation_file[:,k]==0)[0],k] = np.median(np.add(np_scfile[:,scfile_col][imputation_file_row,:],augmented_features[:,scfile_col][imputation_file_row,:]),axis=1) 
            
    df_imputation=pd.DataFrame(imputation_file)
    df_imputation.columns=scfile.columns.tolist()
    df_imputation.index=scfile.index.tolist()
    return df_imputation
