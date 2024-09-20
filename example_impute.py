from mbSparse import Impute
data = pd.read_csv("./data.csv",index_col=0)
data_impute = Impute.imputation(scfile=data, k=2, feature_train_epochs=20, feature_pretrain_epochs=15, cave_train_epochs=1, unnormalized=False)
data_impute.to_csv("./data_impute.csv")
 
