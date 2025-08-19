# mbSparse

## Installation

Before to execute mbSparse, it is necessary to install the following packages:

* torch
* network
* pandas
* joblib
* gensim
* h5py

## Basic Usage

### Example
Run the example script:

```shell
python example_impute.py
```
This script performs imputation on a sample microbiome dataset using the following function call:
```python
data_impute = Impute.imputation(
    scfile=data,
    k=2,
    feature_train_epochs=20,
    feature_pretrain_epochs=15,
    cave_train_epochs=1,
    unnormalized=False,
    samples_are_rows=False
)
```
Function parameters:

* **scfile:** Input data matrix.

* **k:** Number of neighbors used to construct the sample correlation graph.

* **feature_train_epochs:** Number of training epochs for the feature autoencoder.

* **feature_pretrain_epochs:** Number of pretraining epochs for the feature autoencoder.

* **cave_train_epochs:** Number of training epochs for the conditional variational autoencoder.

* **unnormalized:** Whether to normalize the input matrix before imputation.

* **samples_are_rows**: If True, samples are assumed to be in rows; if False, samples are in columns.

**Note:** Experiments in the paper were conducted with samples_are_rows=False.

