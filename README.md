# AlpineFloraLJ
Third party python modules and dataset for crosschecking alpine flora phylogenetics

The modules along with the requrired input and resulting data have been provided as comprehensive dataset by Ding, Wen-Na; Ree, Richard H.; Spicer, Robert A.; Xing, Yao-Wu (2020), Ancient orogenic and monsoon-driven assembly of the world's richest temperate alpine flora, Dryad, Dataset, https://doi.org/10.5061/dryad.k6djh9w3s.

## Low dataset quality

1. The dataset does not match the FAIR principles. While researchers can do without a formal knowledge representation language (I1/I2) in this particular case, the **missing rich metadata in basically all dimensions (F2)** is a poverty certificate for scientific quality control.

1. The quality of the computing code is low. The code as packed in the dataset repository is not functional (see details below). In addition and while there is some documentation in the code, an overall package organisation is missing.

## Required changes to get the modules running

The module *bitstates.py* had to be changed to work with *02_parse_classe.py*. Otherwise, an error is caused by *name = model.areas(a).name* (line 63 in *02_parse_classe.py*)  Therefore, an updated version from of *bitstates.py* from https://github.com/rhr/bitstates/ (commit d8f0725) was used. The change is part of commit 0ff9cb08d9f75e8f55f2cce0f00a643346b6a699 in this repository.

Another error in *02_parse_classe.py* caused by *daughters = set(chain.from_iterable(model.clado_events_by_idx[ancidx]))* (line 87/88 in *02_parse_classe.py*) could not be resolved unitl know. A first attemp to solve it is to recompute the input dataset which is currently done.

