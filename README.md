# polyRBM
## Implementation of a Restricted Boltzmann Machine for the training on polymer conformation data


- the base-level source code for the RBM network, training and evaluation is found in 'analysisUtility.py' and 'polyRBM.py'
    - 'analysisUtility' contains a collection of evaluation functions for calculating observabels of polymer chains
    - 'polyRBM.py' contains the pyTorch implementation of the neural network + transformation and training routines

- 'train.ipynb' sums up training procedure and quick evaluation of results with in-line plots for convenience
- 'results_plot.ipynb' contains the setup for re-creating all types of plots, which were used in the work
