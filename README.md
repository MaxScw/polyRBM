# polyRBM
## Implementation of a Restricted Boltzmann Machine for the training on polymer conformation data


- the base-level source code for the RBM network, training and evaluation is found in <strong>analysisUtility.py</strong> and <strong>polyRBM.py</strong>
    - <strong>analysisUtility</strong> contains a collection of evaluation functions for calculating observabels of polymer chains
    - <strong>polyRBM.py</strong> contains the pyTorch implementation of the neural network + transformation and training routines

- <strong>train.ipynb</strong> sums up training procedure and quick evaluation of results with in-line plots for convenience
- <strong>results_plot.ipynb</strong> contains the setup for re-creating all types of plots, which were used in the work
