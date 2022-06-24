# polyRBM
## Implementation of a Restricted Boltzmann Machine for the training on polymer conformation data


- the base-level source code for the RBM network, training and evaluation is found in 'analysisUtility.py' and 'polyRBM.py'
>    - 'analysisUtility' contains a collection of evaluation functions for calculating observabels of polymer chains
>    - 'polyRBM.py' contains the pyTorch implementation of the neural network + transformation and training routines

- the C++-source code of the slithering-snake Monte-Carlo alogrithm used in the work is contained in '/polysim'
    - insctructions for building the polysim executables with Make/CMake:
     
>>    - on console, go to '/polysim'
         - execute 'cmake .'
         - finish building by executing 'make' in the same directory
     - usage of simulation program: execute 'polysim $N $is2D $isPhantom $nRuns $steps $MCSperRun $isRing $nonReturn'
         - integer variables: N = chain length, nRuns = number of runs to perform, steps = intermediate MCS between two saved conformations,
         MCSperRun = number of MCS to perform for one run
         - boolean variables: is2D = simulation of 2D chains?, isPhantom = simulate excluded volume?, isRing = only sample ring conformations?
         nonReturn = simulate exclusion of back-folding?

- 'train.ipynb' sums up training procedure and quick evaluation of results with in-line plots for convenience
- 'results_plot.ipynb' contains the setup for re-creating all types of plots, which were used in the work
