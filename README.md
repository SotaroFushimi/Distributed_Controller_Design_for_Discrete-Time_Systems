# Numerical Example for "Distributed Controller Design for Discrete-Time Systems Via the Integration of Extended LMI and Clique-Wise Decomposition"

This repository contains the code used for the numerical example in the following paper:

**Paper Title:**  
Distributed Controller Design for Discrete-Time Systems Via the Integration of Extended LMI and Clique-Wise Decomposition  
**Link:** [https://doi.org/10.48550/arXiv.2409.07666](https://doi.org/10.48550/arXiv.2409.07666)

## How to Use

### Stability Check
For stability analysis, run the following script:
`Stbl_script.m`  
This script will return the number of cases where a stabilizer is successfully derived.

### H∞ Performance Analysis
For H∞ performance analysis, run the following script:
`Hinfty_script.m`

### Modifying Parameters
Parameters for the simulations can be modified in the `parameters.m` file.

### Parameters for the paper
Parameters and results used in the paper is stored at `sim_result_paper` folder.

## Continuous-Time Systems
For continuous-time systems, refer to the following paper:  
[https://doi.org/10.48550/arXiv.2404.04576](https://doi.org/10.48550/arXiv.2404.04576)

The code for continuous-time systems is available at:  
[https://github.com/WatanabeYuto/LMI-Based_Distributed_Controller_Design_with_Non-Block-Diagonal_Lyapunov_Functions.git](https://github.com/WatanabeYuto/LMI-Based_Distributed_Controller_Design_with_Non-Block-Diagonal_Lyapunov_Functions.git)
