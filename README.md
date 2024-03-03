This file contains 6 scripts for the simulation study of the paper on hierarchical network autocorrelation models for binary-valued network and weighted-edge network. 
Each script forms the standalone file with data generation and model estimation using Bayesian approach.
HNAM_uniform_1, HNAM_uniform_2 and HNAM_uniform_T are for the base model in the paper with different prior distributions of rho. 
In particular, HNAM_uniform_1 uses a uniform prior of rho with a range from -1 to 1; HNAM_uniform_2 uses a uniform prior of rho with a range from 1/lamda_min to 1/lamda_max; HNAM_uniform_T uses a transformed uniform prior of rho.
Extended_HNAM_uniform_1, Extended_HNAM_uniform_2 and Extended_HNAM_uniform_T are for the extended model in the paper with different prior distributions of rho.
Similarly, Extended_HNAM_uniform_1 uses a uniform prior of rho with a range from -1 to 1; Extended_HNAM_uniform_2 uses a uniform prior of rho with a range from 1/lamda_min to 1/lamda_max; Extended_HNAM_uniform_T uses a transformed uniform prior of rho.
For each script, the output includes the posterior median estimator of coefficients, error terms and rho (and alpha for the extended model). In adiition, it also includes the mean squared error of rho, 95% coverage rate of rho, 95% equal-tailed credible interval of rho and average width of credible intrerval of rho and those of alpha for the extended model.  

Supplementary materials.pdf contains simulation results and derivation of conditional posteriors.
