# Robust-Pleiotropy-Test

This repository contains the code of functions in the paper "A Multivariate Robust Integrated Test for Genetic Pleiotropy" (link: https://dx.doi.org/10.4310/SII.250522054924). 
Here, I show how to use the functions:

1) Loading all functions and dividing the centralized data set into the multivariate Y and the univariate X as a matrix,

2) Use the function MLTS to give the initial values of B0 and Sigma0, 

3) Use Y, X, B0, Sigma0, and given significance level alpha (such as 0.05) to start the MM-estimate (function MM_estimate1X or MM_estimate2X) to get the MM-estimator,
   
4) Use the integrated test (EHT_PFCX or EHT_PFC1X) to detect the existence of Pleiotropy.

One can directly attach to files "EHTC3.R" and "EHTC10.R"; they will automatically source all functions.

Copyright © 2025 Q. Huang. All rights reserved. 

