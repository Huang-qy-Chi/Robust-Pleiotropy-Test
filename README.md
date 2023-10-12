# Robust-Pleiotropy-Test
Copyright © 2023 Q. Huang. All rights reserved.
This repository contains the codes of functions in the paper "A Multivariate Robust Integrated Test for Genetic Pleiotropy", programmed by the author Huang Qiyue. Here I may tell how to use the functions.
After loading all functions, and dividing the centralized data set into the multivariate Y and the univariate X as matrix, please first use the function MLTS to give the initial values of B0 and Sigma0, then use Y, X, B0, Sigma0
and given significance level alpha (such as 0.05) to start the MM-estimate (function MM_estimate1X or MM_estimate2X) to get the MM-estimator, then use the integrated test (EHT_PFCX or EHT_PFC1X).
