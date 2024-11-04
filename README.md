# Empirical Bayes based on genGK to estimate hyperparameters
This repository contains MATLAB files for an empirical Bayes method to 
estimate hyperparameters using an approach based on the generalized 
Golub-Kahan (genGK) bidiagonalization.  The codes accompany the paper: 

"Efficient iterative methods for hyperparameter estimation in large-scale 
linear inverse problems" 
- Hall-Hooper, Saibaba, Chung, and Miller, ACOM, 2024

## Project Description
    We implement an empirical Bayes (EB) method to estimate hyperparameters 
    that maximize the marginal posterior, i.e., the probability density of 
    the hyperparameters conditioned on the data.
    
    For problems where the computation of the square root and inverse of 
    prior covariance matrices are not feasible, we use an approach based on
    the generalized Golub-Kahan bidiagonalization to approximate the 
    marginal posterior and seek hyperparameters that minimize the 
    approximate marginal posterior. 


## Installation 
### Software language

       MATLAB 9.14 (R2023a)
       
### Requirements

# Requirements
These codes require the following packages:
         
         Regularization Tools package: Hansen. Regularization tools: A
             package for analysis and solution of discrete ill-posed 
             problems. Numerical Algorithms, 1994.
             
         genHyBR: generalized hybrid iterative methods
             by Julianne Chung and Arvind K. Saibaba
             https://github.com/juliannechung/genHyBR
         

## How to Use
See Contents.m

### Contributors
        Khalil A. Hall-Hooper
        Department of Mathematics, North Carolina State University

        Arvind K. Saibaba, 
        Department of Mathematics, North Carolina State University
        
        Julianne Chung, 
        Department of Mathematics, Emory University
        
        Scot M. Miller, 
        Department of Environmental Health and Engineering, Johns Hopkins University
        
        
## Licensing

If you use this codes, you *must* cite the original authors:

       [1] Hall-Hooper et al. "Efficient iterative methods for 
            hyperparameter estimation in large-scale linear 
            inverse problems". ACOM, 2024.


[MIT](LICENSE)

## Acknowledgement

This work was partially supported by the National Science Foundation under 
grants DMS-2208294, DMS-2341843, DMS-2026830, and DMS-2026835. 
Any opinions, findings, conclusions or recommendations expressed in this
material are those of the author(s) and do not necessarily reflect the 
views of the National Science Foundation.

