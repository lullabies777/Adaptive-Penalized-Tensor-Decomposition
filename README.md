# Adaptive Penalized Tensor Decomposition

by
Xuanming Zhang,
Fei Huang,
Francis K.C. Hui,
Steven Haberman



## Abstract

Cause-of-death mortality modeling and forecasting is an important topic in demography and actuarial science, as it can provide valuable insights into the risks and factors determining future mortality rates. In this paper, we propose a novel predictive approach for cause-of-death mortality forecasting based on an adaptive penalized tensor decomposition (ADAPT).  The new method jointly models the three dimensions (cause, age, and year) of the data, using adaptively weighted penalty matrices to overcome the computational burden of having to select a large number of tuning parameters when multiple factors are involved. ADAPT can be coupled with a variety of methods for forecasting, and we examine three different ones (random walk with drift, linear extrapolation, and smoothing) for extrapolating the estimated year factors.
Based on an application to US male cause-of-death mortality data, we demonstrate that ADAPT can provide superior out-of-sample predictive performance compared to several existing models (especially) when it comes to mid- and long-term forecasting.

## Software implementation
 This project contains the R source code to implement Adaptive Penalized Tensor Decomposition. However, it is worth noting that, input tensors are all third-order tensors by default and cross-validation are mainly coded for temporal structure tensor. For higher-order tensors or other dataset, the corresponding source code should be modified adaptively.

## Getting the code

You can download a copy of all the files in this repository by cloning the
[git](https://git-scm.com/) repository:

    git clone https://github.com/lullabies777/Adaptive-Penalized-Tensor-Decomposition.git

or [download a zip archive](https://github.com/lullabies777/Adaptive-Penalized-Tensor-Decomposition/archive/refs/heads/main.zip).

# Components

## `APTD-CV.r`
This file contains functions mainly for cross validation of Adaptive Penalized Tensor Decomposition.

## `APTD.r`
This file contains functions for Adaptive Penalized Tensor Decomposition.

## `Basic.r`
This file contains some basic functions for implementing Adaptive Penalized Tensor Decomposition.

## `Lee-Carter.r`
This file contains the functions for implementing the Lee-Carter Model. 

## `MLR.r`
This file contains the functions for implementing the MLR Model. These codes only implement the MLR model for short-term forecast(5 year). For mid-term, long-term or other dataset, they should be modified adaptively.