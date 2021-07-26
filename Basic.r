## --------------------
## --------------------
## --------------------
## This file contains some basic functions for implementing Adaptive Penalized Tensor Decomposition.
## Most function are only suitable for three dimension tensor. For higher-order tensor, codes should be modified adaptively.
## --------------------
## --------------------
## --------------------


## --------------------
## inital values
## --------------------
inits <- function(k) {
    return(t(as.matrix(rnorm(k, 0, 1))))
}

## --------------------
## tensor product
## --------------------
product <- function(X, ind, u, v, w) {
    u <- t(as.matrix(u))
    v <- t(as.matrix(v))
    w <- t(as.matrix(w))


    if (length(ind) == 2) {
        lizt <- list("mat2" = u, "mat3" = v)
        u <- ttl(X, lizt, ms = ind)
        return(as.vector(u@data))
    }

    lizt <- list("mat" = u, "mat2" = v, "mat3" = w)
    u <- ttl(X, lizt, ms = ind)
    return(as.vector(u@data))
}

## --------------------
## generate inital value for tensor decompostion if initial value is null.
## --------------------
gen_startvals <- function(num_factors, tnsr) {
    u <- matrix(inits(dim(tnsr)[1] * num_factors), nrow = num_factors)
    v <- matrix(inits(dim(tnsr)[2] * num_factors), nrow = num_factors)
    w <- matrix(inits(dim(tnsr)[3] * num_factors), nrow = num_factors)
    for (i in 1:num_factors) {
        u[i, ] <- u[i, ] / norm(as.matrix(u[i, ]), "F")
    }
    for (i in 1:num_factors) {
        v[i, ] <- v[i, ] / norm(as.matrix(v[i, ]), "F")
    }
    for (i in 1:num_factors) {
        w[i, ] <- w[i, ] / norm(as.matrix(w[i, ]), "F")
    }

    return(list(u = u, v = v, w = w))
}

## --------------------
## generate mortality tensor.
## --------------------
generate_mortality_tensor <- function(x, log = TRUE, std = TRUE) {
    x <- as.data.table(x)
    len <- ncol(x)
    circu_m <- data.matrix(x[Cause == "Circulatory system"][, (len - 18):len])
    c_m <- data.matrix(x[Cause == "Cancer"][, (len - 18):len])
    r_m <- data.matrix(x[Cause == "Respiratory system"][, (len - 18):len])
    e_m <- data.matrix(x[Cause == "External causes"][, (len - 18):len])
    i_m <- data.matrix(x[Cause == "Infectious and parasitic diseases"][, (len - 18):len])
    o_m <- data.matrix(x[Cause == "Others"][, (len - 18):len])
    ## define the tensor and apply log transformation
    if (log) {
        m_tensor <- array(NA, dim = c(6, 19, (nrow(x) / 6)))
        m_tensor[1, , ] <- t(log(circu_m))
        m_tensor[2, , ] <- t(log(c_m))
        m_tensor[3, , ] <- t(log(r_m))
        m_tensor[4, , ] <- t(log(e_m))
        m_tensor[5, , ] <- t(log(i_m))
        m_tensor[6, , ] <- t(log(o_m))
    } else {
        m_tensor <- array(NA, dim = c(6, 19, (nrow(x) / 6)))
        m_tensor[1, , ] <- t(circu_m)
        m_tensor[2, , ] <- t(c_m)
        m_tensor[3, , ] <- t(r_m)
        m_tensor[4, , ] <- t(e_m)
        m_tensor[5, , ] <- t(i_m)
        m_tensor[6, , ] <- t(o_m)
    }
    ###
    if (std) {
        tensor1 <- (m_tensor - array(mean(m_tensor), dim = dim(m_tensor))) / sd(m_tensor)
        return(as.tensor(tensor1))
    } else {
        tensor1 <- m_tensor
        return(as.tensor(tensor1))
    }
}