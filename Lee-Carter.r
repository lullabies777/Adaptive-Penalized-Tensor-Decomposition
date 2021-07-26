## lee-carter
LC <- function(X, n.pc) {

    ## this function adopt codes from the demorgraphy package, X is the data matrix
    ## row is age group
    ## column is year

    ax <- apply(X, 1, mean, na.rm = TRUE) # ax is mean of X by column

    cX <- sweep(X, 1, ax) ## make sum(kt) = 0

    svd.mx <- svd(cX)



    kt <- matrix(NA, nrow = n.pc, ncol = ncol(X))

    bx <- matrix(NA, nrow = nrow(X), ncol = n.pc)

    for (i in seq(1, n.pc)) {
        sumu <- sum(svd.mx$u[, i])

        ## to make sum(bx)=1

        bx[, i] <- svd.mx$u[, i] / sumu

        kt[i, ] <- svd.mx$d[i] * svd.mx$v[, i] * sumu
    }

    fitted <- sweep(bx %*% kt, 1, -ax)
    svd_d <- svd(X)$d
    return(list(bx = bx, kt = kt, rowmean = ax, fitted = fitted, d = svd_d))
}

# predict the Lee-Carter model, Using RWWD (random walk with drift) model

LC.predict <- function(X, n.pred, n.pc) {
    tmp <- LC(X, n.pc = n.pc)

    kt.pred <- matrix(NA, nrow = n.pc, ncol = n.pred)

    for (i in seq(1, n.pc)) {
        rw.fit <- forecast::rwf(tmp$kt[i, ], drift = TRUE)

        kt.pred[i, ] <- tmp$kt[i, length(tmp$kt[i, ])] + seq(1, n.pred) * rw.fit$model$par$drift
    }

    return(sweep(tmp$bx %*% kt.pred, 1, -tmp$rowmean))
}

LC.forecast <- function(tnsr, n.pred, n.pc, year) {
    cs1 <- LC.predict(tnsr[1, , 1:(dim(tnsr)[3] - year)]@data, n.pred, n.pc)
    c1 <- LC.predict(tnsr[2, , 1:(dim(tnsr)[3] - year)]@data, n.pred, n.pc)
    rs1 <- LC.predict(tnsr[3, , 1:(dim(tnsr)[3] - year)]@data, n.pred, n.pc)
    inf1 <- LC.predict(tnsr[4, , 1:(dim(tnsr)[3] - year)]@data, n.pred, n.pc)
    ext1 <- LC.predict(tnsr[5, , 1:(dim(tnsr)[3] - year)]@data, n.pred, n.pc)
    oth1 <- LC.predict(tnsr[6, , 1:(dim(tnsr)[3] - year)]@data, n.pred, n.pc)
    m_tensor <- array(NA, dim = c(6, 19, n.pred))

    m_tensor[1, , ] <- cs1
    m_tensor[2, , ] <- c1
    m_tensor[3, , ] <- rs1
    m_tensor[4, , ] <- inf1
    m_tensor[5, , ] <- ext1
    m_tensor[6, , ] <- oth1
    ###

    return(as.tensor(m_tensor))
}

LC.fit <- function(tnsr, n.pc, year) {
    cs1 <- LC(tnsr[1, , 1:(dim(tnsr)[3] - year)]@data, n.pc)$fitted
    c1 <- LC(tnsr[2, , 1:(dim(tnsr)[3] - year)]@data, n.pc)$fitted
    rs1 <- LC(tnsr[3, , 1:(dim(tnsr)[3] - year)]@data, n.pc)$fitted
    inf1 <- LC(tnsr[4, , 1:(dim(tnsr)[3] - year)]@data, n.pc)$fitted
    ext1 <- LC(tnsr[5, , 1:(dim(tnsr)[3] - year)]@data, n.pc)$fitted
    oth1 <- LC(tnsr[6, , 1:(dim(tnsr)[3] - year)]@data, n.pc)$fitted
    m_tensor <- array(NA, dim = c(6, 19, dim(tnsr)[3] - year))

    m_tensor[1, , ] <- cs1
    m_tensor[2, , ] <- c1
    m_tensor[3, , ] <- rs1
    m_tensor[4, , ] <- inf1
    m_tensor[5, , ] <- ext1
    m_tensor[6, , ] <- oth1
    ###

    return(as.tensor(m_tensor))
}


################# cross_validation

LC.cv <- function(tnsr, year, k = 5, leave_test = TRUE) {
    if (leave_test) {
        tnsr <- tnsr[, , 1:(dim(tnsr)[3] - year)]
    }
    ### training data
    train_ls <- list()
    for (i in 1:k) {
        train_ls[[i]] <- tnsr[, , 1:(dim(tnsr)[3] - i + 1 - year)]
    }
    ### validating data
    validate_ls <- list()
    for (i in 1:k) {
        validate_ls[[i]] <- tnsr[, , (dim(tnsr)[3] - i + 2 - year):(dim(tnsr)[3] - i + 1)]
    }
    ### record minMSE
    minMSE <- 999999
    mincomponent <- 1
    for (i in 1:19) {
        ### record each component statistics
        temp <- NULL
        for (j in 1:k) {
            temp <- c(temp, fnorm(lcforecast(train_ls[[j]], n.pred = year, n.pc = i, year = 0) - validate_ls[[j]]))
        }
        if (mean(temp) < minMSE) {
            minMSE <- mean(temp)
            mincomponent <- i
        }
    }
    return(list(optimal_component = mincomponent, min_MSE = minMSE))
}