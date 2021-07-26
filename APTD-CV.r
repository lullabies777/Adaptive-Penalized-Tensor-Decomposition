## --------------------
## --------------------
## --------------------
## This file contains functions mainly for cross validation of Adaptive Penalized Tensor Decomposition.
## Most function are only suitable for three dimension tensor. For higher-order tensor, codes should be modified adaptively.
## --------------------
## --------------------
## --------------------


## --------------------
## Initialization for cross validation.
## --------------------
prepare_adaptive_adj <- function(tnsr, num_components, Niter = 1500, tol = 1e-05, year = 5, trace = FALSE, ku, kv, kw, splits = 10) {
    original <- tnsr
    Test <- tnsr[, , (dim(tnsr)[3] - year + 1):dim(tnsr)[3]]
    V1 <- list()
    S1 <- list()
    In1 <- list()
    mean1 <- list()
    sd1 <- list()
    tnsr <- tnsr[, , 1:(dim(tnsr)[3] - year)]
    for (m in 1:splits) {
        message("Training on years: ", paste(1:(dim(tnsr)[3] - splits - year + m), collapse = "+"))
        message("Test on years: ", paste((dim(tnsr)[3] - splits - year + m + 1):((dim(tnsr)[3] - splits - year + m + 1) + year - 1), collapse = "+"))

        temp_s <- tnsr@data[, , 1:(dim(tnsr)[3] - splits - year + m)]
        S1[[m]] <- as.tensor((temp_s - array(mean(temp_s), dim = dim(temp_s))) / sd(temp_s))
        mean1[[m]] <- c(mean(temp_s))
        sd1[[m]] <- c(sd(temp_s))

        V1[[m]] <- as.tensor(tnsr@data[, , (dim(tnsr)[3] - splits - year + m + 1):((dim(tnsr)[3] - splits - year + m + 1) + year - 1)])

        make_init <- gen_startvals(num_components, S1[[m]])

        message("Fitting unpenalized tensor...")
        In1[[m]] <- multiple_tf_simple(
            tnsr = S1[[m]], d_hat = NULL, num_components = num_components,
            u_init = make_init$u, v_init = make_init$v, w_init = make_init$w,
            Niter = Niter, tol = tol,
            c1 = rep(0, num_components), c2 = rep(0, num_components), c3 = rep(0, num_components),
            ku = ku, kv = kv, kw = kw, trace = trace
        )
        message("Fitting completed.")
    }

    ## generate D
    Dut <- list()
    Dvt <- list()
    Dwt <- list()
    for (m in 1:splits) {
        Dut[[m]] <- create_D(In1[[m]]$u, ord = (ku + 1))
        Dvt[[m]] <- create_D(In1[[m]]$v, ord = (kv + 1))
        Dwt[[m]] <- create_D(In1[[m]]$w, ord = (kw + 1))
    }


    return(list(
        Dut = Dut, Dvt = Dvt, Dwt = Dwt, In1 = In1, S1 = S1, V1 = V1, year = year,
        num_components = num_components, mean1 = mean1, sd1 = sd1, Test = Test, original = original
    ))
}


## --------------------
## Function for cross validation.
## --------------------
multiple_cv_adaptive_adj <- function(init, grid, ku = NULL, kv = NULL, kw = NULL, tol = 5e-05, Niter_small = 1500, splits = 10) {

    # trian data
    # input tenor should be training +validation tensor
    # third dimension should be year dimension

    tstart <- Sys.time()
    S1 <- init$S1
    V1 <- init$V1
    Dut <- init$Dut
    Dvt <- init$Dvt
    Dwt <- init$Dwt
    num_components <- init$num_components
    year <- init$year
    In1 <- init$In1
    mean1 <- init$mean1
    sd1 <- init$sd1

    sfInit(parallel = TRUE, cpus = 6)
    # load
    sfLibrary(rTensor)
    sfLibrary(genlasso)
    sfLibrary(forecast)
    sfLibrary(mgcv) # gam
    sfLibrary(Hmisc)
    sfExport(
        "multiple_D_adaptive", "PTD_D", "In1", "tol", "product", "S1", "cal_pmse_gam_adj",
        "cal_pmse_arima_adj", "cal_pmse_linearextra_adj", "num_components", "V1", "year", "Niter_small", "Dut", "Dvt", "Dwt", "mean1", "sd1"
    )

    best <- 0
    meanPMSE_gam <- meanPMSE_linextrap <- meanPMSE_arima <- NULL # record each combination's mean of the prediction RMSEs (averge of 10 )
    allPMSE_gam <- allPMSE_linextrap <- allPMSE_arima <- NULL

    for (i in 1:nrow(grid)) {
        message("Tuning parameters: ", paste(grid[i, ], collapse = "/"))
        sfExport("i", "grid")

        parallel1 <- function(m, num_components, In1, grid, Niter_small, tol, Dut, Dvt, Dwt) {
            temp <- multiple_D_adaptive(
                tnsr = S1[[m]], num_components = num_components, d_hat = In1[[m]]$d,
                u_init = In1[[m]]$u, v_init = In1[[m]]$v, w_init = In1[[m]]$w,
                c1 = rep(grid[i, 1], num_components), c2 = rep(grid[i, 2], num_components), c3 = rep(grid[i, 3], num_components),
                Niter = Niter_small, tol = tol, Du = Dut[[m]], Dv = Dvt[[m]], Dw = Dwt[[m]], trace = FALSE
            )
        }

        tempresult <- sfLapply(1:splits, parallel1, num_components, In1, grid, Niter_small, tol, Dut, Dvt, Dwt)
        sfExport("tempresult")

        temprmse_gam <- as.numeric(sfLapply(1:splits, cal_pmse_gam_adj, year = year))
        temprmse_linextrap <- as.numeric(sfLapply(1:splits, cal_pmse_linearextra_adj, year = year))
        temprmse_arima <- as.numeric(sfLapply(1:splits, cal_pmse_arima_adj, year = year))
        allPMSE_gam <- rbind(allPMSE_gam, temprmse_gam)
        allPMSE_arima <- rbind(allPMSE_arima, temprmse_arima)
        allPMSE_linextrap <- rbind(allPMSE_linextrap, temprmse_linextrap)

        # bias and variance
        meanPMSE_gam <- c(meanPMSE_gam, mean(temprmse_gam))
        meanPMSE_linextrap <- c(meanPMSE_linextrap, mean(temprmse_linextrap))
        meanPMSE_arima <- c(meanPMSE_arima, mean(temprmse_arima))

        # print(meanPMSE_gam)
        # print(meanPMSE_arima)
        # print(meanPMSE_linextrap)
    }

    ## stop parallel
    sfStop()
    tend <- Sys.time()
    message("Fitting completed!")

    colnames(grid) <- c("U", "V", "W")
    performance <- data.frame(grid,
        meanPMSE_gam = meanPMSE_gam, meanPMSE_linextrap = meanPMSE_linextrap, meanPMSE_arima = meanPMSE_arima,
        mean_diffunpenPMSE_gam = meanPMSE_gam - meanPMSE_gam[1], mean_diffunpenPMSE_linextrap = meanPMSE_linextrap - meanPMSE_linextrap[1],
        mean_diffunpenPMSE_arima = meanPMSE_arima - meanPMSE_arima[1]
    )
    performance$sd_diffunpenPMSE_gam <- apply(t(t(allPMSE_gam) - allPMSE_gam[1, ]), 1, sd)
    performance$sd_diffunpenPMSE_arima <- apply(t(t(allPMSE_arima) - allPMSE_arima[1, ]), 1, sd)
    performance$sd_diffunpenPMSE_linextrap <- apply(t(t(allPMSE_linextrap) - allPMSE_linextrap[1, ]), 1, sd)

    out <- list(
        best_parameters_gam = grid[which.min(performance$meanPMSE_gam), ],
        best_parameters_linextrap = grid[which.min(performance$meanPMSE_linextrap), ],
        best_parameters_arima = grid[which.min(performance$meanPMSE_arima), ],
        performance = performance, init = init, time = tend - tstart
    )
    return(out)
}