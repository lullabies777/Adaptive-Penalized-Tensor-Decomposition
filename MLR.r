## --------------------
## --------------------
## --------------------
## These codes implement the MLR for short-term forecast(5 year). For mid-term, long-term or other dataset, these should be modified properly.
## --------------------
## --------------------
## --------------------
rm(list = ls())
library(nnet)
library(tidyverse)

library(reshape2)
library(data.table)
library(rTensor)
set.seed(082020)

## --------------------
## Data attributes
## --------------------
popm <- fread("population_male.csv")
usam <- fread("mortality rate for male with comparability ratio.csv")
colnames(popm)[4:22] <- c(
    "infant", "1-4",
    "5-9", "10-14",
    "15-19", "20-24",
    "25-29", "30-34",
    "35-39", "40-44",
    "45-49", "50-54",
    "55-59", "60-64",
    "65-69", "70-74",
    "75-79", "80-84",
    "85+"
)
colnames(usam)[4:22] <- c(
    "infant", "1-4",
    "5-9", "10-14",
    "15-19", "20-24",
    "25-29", "30-34",
    "35-39", "40-44",
    "45-49", "50-54",
    "55-59", "60-64",
    "65-69", "70-74",
    "75-79", "80-84",
    "85+"
)

###### process population data
popn <- reshape2::melt(popm[, c(2, 4:22)], id.vars = c("Year"))
colnames(popn) <- c("Year", "agegroups", "pop")
popn <- setorder(popn, "Year")
popn$Year <- popn$Year - 1949


###### process mortality rate
cir <- usam[Cause == "Circulatory system"][, c(2, 4:22)]
cir1 <- melt(cir, id.vars = "Year")
setorder(cir1, "Year", "variable", "value")
colnames(cir1) <- c("Year", "agegroups", "Circulatory system")
cir1$Year <- cir1$Year - 1949

can <- usam[Cause == "Cancer"][, c(2, 4:22)]
can1 <- melt(can, id.vars = "Year")
setorder(can1, "Year", "variable", "value")
colnames(can1) <- c("Year", "agegroups", "Cancer")
can1$Year <- can1$Year - 1949

rs <- usam[Cause == "Respiratory system"][, c(2, 4:22)]
rs1 <- melt(rs, id.vars = "Year")
setorder(rs1, "Year", "variable", "value")
colnames(rs1) <- c("Year", "agegroups", "Respiratory system")
rs1$Year <- rs1$Year - 1949

iapd <- usam[Cause == "Infectious and parasitic diseases"][, c(2, 4:22)]
iapd1 <- melt(iapd, id.vars = "Year")
setorder(iapd1, "Year", "variable", "value")
colnames(iapd1) <- c("Year", "agegroups", "Infectious and parasitic diseases")
iapd1$Year <- iapd1$Year - 1949

ec <- usam[Cause == "External causes"][, c(2, 4:22)]
ec1 <- melt(ec, id.vars = "Year")
setorder(ec1, "Year", "variable", "value")
colnames(ec1) <- c("Year", "agegroups", "External causes")
ec1$Year <- ec1$Year - 1949

oth <- usam[Cause == "Others"][, c(2, 4:22)]
oth1 <- melt(oth, id.vars = "Year")
setorder(oth1, "Year", "variable", "value")
colnames(oth1) <- c("Year", "agegroups", "Others")
oth1$Year <- oth1$Year - 1949

probs <- cbind(cir1, can1[, 3], rs1[, 3], iapd1[, 3], ec1[, 3], oth1[, 3])
sur1 <- oth1
colnames(sur1) <- c("Year", "agegroups", "Survival")
sur1$Survival <- 1 - rowSums(probs[, 3:8])
probs <- cbind(sur1, probs[, 3:8])
colnames(probs)[5:9] <- c("Cancer", "Respiratory system", "Infectious and parasitic diseases", "External causes", "Others")

##### counts
resp <- left_join(probs, popn, by = c("Year", "agegroups"))
head(resp)

resp2 <- (resp[, 3:9] * resp$pop) %>%
    round() %>%
    as.matrix()

(rowSums(resp2) - resp$pop) %>%
    summary() ## Sanity check. All of these should be very close to zero
rownames(probs) <- rownames(resp) <- rownames(resp2) <- rownames(popn) <- NULL

## --------------------
## Fit multinomial logistic regression and compare results
## --------------------
years <- (1:58) %>%
    scale(center = FALSE, scale = FALSE)
agegroups <- c(
    "infant", "1-4",
    "5-9", "10-14",
    "15-19", "20-24",
    "25-29", "30-34",
    "35-39", "40-44",
    "45-49", "50-54",
    "55-59", "60-64",
    "65-69", "70-74",
    "75-79", "80-84",
    "85+"
)
dat <- data.frame(years = rep(years, each = length(unique(resp$agegroups))), agegroups = rep(agegroups, length(unique(resp$Year))))
head(dat)

make_offset <- matrix(resp$pop, nrow = nrow(probs), ncol = ncol(resp2))
# simplest
fit_mlr1s <- multinom(resp2[1:1007, ] ~ years + agegroups + offset(make_offset[1:1007, ]), data = dat[1:1007, ], maxit = 1000)
# single
fit_mlr1 <- multinom(resp2[1:1007, ] ~ years + agegroups + years:agegroups + offset(make_offset[1:1007, ]), data = dat[1:1007, ], maxit = 1000)
# quadratic
fit_mlr2 <- multinom(resp2[1:1007, ] ~ years + agegroups + I(years):agegroups + I(years^2) + I(years^2):agegroups + offset(make_offset[1:1007, ]), data = dat[1:1007, ], maxit = 1000)
# cubic
fit_mlr3 <- multinom(resp2[1:1007, ] ~ years + agegroups + I(years^2) + I(years^2):agegroups + I(years):agegroups + I(years^3) + I(years^3):agegroups + offset(make_offset[1:1007, ]), data = dat[1:1007, ], maxit = 5000)

summary(fit_mlr1s)
summary(fit_mlr1)
summary(fit_mlr2)
summary(fit_mlr3)

coefs1s <- coef(fit_mlr1s)
coefs1 <- coef(fit_mlr1)
coefs2 <- coef(fit_mlr2)
coefs3 <- coef(fit_mlr3)

## --------------------
## --------------------
## fit
## --------------------
## --------------------
fit_dat <- data.frame(years = rep(c(1:53), each = length(unique(resp$agegroups))), agegroups = rep(agegroups, 53))
X_fit1s <- model.matrix(~ years + agegroups, data = fit_dat)
X_fit1 <- model.matrix(~ years + agegroups + years:agegroups, data = fit_dat)
X_fit2 <- model.matrix(~ years + agegroups + I(years):agegroups + I(years^2) + I(years^2):agegroups, data = fit_dat)
X_fit3 <- model.matrix(~ years + agegroups + I(years^2) + I(years^2):agegroups + I(years):agegroups + I(years^3) + I(years^3):agegroups, data = fit_dat)


## --------------------
## simplest
## --------------------
fit_probs1s <- tcrossprod(X_fit1s, coefs1s) %>%
    exp()
fit_probs1s <- cbind(1 / (1 + rowSums(fit_probs1s)), fit_probs1s / (1 + rowSums(fit_probs1s)))
colnames(fit_probs1s)[1] <- "Survival"
fit_probs1s <- cbind(fit_dat, fit_probs1s[, 2:7])
colnames(fit_probs1s)[1] <- "Year"

df_wide1s <- as.data.table(fit_probs1s) %>% data.table::dcast(Year ~ agegroups, value.var = colnames(as.data.table(fit_probs1s))[3:8])

fit1s <- df_wide1s %>%
    gather(variable, value, -Year) %>%
    mutate(
        Cause = sub("_.*", "", variable),
        variable = sub(".*_", "", variable)
    ) %>%
    spread(variable, value)

fit1s <- fit1s[, colnames(usam[, c(2:22)])]

## --------------------
## single
## --------------------
fit_probs1 <- tcrossprod(X_fit1, coefs1) %>%
    exp()
fit_probs1 <- cbind(1 / (1 + rowSums(fit_probs1)), fit_probs1 / (1 + rowSums(fit_probs1)))
colnames(fit_probs1)[1] <- "Survival"
fit_probs1 <- cbind(fit_dat, fit_probs1[, 2:7])
colnames(fit_probs1)[1] <- "Year"

df_wide1 <- as.data.table(fit_probs1) %>% data.table::dcast(Year ~ agegroups, value.var = colnames(as.data.table(fit_probs1))[3:8])

fit1 <- df_wide1 %>%
    gather(variable, value, -Year) %>%
    mutate(
        Cause = sub("_.*", "", variable),
        variable = sub(".*_", "", variable)
    ) %>%
    spread(variable, value)

fit1 <- fit1[, colnames(usam[, c(2:22)])]
##### difference for single
fnorm(generate_mortality_tensor(fit1) - generate_mortality_tensor(usam[1:(53 * 6), ]))

## --------------------
## quadratic
## --------------------
fit_probs2 <- tcrossprod(X_fit2, coefs2) %>%
    exp()
fit_probs2 <- cbind(1 / (1 + rowSums(fit_probs2)), fit_probs2 / (1 + rowSums(fit_probs2)))
colnames(fit_probs2)[1] <- "Survival"
fit_probs2 <- cbind(fit_dat, fit_probs2[, 2:7])
colnames(fit_probs2)[1] <- "Year"

df_wide2 <- as.data.table(fit_probs2) %>% data.table::dcast(Year ~ agegroups, value.var = c("Circulatory system", "Cancer", "Respiratory system", "Infectious and parasitic diseases", "External causes", "Others"))

fit2 <- df_wide2 %>%
    gather(variable, value, -Year) %>%
    mutate(
        Cause = sub("_.*", "", variable),
        variable = sub(".*_", "", variable)
    ) %>%
    spread(variable, value)

fit2 <- fit2[, colnames(usam[, c(2:22)])]
##### difference for quadratic
fnorm(generate_mortality_tensor(fit2) - generate_mortality_tensor(usam[1:(53 * 6), ]))


## --------------------
## cubic
## --------------------
fit_probs3 <- tcrossprod(X_fit3, coefs3) %>%
    exp()
fit_probs3 <- cbind(1 / (1 + rowSums(fit_probs3)), fit_probs3 / (1 + rowSums(fit_probs3)))
colnames(fit_probs3)[1] <- "Survival"
fit_probs3 <- cbind(fit_dat, fit_probs3[, 2:7])
colnames(fit_probs3)[1] <- "Year"

df_wide3 <- as.data.table(fit_probs3) %>% data.table::dcast(Year ~ agegroups, value.var = c("Circulatory system", "Cancer", "Respiratory system", "Infectious and parasitic diseases", "External causes", "Others"))

fit3 <- df_wide3 %>%
    gather(variable, value, -Year) %>%
    mutate(
        Cause = sub("_.*", "", variable),
        variable = sub(".*_", "", variable)
    ) %>%
    spread(variable, value)

fit3 <- fit3[, colnames(usam[, c(2:22)])]
##### difference for quadratic
fnorm(generate_mortality_tensor(fit3) - generate_mortality_tensor(usam[1:(53 * 6), ]))


## --------------------
## --------------------
## predict
## --------------------
## --------------------
pre_dat <- data.frame(years = rep(c(54:58), each = length(unique(resp$agegroups))), agegroups = rep(agegroups, 5))
X_pred1s <- model.matrix(~ years + agegroups, data = pre_dat)
X_pred1 <- model.matrix(~ years + agegroups + years:agegroups, data = pre_dat)
X_pred2 <- model.matrix(~ years + agegroups + I(years):agegroups + I(years^2) + I(years^2):agegroups, data = pre_dat)
X_pred3 <- model.matrix(~ years + agegroups + I(years^2) + I(years^2):agegroups + I(years):agegroups + I(years^3) + I(years^3):agegroups, data = pre_dat)

## --------------------
## simplest
## --------------------
pred_probs1s <- tcrossprod(X_pred1s, coefs1s) %>%
    exp()
pred_probs1s <- cbind(1 / (1 + rowSums(pred_probs1s)), pred_probs1s / (1 + rowSums(pred_probs1s)))
colnames(pred_probs1s)[1] <- "Survival"
pred_probs1s <- cbind(pre_dat, pred_probs1s[, 2:7])
colnames(pred_probs1s)[1] <- "Year"

df_wide1ps <- as.data.table(pred_probs1s) %>% data.table::dcast(Year ~ agegroups, value.var = c("Circulatory system", "Cancer", "Respiratory system", "Infectious and parasitic diseases", "External causes", "Others"))

pred1s <- df_wide1ps %>%
    gather(variable, value, -Year) %>%
    mutate(
        Cause = sub("_.*", "", variable),
        variable = sub(".*_", "", variable)
    ) %>%
    spread(variable, value)

pred1s <- pred1s[, colnames(usam[, c(2:22)])]
##### difference
fnorm(generate_mortality_tensor(pred1s) - generate_mortality_tensor(usam[(53 * 6 + 1):(58 * 6), ]))


## --------------------
## single
## --------------------
pred_probs1 <- tcrossprod(X_pred1, coefs1) %>%
    exp()
pred_probs1 <- cbind(1 / (1 + rowSums(pred_probs1)), pred_probs1 / (1 + rowSums(pred_probs1)))
colnames(pred_probs1)[1] <- "Survival"
pred_probs1 <- cbind(pre_dat, pred_probs1[, 2:7])
colnames(pred_probs1)[1] <- "Year"

df_wide1p <- as.data.table(pred_probs1) %>% data.table::dcast(Year ~ agegroups, value.var = c("Circulatory system", "Cancer", "Respiratory system", "Infectious and parasitic diseases", "External causes", "Others"))

pred1 <- df_wide1p %>%
    gather(variable, value, -Year) %>%
    mutate(
        Cause = sub("_.*", "", variable),
        variable = sub(".*_", "", variable)
    ) %>%
    spread(variable, value)

pred1 <- pred1[, colnames(usam[, c(2:22)])]
##### difference
fnorm(generate_mortality_tensor(pred1) - generate_mortality_tensor(usam[(53 * 6 + 1):(58 * 6), ]))



## --------------------
## quadratic
## --------------------
pred_probs2 <- tcrossprod(X_pred2, coefs2) %>%
    exp()
pred_probs2 <- cbind(1 / (1 + rowSums(pred_probs2)), pred_probs2 / (1 + rowSums(pred_probs2)))
colnames(pred_probs2)[1] <- "Survival"
pred_probs2 <- cbind(pre_dat, pred_probs2[, 2:7])
colnames(pred_probs2)[1] <- "Year"

df_wide2p <- as.data.table(pred_probs2) %>% data.table::dcast(Year ~ agegroups, value.var = c("Circulatory system", "Cancer", "Respiratory system", "Infectious and parasitic diseases", "External causes", "Others"))

pred2 <- df_wide2p %>%
    gather(variable, value, -Year) %>%
    mutate(
        Cause = sub("_.*", "", variable),
        variable = sub(".*_", "", variable)
    ) %>%
    spread(variable, value)

pred2 <- pred2[, colnames(usam[, c(2:22)])]
##### difference
fnorm(generate_mortality_tensor(pred2) - generate_mortality_tensor(usam[(53 * 6 + 1):(58 * 6), ]))

## --------------------
## cubic
## --------------------

pred_probs3 <- tcrossprod(X_pred3, coefs3) %>%
    exp()
pred_probs3 <- cbind(1 / (1 + rowSums(pred_probs3)), pred_probs3 / (1 + rowSums(pred_probs3)))
colnames(pred_probs3)[1] <- "Survival"
pred_probs3 <- cbind(pre_dat, pred_probs3[, 2:7])
colnames(pred_probs3)[1] <- "Year"

df_wide3p <- as.data.table(pred_probs3) %>% data.table::dcast(Year ~ agegroups, value.var = c("Circulatory system", "Cancer", "Respiratory system", "Infectious and parasitic diseases", "External causes", "Others"))

pred3 <- df_wide3p %>%
    gather(variable, value, -Year) %>%
    mutate(
        Cause = sub("_.*", "", variable),
        variable = sub(".*_", "", variable)
    ) %>%
    spread(variable, value)

pred3 <- pred3[, colnames(usam[, c(2:22)])]
##### difference
fnorm(generate_mortality_tensor(pred3) - generate_mortality_tensor(usam[(53 * 6 + 1):(58 * 6), ]))