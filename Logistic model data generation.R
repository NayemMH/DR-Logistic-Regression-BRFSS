#GLM 
set.seed(123)
mc_simulation <- function()
{
  #set.seed(123)
  nobs <- 5000000
  x1 <- runif(nobs)
  x2 <- runif(nobs)
  x3 <- runif(nobs)
  x4 <- runif(nobs)
  x5 <- runif(nobs)
  x6 <- runif(nobs)
  x7 <- runif(nobs)
  x8 <- runif(nobs)
  x9 <- runif(nobs)
  x10 <- runif(nobs)
  
  b <- c(2, .75, -1.25, .5, .6, 1.45, -.4, 1.95, .55, 1.10, -.80)
  X <- cbind(1, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  eta <- X %*% b
  inv.logit = function(p)
  {
    return(exp(p)/(1+exp(p)))
  }
  
  y <- rbinom(nobs, 1, inv.logit(eta))
  
  
  model_glm <- glm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10, family = binomial)
  
  pr <- sum(residuals(model_glm, type = "pearson")^2)
  prdisp <- pr / model_glm$df.residual
  beta <- model_glm$coef
  se <- sqrt(diag(vcov(model_glm)))
  
  list(beta = beta, se = se, prdisp = prdisp)
}

B1_logistic <- replicate(1000, mc_simulation(), simplify = FALSE)


# DRGLM 100

set.seed(123)
mc_simulation <- function()
{
  # set.seed(123)
  nobs <- 5000000
  x1 <- runif(nobs)
  x2 <- runif(nobs)
  x3 <- runif(nobs)
  x4 <- runif(nobs)
  x5 <- runif(nobs)
  x6 <- runif(nobs)
  x7 <- runif(nobs)
  x8 <- runif(nobs)
  x9 <- runif(nobs)
  x10 <- runif(nobs)
  
  b <- c(2, .75, -1.25, .5, .6, 1.45, -.4, 1.95, .55, 1.10, -.80)
  X <- cbind(1, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  eta <- X %*% b
  inv.logit = function(p)
  {
    return(exp(p)/(1+exp(p)))
  }
  
  y <- rbinom(nobs, 1, inv.logit(eta))
  
  data= data.frame(y,X)
  model_drglm<-  drglm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10, 
                       family="binomial", fitfunction = "glm", k=100, data)
  
  
  
  pr<-  sum(drglm_residuals(model_drglm, type="response")^2)
  prdisp<-  pr/model_drglm$df.residual
  beta<-  model_drglm$coef
  se <-  model_drglm$Estimates[,3]
  
  
  list(beta = beta, se = se, prdisp = prdisp)
}

B2_logistic <- replicate(1000, mc_simulation(), simplify = FALSE)
save(B2_logistic, file="B2_logistic.RData")



# DRGLM 50

set.seed(123)
mc_simulation <- function()
{
  nobs <- 5000000
  x1 <- runif(nobs)
  x2 <- runif(nobs)
  x3 <- runif(nobs)
  x4 <- runif(nobs)
  x5 <- runif(nobs)
  x6 <- runif(nobs)
  x7 <- runif(nobs)
  x8 <- runif(nobs)
  x9 <- runif(nobs)
  x10 <- runif(nobs)
  
  b <- c(2, .75, -1.25, .5, .6, 1.45, -.4, 1.95, .55, 1.10, -.80)
  X <- cbind(1, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  eta <- X %*% b
  inv.logit = function(p)
  {
    return(exp(p)/(1+exp(p)))
  }
  y <- rbinom(nobs, 1, inv.logit(eta))
  
  data= data.frame(y,X)
  model_drglm<-  drglm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10, 
                       family="binomial", fitfunction = "glm", k=50, data)
  
  
  
  pr<-  sum(drglm_residuals(model_drglm, type="response")^2)
  prdisp<-  pr/model_drglm$df.residual
  beta<-  model_drglm$coef
  se <-  model_drglm$Estimates[,3]
  
  
  list(beta = beta, se = se, prdisp = prdisp)
}

B3_logistic <- replicate(1000, mc_simulation(), simplify = FALSE)

save(B3_logistic, file="B3_logistic.RData")



# DRGLM 25

set.seed(123)
mc_simulation <- function()
{
  nobs <- 5000000
  x1 <- runif(nobs)
  x2 <- runif(nobs)
  x3 <- runif(nobs)
  x4 <- runif(nobs)
  x5 <- runif(nobs)
  x6 <- runif(nobs)
  x7 <- runif(nobs)
  x8 <- runif(nobs)
  x9 <- runif(nobs)
  x10 <- runif(nobs)
  
  b <- c(2, .75, -1.25, .5, .6, 1.45, -.4, 1.95, .55, 1.10, -.80)
  X <- cbind(1, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  eta <- X %*% b
  inv.logit = function(p)
  {
    return(exp(p)/(1+exp(p)))
  }
  y <- rbinom(nobs, 1, inv.logit(eta))
  
  data= data.frame(y,X)
  model_drglm<-  drglm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10, 
                       family="binomial", fitfunction = "glm", k=25, data)
  
  
  pr<-  sum(drglm_residuals(model_drglm, type="response")^2)
  prdisp<-  pr/model_drglm$df.residual
  beta<-  model_drglm$coef
  se <-  model_drglm$Estimates[,3]
  
  
  list(beta = beta, se = se, prdisp = prdisp)
}

B4_logistic <- replicate(1000, mc_simulation(), simplify = FALSE)

save(B4_logistic, file="B4_logistic.RData")