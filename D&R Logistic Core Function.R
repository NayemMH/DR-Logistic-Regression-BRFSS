drglm<-function(formula,family,data,k,fitfunction)
{
  result <- list()
  class(result) <- "drglm"
  result$call <- match.call()
  result$formula <- formula
  result$data <- data
  result$family <- family
  
  if(family=="binomial" & fitfunction=="glm")
  {
    n = nrow(data)
    rows_per_chunk = ceiling(n / k)
    split_data = lapply(1:k, function(i) {
      start_row = (i-1) * rows_per_chunk + 1
      end_row = min(i * rows_per_chunk, n)
      return(data[start_row:end_row, ])
    })
    
    model <- list()
    
    for (i in 1:length(split_data)) {
      
      # fit the glm() function for each element
      model[[i]] <- glm(formula,split_data[[i]], family = binomial)
    }
    
    vcov <- list()
    
    for (i in 1:length(model)) {
      # calculate the vcov() function for each element
      vcov[[i]] <- vcov(model[[i]])
    }
    
    v <- list()
    for (i in 1:length(vcov))
    {
      v[[i]] <- as.matrix(diag(vcov[[i]]))/k
    }
    
    v_com <- Reduce("+", v)/k
    
    se_com <- sqrt(v_com) # Take the square root of v_com
    
    H <- list()
    
    for (i in 1:length(vcov))
    {
      H[[i]] <- as.matrix(solve(vcov[[i]]))
    }
    
    # Beta's of models
    b <- list()
    for (i in 1:length(model))
    {
      b[[i]] <- as.matrix(coef(model[[i]]))
    }
    
    HE <- Reduce("+", H)
    
    IH=solve(HE)
    
    HB <- lapply(1:k, function(i) H[[i]]%*%b[[i]])
    # Sum the elements of the list
    HB <- Reduce("+", HB)
    
    B=(IH%*%HB)
    
    OR=exp(B)
    
    alpha <- 0.05
    
    z<-qnorm(1-alpha/2)
    
    t<-qt(1-alpha/2,df = length(B))
    
    l_normal= B-z*se_com
    u_normal=B+z*se_com
    
    lower_t=B-t*se_com
    upper_t=B+t*se_com
    
    Z=B/se_com
    
    p_value=2*(1- pnorm(abs(Z)))
    
    # creating a data frame with four columns
    table <- data.frame("Estimate"=B,
                        "Odds Ratio"=OR,
                        "standard error"=se_com ,
                        "z value"= Z,
                        "Pr(>|z|)"=p_value,
                        "95% CI" = paste("[", round(l_normal, 2), ",", round(u_normal, 2), "]"),
                        check.names = FALSE)
    
    
    result$coefficients <- B
    result$Estimates <- table
    result$models <- model
    result$df.residual <- nrow(data) - nrow(B)
    
  }
  else if (family=="binomial" &  fitfunction=="speedglm")
  {
    n = nrow(data)
    rows_per_chunk = ceiling(n / k)
    split_data = lapply(1:k, function(i) {
      start_row = (i-1) * rows_per_chunk + 1
      end_row = min(i * rows_per_chunk, n)
      return(data[start_row:end_row, ])
    })
    
    model <- list()
    
    for (i in 1:length(split_data)) {
      
      model[[i]] <- speedglm::speedglm(formula,split_data[[i]], family = binomial())
    }
    
    vcov <- list()
    
    for (i in 1:length(model)) {
      # calculate the vcov() function for each element
      vcov[[i]] <- vcov(model[[i]])
    }
    
    v <- list()
    for (i in 1:length(vcov))
    {
      v[[i]] <- as.matrix(diag(vcov[[i]]))/k
    }
    
    v_com <- Reduce("+", v)/k
    
    se_com <- sqrt(v_com) # Take the square root of v_com
    
    H <- list()
    
    for (i in 1:length(vcov))
    {
      H[[i]] <- as.matrix(solve(vcov[[i]]))
    }
    
    b <- list()
    for (i in 1:length(model))
    {
      b[[i]] <- as.matrix(coef(model[[i]]))
    }
    
    HE <- Reduce("+", H)
    IH=solve(HE)
    
    HB <- lapply(1:k, function(i) H[[i]]%*%b[[i]])
    # Sum the elements of the list
    HB <- Reduce("+", HB)
    
    B=(IH%*%HB)
    OR=exp(B)
    
    alpha <- 0.05
    
    z<-qnorm(1-alpha/2)
    
    t<-qt(1-alpha/2,df = length(B))
    
    l_normal= B-z*se_com
    u_normal=B+z*se_com
    
    lower_t=B-t*se_com
    upper_t=B+t*se_com
    
    Z=B/se_com
    
    p_value=2*(1- pnorm(abs(Z)))
    
    table <- data.frame("Estimate"=B,
                        "Odds Ratio"=OR,
                        "standard error"=se_com ,
                        "z value"= Z,
                        "Pr(>|z|)"=p_value,
                        "95% CI" = paste("[", round(l_normal, 2), ",", round(u_normal, 2), "]"),
                        check.names = FALSE)
    
    result$coefficients <- B
    result$Estimates <- table
    result$models <- model
    result$df.residual <- nrow(data) - nrow(B)
  }
  else
  {
    stop("Unsupported family")
  }
  
  class(result) <- "drglm"
  return(result)
}


print.drglm <- function(x, ...) {
  cat("Generalized Linear Model in Divide-and-Recombine Approach\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nFamily:", x$family, "\n")
  cat("\nCoefficients:\n")
  print(x$Estimates)
  invisible(x)
}


drglm_residuals <- function(model, type = "response") {
  B <- as.matrix(model$Estimates[, "Estimate"])
  mf <- model.frame(model$formula, data = model$data)
  y <- model.response(mf)
  X <- model.matrix(model$formula, data = model$data)
  eta <- as.vector(X %*% B)
  mu <- switch(model$family,
               "gaussian" = eta,
               "binomial" = 1 / (1 + exp(-eta)),
               "poisson" = exp(eta),
               stop("Unsupported family for residuals"))
  if (type == "response") {
    return(y - mu)
  } else if (type == "deviance") {
    res <- switch(model$family,
                  "gaussian" = y - mu,
                  "binomial" = {
                    eps <- .Machine$double.eps
                    y <- pmin(pmax(y, eps), 1 - eps)
                    mu <- pmin(pmax(mu, eps), 1 - eps)
                    sign(y - mu) * sqrt(2 * (y * log(y / mu) + (1 - y) * log((1 - y) / (1 - mu))))
                  },
                  "poisson" = {
                    eps <- .Machine$double.eps
                    y <- pmax(y, eps)
                    mu <- pmax(mu, eps)
                    sign(y - mu) * sqrt(2 * (y * log(y / mu) - (y - mu)))
                  },
                  stop("Unsupported family for deviance residuals"))
    return(res)
  } else if (type == "pearson") {
    var_fun <- switch(model$family,
                      "gaussian" = rep(1, length(mu)),
                      "binomial" = mu * (1 - mu),
                      "poisson" = mu,
                      stop("Unsupported family for Pearson residuals"))
    return((y - mu) / sqrt(var_fun))
  } else {
    stop("Unsupported residual type. Use 'response', 'deviance', or 'pearson'.")
  }
}
