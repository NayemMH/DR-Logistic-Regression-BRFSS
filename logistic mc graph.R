#GLM

beta_hat = do.call(rbind, lapply(B1_logistic, function(x) x$beta))

par(mfrow = c(4, 3), mar = c(4.5, 4, 1.5, 1), las = 1)

for(i in 0:10) {
  
  x_vals = beta_hat[, i + 1]
  mean_x = mean(x_vals)
  sd_x = sd(x_vals)
  
  hist(x_vals, breaks = 50, freq = FALSE, main = "",
       xlab = bquote(paste("simulated ", beta[.(i)])),
       xlim = c(mean_x - 4*sd_x, mean_x + 4*sd_x),
       col = "gray90", border = "gray40")
  
  curve(dnorm(x, mean = mean_x, sd = sd_x),
        col = "gray40", lwd = 2, add = TRUE)
  
  abline(v = c(2, .75, -1.25, .5, .6, 1.45, -.4, 1.95, .55, 1.10, -.80),
         lwd = 2, col = "gray40")
}

plot.new()  
box()       
text(0.5, 0.90, "Model: Logistic Model", cex = 1.2)
text(0.5, 0.60,  expression(n == "5,000,000"), cex = 1.2)
text(0.5, 0.35, expression(k == 1), cex = 1.2)
text(0.5, 0.10, "fitfunction = glm()", cex = 1.2)




#DRGLM 100
beta_hat = do.call(rbind, lapply(B2_logistic, function(x) t(x$beta)))

par(mfrow = c(4, 3), mar = c(4.5, 4, 1.5, 1), las = 1)

for(i in 0:10) {
  
  x_vals = beta_hat[, i + 1]
  mean_x = mean(x_vals)
  sd_x = sd(x_vals)
  
  hist(x_vals, breaks = 50, freq = FALSE, main = "",
       xlab = bquote(paste("simulated ", beta[.(i)])),
       xlim = c(mean_x - 4*sd_x, mean_x + 4*sd_x),
       col = "gray90", border = "gray40")
  
  curve(dnorm(x, mean = mean_x, sd = sd_x),
        col = "gray40", lwd = 2, add = TRUE)
  
  abline(v = c(2, .75, -1.25, .5, .6, 1.45, -.4, 1.95, .55, 1.10, -.80),
         lwd = 2, col = "gray40")
}

plot.new()  
box()       
text(0.5, 0.90, "Model: Logistic Model", cex = 1.2)
text(0.5, 0.60,  expression(n == "5,000,000"), cex = 1.2)
text(0.5, 0.35, expression(k == 100), cex = 1.2)
text(0.5, 0.10, "fitfunction = D&R", cex = 1.2)


#DRGLM 50
beta_hat = do.call(rbind, lapply(B3_logistic, function(x) t(x$beta)))

par(mfrow = c(4, 3), mar = c(4.5, 4, 1.5, 1), las = 1)

for(i in 0:10) {
  
  x_vals = beta_hat[, i + 1]
  mean_x = mean(x_vals)
  sd_x = sd(x_vals)
  
  hist(x_vals, breaks = 50, freq = FALSE, main = "",
       xlab = bquote(paste("simulated ", beta[.(i)])),
       xlim = c(mean_x - 4*sd_x, mean_x + 4*sd_x),
       col = "gray90", border = "gray40")
  
  curve(dnorm(x, mean = mean_x, sd = sd_x),
        col = "gray40", lwd = 2, add = TRUE)
  
  abline(v = c(2, .75, -1.25, .5, .6, 1.45, -.4, 1.95, .55, 1.10, -.80),
         lwd = 2, col = "gray40")
}

plot.new()  
box()       
text(0.5, 0.90, "Model: Logistic Model", cex = 1.2)
text(0.5, 0.60,  expression(n == "5,000,000"), cex = 1.2)
text(0.5, 0.35, expression(k == 50), cex = 1.2)
text(0.5, 0.10, "fitfunction = D&R", cex = 1.2)

#DRGLM 25

beta_hat = do.call(rbind, lapply(B4, function(x) t(x$beta)))

par(mfrow = c(4, 3), mar = c(4.5, 4, 1.5, 1), las = 1)

for(i in 0:10) {
  
  x_vals = beta_hat[, i + 1]
  mean_x = mean(x_vals)
  sd_x = sd(x_vals)
  
  hist(x_vals, breaks = 50, freq = FALSE, main = "",
       xlab = bquote(paste("simulated ", beta[.(i)])),
       xlim = c(mean_x - 4*sd_x, mean_x + 4*sd_x),
       col = "gray90", border = "gray40")
  
  curve(dnorm(x, mean = mean_x, sd = sd_x),
        col = "gray40", lwd = 2, add = TRUE)
  
  abline(v = c(2, .75, -1.25, .5, .6, 1.45, -.4, 1.95, .55, 1.10, -.80),
         lwd = 2, col = "gray40")
}

plot.new()  
box()       
text(0.5, 0.90, "Model: Logistic Model", cex = 1.2)
text(0.5, 0.60,  expression(n == "5,000,000"), cex = 1.2)
text(0.5, 0.35, expression(k == 25), cex = 1.2)
text(0.5, 0.10, "fitfunction = D&R", cex = 1.2)

