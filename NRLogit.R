set.seed(1000)

x = rnorm(100)
y = rbinom(100, 1,0.4)

x0 = rep(1,100) #b_0對應的x_i0=1，可以算出截距
X = cbind(x0,x)
colnames(X) <- c("Int", "x")

beta = rep(0,ncol(X))#先從0開始猜

log_reg_p<- function(X, beta) {
  # 先算機率，機率是透過beta轉置矩陣與X算出
  X <- as.matrix(X)
  beta <- as.vector(beta)
  p <- exp(X %*% beta) / (1 + exp(X %*% beta))
  return(p)
}

log_reg_like <- function(y, p){
  #計算log-likelihood，會利用log-likelihood做數值方法的控制
  l <- t(y) %*% log(p) + t(1 - y) %*% log(1 - p)
  return(l)
}

log_reg_newton <- function(X, y, beta, thres1 = 1e-10, thres2 = 1e-10, iter_count = 50){

  # thres1 = 考慮beta的容忍度大小
  # thres2 = 考慮likelihood的容忍度大小
  # iter_count = 最大迭代次數

beta_1 <- rep(-Inf, length(beta)) # logistic有global optima，安心地把beta_1設定在負無限大沒問題！
diff.beta <- sqrt(sum((beta - beta_1)^2))
log_like.1 <- log_reg_like(y,log_reg_p(X, beta)) # 更新 log-likelihood
log_like.2 <- log_reg_like(y,log_reg_p(X, beta_1))
diff.like <- abs(log_like.1 - log_like.2) 
if (is.nan(diff.like)) { diff.like <- 1e9 }

i <- 1 #迭代次數

alpha.step <- seq(-1, 2, by = 0.1)[-11] # line search step sizes, excluding 0
NR.hist <- data.frame(i, diff.beta, diff.like, log_like.1, step.size = 1) # iteration history
beta.hist <- matrix(beta, nrow = 1)
while ((i <= iter_count) & (diff.beta > thres1) & (diff.like > thres2)){
  i <- i + 1 # increment iteration
  # update beta
  beta_1 <- beta # old guess is current guess
  mu.2 <- log_reg_p(X, beta_1) # m * p is mean
  # variance matrix
  v.2 <- diag(as.vector(log_reg_p(X, beta_1) * (1 - log_reg_p(X, beta_1))))
  score.2 <- t(X) %*% (y - mu.2) # score function
  # this increment version inverts the information matrix
  # Iinv.2 <- solve(t(X) %*% v.2 %*% X) # Inverse information matrix
  # increm <- Iinv.2 %*% score.2 # increment, solve() is inverse
  # this increment version solves for (beta_1-beta.1) without inverting Information
  increm <- solve(t(X) %*% v.2 %*% X, score.2) # solve for increment
  # line search for improved step size
  log_like.alpha.step <- rep(NA, length(alpha.step)) # init log_like for line search
  for (i.alpha.step in 1:length(alpha.step)) {
  log_like.alpha.step[i.alpha.step] <- log_reg_like(y, log_reg_p(X, beta_1 + alpha.step[i.alpha.step] * increm))
  }
  # step size index for max increase in log-likelihood (if tie, [1] takes first)
  ind.max.alpha.step <- which(log_like.alpha.step == max(log_like.alpha.step))[1]
  beta <- beta_1 + alpha.step[ind.max.alpha.step] * increm # update beta
  diff.beta <- sqrt(sum((beta - beta_1)^2)) # Euclidean distance
  log_like.2 <- log_like.1 # age likelihood value
  log_like.1 <- log_reg_like(y, log_reg_p(X, beta)) # update loglikelihood
  diff.like <- abs(log_like.1 - log_like.2) # diff
  # iteration history
  NR.hist <- rbind(NR.hist, c(i, diff.beta, diff.like, log_like.1, alpha.step[ind.max.alpha.step]))
  beta.hist <- rbind(beta.hist, matrix(beta, nrow = 1))
}
# prepare output
out <- list()
out$beta.MLE <- beta
out$iter <- i - 1
out$NR.hist <- NR.hist
out$beta.hist <- beta.hist
v.1 <- diag(as.vector(log_reg_p(X, beta) * (1 - log_reg_p(X, beta))))
Iinv.1 <- solve(t(X) %*% v.1 %*% X) # Inverse information matrix
out$beta.cov <- Iinv.1
if (!(diff.beta > thres1) & !(diff.like > thres2)) {
out$note <- paste("Absolute convergence of", thres1, "for betas and"
                  , thres2, "for log-likelihood satisfied")
}
if (i > iter_count) {
out$note <- paste("Exceeded max iterations of ", iter_count)
}
return(out)
}
out <- log_reg_newton(X, y, beta)
out

options(digits = 8)
GLM_1= glm(y~x, family = "binomial")
Z<-summary(GLM_1,corr=T)
Z$cov.unscaled

