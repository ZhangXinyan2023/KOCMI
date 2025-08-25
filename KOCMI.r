library('MASS')
library('Matrix')
#---- function: cause inference ----
KOCMI.net=function(expr, k=3, M, tf=colnames(expr), pcc = 0){
  expr <- as.data.frame(expr)
  expr <- expr[,tf]
  weightNet <- transMatrix(expr)
  if(pcc == 0){
    df.ko0 <- x_knockoff(expr,M)
    a <- apply(weightNet,1,function(x){KOCMI(expr,x[1],x[2],k,M,df.ko0)})
  }else{
    sigma <- cor(expr)
    a <- apply(weightNet,1,function(x){KOCMI(expr,x[1],x[2],k,M,sigma=sigma,pcc=pcc)})
  }
  p.value <- lapply(a,function(x){x$pvalue})
  t1 <- lapply(a,function(x){x$t1})
  weightNet$pvalue <- as.numeric(p.value)
  weightNet$p_adj  <- as.numeric(p.adjust(p.value,method = 'BH'))
  weightNet$t1 <- as.numeric(t1)
  weightNet$cs <- abs(as.numeric(t1))
  weightNet
  return(weightNet)
}

transMatrix = function(expr){
  tf <- colnames(expr)
  n <- length(tf)
  re <- matrix(0, n*(n-1), 4, dimnames = list(c(), c('regulator', 'target', 'cs', 'pvalue')))
  re <- as.data.frame(re)
  re$regulator <- rep(tf, each = (n-1))
  target <- c()
  for(i in 1:n){
    target <- c(target, tf[-i])
  }
  re$target <- target
  return(re)
}

#---- function: cause inference ----
KOCMI = function(data, cause, effect, k=3, M, df.ko0 = NULL, sigma = NULL, pcc = 0){
  if(pcc != 0 & !is.null(sigma)){
    sigma <- sigma[, colnames(sigma) != effect]
    pccnode <- names(which(abs(sigma[cause,]) >= pcc))
    pccnode <- c(pccnode, effect)
    if(length(pccnode) <= 2){
      pccnode <- names(sort(abs(sigma[cause,]), decreasing = T)[1:2])
      pccnode <- c(pccnode, effect)
    }
    condition <- setdiff(pccnode, c(cause, effect))
    
    data <- data[, pccnode]
  }else{
    condition <- setdiff(colnames(data), c(cause, effect))
  }
  
  df.cause <- data[, cause]
  df.effect <- data[, effect]
  df.condition <- data[, condition]
  
  if(is.null(df.ko0)){df.ko0 <- x_knockoff(data, M)}
  
  cmi0 <- KNN.CMI(data, cause, effect, k)
  
  cmi0.knockoff <- apply(df.ko0, 3, function(x){KNN.CMI(cbind(cause = x[,cause],effect = df.effect, df.condition), 'cause', 'effect', k)})
  
  df.knockoff <- x_knockoff(cbind(df.cause, df.condition), M)
  
  cmi.knockoff <- apply(df.knockoff, 3, function(x){KNN.CMI(cbind(cause = x[,1], effect = df.effect, df.condition), 'cause', 'effect', k)})

  cmi.knockoff <- round(cmi.knockoff, digits = 8)
  cmi0.knockoff <- round(cmi0.knockoff, digits = 8)

  D = cmi0.knockoff-cmi.knockoff
  p = permutation_test_mean(D)$p_value
  t1 = mean(D)/sd(D)
  #t3 = mean(cmi0-cmi.knockoff)/sd(cmi0-cmi.knockoff)
  #t2 = (mean(cmi0.knockoff) - mean(cmi.knockoff)) / sqrt(var(cmi0.knockoff)/M + var(cmi0.knockoff)/M)
  
  return(list(pvalue=p,t1=t1,cmi0=cmi0,
              cmi0.knockoff=cmi0.knockoff,cmi.knockoff=cmi.knockoff))
}

# KNN CMI
KNN.CMI = function(data, cause = 'x', effect = 'y', k = 3){
  
  condition <- setdiff(colnames(data), c(cause, effect))
  df_xz <- data[,c(cause, condition)]
  df_yz <- data[,c(effect, condition)]
  df_z <- as.matrix(data[,condition])
  
  D_all <- as.matrix(dist(data, method = 'maximum'))
  D_z <- as.matrix(dist(df_z, method = 'maximum'))
  D_xz <- as.matrix(dist(df_xz, method = 'maximum'))
  D_yz <- as.matrix(dist(df_yz, method = 'maximum'))
  
  N <- nrow(data)
  n_xz <- c(); n_yz <- c(); n_z <- c()
  
  for(i in 1:N){
    epsilon <- sort(D_all[i,])[k+1]
    n_xz <- c(n_xz, length(which(D_xz[i,] < epsilon)))
    n_yz <- c(n_yz, length(which(D_yz[i,] < epsilon)))
    n_z <- c(n_z, length(which(D_z[i,] < epsilon)))
  }
  
  c <- digamma(k) - mean(digamma(n_xz)) - mean(digamma(n_yz)) + mean(digamma(n_z))
  return(c)
}

permutation_test_mean <- function(x, n_perm = 10000) {
  # 观测统计量（绝对值）
  observed_stat <- abs(mean(x))
  
  # 置换生成的统计量
  permuted_stats <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    # 随机翻转符号
    signs <- sample(c(-1, 1), size = length(x), replace = TRUE)
    permuted_x <- x * signs
    permuted_stats[i] <- abs(mean(permuted_x))
  }
  
  # 计算 p 值：大于等于观测值的比例
  p_value <- mean(permuted_stats >= observed_stat)
  
  return(list(
    observed_statistic = observed_stat,
    p_value = p_value
  ))
}

#---- function: knockoff ----
#library('GhostKnockoff')
#cite from package "Ghostknockoff"

#GhostKnockoff.prelim()
# getAnywhere(GhostKnockoff.prelim) 查看源代码
is_posdef = function (A, tol = 1e-09) 
{
  p = nrow(matrix(A))
  if (p < 500) {
    lambda_min = min(eigen(A)$values)
  }
  else {
    oldw <- getOption("warn")
    lambda_min = suppressWarnings(RSpectra::eigs(A, 1, which = "SM", 
                                                 opts = list(retvec = FALSE, maxitr = 100, tol))$values)
    options(warn = oldw)
    if (length(lambda_min) == 0) {
      lambda_min = min(eigen(A)$values)
    }
  }
  return(lambda_min > tol * 10)
}

GhostKnockoff.prelim = function (cor.G, M = 5, method = "asdp", max.size = 500) 
{
  temp.index <- 1:nrow(cor.G)
  n.G <- nrow(cor.G)
  permute.index <- rep(0, n.G)
  permute.index[-temp.index] <- 1
  Normal_50Studies <- matrix(rnorm(n.G * M * 50), n.G * M, 
                             50)
  P.each <- matrix(0, n.G, n.G)
  if (length(temp.index) != 0) {
    Sigma <- cor.G[temp.index, temp.index, drop = F]
    SigmaInv <- ginv(Sigma) #solve(Sigma) ginv
    if (method == "sdp") {
      temp.s <- create.solve_sdp_M(Sigma, M = M)
    }
    if (method == "asdp") {
      temp.s <- create.solve_asdp_M(Sigma, M = M, max.size = max.size)
    }
    s <- temp.s
    diag_s <- diag(s, length(s))
    if (sum(s) == 0) {
      V.left <- matrix(0, length(temp.index) * M, length(temp.index) * 
                         M)
    }
    else {
      Sigma_k <- 2 * diag_s - s * t(s * SigmaInv)
      V.each <- Matrix(forceSymmetric(Sigma_k - diag_s))
      V <- matrix(1, M, M) %x% V.each
      diag(V) <- diag(V) + rep(s, M)
      V.left <- try(t(chol(V)), silent = T)
      if (class(V.left) == "try-error") {
        svd.fit <- svd(V)
        u <- svd.fit$u
        svd.fit$d[is.na(svd.fit$d)] <- 0
        cump <- cumsum(svd.fit$d)/sum(svd.fit$d)
        n.svd <- which(cump >= 0.999)[1]
        if (is.na(n.svd)) {
          n.svd <- nrow(V)
        }
        svd.index <- intersect(1:n.svd, which(svd.fit$d != 
                                                0))
        V.left <- t(sqrt(svd.fit$d[svd.index]) * t(u[, 
                                                     svd.index, drop = F]))
      }
    }
    P.each[temp.index, temp.index] <- diag(1, length(s)) - 
      s * SigmaInv
    V.index <- rep(temp.index, M) + rep(0:(M - 1), each = length(temp.index)) * 
      n.G
    Normal_50Studies[V.index, ] <- as.matrix(V.left %*% matrix(rnorm(ncol(V.left) * 
                                                                       50), ncol(V.left), 50))
    permute.index[temp.index[s == 0]] <- 1
  }
  permute.V.index <- rep(permute.index, M)
  P.each[permute.index == 1, ] <- 0
  Normal_50Studies[permute.V.index == 1, ] <- matrix(rnorm(sum(permute.index) * 
                                                             M * 50), sum(permute.index) * M, 50)
  return(list(P.each = as.matrix(P.each), V.left = V.left, 
              Normal_50Studies = as.matrix(Normal_50Studies), permute.index = permute.index, 
              M = M))
}

create.solve_sdp_M = function(Sigma, M = 1, gaptol = 1e-06, maxit = 1000, verbose = FALSE) 
{
  stopifnot(isSymmetric(Sigma))
  G = cov2cor(Sigma)
  p = dim(G)[1]
  if (!is_posdef(G)) {
    warning("The covariance matrix is not positive-definite: knockoffs may not have power.", 
            immediate. = T)
  }
  Cl1 = rep(0, p)
  Al1 = -Matrix::Diagonal(p)
  Cl2 = rep(1, p)
  Al2 = Matrix::Diagonal(p)
  d_As = c(diag(p))
  As = Matrix::Diagonal(length(d_As), x = d_As)
  As = As[which(Matrix::rowSums(As) > 0), ]
  Cs = c((M + 1)/M * G)
  A = cbind(Al1, Al2, As)
  C = matrix(c(Cl1, Cl2, Cs), 1)
  K = NULL
  K$s = p
  K$l = 2 * p
  b = rep(1, p)
  OPTIONS = NULL
  OPTIONS$gaptol = gaptol
  OPTIONS$maxit = maxit
  OPTIONS$logsummary = 0
  OPTIONS$outputstats = 0
  OPTIONS$print = 0
  if (verbose) 
    cat("Solving SDP ... ")
  sol = Rdsdp::dsdp(A, b, C, K, OPTIONS)
  if (verbose) 
    cat("done. \n")
  if (!identical(sol$STATS$stype, "PDFeasible")) {
    warning("The SDP solver returned a non-feasible solution. Knockoffs may lose power.")
  }
  s = sol$y
  s[s < 0] = 0
  s[s > 1] = 1
  if (verbose) 
    cat("Verifying that the solution is correct ... ")
  psd = 0
  s_eps = 1e-08
  while ((psd == 0) & (s_eps <= 0.1)) {
    if (is_posdef((M + 1)/M * G - diag(s * (1 - s_eps), length(s)), 
                  tol = 1e-09)) {
      psd = 1
    }
    else {
      s_eps = s_eps * 10
    }
  }
  s = s * (1 - s_eps)
  s[s < 0] = 0
  if (verbose) 
    cat("done. \n")
  if (all(s == 0)) {
    warning("In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.", 
            immediate. = T)
  }
  return(s * diag(Sigma))
}


create.solve_asdp_M = function (Sigma, M = 1, max.size = 500, gaptol = 1e-06, maxit = 1000, 
                                verbose = FALSE) 
{
  stopifnot(isSymmetric(Sigma))
  if (ncol(Sigma) <= max.size) 
    return(create.solve_sdp_M(Sigma, M = M, gaptol = gaptol, 
                              maxit = maxit, verbose = verbose))
  if (verbose) 
    cat(sprintf("Dividing the problem into subproblems of size <= %s ... ", 
                max.size))
  cluster_sol = divide.sdp(Sigma, max.size = max.size)
  n.blocks = max(cluster_sol$clusters)
  if (verbose) 
    cat("done. \n")
  if (verbose) 
    cat(sprintf("Solving %s smaller SDPs ... \n", n.blocks))
  s_asdp_list = list()
  if (verbose) 
    pb <- utils::txtProgressBar(min = 0, max = n.blocks, 
                                style = 3)
  for (k in 1:n.blocks) {
    s_asdp_list[[k]] = create.solve_sdp_M(as.matrix(cluster_sol$subSigma[[k]]), 
                                          M = M, gaptol = gaptol, maxit = maxit)
    if (verbose) 
      utils::setTxtProgressBar(pb, k)
  }
  if (verbose) 
    cat("\n")
  p = dim(Sigma)[1]
  idx_count = rep(1, n.blocks)
  s_asdp = rep(0, p)
  for (j in 1:p) {
    cluster_j = cluster_sol$clusters[j]
    s_asdp[j] = s_asdp_list[[cluster_j]][idx_count[cluster_j]]
    idx_count[cluster_j] = idx_count[cluster_j] + 1
  }
  if (verbose) 
    cat(sprintf("Combinining the solutions of the %s smaller SDPs ... ", 
                n.blocks))
  tol = 1e-09
  maxitr = 1000
  gamma_range = c(seq(0, 0.1, len = 11)[-11], seq(0.1, 1, len = 10))
  gamma_opt = gtools::binsearch(function(i) {
    G = (M + 1)/M * Sigma - gamma_range[i] * diag(s_asdp)
    lambda_min = suppressWarnings(RSpectra::eigs(G, 1, which = "SR", 
                                                 opts = list(retvec = FALSE, maxitr = maxitr, tol = tol))$values)
    if (length(lambda_min) == 0) {
      lambda_min = min(eigen(G)$values)
    }
    lambda_min
  }, range = c(1, length(gamma_range)))
  s_asdp_scaled = gamma_range[min(gamma_opt$where)] * s_asdp
  options(warn = 0)
  if (verbose) 
    cat("done. \n")
  if (verbose) 
    cat("Verifying that the solution is correct ... ")
  if (!is_posdef((M + 1)/M * Sigma - diag(s_asdp_scaled, length(s_asdp_scaled)))) {
    warning("In creation of approximate SDP knockoffs, procedure failed. Knockoffs will have no power.", 
            immediate. = T)
    s_asdp_scaled = 0 * s_asdp_scaled
  }
  if (verbose) 
    cat("done. \n")
  s_asdp_scaled
}

x_knockoff <- function(x, M = 50){
  cor.G <- cor(x)
  n.sample <- nrow(x)
  n.G <- ncol(x)
  x.knockoff <- array(NA, dim = c(n.sample, n.G, M), dimnames = list(rownames(x), colnames(x), 1:M))  # 创建X.knockoff空矩阵
  
  knockoff <- GhostKnockoff.prelim(cor.G, M = M, method = 'asdp')
  P.each <- knockoff$P.each
  V.left <- as.matrix(knockoff$V.left)
  permute.index <- knockoff$permute.index
  permute.V.index <- rep(permute.index, M)
  
  set.seed('36336')
  Normal_50Studies <- as.matrix(V.left %*% matrix(rnorm(ncol(V.left) * n.sample), ncol(V.left), n.sample))  # 50 -> sample size
  Normal_50Studies[permute.V.index == 1, ] <- matrix(rnorm(sum(permute.index) * 
                                                             M * n.sample), sum(permute.index) * M, n.sample)
  
  
  for(i in 1:n.sample){
    Normal_k <- matrix(Normal_50Studies[, i], nrow = n.G)
    
    x_ik <- as.vector(P.each %*% t(x[i, , drop = F])) + Normal_k  # pxM, p:gene size, M:knockoff size 第i个样本的M个knockoff
    
    for(j in 1:M){
      x.knockoff[i,,j] <- x_ik[,j]
    }
  }
  return(x.knockoff)
}




#---- example ----
set.seed(10824)
beta <- 0.01
z <- rnorm(100)
x <- sin(z) + beta * rnorm(100)
y <- sin(x) + exp(z) + beta * rnorm(100)
data <- cbind(x, y, z)

KOCMI.net(data, k = 3, M = 40)
