#generate data
# set.seed(1234)
setwd("/mbSparse_Impute/Simulation")
control_real_data <- read.table("control.txt", header = TRUE, sep = "\t", row.names=1)
control_real_data[is.na(control_real_data)] <- 0
CA_real_data <- read.table("CA.txt", header = TRUE, sep = "\t", row.names=1)
CA_real_data[is.na(CA_real_data)] <- 0

union_rows <- union(row.names(control_real_data), row.names(CA_real_data))
data <- CA_real_data
existing_row_names <- rownames(data)
new_row_names <- union_rows[!(union_rows %in% existing_row_names)]
for (i in 1:length(new_row_names)) {
    data[nrow(data) + 1,] <- rep(0, ncol(data))
    rownames(data)[nrow(data)] <- new_row_names[i]
}
CA_real_data <- data[union_rows,]

data <- control_real_data
existing_row_names <- rownames(data)
new_row_names <- union_rows[!(union_rows %in% existing_row_names)]
for (i in 1:length(new_row_names)) {
    data[nrow(data) + 1,] <- rep(0, ncol(data))
    rownames(data)[nrow(data)] <- new_row_names[i]
}
control_real_data <- data[union_rows,]

control_real_data <- t(control_real_data)
CA_real_data <- t(CA_real_data)
real_data <- rbind(control_real_data, CA_real_data)
row.names(real_data) <- sub("\\..*$", "", row.names(real_data))
real_data <- log10(real_data + 1.01)

library(xlsx)
meta_real_data <- read.xlsx("metadata.xlsx", sheetName = "Sheet1", index_col=0)
row.names(meta_real_data) <- meta_real_data[,1]
meta_real_data <- meta_real_data[, 2:(ncol(meta_real_data)-1)]
meta_real_data <- meta_real_data[order(match(row.names(meta_real_data), row.names(real_data))), ]
meta_real_data <- meta_real_data[row.names(real_data),]



#droupt zero
real_data_zi_rate <- apply(real_data, 2, FUN = function(x){
  sum(x < (log10(1.01) + 0.01))/length(x)
})
chosen_taxa <- which(real_data_zi_rate < 0.4)
simulated2 <- real_data[,chosen_taxa]
dim(simulated2)

# introduce zeros based on real data.
gamma_norm_mix <- function(y, X){
  loglik <- function(p, alpha, beta, cov_par, var1, X, y){
    n = length(y)
    lkval <- 0
    fgam <- dgamma(y, shape = alpha, rate = beta)
    for(i in 1:n){
      if(!is.vector(X)){
        lkval <- lkval + log10( p*fgam[i] + (1-p)*dnorm(y[i], mean = X[i,] %*% cov_par, sd = sqrt(var1)))
      }else{
        lkval <- lkval + log10( p*fgam[i] + (1-p)*dnorm(y[i], mean = X[i] * cov_par, sd = sqrt(var1)))
      }
      
    }
    return(lkval)
  }
  n = length(y)
  alpha_init <- 1
  beta_init <- 10
  p_init <- 0.5
  cov_par_init <- solve(t(X) %*% X) %*% t(X) %*% y
  var_init <- t(y - X %*% cov_par_init) %*% (y - X %*% cov_par_init) / n
  
  alpha_t <- alpha_init
  beta_t <- beta_init
  cov_par_t <- cov_par_init
  var_t <- var_init
  p_t <- p_init
  
  #update gamam param
  #Wei's Method
  ### root-finding equation
  fn = function(alpha, target){
    log(alpha) - digamma(alpha) - target
  }
  update_gmm_pars = function(x, wt){
    if(max(wt) > 0.00001){
      tp_s = sum(wt)
      tp_t = sum(wt * x)
      tp_u = sum(wt * log(x))
      tp_v = -tp_u / tp_s - log(tp_s / tp_t)
      if (tp_v <= 0){
        alpha = 20
      }else{
        alpha0 = (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / 12 / tp_v
        if (alpha0 >= 20){alpha = 20
        }else{
          alpha = uniroot(fn, c(0.9, 1.1) * alpha0, target = tp_v,
                          extendInt = "yes")$root
        }
      }
      ## need to solve log(x) - digamma(x) = tp_v
      ## We use this approximation to compute the initial value
      beta = tp_s / tp_t * alpha
    }else{
      alpha = 0.001
      beta = 1000
    }
    
    return(c(alpha, beta))
  }
  #convergence criteria
  flag = TRUE
  maxitr = 300
  itr = 0
  while(flag){
    #E_step
    mean_t <- X %*% cov_par_t
    n = length(y)
    a_hat_t <- rep(0, n)
    dg_t <- dgamma(y, shape = alpha_t, rate = beta_t)
    for(i in 1:n){
      if(dg_t[i] == 0){
        a_hat_t[i] = 0
      }else{
        a_hat_t[i] <- p_t * dg_t[i]/(p_t * dg_t[i]+ (1-p_t)*dnorm(y[i], mean = mean_t[i], sd = sqrt(var_t)))
      }
    }
    #maximization
    #fit p
    p_t1 <- sum(a_hat_t)/n
    X_tilta <- sqrt(1-a_hat_t) * X
    y_tilta <- sqrt(1-a_hat_t) * y
    #fit normal
    out <- tryCatch(
      {
        # Just to highlight: if you want to use more than one
        # R expression in the "try" part then you'll have to
        # use curly brackets.
        # 'tryCatch()' will return the last evaluated expression
        # in case the "try" part was completed successfully
        cov_par_t1 <- solve(t(X_tilta) %*% X_tilta) %*% t(X_tilta) %*% y_tilta
        TRUE
      },
      error=function(cond) {
        FALSE
      }
    )
    if(!out){
      return(list("d" = y < log10(1.01) + 10^(-3) ))
    }
    cov_par_t1 <- solve(t(X_tilta) %*% X_tilta) %*% t(X_tilta) %*% y_tilta
    var_t1 <- sum((1 - a_hat_t) * (y - X %*% cov_par_t)^2) / sum(1-a_hat_t)
    #fit gamma
    par_gamm <- update_gmm_pars(x = y, wt = a_hat_t)
    alpha_t1 = par_gamm[1]
    beta_t1 <- par_gamm[2]
    loglik1 <- loglik(p = p_t, alpha = alpha_t, beta = beta_t, cov_par = cov_par_t, var1 = var_t, X = X, y = y)
    loglik2 <- loglik(p = p_t1, alpha = alpha_t1, beta = beta_t1, cov_par = cov_par_t1, var1 = var_t1, X = X, y = y)
    if((abs(loglik1 - loglik2)) < 0.05 || itr > maxitr){
      flag = FALSE
    }else{
      alpha_t <- alpha_t1
      beta_t <- beta_t1
      cov_par_t <- cov_par_t1
      var_t <- var_t1
      p_t <- p_t1
      itr = itr + 1
    }
  }
  #fit a normal curve
  eta_hat <- cov_par_init
  omega_hat <- var_init
  #calculate new likelihood
  norm_log_lik = 0
  for(i in 1:n){
    if(!is.vector(X)){
      norm_log_lik = norm_log_lik + log10(dnorm(y[i], mean = X[i,] %*% eta_hat, sd = sqrt(omega_hat)))
    }else{
      norm_log_lik = norm_log_lik + log10(dnorm(y[i], mean = X[i] * eta_hat, sd = sqrt(omega_hat)))
    }
  }
  Dev = -2 * norm_log_lik - (-2 * loglik2)
  judgement = pchisq(Dev, df = 3, lower.tail = FALSE) < 0.05
  
  if(!judgement || (alpha_t/beta_t > 1)){
    p_t = 0
    a_hat_t <- rep(0,n)
  }
  return(list("p" = p_t, "alpha" = alpha_t, "beta" = beta_t, "cov_par" = cov_par_t, "var" = var_t, "d" = a_hat_t, "eta" = eta_hat, "omega" = omega_hat, "Deviance" = Dev))
}

otable <- real_data
meta_tab <- meta_real_data
dim(meta_tab)
dim(real_data)
# read in the gamma_norm_mix function
mean_record <- c()
percentage_record <- c()
D_vale_record <- c()
beta_record <- list()
for(j in 1:dim(otable)[2]){
  result <- gamma_norm_mix(otable[,j], data.matrix(meta_tab))
  beta_record[[j]] = result$cov_par
  mean_record <- c(mean_record, mean(otable[which(result$d < 0.5),j]))
  percentage_record <- c(percentage_record, sum(result$d > 0.5)/dim(otable)[1])
  D_vale_record <- c(D_vale_record, result$d)
}
#filter <- which(is.nan(mean_record) | mean_record < log10(1.01) + 0.001)

# filter out the nan values.
# plot(mean_record, percentage_record)

# build a map between the mean and percentage of missing
missing_rate <- function(mean, emp_mean, emp_miss){
  win_len <- (range(emp_mean)[2] - range(emp_mean)[1])/5
  mean_up <- mean + win_len 
  mean_lo <- mean - win_len
  sample(emp_miss[emp_mean > mean_lo & emp_mean < mean_up], 1)
}
# missing_rate(1, mean_record, percentage_record)

y_sim <- simulated2
y_preserve <- y_sim
col_mean <- colMeans(y_sim)
zero_rate <- unlist(lapply(col_mean, FUN = function(x){
#   print(x)
  return(missing_rate(x, mean_record, percentage_record))
}))
#zero_rate[zero_rate > 0.9] = 0.9
zero_rate[zero_rate > 0.8] = 0.8
n = dim(y_sim)[1]
m = dim(y_sim)[2]
zero_mat <- matrix(NA, nrow = n, ncol = m)
for(i in 1:m){
  zero_mat[,i] = rbinom(n,1, 1-zero_rate[i])
}
sim_tab_zi = y_sim * zero_mat
sim_tab_zi[sim_tab_zi < log10(1.01)+1e-6] = log10(1.01)

sim_tab_zi <- 10^(sim_tab_zi) - 1.01
sim_tab_zi_trans <- t(sim_tab_zi)
write.csv(sim_tab_zi_trans, "Health_CA_zero_2.csv")
simulated2_trans <-t(10^(simulated2) - 1.01)
write.csv(simulated2_trans, "Health_CA_full_2.csv")
meta_real_data$Age <- c(rep("Health", 135), rep("Disease", nrow(meta_real_data)-135))
colnames(meta_real_data)[1] <- "Gruop"
meta_real_data <- meta_real_data[, 1, drop = FALSE]
meta_real_data <- t(meta_real_data)
write.csv(meta_real_data, "Health_CA_Metadata.csv")



#mbimpute
library(mbImpute)
library(glmnet)
library(Matrix)
data_zero <- read.csv("Health_CA_zero_2.csv", row.names = "X")
data_zero <- t(data_zero)
mbimpute_mat <- mbImpute(otu_tab = data_zero)
mbimpute_mat = mbimpute_mat$imp_count_mat_lognorm
data_full <- read.csv("Health_CA_full_2.csv", row.names = "X")
data_full_eval <- t(data_full)
scale <- rowSums(data_full_eval) / 10^6
data_full_eval <- data_full_eval / scale
data_full_eval <- log10(data_full_eval + 1)
for(j in 1:dim(mbimpute_mat)[2]){
  mbimpute_mat[,j] <- mbimpute_mat[,j] * max(data_full_eval[,j]) / max(mbimpute_mat[,j])
}
size <- dim(data_full_eval)[1] * dim(data_full_eval)[2] 
mse <- sum((data_full_eval- mbimpute_mat)^2) / size

pearson_coefficient_sum <- 0
for (i in 1:dim(data_full_eval)[1]) {
    pearson_coefficient_sum <- pearson_coefficient_sum + cor(data_full_eval[i,], mbimpute_mat[i,])
}
pearson_coefficient_mean <- pearson_coefficient_sum / dim(data_full_eval)[1]
pearson_coefficient_mean
cat(sprintf("mse=%f \n", mse))
cat(sprintf("pearson=%f \n", pearson_coefficient_mean))


#orginal
data_zero <- read.csv("Health_CA_zero_2.csv", row.names = "X")
data_zero <- t(data_zero)
scale <- rowSums(data_zero) / 10^6
data_zero <- data_zero / scale
data_zero <- log10(data_zero + 1)
data_full <- read.csv("Health_CA_full_2.csv", row.names = "X")
data_full_eval <- t(data_full)
scale <- rowSums(data_full_eval) / 10^6
data_full_eval <- data_full_eval / scale
data_full_eval <- log10(data_full_eval + 1)
size <- dim(data_full_eval)[1] * dim(data_full_eval)[2] 
mse <- sum((data_full_eval- data_zero)^2) / size
pearson_coefficient_sum <- 0
for (i in 1:dim(data_full_eval)[1]) {
    pearson_coefficient_sum <- pearson_coefficient_sum + cor(data_full_eval[i,], data_zero[i,])
}
pearson_coefficient_mean <- pearson_coefficient_sum / dim(data_full_eval)[1]
cat(sprintf("mse=%f \n", mse))
cat(sprintf("pearson=%f \n", pearson_coefficient_mean))



#SAVER
library(SAVER)
data_full <- read.csv("Health_CA_full_2.csv", row.names = "X")
data_full_eval <- t(data_full)
scale <- rowSums(data_full_eval) / 10^6
data_full_eval <- data_full_eval / scale
data_full_eval <- log10(data_full_eval + 1)
data_zero <- read.csv("Health_CA_zero_2.csv", row.names = "X")
saver_mat <- saver(data_zero, ncores = 1, estimates.only = TRUE)
saver_mat_eval <- log10(saver_mat + 1)
saver_mat_eval <- t(saver_mat_eval)
for(j in 1:dim(saver_mat_eval)[2]){
  saver_mat_eval[,j] <- saver_mat_eval[,j] * max(data_full_eval[,j]) / max(saver_mat_eval[,j])
}
size <- dim(data_full_eval)[1] * dim(data_full_eval)[2] 
mse <- sum((data_full_eval- saver_mat_eval)^2) / size

pearson_coefficient_sum <- 0
for (i in 1:dim(data_full_eval)[1]) {
    pearson_coefficient_sum <- pearson_coefficient_sum + cor(data_full_eval[i,], saver_mat_eval[i,])
}
pearson_coefficient_mean <- pearson_coefficient_sum / dim(data_full_eval)[1]
cat(sprintf("mse=%f \n", mse))
cat(sprintf("pearson=%f \n", pearson_coefficient_mean))


#GE-Impute
data_full <- read.csv("Health_CA_full_2.csv", row.names = "X")
data_full_eval <- t(data_full)
scale <- rowSums(data_full_eval) / 10^6
data_full_eval <- data_full_eval / scale
data_full_eval <- log10(data_full_eval + 1)
GEImpute <- t(read.csv("geimpute_Health_CA_impute_2.csv", row.names = "X"))
GEImpute <- log10(GEImpute + 1)
for(j in 1:dim(GEImpute)[2]){
  GEImpute[,j] <- GEImpute[,j] * max(data_full_eval[,j]) / max(GEImpute[,j])
}
size <- dim(data_full_eval)[1] * dim(data_full_eval)[2] 
mse <- sum((data_full_eval- GEImpute)^2) / size
pearson_coefficient_sum <- 0                                  
for (i in 1:dim(data_full_eval)[1]) {
    pearson_coefficient_sum <- pearson_coefficient_sum + cor(data_full_eval[i,], GEImpute[i,])
}
pearson_coefficient_mean <- pearson_coefficient_sum / dim(data_full_eval)[1]
cat(sprintf("mse=%f \n", mse))
cat(sprintf("pearson=%f \n", pearson_coefficient_mean))


#mbSparse-Impute
data_full <- read.csv("Health_CA_full_2.csv", row.names = "X")
data_full_eval <- t(data_full)
scale <- rowSums(data_full_eval) / 10^6
data_full_eval <- data_full_eval / scale
data_full_eval <- log10(data_full_eval + 1)
mbSparseImpute <- t(read.csv("Health_CA_impute_normalize_2.csv", row.names = "X"))
size <- dim(data_full_eval)[1] * dim(data_full_eval)[2] 
mse <- sum((data_full_eval- mbSparseImpute)^2) / size
pearson_coefficient_sum <- 0                                  
for (i in 1:dim(data_full_eval)[1]) {
    pearson_coefficient_sum <- pearson_coefficient_sum + cor(data_full_eval[i,], mbSparseImpute[i,])
}
pearson_coefficient_mean <- pearson_coefficient_sum / dim(data_full_eval)[1]
cat(sprintf("mse=%f \n", mse))
cat(sprintf("pearson=%f \n", pearson_coefficient_mean))

#softImpute
library(softImpute)
si_impute <- function(sim_tab){
  set.seed(12345)
  y_sim = sim_tab
  # add filter
  filter <- lapply(X = 1:ncol(y_sim), FUN = function(col_i){
    y = y_sim[,col_i]
    n = length(y)
    nz <- sum(y <= (log10(1) + 1e-6))
    pz = 1 - nz/n
    test = pz - 1.96 * sqrt(pz * (1-pz)/n)
    if(nz == n || test <= 0){
      return(0)
    }else{
      return(1)
    }
  })
  y_imp <- y_sim
  #perform imputation on the rest
  filter_vec <- which(unlist(filter) == 1)
  y_sim = y_sim[, filter_vec]
  
  na_idx <- y_sim < log10(1) + 0.001
  y_sim[y_sim < log10(1) + 0.001] = NA
  
  # y_sim_cv <- unlist(y_sim)
  # y_sim_validate <- matrix(y_sim_cv, nrow(y_sim), ncol = ncol(y_sim))
  # identical(y_sim_validate, y_sim)
  y_sim_cv <- unlist(y_sim)
  na_intro <- sample(which(!is.na(y_sim_cv)), floor(sum(!is.na(y_sim_cv))/10))
  y_sim_cv_intro <- y_sim_cv
  y_sim_cv_intro[na_intro] = NA
  y_sim_cv_intro <- matrix(y_sim_cv_intro, nrow = nrow(y_sim), ncol = ncol(y_sim))
  
  j = 1
  se = 1e10
  for(i in 1:5){
    si_cv_1 <- softImpute(y_sim_cv_intro, rank.max = i, lambda = 0)
    y_imp_cv <- complete(y_sim_cv_intro, si_cv_1)
    y_sim_vali <- as.vector(y_imp_cv)
    se2 <- sum((y_sim_cv[na_intro] - y_sim_vali[na_intro])^2)
    print(se2)
    if(se2 < se){
      se = se2
      j = i
    }
  }
  print(j)
  si1 <- softImpute(y_sim_cv_intro, rank.max = 10, lambda = 0, trace.it = TRUE)
  impute_mat <- complete(y_sim_cv_intro, si1)
  
  y_imp[, filter_vec] = impute_mat
  
  return(y_imp)
}

data_full <- read.csv("Health_CA_full_2.csv", row.names = "X")
data_full_eval <- t(data_full)
scale <- rowSums(data_full_eval) / 10^6
data_full_eval <- data_full_eval / scale
data_full_eval <- log10(data_full_eval + 1)
data_zero <- read.csv("Health_CA_zero_2_2.csv", row.names = "X")
softImpute <- si_impute(data_zero)
softImpute[softImpute < 0] = 0
softImpute <- t(softImpute)
softImpute <- log10(softImpute+1)
for(j in 1:dim(softImpute)[2]){
  softImpute[,j] <- softImpute[,j] * max(data_full_eval[,j]) / max(softImpute[,j])
}
size <- dim(data_full_eval)[1] * dim(data_full_eval)[2] 
mse <- sum((data_full_eval- softImpute)^2) / size
pearson_coefficient_sum <- 0                                  
for (i in 1:dim(data_full_eval)[1]) {
    pearson_coefficient_sum <- pearson_coefficient_sum + cor(data_full_eval[i,], softImpute[i,])
}
pearson_coefficient_mean <- pearson_coefficient_sum / dim(data_full_eval)[1]
cat(sprintf("mse=%f \n", mse))
cat(sprintf("pearson=%f \n", pearson_coefficient_mean))