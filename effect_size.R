#-----------------------------------------------------------------------
# Packages

library(SimCorMultRes)
library(tidyverse)
library(furrr)


#-----------------------------------------------------------------------
# beta_calc determines the value of the betas vector parameters

# j:input variables 
# ef_size: represents the effect of the input variables on the response variable (How much they explain)
# imp: predictive power of the input variables (How do they explain)
 
# A density function of the geometric distribution was used to determine the value of the betas
# For the geometric function is passed from the parameters: k and imp (prob of function)

# The effect control the accuracy
# The predictive power (imp) controls how input variables explain the answer
 
# The function returns a numerical vector with k elements which are the values
# parametric of the effects of input variables. The sum of the elements of the
# vector equals ef_size

beta_calc <- function(j = 10, ef_size = 0, imp = 0.5) {
  beta <- dgeom(0:(j - 1), prob = imp)
  beta <- ef_size * beta/sum(beta)
  return(beta)
}


#-----------------------------------------------------------------------

# Returns the matrix of input variables

# n: sample size.
# j: input variables.
# dist: distribution of input variables.
# rho: correlation between input variables in the Toeplitz structure.

# The function returns the matrix of input variables with n rows and j columns.
 
generate_predictors <- function(n = 100,
                                j = 10,
                                dist = "normal",
                                rho = 0) {
  correl <- toeplitz(x = c(1, rep(rho, j - 1)))
  X <- switch(
    EXPR = dist,
    "normal" = {
      qpar <- replicate(j,
                        list(mean = 0, sd = 1),
                        simplify = FALSE)
      rnorta(n,
             cor.matrix = correl,
             distr = rep("qnorm", j),
             qparameters = qpar)x
    },
    "binomial" = {
      qpar <- replicate(j,
                        list(size = 1, prob = 0.5),
                        simplify = FALSE)
      (rnorta(n,
              cor.matrix = correl,
              distr = rep("qbinom", j),
              qparameters = qpar) - 0.5)/0.5
    },
    "poisson" = {
      # (matrix(rpois(n * k, lambda = 10), ncol = k) - 10)/3.162278
      qpar <- replicate(j,
                        list(lambda = 10),
                        simplify = FALSE)
      (rnorta(n,
              cor.matrix = correl,
              distr = rep("qpois", j),
              qparameters = qpar) - 10)/3.162278
    })
  return(X)
}


#-----------------------------------------------------------------------

# Determines response variable, applies model and returns accuracy

# beta: effect vector of the input varible determined by calling the function beta_calc()
# X: matrix representing the input variable used to simulate the response variable

# The function returns the accuracy with simulated data according to the specifications.

simul_acc <- function(X = generate_predictors(n = 10,
                                              j = 3,
                                              dist = "normal",
                                              rho = 0),
                      beta = beta_calc(j = 3)) {
  eta <- X %*% beta
  y <- rbinom(nrow(X), size = 1, prob = binomial()$linkinv(eta))
  m0 <- glm.fit(x = cbind(1, X), y = y, family = binomial())
  acc <- sum((m0$fitted.values > 0.5) == y)/nrow(X)
  return(acc)
}


#-----------------------------------------------------------------------
# Function that resolves the effect size for a fixed input data matrix in all steps.
# Get the effect size value given simulation conditions.

# acc: fixed accuracy to determine the ef_size
# n: sample size (rows)
# j: input variables
# imp: predictive power of the input variables
# dist: distributions of input variables
# rho: correlation between the input variables
# n_repl: Number of independent replications


solve_ef_size <- function(acc = 0.9, 
                          n = 1000,  
                          j = 10,    
                          imp = 0.5, 
                          dist = "normal", 
                          rho = 0,         
                          n_repl = 10,    
                          verbose = TRUE) {

  f_obj <- function(ef_size, acc, j, imp, X) {
    beta_pars <- beta_calc(j = j,
                           ef_size = ef_size,
                           imp = imp)

    simul_acc <- replicate(n = n_repl,
                           simul_acc(X = X, beta = beta_pars))
    cat(".")
    acc - mean(simul_acc)
  }

  X <- generate_predictors(n = n,
                           j = j,
                           dist = dist,
                           rho = rho)

  root <- try(uniroot(f_obj,
                      interval = c(0, 30),
                      acc = acc,
                      j = j,
                      imp = imp,
                      X = X))
  
  acc_opt <- replicate(n = n_repl,
                       simul_acc(X = X,
                                 beta = beta_calc(j = j,
                                                  ef_size = root$root,
                                                  imp = imp)))
  acc_opt <- mean(acc_opt)
  if (verbose) {
    cat("\n")
    fmt <- paste("target acc: %0.3f",
                 "k: %d",
                 "imp: %0.3f",
                 "dist: %s",
                 "rho: %0.3f",
                 "EF_SIZE: %0.3f",
                 "MEAN ACC: %0.3f.\n",
                 sep = "; ")
    cat(sprintf(fmt = fmt,
                acc, j, imp, dist, rho, root$root, acc_opt))
  }
  if (inherits(root, "try-error")) {

    cat(sprintf("target acc: %0.3f; j: %d; imp: %0.3f; dist: %s; rho: %0.3f\n",
                acc, j, imp, dist, rho))
    return(data.frame(root = NA, f.root = NA, iter = NA,
                      init.it = NA, estim.prec = NA, acc = NA))
  } else {
    
    cbind(as.data.frame(root), acc_opt = acc_opt)
  }
}


#-----------------------------------------------------------------------
# Function that solves the effect size for various input data matrices.
# Get the effect size value given simulation conditions.
# acc: fixed accuracy to determine the ef_size
# n: sample size (row)
# j:  input variables 
# imp: predictive power of the input variables
# dist: distribution of input variables
# rho: correlation between input variables
# n_repl: Number of independent replications


solve_ef_size_final <- function(acc = 0.9, 
                                n = 1000,  
                                j = 10,    
                                imp = 0.5, 
                                dist = "normal", 
                                rho = 0,        
                                n_repl = 10,     
                                verbose = TRUE) {
  
  
  res <- replicate(n = 15,
                   simplify = FALSE,
                   solve_ef_size(acc = acc,
                                 n = n,
                                 j = j,
                                 imp = imp,
                                 dist = dist,
                                 rho = rho,
                                 n_repl = n_repl,
                                 verbose = TRUE))
  
  tb <- do.call(rbind, res)
  
  # Estimativa de `ef_size` pela média.
  root_m <- round(mean(tb$root), 2)
  my_tb <- cbind(acc = acc, n = n, j = j, imp = imp, dist = dist, rho = rho, root_m = as.data.frame(root_m))
  
  write.table(my_tb, file = "teste.txt", append = TRUE, sep = "\t",col.names = FALSE)
  return(my_tb)
  
  
}

#-----------------------------------------------------------------------
# All combinations of scenario

grid_ef_size <- as.data.frame(crossing(j = c(10, 50, 200),
                         n = c(500, 1000, 10000),
                         acc = c(0.9),
                         imp = c(0.2, 0.5, 0.8),
                         rho = c(0, 0.5, 0.8),
                         dist = c("normal", "binomial", "poisson")))


#-----------------------------------------------------------------------
#Parallel cod 
require(doParallel)
ncl <- detectCores() # Checa quantos núcleos existem na máquina
cl <- makeCluster(ncl)
registerDoParallel(cl) # Registra os clusters a serem utilizados

system.time({ 
  foreach(w=1:nrow(grid_ef_size),.packages=c("SimCorMultRes","furrr", "tidyverse")) %dopar% solve_ef_size_final(acc = grid_ef_size[w,3],
                                                                                                                n = grid_ef_size[w,2],
                                                                                                                j = grid_ef_size[w,1],
                                                                                                                imp = grid_ef_size[w,4],
                                                                                                                dist = grid_ef_size[w,6],
                                                                                                                rho = grid_ef_size[w,5])
})

stopCluster(cl)
