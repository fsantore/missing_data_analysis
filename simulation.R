#--------------------------------------------------------------------------------
#Packages
require(SimCorMultRes) #rnorta
require(caret) #create partition 
require(bindata)
require(ggplot2)
library(tidyverse)

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


#--------------------------------------------------------------------------------

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
             qparameters = qpar)
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

# The function simul_model(train, test) adjusts a glm to the training data 
# and performs the prediction on the test data
# returning the accuracy amd some others measurement of the model 

# Metrics : Accuracy, Specificity, Precision, Recall, F1

simul_model <- function(train, test){

  
  fit <- glm(Y~. ,family = binomial(), data = train)
  pred <- predict(fit, test, type = "response")
  
  cm <- confusionMatrix(data = as.factor(as.numeric(pred>0.5)), reference = as.factor(test$Y))

  ac <- cm$overall['Accuracy'] 
  spe <- cm$byClass['Specificity']
  prec <- cm$byClass['Precision']
  rec <- cm$byClass['Recall']
  f1 <- cm$byClass['F1']

  metricas<-c(ac,spe,prec,rec,f1) 
  return(metricas)
  
}


#-----------------------------------------------------------------------

#The missing function inserts artificially missing values into the original data,
#by different mechanisms and in different amount 

# q: amount of missing
# t: mechanism generator of missing (COM (complete), MCAR, MAR, MNAR)
# r: correlation parameter fized to multivariate Bernoulli distribution


missing<-function(data,q,t,r){
  
  l <- nrow(data)
  c <- ncol(data)
  
  if (t == "MCAR"){
    
    COR <- toeplitz(c(1, rep(r, (c-1)))) #Matriz de correlação
    
    if (r == 0){
      
      m <- matrix(rep(1,l*c),nrow = l, ncol = c)
      p_missing <- rbinom(l, p = 0.1, size = 1)
      p_sum <- sum(p_missing)
      
      index <- which(p_missing == 1)
      
      for (i in index) {
        
        m[i,] <- rbinom(c, p = (0.9/p_sum), size = 1)
      }
   
      
    }else if( r != 0 ){
      
      m <- rmvbin(l, margprob = c(rep(q,c)), bincorr = COR)
      
    }
    
    m[m == 0] <- NA
    
    
    b <- which(is.na(m),arr.ind =TRUE)
    e <- nrow(b)
    z <- 1
    
    while(z<=e){
      data[b[z,1],b[z,2]] <- NA
      z <- z+1
    }
    
  } else if (t == "MAR"){
  
    COR <- toeplitz(c(1, rep(r, (c-1)))) #Matriz de correlação
    p_missing <- sort(runif(l, min = 0.5, max = 1))
    p_mean <- mean(p_missing)
    p_target <- 1-q
    fc <- p_mean/p_target
    p_missing <- p_missing/fc
    m <- c()
    
    for (w in p_missing){
      
      mis <- rmvbin(1, margprob = c(rep(1-w,c)), bincorr = COR)
      m <- rbind(m,mis)
      
    }
    
    m[m == 0] <- NA
    
    b <- which(is.na(m),arr.ind =TRUE)
    e <- nrow(b)
    z <- 1
    
    data[] <- data[order(data[,1]), ]
    
    while(z<=e){
      data[b[z,1],b[z,2]] <- NA
      z <- z+1
    }
    
  }else if(t == "MNAR") {
    
    COR <- toeplitz(c(1, rep(0.8, ceiling((c-1)/2)),rep(0.5, floor((c-1)/2)))) #Matriz de correlação
    m <- rmvbin(l, margprob = c(rep(q,c)), bincorr = COR)
    
    m[m == 0] <- NA
    
    
    b <- which(is.na(m),arr.ind =TRUE)
    e <- nrow(b)
    z <- 1
    
    while(z<=e){
      data[b[z,1],b[z,2]] <- NA
      z <- z+1
    }
    
  }

  
  return(data)
  
}


#--------------------------------------------------------------------------------

# The function solve_simulation get the accuracy value given simulation conditions.

# n: sample size (row)
# j:  input variables 
# i: predictive power of the input variables
# dist: distribution of input variables
# p: correlation between input variables
# ef: effect of the input variables on the response variable (How much they explain)

solve_simulation <- function(j, n, dist, p, i, ef){
  
  type <- c("MCAR","MAR","MNAR")

  rho <- c(0.5,0.8)
  amount <- c(0.9,0.7)
  
  results <- c()
  resultsFinal <- c()
  s <- 0
  u <- 0
  while(u < 10){
    
    X <- generate_predictors( n = n , j = j, dist = dist, rho = p)
    beta <- beta_calc (j = j, ef_size = ef, imp = i)
    eta <- X %*% beta
    
    while(s < 15){
      
      
      Y <- rbinom(nrow(X), size = 1, prob = binomial()$linkinv(eta))
      data <- data.frame(cbind(X, Y))
      
      intrain <- createDataPartition(y = data$Y, p= 0.7, list = FALSE)
      
      train <- data[intrain,]
      test <- data[-intrain,]
      met <- simul_model(train, test)
      
      res <- c(j, n, dist, p, i, ef,"COM", 1, 1, met)
      results <- rbind(results,res)
      
      for(t in type){
        for( r in rho){
          for(q in amount){
            
            M <- missing(data[,-ncol(data)],q,t,r)
            M$Y <- data$Y
            trainM <- M[intrain,]
            testM <- M[-intrain,]
            met <- simul_model(trainM, testM)
            
            res <- c(j, n, dist, p, i, ef, t, r, q, met)
            results <- rbind(results,res)
          }
        }
      }
      
      s <- s + 1    
    }
    resultsFinal <- rbind(resultsFinal, results)
    u <- u + 1
  
  }
  
  str(resultsFinal)

  
  colnames(resultsFinal)[1:9] <- c("col", "row", "distribution", "corr", "power","effect", "type", "corrT", "amount")

  name<-paste("~/",j, n, dist, p, i, ef,".RData") 
  
  saveRDS(resultsFinal, file= name )#salvando 
  
}


#--------------------------------------------------------------------------------
#Parallel cod 
# The experimentos data frame contain the results of calibration effect size
# given simulation conditions.

require(doParallel)
ncl <- detectCores() 
cl <- makeCluster(ncl)
registerDoParallel(cl)

system.time({ 
  foreach(w=1:nrow(experimentos),.packages=c("SimCorMultRes","caret","bindata")) %dopar% solve_simulation (j = experimentos[w,1], n = experimentos[w,2], dist = experimentos[w,3], p = experimentos[w,4], i = experimentos[w,5], ef = experimentos [w,6])
})

stopCluster(cl)

