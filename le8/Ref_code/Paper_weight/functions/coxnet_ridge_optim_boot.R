####################################################################
#### Script Name: coxnet_optim
#### Author: Julia Carrasco Zanini Sanchez
#### Description: Function to run coxnet optimization and testing  
#### Date : 17/07/2022

## Citation:
## From: Carrasco-Zanini, J. et al. Proteomic signatures improve risk prediction for common and rare diseases. Nat Med. 30, 2489–2498 (2024).

coxnet.optim.r <- function(Train.data, pred.vec, Train.surv.data, times, Test.data, Test.surv.data){
  ## Train.data : data set containing predictor variables used for model optimization 
  ## pred.vec : vector containing names of predictor variable 
  ## Train.surv.data : survival object from optimization set
  ## times : number of boostrap samples to compute confidence intervals around the C-index
  ## Test.data : data set containing predictor variables for the test set
  ## Test.surv.data : survival object from test set
  
  if(sum(!complete.cases(Train.data[,pred.vec]))!=0)
    stop("Missing training values for predictors")
  if(sum(!complete.cases(Test.data[,pred.vec]))!=0)
    stop("Missing testing values for predictors")
  if(sum(is.na(Train.surv.data))!=0)
    stop("Missing training survival values")
  if(sum(is.na(Test.surv.data))!=0)
    stop("Missing testing survival values")
  
  require(caret)
  require(glmnet)
  require(ROSE)
  
  ## run lasso, only if two or more predictors
  if(length(pred.vec) > 1){
    ## optimize
    set.seed(354)
    las.morb <- cv.glmnet(as.matrix(Train.data[,pred.vec, drop=F]),
                          as.matrix(Train.surv.data),
                          family="cox",
                          alpha=0,
                          lambda = 10^-seq(10,.25,-.25),
                          nfolds = 5,
                          sampling = "rose",
                          parallel = T)
    
    ## parameters for the best model
    coxnet.opt  <- las.morb$glmnet.fit
    s.opt       <- las.morb$lambda.min
    opt.coef    <- coef(coxnet.opt, s = s.opt)
  }else{
    ## generate the data to do so
    tmp        <- cbind(Train.data[,pred.vec, drop=F], as.matrix(Train.surv.data))
    ## assign names
    names(tmp) <- c(pred.vec, "time", "bin")
    ## compute simple Cox model
    coxnet.opt <- coxph(as.formula(paste0("Surv(time, bin) ~ ", pred.vec)), tmp, ties = "breslow")
    s.opt      <- NA
    opt.coef   <- coef(coxnet.opt)
  }
  
  ## create samples for boot-strapping
  set.seed(42)
  boot.samp   <- lapply(1:times, function(x) sample(1:nrow(Test.data), nrow(Test.data), replace = T))
  
  ## Testing - C-index by boostraping
  boot.cindex <- mclapply(boot.samp, function(x){
    ## compute the prediction
    if(class(coxnet.opt)[1] == "coxph"){
      print("single predictor")
      ## compute predictions (https://stats.stackexchange.com/questions/48298/computing-c-index-for-an-external-validation-of-a-cox-ph-model-with-r)
      tmp.pred <- as.numeric(predict(coxnet.opt, type =  "lp", 
                                     newdata = Test.data[ x , pred.vec, drop=F]))
      tmp.pred <- rcorr.cens(tmp.pred, Test.surv.data[x])[2]
      return((-tmp.pred+1)/2)
    }else{
      ## compute predictions
      tmp.pred <- as.numeric(predict(coxnet.opt, type = "response", 
                                     newx = as.matrix(Test.data[ x , pred.vec, drop=F]),
                                     s = s.opt))
      ## compute index
      return(glmnet::Cindex(tmp.pred, as.matrix(Test.surv.data[x, ])))
    }
  }, mc.cores=10)
  ## unlist
  boot.cindex <- unlist(boot.cindex)
  print(head(boot.cindex))

  
  if(class(coxnet.opt)[1] != "coxph"){
    ## linear predictor for the entire test set 
    pred.lp <- as.numeric(predict(coxnet.opt, 
                                  type = "link", 
                                  newx = as.matrix(Test.data[, pred.vec, drop=F]),
                                  s = s.opt))
    
    ## relative risk estimates for the test set 
    pred.risk <- as.numeric(predict(coxnet.opt, 
                                    type = "response", 
                                    newx = as.matrix(Test.data[, pred.vec, drop=F]),
                                    s = s.opt))
  }else{
    ## linear predictor for the entire test set 
    pred.lp <- as.numeric(predict(coxnet.opt, 
                                  type = "lp", 
                                  newx = as.matrix(Test.data[, pred.vec, drop=F])))
    
    ## relative risk estimates for the test set 
    pred.risk <- as.numeric(predict(coxnet.opt, 
                                    type = "risk", 
                                    newx = as.matrix(Test.data[, pred.vec, drop=F])))
  }

  ## aggregate all results
  coxnet.res <- list(
    glmnet.opt = coxnet.opt,
    lambda.opt = s.opt,
    opt.coefficients = opt.coef,
    Cindex.vec = boot.cindex,
    Cindex.mean = mean(boot.cindex,na.rm=T),
    ci.low = quantile(boot.cindex, 0.025,na.rm=T),
    ci.upp = quantile(boot.cindex, 0.975,na.rm=T),
    linear.predictor = pred.lp,
    relative.risk = pred.risk
  )
  return(coxnet.res)
}
