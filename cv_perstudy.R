cur.mod <- 6

library(doParallel)
cl <- makePSOCKcluster(4)
registerDoParallel(cl)

result.cv.perstudy <- expand.grid(study   = trials,
                                  model   = models[cur.mod],
                                  outcome = c("unfav", "mort"),
                                  measure = c("a", "b", "c", "brier"),
                                  est     = 0,
                                  lo      = 0,
                                  hi      = 0,
                                  nice    = 0)

pdf(paste("Output/results_cv/calibration_plots_studycv",models[cur.mod],".pdf", sep="_"))
set.seed(12)
for (i in models[cur.mod]){
  for (j in 1:length(trials)){
    
    ####TRAIN MODELS
    if(i == "LR"){
      mort_fit <- glm(mort ~ . , 
                      data = bothi[!bothi$trial%in%c(trials[j],"center"),vars_mort],
                      family = "binomial")
      unfav_fit <- glm(unfav ~ . , 
                       data = bothi[!bothi$trial%in%c(trials[j],"center"),vars_unfav],
                       family = "binomial")
    }else{
      if(i=="lasso"){
        dat.mort <- bothi[!bothi$trial%in%c(trials[j],"center"),vars_mort]
        
        y <- as.numeric(dat.mort$mort)
        y[y==1] <- 0
        y[y==2] <- 1
        mort_fit<-cv.glmnet(x = data.matrix(dat.mort[,-ncol(dat.mort)]),
                            y=  y, 
                            family="binomial",
                            alpha=1,
                            type.measure="deviance")
        
        dat.unfav <- bothi[!bothi$trial%in%c(trials[j],"center"),vars_unfav]
        
        y <- as.numeric(dat.unfav$unfav)
        y[y==1] <- 0
        y[y==2] <- 1
        unfav_fit<-cv.glmnet(x = data.matrix(dat.unfav[,-ncol(dat.unfav)]),
                             y=  y,
                             alpha=1, 
                             family="binomial",
                             type.measure="deviance")
      }else{
        if(i=="ridge"){
          dat.mort <- bothi[!bothi$trial%in%c(trials[j],"center"),vars_mort]
          y <- as.numeric(dat.mort$mort)
          y[y==1] <- 0
          y[y==2] <- 1
          mort_fit<-cv.glmnet(x = data.matrix(dat.mort[,-ncol(dat.mort)]),
                              y=  y, 
                              family="binomial",
                              alpha=0,
                              type.measure="deviance")
          
          dat.unfav <- bothi[!bothi$trial%in%c(trials[j],"center"),vars_unfav]
          ### no/yes into 0/1
          y <- as.numeric(dat.unfav$unfav)
          y[y==1] <- 0
          y[y==2] <- 1
          unfav_fit<-cv.glmnet(x = data.matrix(dat.unfav[,-ncol(dat.unfav)]),
                               y= y , 
                               family="binomial",
                               alpha=0,
                               type.measure="deviance")
        }else{
          if(i=="ranger"){
            tune_grid <- expand.grid(mtry=c(2,5, 10, 18),
                                     splitrule="gini",
                                     min.node.size=10)
          }
          if(i=="nnet"){
            tune_grid <- expand.grid(size=c(1,3,5,10,15),
                                     decay=c(0,0.001,0.01,0.1))
          }
          if(i =="svmLinear"){
            tune_grid <- expand.grid(C=c(0.5,0.6,0.7,0.8,0.9,1))
          }
          if(i =="gbm"){
            tune_grid <- expand.grid(n.trees=c(50,100,150,200,300),
                                     interaction.depth=c(1,2,3),
                                     shrinkage=c(0,0.01,0.1),
                                     n.minobsinnode=c(50,100,150))
          }
          mort_fit <- train(mort ~ . , 
                            data = bothi[!bothi$trial%in%c(trials[j],"center"),vars_mort], 
                            method = i, 
                            trControl = fitControl,
                            ## This last option is actually one
                            ## for gbm() that passes through
                            metric = "LogLik",
                            tuneGrid=tune_grid)
          unfav_fit <- train(unfav ~ . , 
                             data = bothi[!bothi$trial%in%c(trials[j],"center"),vars_unfav], 
                             method = i, 
                             trControl = fitControl,
                             ## This last option is actually one
                             ## for gbm() that passes through
                             metric = "LogLik",
                             tuneGrid=tune_grid)
          
        }#close no ridge
      }#close no lasso
    }#close no lr
    
    #####validate
    
    calc.val_abc <- function(data, measure=NULL,outcome=NULL){
      ##calculate predictions 
      if(outcome=="mort"){
        if(i=="LR"){
          pred <- predict(mort_fit, newdata=data,type="response")
        }else{
          if(i%in%c("lasso", "ridge")){
            data <- data[,colnames(dat.mort)]
            
            pred <- c(predict(mort_fit, 
                              newx=data.matrix(data[,-ncol(data)]),
                              s="lambda.min",type="response"))
          }else{
            pred <- predict(mort_fit, newdata=data,type="prob") 
          }
        }
      }else{
        if(i=="LR"){
          pred <- predict(unfav_fit, newdata=data,type="response")
        }else{
          if(i%in%c("lasso", "ridge")){
            data <- data[,colnames(dat.unfav)]
            
            pred <- c(predict(unfav_fit, 
                            newx=data.matrix(data[,-ncol(data)]),
                            s="lambda.min",type="response"))
          }else{
            pred <- predict(unfav_fit, newdata=data,type="prob") 
          }
        }
      }
      
      ##OBTAIN PREDICTIONS
      if(i%in%c(models[2:5])){ #for svm/gbm/nnet and rf
        if(outcome=="mort"){
          pred <- pred$dead
        }else{
          pred <- pred$unfav
        }
      }
      #IF predictions are 1 or 0, make them almost that so that LR works 
      pred[pred==1] <- 0.999
      pred[pred==0] <- 0.001
      
      ##FIT CALLIBRATION MODEL
      
      
      if(measure=="a"){
        vp <- rms::lrm(data[,outcome]~offset(qlogis(pred)))
        beta <- vp$coefficients[1]
        se <- sqrt(diag(vp$var))[1] 
        return(c(beta, beta-1.96*se, beta+1.96*se))  
      }
      if(measure=="b"){
        vp <- rms::lrm(data[,outcome]~qlogis(pred))
        beta <-  vp$coefficients[2]
        se <- sqrt(diag(vp$var))[2] 
        return(c(beta, beta-1.96*se, beta+1.96*se))
      }
      if(measure=="c"){
        if(outcome=="mort"){
          obs <- data$mort
        }else{
          obs <- data$unfav
        }
        res <- roc(obs~pred)
        
        
        #calibration plot
        obs <- as.numeric(obs)
        obs[obs==1] <- 0
        obs[obs==2] <- 1
        val.prob.ci.2(p = pred, y = obs, logit = "p",
                      main=paste(i, trials[j], outcome))
        return(c(pROC::ci.auc(res))[c(2,1,3)])
      }
      if(measure=="brier"){
        if(outcome=="mort"){
          obs <- as.numeric(data$mort=="dead")
        }else{
          obs <- as.numeric(data$unfav=="unfav")
        }
        calc_brier <- function(data, indices){
          d <- data[indices,]
          return(with(d, mean((as.numeric(pred)-obs)^2)))
        }
        brier_res <- boot::boot(data = data.frame(pred,obs),
                                statistic = calc_brier,
                                R = 500)
        return(quantile(brier_res$t, probs=c(0.5,0.025,0.975)))
      }
    }
    
    ###SAVE RESULTS
    for(msr in c("a","b","c","brier")){
      for(outcm in c("mort", "unfav")){
        result.cv.perstudy[result.cv.perstudy$study==trials[j]&
                             result.cv.perstudy$model==i&
                             result.cv.perstudy$measure==msr&
                             result.cv.perstudy$outcome==outcm,
                         c("est", "lo","hi")] <- calc.val_abc(
                           data = bothi[bothi$trial==trials[j],],
                           measure = msr,
                           outcome = outcm)
      }
    }
    
    ## add current results to the saved version
    save(result.cv.perstudy, file=paste("Output/results_cv/results.cs.perstudy",models[cur.mod],".RData", sep="_"))
  }
}
dev.off()
stopCluster(cl)
