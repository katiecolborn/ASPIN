#install.packages(c("lubridate", "glmnet", "caret", "DMwR", "pROC", "sas7bdat", "rmarkdown"))
library(lubridate)
library(glmnet)
library(caret)
library(DMwR)
library(pROC)
library(sas7bdat)
library(glmnet)
#install.packages("knockoff")
library(knockoff)
library(knitr)
library(doMC)

Generate_statistics_res_SSI <- function(nfolds=10,train_date_end="2017-03-01",column_sums=c(5,10,15,20),fdr_col=c(0,1,0.12,0.15,0.2),offset=1,parallelcores ="max",data_type="Rdata",
                                        data_path="//data.ucdenver.pvt/dept/SOM/ACCORDS/PiFolders/PI_Colborn/UCHealthSOR/Data/Data_20210113/final_analytical_datasets_0405/SSI/",
                                        fileName ="final_data_labpanel_SSIUTIPNEU",save_res=TRUE,save_dir="C:/Users/ASRock/OneDrive - The University of Colorado Denver/Katie Project/2021_1109/SSI/res/")
{
  if (data_type=="CSV"| data_type=="csv"| data_type=="Csv")
  {AnyComp <-paste(data_path,fileName,".csv",sep = "")
  AnyComp <- AnyComp[,-1]}
  if (data_type=="Rdata"| data_type=="rdata"| data_type=="RDATA")
  {load(paste(data_path,fileName,".Rdata",sep = ""))
    AnyComp <- final_data_loinc_v3}
  
  detectedCores <- detectCores()
  numCores <- parallelcores
  
  if (parallelcores =="max") {
    warning("Parallel computing use all cores")
    numCores <- detectCores()
  }
  
  if (numCores > detectedCores)
  {
    warning("The # of cores you entered is larger than you have on your computer, max # of cores being used")
    numCores <- detectedCores
  }
  #colnames(AnyComp)
  #dim(AnyComp)
  
  #table(AnyComp$score18)
  #sum(table(AnyComp$score18)) ### no missing values ###
  #AnyComp$Score18Rate <- ifelse(is.na(AnyComp$Score18Rate)==T, 0, AnyComp$Score18Rate)
  #table(AnyComp$E11.9)
  combres <- list()
  comb_res_fdr_colnum <- list()
  store_save_coef_list <- list()
  store_coef_comb <- list()
  ### Read in AMA (not used for now) ###
  #a$data = substr(a$data,1,nchar(a$data)-3)
  AMA_path <- substr(data_path,1,nchar(data_path)-49)
  print("loading AMA dataset")
  #AMA <- read.sas7bdat(file = "//data.ucdenver.pvt/dept/SOM/ACCORDS/PiFolders/PI_Colborn/UCHealthSOR/Data/AMA/cptcodes2020.sas7bdat")
  AMA <- read.sas7bdat(file = paste(AMA_path,"AMA/cptcodes2020.sas7bdat",sep = ""))
  #dim(AMA)
  #colnames(AMA)
  
  AMA_merge <- AMA[!duplicated(AMA$CPT),]
  
  #table(AnyComp$asaclas)
  #sum(table(AnyComp$asaclas)) ### no missing values ###
  
  ### check asa class and score 18 ###
  #prop.table(table(AnyComp$asaclas, AnyComp$score18),1)
  
  AnyComp$asa_bin <- ifelse(AnyComp$asaclas %in% c("4-Life Threat", "4-LifeThreat","5-Moribund"), 1, 0)
  
  #class(AnyComp$as_date_operation)
  AnyComp$as_date_operation <- as.Date(AnyComp$as_date_operation)
  #max(AnyComp$as_date_operation)
  #min(AnyComp$as_date_operation)
  
  ### Use operations before 2018-10-01 as tranning and after 2018-10-01 as testing ###
  #AnyComp$train[AnyComp$as_date_operation <= as.Date("2018-10-01")] <- 1
  #AnyComp$train[AnyComp$as_date_operation > as.Date("2018-10-01")] <- 0
  #table(AnyComp$train)
  #sum(table(AnyComp$train))
  
  #### try a larger test set #####
  AnyComp$train[AnyComp$as_date_operation <= as.Date(train_date_end)] <- 1
  AnyComp$train[AnyComp$as_date_operation > as.Date(train_date_end)] <- 0
  #table(AnyComp$train)
  #sum(table(AnyComp$train))
  
  training_set <- AnyComp[AnyComp$train==1,]
  testing_set <- AnyComp[AnyComp$train==0,]
  #dim(training_set)
  #dim(testing_set)
  
  #colnames(training_set)
  #colnames(training_set)[300]
  #colnames(training_set)[300:400]
  #colnames(training_set)[290]
#  colnames(training_set)[285]
#  dim(training_set)
  #colnames(training_set)[289]
  #colnames(training_set)[ncol(training_set)-1]
  #colnames(training_set)[ncol(training_set)-2]
  
  #dim(training_set)
  #which(colnames(training_set)=="train")
  
  X <- data.matrix(training_set[, c(285:(ncol(training_set)-2))]) 
  #dim(X)
  #column_sums <- c(200,300)
  pb <- txtProgressBar(min = 0, max = length(column_sums), style = 3)
  for (i in 1:length(column_sums)) {
    col_cutoff <- column_sums[i]
    print(paste("Selecting cutoff: # of patients > ",column_sums[i],", Progress percentage: ",round(((i-1)/length(column_sums))*100,2),"%",sep=""))
    
    l <- data.matrix(X[, colSums(X) > col_cutoff]) ### might need to be using percentage
    #l <- data.matrix(X[, colSums(X) > 100]) ### might need to be using percentage 
    
    #dim(l)
    
    set.seed(1234)
    X_k = create.second_order(l, method = c("asdp"), shrink = F)
    #dim(X_k)
    ### Infectious complications (at least one). ASA class included (CPT event-rate excluded). ###
    #numCores <- 16
    #nfolds <- 10
    set.seed(1234)
    W = stat.glmnet_coefdiff(l, X_k, y=training_set$SSIinfBin, nfolds=nfolds, family="binomial",cores = numCores)
    #hist(W,breaks = 100)
    #hist(W)
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, i)
    #fdr_col <- c(0.2,0.25)
    #fdr_col <- c(0.05,0.1,0.2)
    
    for (j in 1:length(fdr_col)) {
      false_discovery_rate <- fdr_col[j]
      print(paste("Selecting fdr ",fdr_col[j],", Progress percentage: ",round(((j-1)/length(fdr_col))*100,2),"%",sep=""))
      set.seed(1234)
      #false_discovery_rate <- 0.2
      thres = knockoff.threshold(W, fdr=false_discovery_rate, offset=offset)
      #thres = knockoff.threshold(W, fdr=0.15, offset=1)
      #thres = knockoff.threshold(W, fdr=0.12, offset=1)
      #thres = knockoff.threshold(W, fdr=0.1, offset=1)
      #thres = knockoff.threshold(W, fdr=0.15, offset=1)
      #thres2 = knockoff.threshold(W, fdr=0.05, offset=1)
      #thres1_5 = knockoff.threshold(W, fdr=0.2, offset=1)
      #thres = knockoff.threshold(W, fdr=0.1, offset=1)
      #?knockoff.threshold()
      
      selected = which(W >= thres)
      if (length(selected)==0)
      {
        warning(paste("A fdr of ",fdr_col[j]," is too low for colSum cutoff ",column_sums[i]," ,thus no predictor selected",sep = ""))
        print(paste("A fdr of ",fdr_col[j]," is too low for colSum cutoff ",column_sums[i]," ,thus no predictor selected",sep = ""))
        next()
      }
      #selected = which(W >= thres3)
      #selected = which(W >= thres1_5)
      
      #print(selected)
      #length(print(selected))
      
      betas <- W[W >= thres]
      
      #cbind(selected, betas)
      n <- as.data.frame(l)
      m <- n[c(selected)]
      
      train2 <- training_set[c("SSIinfBin", names(m))]
      glm1 <- glm(SSIinfBin ~ ., data = train2, family="binomial")
      
      #names(m)
      selected <- c(0,selected)
      
      save_coef <- cbind(names(glm1$coefficients),selected, round(glm1$coefficients,4))
      #rownames(save_coef) <- names(m)
      store_rownames <- rep(NA,nrow(save_coef))
      for (abc in 1:nrow(save_coef)) {
        store_rownames[abc] <- paste("cutoff",column_sums[i],"fdr",fdr_col[j],abc,sep = "_")
        
      }
      rownames(save_coef) <- store_rownames
      #rownames(save_coef) <- paste("cutoff",column_sums[i],"fdr",fdr_col[j],sep = "")
      
      colnames(save_coef) <- c("Predictor_names","col_num_selected","betas")
      
      fit1 <- predict(glm1, training_set, type="response")
      #summary(fit1)
      
      cmFull.1 <- as.data.frame(cbind(training_set$SSIinfBin, fit1))
      colnames(cmFull.1) <- c("SSIinfBin", "PredSSI")
      
      ### Training set fit
      
      roc_lasso1.1 <- roc(cmFull.1$SSIinfBin, cmFull.1$PredSSI)
      #roc_lasso1.1
      youd_lasso1.1 <- coords(roc_lasso1.1, x = "best", best.method = "youden", ret = c("threshold", "specificity", "sensitivity","accuracy", "npv", "ppv", "fn", "fp"))
      #round(youd_lasso1.1,2)
      #plot(roc_lasso1.1)
      
      fit2 <- predict(glm1, testing_set, type="response")
      #summary(fit2)
      
      cmFull.2 <- as.data.frame(cbind(testing_set$SSIinfBin, fit2))
      colnames(cmFull.2) <- c("SSIinfBin", "PredSSI")
      
      ### test set fit
      roc_lasso1.2 <- roc(cmFull.2$SSIinfBin, cmFull.2$PredSSI)
      #roc_lasso1.2
      youd_lasso1.2 <- coords(roc_lasso1.2, x = as.numeric(youd_lasso1.1$threshold), ret = c("threshold", "specificity", "sensitivity","accuracy", "npv", "ppv", "fn", "fp"))
      #round(youd_lasso1.2,2)
      #plot(roc_lasso1.2)
      #ci <- round(ci.coords(roc_lasso1.2, x=.04, ret = c("specificity", "sensitivity","accuracy", "npv", "ppv", "fn", "fp"))*100,2)
      
      #ci <- round(ci.coords(roc_lasso1.2, x=.04, ret = c("specificity", "sensitivity","accuracy", "npv", "ppv", "fn", "fp"))*100,2)
      #class(ci)
      
      #library(knitr)
      #kable(ci)
      
      
      ##### results #####
      #round(youd_lasso1.1,2)
      #round(youd_lasso1.2,2)
      #class(youd_lasso1.1)
      #length(selected)
      #col_cutoff <- 100
      comb_res_train <- round(cbind(col_cutoff,ncol(l),false_discovery_rate,(length(selected)-1),youd_lasso1.1,as.numeric(auc(roc_lasso1.1))),3)
      comb_res_train$training_or_testing <- "training"
      names(comb_res_train)[13] <- "AUC"
      comb_res_test <- round(cbind(col_cutoff,ncol(l),false_discovery_rate,(length(selected)-1),youd_lasso1.2,as.numeric(auc(roc_lasso1.2))),3)
      comb_res_test$training_or_testing <- "testing"
      names(comb_res_test)[13] <- "AUC"
      res1 <- rbind(comb_res_train,comb_res_test)
      colnames(res1) <- c("cutoff: # of patients > "," # of candidate parameters","fdr","# of predictors selected","threshold", 
                          "specificity", "sensitivity","accuracy", "npv", "ppv", "fn", "fp","auc","training_or_testing")
      #rownames(res1) <- c("training","testing")
      #rownames(res1) <- c("training","testing")
      rownames(res1) <- c(paste("training","cutoff",column_sums[i],"fdr",fdr_col[j],"_",sep = "_"),paste("testing","cutoff",column_sums[i],"fdr",fdr_col[j],"_",sep = "_"))
      res1 <- as.data.frame(res1)
      combres[[j]] <- res1
      
      save_coef <- as.data.frame(save_coef)
      #save_coef$cutoff <- t(rep(column_sums[i],nrow(save_coef)))
      save_coef$cutoff <- column_sums[i]
      save_coef$fdr <- fdr_col[j]
      
      store_save_coef_list[[j]] <- save_coef
      #store_save_coef_list[[j]] <- save_coef
      #store_save_coef_list[[j]] <- list(paste("cutoff",column_sums[i],"fdr",fdr_col[j],sep = "")=save_coef)
      #names(store_save_coef_list[[j]]) <- paste("cutoff",column_sums[i],"fdr",fdr_col[j],sep = "")
      #list1 <- list(for (k in 1:length(fdr_col)) { paste("cutoff",column_sums[i],"fdr",fdr_col[k],sep = "")= store_save_coef_list[[j]]
      #  })
      
    }
    
    #list1 <- list(for (k in 1:length(fdr_col)) { paste("cutoff",column_sums[i],"fdr",fdr_col[k],sep = "")= store_save_coef_list[[k]]
    #})
    df1 <- do.call(rbind,store_save_coef_list)
    res_comb_fdr <- do.call(rbind,combres)
    #rownames(res_comb_fdr) <- rep(c("training","testing"),(length(fdr_col)/2))
    #res_comb_fdr
    #res_comb_fdr <- as.data.frame(res_comb_fdr)
    comb_res_fdr_colnum[[i]] <- res_comb_fdr
    store_coef_comb[[i]] <- df1
    
  }
  coef_matrix <- do.call(rbind,store_coef_comb)
  result_matrix <- do.call(rbind,comb_res_fdr_colnum)
  result_matrix$record_index <- rep(1:(nrow(result_matrix)/2), each = 2)
  #result_matrix
  close(pb)
  train_test_info <- c(nrow(training_set),round((nrow(training_set)/nrow(AnyComp))*100,2),nrow(testing_set),round((nrow(testing_set)/nrow(AnyComp))*100,2),
                       round(nrow(training_set)/nrow(testing_set),2))
  names(train_test_info) <- c("# of patients in training set", "% of training set","# of patients in testing set","% of testing set","training/testing ratio")
  store_saved_res_coef_statistics <- list("result"=result_matrix,"coef"=coef_matrix,"train_test_info"=train_test_info)
  
  if (save_res==TRUE){
    if (is.na(save_dir)) {
      warning("When save_res is true, you need a file directory to save results")}
    #  modname <- paste("colSum_cutoff",column_sums[1],"to",column_sums[length(column_sums)],"fdr",fdr_col[1],"to",fdr_col[length(fdr_col)],
    #  ", rho =",rho,", G =",G,", N =",N,", M =",M,", # of perms = ",permnum,", alpha = ",alpha,", type = ",type,", adj_power = ",adjpower)
    modname <- paste(fileName,"train_test_ratio",round(nrow(training_set)/nrow(testing_set),2),"n_of_cutoff",length(column_sums),"from",
                     column_sums[1],"to",column_sums[length(column_sums)],"n_of_fdr",length(fdr_col),"from",fdr_col[1],"to",fdr_col[length(fdr_col)],sep = "_")
    
    if(file.exists(paste(save_dir,modname,"store_saved_res_coef_statistics.Rdata",sep = "_"))==TRUE){
      warning("A previous result has already being saved, press enter to continue")
    }
    #  modname <- paste(modtype,", rho =",rho,", G =",G,", N =",N,", M =",M,", # of perms = ",permnum,", alpha = ",alpha,", type = ",type,", adj_power = ",adjpower)
    #  store_saved_res_coef_statistics <- list("result"=result_matrix,"coef"=coef_matrix,"train_test_info"=train_test_info)
    save(store_saved_res_coef_statistics,file = paste(save_dir,modname,"store_saved_res_coef_statistics.Rdata",sep = "_"))
  }
  invisible(return(store_saved_res_coef_statistics))
}
