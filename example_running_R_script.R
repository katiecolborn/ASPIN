source("C:/Users/ASRock/OneDrive - The University of Colorado Denver/Katie Project/2021_0425/SSI/functions/Generate_statistics_source_function_0429_addCoef_SSI.R")
source("C:/Users/ASRock/OneDrive - The University of Colorado Denver/Katie Project/2021_0425/UTI/functions/Generate_statistics_source_function_0429_addCoef_UTI.R")
source("C:/Users/ASRock/OneDrive - The University of Colorado Denver/Katie Project/2021_0425/SYSEP/functions/Generate_statistics_source_function_0429_addCoef_SYSEP.R")
source("C:/Users/ASRock/OneDrive - The University of Colorado Denver/Katie Project/2021_0425/Pneumonia/functions/Generate_statistics_source_function_0429_addCoef_Pneumonia.R")

start8 <- Sys.time()
trial1 <- Generate_statistics_res_SSI(column_sums=c(500,600),fdr_col=c(0.225,0.25))
end8 <- Sys.time()
spendtime8 <- end8-start8
spendtime8

trial1$result
trial1$coef
trial1$train_test_info

start8 <- Sys.time()
trial2 <- Generate_statistics_res_UTI(column_sums=c(500,600),fdr_col=c(0.225,0.25))
end8 <- Sys.time()
spendtime8 <- end8-start8
spendtime8

trial2$result
trial2$coef
trial2$train_test_info

start8 <- Sys.time()
trial3 <- Generate_statistics_res_SYSEP(column_sums=c(500,600),fdr_col=c(0.225,0.25))
end8 <- Sys.time()
spendtime8 <- end8-start8
spendtime8

trial3$result
trial3$coef
trial3$train_test_info

start8 <- Sys.time()
trial4 <- Generate_statistics_res_Pneu(column_sums=c(500,600),fdr_col=c(0.225,0.25))
end8 <- Sys.time()
spendtime8 <- end8-start8
spendtime8

trial4$result
trial4$coef
trial4$train_test_info


##### SSI #####

#### trial 8 #### change the train test ratio ### Run date:0429

start8 <- Sys.time()
trial8 <- Generate_statistics_res_SSI(column_sums=c(5,10,15,20,25,30,35,40,45,50),fdr_col=c(0.1,0.125,0.15,0.175,0.2,0.225,0.25))
end8 <- Sys.time()
spendtime8 <- end8-start8
spendtime8

trial8$result
trial8$coef
trial8$train_test_info


##### UTI #####

#### trial 8 #### change the train test ratio

start8 <- Sys.time()
trial8 <- Generate_statistics_res_UTI(column_sums=c(5,10,15,20,25,30,35,40,45,50),fdr_col=c(0.1,0.125,0.15,0.175,0.2,0.225,0.25))
end8 <- Sys.time()
spendtime8 <- end8-start8
spendtime8

trial8$result
trial8$coef
trial8$train_test_info

##### SYSEP #####

#### trial 8 #### change the train test ratio

start8 <- Sys.time()
trial8 <- Generate_statistics_res_SYSEP(column_sums=c(5,10,15,20,25,30,35,40,45,50),fdr_col=c(0.1,0.125,0.15,0.175,0.2,0.225,0.25))
end8 <- Sys.time()
spendtime8 <- end8-start8
spendtime8

trial8$result
trial8$coef
trial8$train_test_info

##### PNEU #####

#### trial 8 #### change the train test ratio

start8 <- Sys.time()
trial8 <- Generate_statistics_res_Pneu(column_sums=c(5,10,15,20,25,30,35,40,45,50),fdr_col=c(0.1,0.125,0.15,0.175,0.2,0.225,0.25))
end8 <- Sys.time()
spendtime8 <- end8-start8
spendtime8

trial8$result
trial8$coef
trial8$train_test_info

