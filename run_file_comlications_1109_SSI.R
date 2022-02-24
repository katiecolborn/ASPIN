#source("C:/Users/ASRock/OneDrive - The University of Colorado Denver/Katie Project/2021_0503/bleeding/functions/Generate_statistics_source_function_0504_addCoef_bleeding.R")
source("C:/Users/ASRock/OneDrive - The University of Colorado Denver/Katie Project/2021_1109/SSI/Generate_statistics_source_function_1109_SSI.R")

start8 <- Sys.time()
trial1 <- Generate_statistics_res_SSI(column_sums=c(500,600),fdr_col=c(0.225,0.25))
end8 <- Sys.time()
spendtime8 <- end8-start8
spendtime8

trial1$result
trial1$coef
trial1$train_test_info

##### pre_op SSI #####

#### trial 8 #### change the train test ratio ### Run date:0503

start8 <- Sys.time()
trial8 <- Generate_statistics_res_SSI(column_sums=c(5,10,15,20,25,30,35,40,45,50),fdr_col=c(0.1,0.125,0.15,0.175,0.2,0.225,0.25))
end8 <- Sys.time()
spendtime8 <- end8-start8
spendtime8

trial8$result
trial8$coef
trial8$train_test_info

write.csv(trial8$result,file = "C:/Users/ASRock/OneDrive - The University of Colorado Denver/Katie Project/2021_1109/SSI/res/SSI_PheCode.csv")