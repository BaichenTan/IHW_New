#This files tests the optimal_ts function
source("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UchicagoStudy/Nikolaos Research/IHW_New/functions_new_methods.R")

#load the dataset
gren_objects <- readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UchicagoStudy/Nikolaos Research/gren_objects.Rds")

#load parameters
train_gre <- gren_objects[["grenander_list"]]
lambdas <- c(Inf, 1.0 / 1.2^(1:7), 1/ 2^(2:13), 0)
m_groups <- gren_objects[["ms"]]

#run optimal_ts
#time1 <- system.time(results_1 <- IHW:::optimal_ts(train_gre, m_groups, 0.1, lambdas))[1]
time1 <- system.time(results_1 <- IHW:::optimal_ts(train_gre, m_groups, 0.1, c(1.0 / 1.2^(1:7), 1/ 2^(2:13))))[1]
time2 <- system.time(results_2 <- optimal_ts_new(train_gre, m_groups, 0.1, c(1.0 / 1.2^(1:7), 1/ 2^(2:13))))[1]
#source("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UchicagoStudy/Nikolaos Research/New IHW/R/ihw_convex.R")
#time3 <- system.time(results_3 <- optimal_cvxr(train_gre, m_groups, 0.1, c(1.0 / 1.2^(1:7), 1/ 2^(2:13))))[1]

print(abs(results_1 - results_2))
print(c(time1, time2))