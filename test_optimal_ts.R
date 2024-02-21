#This files tests the optimal_ts function

#load the dataset
gren_objects <- readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UchicagoStudy/Nikolaos Research/gren_objects.Rds")

#load parameters
train_gre <- gren_objects[["grenander_list"]]
lambdas <- c(1.0 / 2 ^ (1:13))
m_groups <- gren_objects[["ms"]]

#run optimal_ts
time1 <- system.time(results_1 <- IHW:::optimal_ts(train_gre, m_groups, 0.1, lambdas))[1]
time1 <- system.time(results_1 <- IHW:::optimal_ts(train_gre, m_groups, 0.1, lambdas))[1]
#source("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UchicagoStudy/Nikolaos Research/New IHW/R/ihw_convex.R")
time2 <- system.time(results_2 <- optimal_cvxr(train_gre, m_groups, 0.1, lambdas))[1]

print(c(time1, time2))