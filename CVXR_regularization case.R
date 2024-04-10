#This file uses CVXR to testify the implementation of the New IHW method
library(CVXR)
#load the dataset
gren_objects <- readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UchicagoStudy/Nikolaos Research/gren_objects.Rds")

#function that calculates F_g = b_j^g*t
calculate_F_g <- function(group){
  #get the slopes in group g: F_g = min b_j*t
  return(gren_objects[["grenander_list"]][[group]][["Fs"]][-1])
}

calculate_f_g <- function(group){
  #get the slopes in group g: f_g = min b_j^g
  if(length(gren_objects[["grenander_list"]][[group]][["fs"]]) < 2){
    return(gren_objects[["grenander_list"]][[group]][["fs"]])
  }else{
    return(gren_objects[["grenander_list"]][[group]][["fs"]])
  }
  #return(dist_f_g)  
}

calculate_locs <- function(group){
  return(gren_objects[["grenander_list"]][[group]][["locs"]][-1])
}

#we first test the case when we do not do the regularization part: lambda = 0

alpha <- 0.1 

#G is the number of groups
G <- length(gren_objects[["grenander_list"]])

#rejection threshold for each group, in total 47 groups
ts <- Variable(G)
F_bar <- Variable(G)

#m is the total number of hypothesis testings
m <-sum(gren_objects[["ms"]])
#m_g is the number of hypothesidis in each group
m_g <- gren_objects[["ms"]]

#q_g is the ratio of m_g/m
q_g <- as.numeric(m_g/m)

F_g_list <- sapply(c(1:47), calculate_F_g)
f_g_list <- sapply(c(1:47), calculate_f_g)
locs_list <- sapply(c(1:47), calculate_locs)

#obtain a_j from F_g - b_j^g*t
a_j <- list()

for(i in 1:47){
  temp <- F_g_list[[i]] - locs_list[[i]]*f_g_list[[i]]
  a_j[[i]] <- temp
}

#construct 47 lists of constraints on F_bar
F_bar_constraints <- list()
for(i in 1:47){
  temp <- a_j[[i]] + f_g_list[[i]]*ts[i] 
  F_bar_constraints [[i]] <- F_bar[i] <= temp
}

# Create an empty data frame
results <- matrix(nrow = 47, ncol = 15)

#Now we test the regularization case when lambda != 0.
#iterate through each lambda
for(i in 1: length(gren_objects[["lambdas"]])){
  lambda <- gren_objects[["lambdas"]][i]

  #construct the objective function
  objective <- Maximize( sum(q_g*F_bar) - lambda*sum(abs(diff(ts))))
  constraints <- list(ts >= 0, ts <= 1, sum(q_g*(ts - alpha*F_bar)) <= 0, F_bar >= 0)
  constraints <- c(constraints, F_bar_constraints)
  
  #obtain the solution of the no regularization case
  solution <- solve(Problem(objective, constraints), solver="GUROBI",
                    verbose=TRUE)
  
  solution$status
  solution$value
  results <- cbind(results, solution$getValue(ts))
}

print(results)
results <- results[, c(16:30)]
save(results, file = "regularization.RData")
error <- abs(results - gren_objects[["ts"]])
save(error, file = "error.RData")
#solution$getValue(F_bar)


## This part test the case of lambda = Inf, when all the tg are forced to be the same
#construct the objective function

#Trial 1: Set a very large lambda = 1e100
lambda <- 1e10

#construct the objective function
objective <- Maximize( sum(q_g*F_bar) - lambda*p_norm(diff(x = ts), 1) )
constraints <- list(ts >= 0, ts <= 1, sum(q_g*(ts - alpha*F_bar)) <= 0, F_bar >= 0)
constraints <- c(constraints, F_bar_constraints)

#obtain the solution of the no regularization case
solution <- solve(Problem(objective, constraints), solver="GUROBI",
                  verbose=TRUE)

solution$status
solution$getValue(ts)


#Trial 2: the constraint that all tg are equal
ts_constraints <- list()
for (i in 2:G) {
  ts_constraints[[i-1]] <- ts[i] == ts[i-1]
}


objective <- Maximize( sum(q_g*F_bar))
constraints <- list(ts >= 0, ts <= 1, sum(q_g*(ts - alpha*F_bar)) <= 0, F_bar >= 0)
constraints <- c(constraints, F_bar_constraints, ts_constraints)

#obtain the solution of the no regularization case
solution <- solve(Problem(objective, constraints), solver="GUROBI",
                  verbose=TRUE)

solution$status

solution$getValue(ts)

