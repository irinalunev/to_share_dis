##### LUNEVA, 2022 #####

# Packages
library(parallel)
library(doParallel)
library(foreach)
library(compiler)
library(future)
library(dplyr)
library(data.table)

# Cluster for parallel computing
#cl_fullcode <- parallel::makeCluster(detectCores()-2)
#plan(cluster, workers=cl_fullcode)
n.cores <- availableCores()-1
cl1 <- makeCluster(n.cores)
#register it to be used by %dopar%
#doParallel::registerDoParallel(cl = cl1)
#check if it is registered (optional)
#foreach::getDoParRegistered()
#require(compiler)
#enableJIT(3)

# General functions
norm_pdf <- function(x, mu, sigma){
  return(dnorm(x, mean = mu, sd = sigma))
}

# Hyperparametes
#mus <- 0
#mud <- 0
sigmas <- 0.25 #hard info
sigmaepsilons <- 0.05 #noise in hard info
sigmad <- 0.25 #soft info
sigmaepsilond <- 0.001 #noise in soft info
mv <- 0.001
mnv <- 0.003
v <- 0.1
k <- 0.5
L1 <- 0
L2 <- -0.1
I <- -0.2
x_grid <- seq(-3, 3, 0.3)
x_double_grid <- as.matrix(expand.grid(x_grid, x_grid))
names(x_double_grid) <- c("z", "K")
x_grid_small <- seq(-1,1,0.3)
x_double_grid_small <- as.matrix(expand.grid(x_grid_small, x_grid_small))
names(x_double_grid_small) <- c("z", "K")
s_grid <- seq(-1,1,0.3)
d_grid <- seq(-1,1,0.3)
sd_grid <- as.matrix(expand.grid(s_grid, d_grid))
names(sd_grid) <- c("s", "d")

# Vectors of actions
mons <- c(NA,NA)
viols <- c(NA,NA)
liquids <- c(NA, NA)

# Vectors of informations
ss <- c(s_grid[4],s_grid[3])
ds <- c(d_grid[5],d_grid[5])


### Pre-terminal period, t=5 ####
solve_period_2 <- function(mon1, mon2, viol1,viol2, liquids, z1, K1, s1,s2, d1,d2)
{
  ss <- c(s1,s2)
  ds <- c(d1,d2)
  mons <- c(mon1,mon2)
  viols <- c(viol1, viol2)
  liquids <- liquids
  z1 <- z1
  K1 <- K1
  
  
  # Did violation happen?
  ifelse(ss[2]>=z1, viols[2] <- 0, viols[2] <- 1)
  
  if(viols[2]==0) {
    # Monitor or not? 
    # In the terminal period, since the bank can not obtain higher than K4 but monitoring is costly, the bank NEVER monitors
    mons[2] <- 0
    K2 <- K1
    liquids[2] <- 0
    Ex2_nv_nm <- (sum(ss[1:2])*sigmas^2/(sigmas^2+sigmaepsilons^2)
                  +as.numeric(sum(mons, na.rm = T)>0)*sum(ds, na.rm = T)*sigmad^2/(sigmad^2+sigmaepsilond^2))
    Varx2_nv_nm <- (3*sigmas^2-2*sigmas^4/(sigmas^2+sigmaepsilons^2)
                    +3*sigmad^2-(sum(mons, na.rm = T))*sigmad^4/(sigmad^2+sigmaepsilond^2))
  }
  
  if (viols[2]==1) {
    # Violation, monitoring 
    Ex2_v_m <- (sum(ss[1:2])*sigmas^2/(sigmas^2+sigmaepsilons^2)
                +sum(ds[1:2])*sigmad^2/(sigmad^2+sigmaepsilond^2))
    Varx2_v_m <- (3*sigmas^2-2*sigmas^4/(sigmas^2+sigmaepsilons^2)
                  +3*sigmad^2-(sum(mons, na.rm = T)+1)*sigmad^4/(sigmad^2+sigmaepsilond^2))
    
    find_K2_v_m <- function(viols, Ex2_v_m, Varx2_v_m)
    {
      # Do we liquidate or not?
      ifelse(L2>Ex2_v_m-sum(viols)*v, return(c("liquidate", L2)),
             {
               func_for_K2_best_v_m <- function(K2_candidate)
               {
                 toint_equation1 <- function(x) {(x-K2_candidate-sum(viols)*v)*norm_pdf(x, Ex2_v_m, sqrt(Varx2_v_m))}
                 equation <- abs(integrate(toint_equation1, K2_candidate+sum(viols)*v, Inf)$value - k*(Ex2_v_m-sum(viols)*v-L2))
                 return(equation)
               }
               func_for_K2_best_v_m_vect <- Vectorize(func_for_K2_best_v_m, vectorize.args = "K2_candidate")
               equations_v_m <- func_for_K2_best_v_m_vect(x_grid_small)
               
               K2_best_v_m <- x_grid_small[which(equations_v_m==min(equations_v_m))]
               
               toint_bankspayoff1_v_m <- function(x) {(x-sum(viols)*v)*norm_pdf(x, Ex2_v_m, sqrt(Varx2_v_m))}
               toint_bankspayoff2_v_m <- function(x) {K2_best_v_m*norm_pdf(x, Ex2_v_m, sqrt(Varx2_v_m))}
               
               banks_payoff_v_m <- integrate(toint_bankspayoff1_v_m, -Inf, K2_best_v_m+sum(viols)*v)$value + integrate(toint_bankspayoff2_v_m, K2_best_v_m+sum(viols)*v, Inf)$value
               return(c(K2_best_v_m, banks_payoff_v_m))
             })
    }
    
    # Violation, no monitoring
    Ex2_v_nm <- (sum(ss[1:2])*sigmas^2/(sigmas^2+sigmaepsilons^2)
                 +as.numeric(sum(mons, na.rm = T)>0)*sum(ds,na.rm=T)*sigmad^2/(sigmad^2+sigmaepsilond^2))
    Varx2_v_nm <- (3*sigmas^2-2*sigmas^4/(sigmas^2+sigmaepsilons^2)
                   +3*sigmad^2-(sum(mons, na.rm = T))*sigmad^4/(sigmad^2+sigmaepsilond^2))
    
    find_K2_v_nm <- function(viols, Ex2_v_nm, Varx2_v_nm)
    {
      # Do we liquidate or not?
      ifelse(L2>Ex2_v_nm-sum(viols)*v, return(c("liquidate", L2)),
             {
               func_for_K2_best_v_nm <- function(K2_candidate)
               {
                 toint_equation1_v_nm <- function(x) {(x-K2_candidate-sum(viols)*v)*norm_pdf(x, Ex2_v_nm, sqrt(Varx2_v_nm))}
                 equation <- abs(integrate(toint_equation1_v_nm, K2_candidate+sum(viols)*v, Inf)$value - k*(Ex2_v_nm-sum(viols)*v-L2))
                 return(equation)
               }
               func_for_K2_best_v_nm_vect <- Vectorize(func_for_K2_best_v_nm, vectorize.args = "K2_candidate")
               equations <- func_for_K2_best_v_nm_vect(x_grid_small)
               K2_best_v_nm <- x_grid_small[which(equations==min(equations))]
               
               toint_bankspayoff1_v_nm <- function(x) {(x-sum(viols)*v)*norm_pdf(x, Ex2_v_nm, sqrt(Varx2_v_nm))}
               toint_bankspayoff2_v_nm <- function(x) {K2_best_v_nm*norm_pdf(x, Ex2_v_nm, sqrt(Varx2_v_nm))}
               
               banks_payoff_v_nm <- integrate(toint_bankspayoff1_v_nm, -Inf, K2_best_v_nm+sum(viols)*v)$value + integrate(toint_bankspayoff2_v_nm, K2_best_v_nm+sum(viols)*v, Inf)$value
               
               return(c(K2_best_v_nm, banks_payoff_v_nm))
             })
    }
    
    # Monitor or not? 
    # Expected payoff if monitor
    expected_payoff_ifmonitor_2_v <- function(viols)
    {
      
      forapply1_v <- function(d2)
      { 
        Ex2_v_m <- (sum(ss[1:2])*sigmas^2/(sigmas^2+sigmaepsilons^2)
                    +sum(c(ds[1],d2))*sigmad^2/(sigmad^2+sigmaepsilond^2))
        Varx2_v_m <- (3*sigmas^2-2*sigmas^4/(sigmas^2+sigmaepsilons^2)
                      +3*sigmad^2-(sum(mons, na.rm = T)+1)*sigmad^4/(sigmad^2+sigmaepsilond^2))
        
        return(as.numeric(find_K2_v_m(viols, Ex2_v_m, Varx2_v_m)[2])*norm_pdf(d2, 0, sqrt(sigmad^2+sigmaepsilond^2))*(x_grid[2]-x_grid[1]))
      }
      forapply1_v_vect <- Vectorize(forapply1_v, vectorize.args = "d2")
      payoffs <- forapply1_v_vect(d_grid)
      
      exp_payoff_v <- sum(payoffs)
      
      return(exp_payoff_v)
    }
    
    ifmonitor_v <- expected_payoff_ifmonitor_2_v(viols)-mv
    ifnotmonitor_v <- find_K2_v_nm(viols, Ex2_v_nm, Varx2_v_nm)[2]
    
    ifelse(ifmonitor_v>=ifnotmonitor_v, mons[2] <- 1, mons[2] <- 0)
    ifelse(ifmonitor_v>=ifnotmonitor_v, K2 <- find_K2_v_m(viols, Ex2_v_m, Varx2_v_m)[1], K2 <- find_K2_v_nm(viols, Ex2_v_nm, Varx2_v_nm)[1])
    ifelse(K2!="liquidate", liquids[2] <- 0, liquids[2] <- 1)
  }
  
  #manager's and bank's payoffs
  if (viols[2]==0) 
  { toint_for_man_payoff_nv <- function(x) {(x-sum(viols)*v-K2)*norm_pdf(x, Ex2_nv_nm, sqrt(Varx2_nv_nm))}
  man_payoff <- integrate(toint_for_man_payoff_nv, K2+sum(viols)*v, Inf)$value
  
  toint1_for_bank_payoff_nv <- function(x) {(x-sum(viols)*v)*norm_pdf(x, Ex2_nv_nm, sqrt(Varx2_nv_nm))}
  toint2_for_bank_payoff_nv <- function(x) {K2*norm_pdf(x, Ex2_nv_nm, sqrt(Varx2_nv_nm))}
  bank_payoff <- integrate(toint1_for_bank_payoff_nv, -Inf, K2+sum(viols)*v)$value+integrate(toint2_for_bank_payoff_nv, K2+sum(viols)*v, Inf)$value}
  
  if (viols[2]==1)
  {
    if(liquids[2]==1)
    {
      man_payoff <- 0
      bank_payoff <- L2-mv*as.numeric(mons[2]==1)
    }
    if (liquids[2]==0)
    {
      if(mons[2]==0){
        toint_for_man_payoff_v <- function(x) {(x-sum(viols)*v-K2)*norm_pdf(x, Ex2_v_nm, sqrt(Varx2_v_nm))}
        man_payoff <- integrate(toint_for_man_payoff_v, K2+sum(viols)*v, Inf)$value
        
        toint1_for_bank_payoff_v <- function(x) {(x-sum(viols)*v)*norm_pdf(x, Ex2_v_nm, sqrt(Varx2_v_nm))}
        toint2_for_bank_payoff_v <- function(x) {K2*norm_pdf(x, Ex2_v_nm, sqrt(Varx2_v_nm))}
        bank_payoff <- integrate(toint1_for_bank_payoff_v, -Inf, K2+sum(viols)*v)$value+integrate(toint2_for_bank_payoff_v, K2+sum(viols)*v, Inf)$value
      }
      if(mons[2]==1){
        toint_for_man_payoff_v <- function(x) {(x-sum(viols)*v-K2)*norm_pdf(x, Ex2_v_m, sqrt(Varx2_v_m))}
        man_payoff <- integrate(toint_for_man_payoff_v, K2+sum(viols)*v, Inf)$value
        
        toint1_for_bank_payoff_v <- function(x) {(x-sum(viols)*v)*norm_pdf(x, Ex2_v_m, sqrt(Varx2_v_m))}
        toint2_for_bank_payoff_v <- function(x) {K2*norm_pdf(x, Ex2_v_m, sqrt(Varx2_v_m))}
        bank_payoff <- integrate(toint1_for_bank_payoff_v, -Inf, K2+sum(viols)*v)$value+integrate(toint2_for_bank_payoff_v, K2+sum(viols)*v, Inf)$value-mv
      }
    }
  }
  
  return(list(viols=viols, mons=mons, liquids=liquids, K2=K2, man_payoff=man_payoff, bank_payoff=bank_payoff))
}
solve_period_2_comp <- cmpfun(solve_period_2)

solve_period_2_new <- function(x)
{ z1 <- x[1]
K1 <- x[2]
d2 <- x[3]
s2 <- x[4]
d1 <- x[5]
mon1 <- x[6]
viol1 <- x[7]
solve_period_2_comp(mon1,NA, viol1, NA, liquids=liquids,z1=z1,K1=K1,s1=ss[1],s2=s2,d1=d1,d2=d2)}
solve_period_2_new_comp <- cmpfun(solve_period_2_new)
total_grid_forper2 <- expand.grid(x_grid_small, x_grid_small, d_grid, s_grid, d_grid, c(1,0), c(1,0))
names(total_grid_forper2) <- c("z1","K1","d2","s2","d1","mon1","viol1")
#f_per5_mon40_new <- apply(total_grid_forper5_new, 1, solve_period_5_new)
clusterExport(cl1, list("x_double_grid_small","d_grid","ds","I","k","L1","L2","liquids","mons","mnv","mv","s_grid","sigmad","sigmaepsilond","sigmas","sigmaepsilons","ss","v","viols","x_grid","x_grid_small","norm_pdf","solve_period_2","solve_period_2_new","solve_period_2_comp"))
f_per2 <- parApply(cl1, total_grid_forper2, 1, solve_period_2_new_comp)
f_per2 <- do.call("cbind",f_per2)
solution_period_2 <- cbind(total_grid_forper2, t(f_per2))
solution_period_2_dt <- as.data.table(solution_period_2)
solution_period_2_dt_v_m <- solution_period_2_dt[solution_period_2_dt$mon1==1 & solution_period_2_dt$viol1==1,]
solution_period_2_dt_v_nm <- solution_period_2_dt[solution_period_2_dt$mon1==0 & solution_period_2_dt$viol1==1,]
solution_period_2_dt_nv_m <- solution_period_2_dt[solution_period_2_dt$mon1==1 & solution_period_2_dt$viol1==0,]
solution_period_2_dt_nv_nm <- solution_period_2_dt[solution_period_2_dt$mon1==0 & solution_period_2_dt$viol1==0,]
#stopCluster(cl1)


#f_per5_mon40 <- future(solve_period_5_vect(c(mons[1:3],0, NA),viols,liquids,total_grid_forper5[,1],total_grid_forper5[,2],ss[1],ss[2],ss[3],ss[4],total_grid_forper5[,4],ds[1],ds[2],ds[3],total_grid_forper5[,5],total_grid_forper5[,3]))
#value(f_per5_mon40)
#solution_period_5_mon40 <- cbind(total_grid_forper5, t(value(f_per5_mon40)))
#f_per5_mon41 <- future(solve_period_5_vect(c(mons[1:3],1, NA),viols,liquids,total_grid_forper5[,1],total_grid_forper5[,2],ss[1],ss[2],ss[3],ss[4],total_grid_forper5[,4],ds[1],ds[2],ds[3],total_grid_forper5[,5],total_grid_forper5[,3]))
#value(f_per5_mon41)
#solution_period_5_mon41 <- cbind(total_grid_forper5, t(value(f_per5_mon41)))



### Regular period, t=1 ####
solve_period_1 <- function(mons, viols, liquids, z0, K0, s1,s2,d1,d2)
{
  ss <- c(s1,s2)
  ds <- c(d1,d2)
  
  # Did violation happen?
  ifelse(ss[1]>=z0, viols[1] <- 0, viols[1] <- 1)
  
  if(viols[1]==0)
  {
    ###liquids[4] <- 0
    # No violation, monitoring
    find_bestz1_K1_nv_mon <- function(for_d1)
    {
      solution_period_2_ford1_nv_m <- solution_period_2_dt_nv_m[solution_period_2_dt_nv_m$d1==for_d1,]
      func_payoffs_nv_m <- function(for_z1, for_K1)
      {
        solution_period_2_small_nv_m <- solution_period_2_ford1_nv_m[solution_period_2_ford1_nv_m$z1==for_z1 & solution_period_2_ford1_nv_m$K1==for_K1,]
        manager_total_payoff <- (sum(unlist(solution_period_2_small_nv_m$man_payoff)
                                     *(s_grid[2]-s_grid[1])*norm_pdf(unlist(solution_period_2_small_nv_m$s2),0,sqrt(sigmas^2+sigmaepsilons^2))
                                     *(d_grid[2]-d_grid[1])*norm_pdf(unlist(solution_period_2_small_nv_m$d2),0,sqrt(sigmad^2+sigmaepsilond^2))))
        bank_total_payoff <- (sum(unlist(solution_period_2_small_nv_m$bank_payoff)
                                  *(s_grid[2]-s_grid[1])*norm_pdf(unlist(solution_period_2_small_nv_m$s2),0,sqrt(sigmas^2+sigmaepsilons^2))
                                  *(d_grid[2]-d_grid[1])*norm_pdf(unlist(solution_period_2_small_nv_m$d2),0,sqrt(sigmad^2+sigmaepsilond^2))))
        
        solution_period_2_small_alt_nv_m <- solution_period_2_ford1_nv_m[solution_period_2_ford1_nv_m$z1==z0 & solution_period_2_ford1_nv_m$K1==K0,]
        bank_total_payoff_alt <- (sum(unlist(solution_period_2_small_alt_nv_m$bank_payoff)
                                      *(s_grid[2]-s_grid[1])*norm_pdf(unlist(solution_period_2_small_alt_nv_m$s2),0,sqrt(sigmas^2+sigmaepsilons^2))
                                      *(d_grid[2]-d_grid[1])*norm_pdf(unlist(solution_period_2_small_alt_nv_m$d2),0,sqrt(sigmad^2+sigmaepsilond^2))))
        
        return(c(manager_total_payoff, bank_total_payoff, bank_total_payoff_alt))
      }
      func_payoffs_nv_m_vect <- Vectorize(func_payoffs_nv_m, vectorize.args = c("for_z1", "for_K1"))
      payoffs_nv_m <- t(func_payoffs_nv_m_vect(x_double_grid_small[,1],x_double_grid_small[,2]))
      payoffs_nv_m <- data.frame(manager_total_payoff=payoffs_nv_m[,1],bank_total_payoff=payoffs_nv_m[,2],bank_total_payoff_alt=payoffs_nv_m[,3])
      
      payoffs_nv_m$to_choose <- as.numeric(payoffs_nv_m$bank_total_payoff>=payoffs_nv_m$bank_total_payoff_alt)
      payoffs_nv_m$to_choose <- payoffs_nv_m$to_choose*payoffs_nv_m$manager_total_payoff
      
      
      z1_best_nv_m <- x_double_grid_small[which(payoffs_nv_m$to_choose==max(payoffs_nv_m$to_choose)),1]
      K1_best_nv_m <- x_double_grid_small[which(payoffs_nv_m$to_choose==max(payoffs_nv_m$to_choose)),2]
      bank_payoff_for_bestz1K1_nv_m <- payoffs_nv_m$bank_total_payoff[which(payoffs_nv_m$to_choose==max(payoffs_nv_m$to_choose))]
      manager_payoff_for_bestz1K1_nv_m <- payoffs_nv_m$manager_total_payoff[which(payoffs_nv_m$to_choose==max(payoffs_nv_m$to_choose))]
      
      #bank_payoff <- payoffs[which(payoffs[,1]==max(payoffs[,1])),3]
      #bank_alt_payoff <- payoffs[which(payoffs[,1]==max(payoffs[,1])),4]
      #man_payoff <- payoffs[which(payoffs[,1]==max(payoffs[,1])),2]
      return(c(z1_best_nv_m, K1_best_nv_m, bank_payoff_for_bestz1K1_nv_m, manager_payoff_for_bestz1K1_nv_m))
    }
    #find_bestz4_K4_mon_comp <- cmpfun(find_bestz4_K4_mon)
    find_bestz1_K1_nv_mon_vect <- Vectorize(find_bestz1_K1_nv_mon, vectorize.args = "for_d1")
    #find_bestz4_K4_mon_vect_comp <- cmpfun(find_bestz4_K4_mon_vect)
    
    # No violation, no monitoring
    find_bestz1_K1_nv_nomon <- function()
    {
      solution_period_2_mini_nv_nm <- solution_period_2_dt_nv_nm[solution_period_2_dt_nv_nm$d1==d_grid[1],]
      func_payoffs_nv_nm <- function(for_z1, for_K1)
      {
        solution_period_2_small_nv_nm <- solution_period_2_mini_nv_nm[solution_period_2_mini_nv_nm$z1==for_z1 & solution_period_2_mini_nv_nm$K1==for_K1,]
        manager_total_payoff <- (sum(unlist(solution_period_2_small_nv_nm$man_payoff)
                                     *(s_grid[2]-s_grid[1])*norm_pdf(unlist(solution_period_2_small_nv_nm$s2),0,sqrt(sigmas^2+sigmaepsilons^2))
                                     *(d_grid[2]-d_grid[1])*norm_pdf(unlist(solution_period_2_small_nv_nm$d2),0,sqrt(sigmad^2+sigmaepsilond^2))))
        bank_total_payoff <- (sum(unlist(solution_period_2_small_nv_nm$bank_payoff)
                                  *(s_grid[2]-s_grid[1])*norm_pdf(unlist(solution_period_2_small_nv_nm$s2),0,sqrt(sigmas^2+sigmaepsilons^2))
                                  *(d_grid[2]-d_grid[1])*norm_pdf(unlist(solution_period_2_small_nv_nm$d2),0,sqrt(sigmad^2+sigmaepsilond^2))))
        
        solution_period_2_small_alt_nv_nm <- solution_period_2_mini_nv_nm[solution_period_2_mini_nv_nm$z1==z0 & solution_period_2_mini_nv_nm$K1==K0,]
        bank_total_payoff_alt <- (sum(unlist(solution_period_2_small_alt_nv_nm$bank_payoff)
                                      *(s_grid[2]-s_grid[1])*norm_pdf(unlist(solution_period_2_small_alt_nv_nm$s2),0,sqrt(sigmas^2+sigmaepsilons^2))
                                      *(d_grid[2]-d_grid[1])*norm_pdf(unlist(solution_period_2_small_alt_nv_nm$d2),0,sqrt(sigmad^2+sigmaepsilond^2))))
        
        return(c(manager_total_payoff, bank_total_payoff, bank_total_payoff_alt))
      }
      func_payoffs_nv_nm_vect <- Vectorize(func_payoffs_nv_nm, vectorize.args = c("for_z1", "for_K1"))
      payoffs_nv_nm <- t(func_payoffs_nv_nm_vect(x_double_grid_small[,1],x_double_grid_small[,2]))
      payoffs_nv_nm <- data.frame(manager_total_payoff=payoffs_nv_nm[,1],bank_total_payoff=payoffs_nv_nm[,2],bank_total_payoff_alt=payoffs_nv_nm[,3])
      
      payoffs_nv_nm$to_choose <- as.numeric(payoffs_nv_nm$bank_total_payoff>=payoffs_nv_nm$bank_total_payoff_alt)
      payoffs_nv_nm$to_choose <- payoffs_nv_nm$to_choose*payoffs_nv_nm$manager_total_payoff
      
      
      z1_best_nv_nm <- x_double_grid_small[which(payoffs_nv_nm$to_choose==max(payoffs_nv_nm$to_choose)),1]
      K1_best_nv_nm <- x_double_grid_small[which(payoffs_nv_nm$to_choose==max(payoffs_nv_nm$to_choose)),2]
      bank_payoff_for_bestz1K1_nv_nm <- payoffs_nv_nm$bank_total_payoff[which(payoffs_nv_nm$to_choose==max(payoffs_nv_nm$to_choose))]
      manager_payoff_for_bestz1K1_nv_nm <- payoffs_nv_nm$manager_total_payoff[which(payoffs_nv_nm$to_choose==max(payoffs_nv_nm$to_choose))]
      
      #bank_payoff <- payoffs[which(payoffs[,1]==max(payoffs[,1])),3]
      #bank_alt_payoff <- payoffs[which(payoffs[,1]==max(payoffs[,1])),4]
      #man_payoff <- payoffs[which(payoffs[,1]==max(payoffs[,1])),2]
      return(c(z1_best_nv_nm, K1_best_nv_nm, bank_payoff_for_bestz1K1_nv_nm, manager_payoff_for_bestz1K1_nv_nm))
    }
    #find_bestz4_K4_nomon_comp <- cmpfun(find_bestz4_K4_nomon)
    
    # Monitor or not?
    # Expected payoff if monitor
    for_expected_payoffs_nv <- find_bestz1_K1_nv_mon_vect(d_grid)
    bank_expected_payoff_ifmon_nv <- sum(for_expected_payoffs_nv[3,]*norm_pdf(d_grid,0,sqrt(sigmad^2+sigmaepsilond^2))*(d_grid[2]-d_grid[1]))
    sol_ifmon_nv <- find_bestz1_K1_nv_mon(for_d1=ds[1])
    sol_ifnomon_nv <- find_bestz1_K1_nv_nomon()
    ifelse(bank_expected_payoff_ifmon_nv>=sol_ifnomon_nv[3]+mnv, mons[1] <- 1, mons[1] <- 0)
    ifelse(bank_expected_payoff_ifmon_nv>=sol_ifnomon_nv[3]+mnv, z1 <- sol_ifmon_nv[1], z1 <- sol_ifnomon_nv[1])
    ifelse(bank_expected_payoff_ifmon_nv>=sol_ifnomon_nv[3]+mnv, K1 <- sol_ifmon_nv[2], K1 <- sol_ifnomon_nv[2])
    ifelse(bank_expected_payoff_ifmon_nv>=sol_ifnomon_nv[3]+mnv, bank_payoff <- sol_ifmon_nv[3]-mnv, bank_payoff <- sol_ifnomon_nv[3])
    ifelse(bank_expected_payoff_ifmon_nv>=sol_ifnomon_nv[3]+mnv, man_payoff <- sol_ifmon_nv[4], man_payoff <- sol_ifnomon_nv[4])
  }
  
  if(viols[1]==1)
  {
    # Violation, monitoring
    solve_viol_mon_per1 <- function(for_d1)
    {
      ds <- c(for_d1,ds[2])
      Ex1_v_m <- (s1*sigmas^2/(sigmas^2+sigmaepsilons^2)
                  +d1*sigmad^2/(sigmad^2+sigmaepsilond^2))
      Varx1_v_m <- (3*sigmas^2-sigmas^4/(sigmas^2+sigmaepsilons^2)
                    +3*sigmad^2-sigmad^4/(sigmad^2+sigmaepsilond^2))
      
      # Do we liquidate or not?
      if(L1>Ex1_v_m-viols[1]*v)
      {
        ###liquids[4] <- 1
        manager_payoff1_viol_mon <- 0
        bank_payoff1_viol_mon <- L1
        z1_best <- "liquidate"
        K1_best <- "liquidate"
      }
      else
      {
        ###liquids[4] <- 0
        # Find optimal z4
        solution_period_2_ford1_v_m <- solution_period_2_dt_v_m[solution_period_2_dt_v_m$d1==for_d1,]
        to_find_best_z1_viol_mon <- function(for_z1)
        {
          sol_per2_small_v_m <- solution_period_2_ford1_v_m[solution_period_2_ford1_v_m$z1==for_z1 & solution_period_2_ford1_v_m$K1==x_grid_small[2],]
          social_payoff_v_m <- sum((unlist(sol_per2_small_v_m$bank_payoff)+unlist(sol_per2_small_v_m$man_payoff))
                               *norm_pdf(unlist(sol_per2_small_v_m$d2),0,sqrt(sigmad^2+sigmaepsilond^2))*norm_pdf(unlist(sol_per2_small_v_m$s2),0,sqrt(sigmas^2+sigmaepsilons^2))
                               *(d_grid[2]-d_grid[1])*(s_grid[2]-s_grid[1]))
          return(social_payoff_v_m)
        }
        to_find_best_z1_viol_mon_vect <- Vectorize(to_find_best_z1_viol_mon, vectorize.args = "for_z1")
        social_payoffs_v_m <- to_find_best_z1_viol_mon_vect(x_grid_small)
        z1_best <- x_grid_small[which(social_payoffs_v_m==max(social_payoffs_v_m))]
        
        # Find optimal K4
        solution_period_2_ford1_withoptz1_v_m <- solution_period_2_ford1_v_m[solution_period_2_ford1_v_m$z1==z1_best,]
        to_find_best_K1_viol_mon <- function(for_K1)
        {
          sol_per2_small_v_m <- solution_period_2_ford1_withoptz1_v_m[solution_period_2_ford1_withoptz1_v_m$K1==for_K1,]
          man_payoff <- sum(unlist(sol_per2_small_v_m$man_payoff)
                            *norm_pdf(unlist(sol_per2_small_v_m$s2),0,sqrt(sigmas^2+sigmaepsilons^2))*norm_pdf(unlist(sol_per2_small_v_m$d2),0,sqrt(sigmad^2+sigmaepsilond^2))
                            *(s_grid[2]-s_grid[1])*(d_grid[2]-d_grid[1]))
          equation <- abs(man_payoff-k*(Ex1_v_m-viols[1]*v-L1))
          bank_payoff <- sum(unlist(sol_per2_small_v_m$bank_payoff)
                             *norm_pdf(unlist(sol_per2_small_v_m$s2),0,sqrt(sigmas^2+sigmaepsilons^2))*norm_pdf(unlist(sol_per2_small_v_m$d2),0,sqrt(sigmad^2+sigmaepsilond^2))
                             *(s_grid[2]-s_grid[1])*(d_grid[2]-d_grid[1]))
          return(c(equation, man_payoff, bank_payoff))
        }
        to_find_best_K1_viol_mon_vect <- Vectorize(to_find_best_K1_viol_mon, vectorize.args = "for_K1")
        equations_v_m <- t(to_find_best_K1_viol_mon_vect(x_grid_small))
        K1_best <- x_grid_small[which(equations_v_m[,1]==min(equations_v_m[,1]))]
        manager_payoff1_viol_mon <- equations_v_m[,2][which(equations_v_m[,1]==min(equations_v_m[,1]))]
        bank_payoff1_viol_mon <- equations_v_m[,3][which(equations_v_m[,1]==min(equations_v_m[,1]))]
      }
      
      return(c(z1_best,K1_best,manager_payoff1_viol_mon,bank_payoff1_viol_mon))
    }
    solve_viol_mon_per1_vect <- Vectorize(solve_viol_mon_per1, vectorize.args = "for_d1")
    
    # Violation, no monitoring
    solve_viol_nomon_per1 <- function()
    {
      Ex1_v_nm <- (s1*sigmas^2/(sigmas^2+sigmaepsilons^2))
        
      Varx1_v_nm <- (3*sigmas^2-sigmas^4/(sigmas^2+sigmaepsilons^2)+3*sigmad^2)
      
      # Do we liquidate or not?
      if(L1>Ex1_v_nm-viols[1]*v)
      {
        ###liquids[4] <- 1
        manager_payoff1_viol_nomon <- 0
        bank_payoff1_viol_nomon <- L1
        z1_best <- "liquidate"
        K1_best <- "liquidate"
      }
      else
      {
        ###liquids[4] <- 0
        # Find optimal z4
        solution_period_2_mini_v_nm <- solution_period_2_dt_v_nm[solution_period_2_dt_v_nm$d1==d_grid[1],]
        to_find_best_z1_viol_nomon <- function(for_z1)
        {
          sol_per2_small_v_nm <- solution_period_2_mini_v_nm[solution_period_2_mini_v_nm$z1==for_z1 & solution_period_2_mini_v_nm$K1==x_grid_small[2],]
          social_payoff_v_nm <- sum((unlist(sol_per2_small_v_nm$bank_payoff)+unlist(sol_per2_small_v_nm$man_payoff))
                               *norm_pdf(unlist(sol_per2_small_v_nm$d2),0,sqrt(sigmad^2+sigmaepsilond^2))*norm_pdf(unlist(sol_per2_small_v_nm$s2),0,sqrt(sigmas^2+sigmaepsilons^2))
                               *(d_grid[2]-d_grid[1])*(s_grid[2]-s_grid[1]))
          return(social_payoff_v_nm)
        }
        to_find_best_z1_viol_nomon_vect <- Vectorize(to_find_best_z1_viol_nomon, vectorize.args = "for_z1")
        social_payoffs_v_nm <- to_find_best_z1_viol_nomon_vect(x_grid_small)
        z1_best <- x_grid_small[which(social_payoffs_v_nm==max(social_payoffs_v_nm))]
        
        # Find optimal K4
        solution_period_2_mini_withoptz1_v_nm <- solution_period_2_mini_v_nm[solution_period_2_mini_v_nm$z1==z1_best,]
        to_find_best_K1_viol_nomon <- function(for_K1)
        {
          sol_per2_small_v_nm <- solution_period_2_mini_withoptz1_v_nm[solution_period_2_mini_withoptz1_v_nm$K1==for_K1,]
          man_payoff <- sum(unlist(sol_per2_small_v_nm$man_payoff)
                            *norm_pdf(unlist(sol_per2_small_v_nm$s2),0,sqrt(sigmas^2+sigmaepsilons^2))*norm_pdf(unlist(sol_per2_small_v_nm$d2),0,sqrt(sigmad^2+sigmaepsilond^2))
                            *(s_grid[2]-s_grid[1])*(d_grid[2]-d_grid[1]))
          equation <- abs(man_payoff-k*(Ex1_v_nm-viols[1]*v-L1))
          bank_payoff <- sum(unlist(sol_per2_small_v_nm$bank_payoff)
                             *norm_pdf(unlist(sol_per2_small_v_nm$s2),0,sqrt(sigmas^2+sigmaepsilons^2))*norm_pdf(unlist(sol_per2_small_v_nm$d2),0,sqrt(sigmad^2+sigmaepsilond^2))
                             *(s_grid[2]-s_grid[1])*(d_grid[2]-d_grid[1]))
          return(c(equation, man_payoff, bank_payoff))
        }
        to_find_best_K1_viol_nomon_vect <- Vectorize(to_find_best_K1_viol_nomon, vectorize.args = "for_K1")
        equations_v_nm <- t(to_find_best_K1_viol_nomon_vect(x_grid_small))
        K1_best <- x_grid_small[which(equations_v_nm[,1]==min(equations_v_nm[,1]))]
        manager_payoff1_viol_nomon <- equations_v_nm[,2][which(equations_v_nm[,1]==min(equations_v_nm[,1]))]
        bank_payoff1_viol_nomon <- equations_v_nm[,3][which(equations_v_nm[,1]==min(equations_v_nm[,1]))]
      }
      
      return(c(z1_best,K1_best,manager_payoff1_viol_nomon,bank_payoff1_viol_nomon))
    }
    
    # Monitor or not?
    # Expected payoff if monitor
    for_expected_payoffs_v <- solve_viol_mon_per1_vect(d_grid)
    bank_expected_payoff_ifmon <- sum(as.numeric(for_expected_payoffs_v[4,])*norm_pdf(d_grid,0,sqrt(sigmad^2+sigmaepsilond^2))*(d_grid[2]-d_grid[1]))
    sol_ifmon_v <- solve_viol_mon_per1(for_d1=ds[1])
    sol_ifnomon_v <- solve_viol_nomon_per1()
    ifelse(bank_expected_payoff_ifmon>=as.numeric(sol_ifnomon_v[4])+mv, mons[1] <- 1, mons[1] <- 0)
    ifelse(bank_expected_payoff_ifmon>=as.numeric(sol_ifnomon_v[4])+mv, z1 <- sol_ifmon_v[1], z1 <- sol_ifnomon_v[1])
    ifelse(bank_expected_payoff_ifmon>=as.numeric(sol_ifnomon_v[4])+mv, K1 <- sol_ifmon_v[2], K1 <- sol_ifnomon_v[2])
    ifelse(bank_expected_payoff_ifmon>=as.numeric(sol_ifnomon_v[4])+mv, bank_payoff <- as.numeric(sol_ifmon_v[4])-mv, bank_payoff <- sol_ifnomon_v[4])
    ifelse(bank_expected_payoff_ifmon>=as.numeric(sol_ifnomon_v[4])+mv, man_payoff <- sol_ifmon_v[3],  man_payoff <- sol_ifnomon_v[3])
    
  }
  
  ifelse(z1=="liquidate", liquids[1] <- 1, liquids[1] <- 0)
  
  return(list(viols=viols, mons=mons, liquids=liquids, z1=z1, K1=K1, man_payoff=man_payoff, bank_payoff=bank_payoff))
}
solve_period_1_comp <- cmpfun(solve_period_1)

solve_period_1_new <- function(x)
{ z0 <- x[1]
K0 <- x[2]
s1 <- x[3]
d1 <- x[4]
solve_period_1_comp(mons=c(NA,NA), viols=c(NA,NA), liquids=liquids,z0=z0,K0=K0,s1=ss[1],s2=ss[2],d1=ds[1],d2=ds[2])}
solve_period_1_new_comp <- cmpfun(solve_period_1_new)
total_grid_forper1 <- expand.grid(x_grid_small, x_grid_small, s_grid, d_grid)
names(total_grid_forper1) <- c("z0","K0","s1", "d1")
clusterExport(cl1, list("solve_period_1","%>%","solution_period_2_dt_v_m", "solve_period_1_comp", "solution_period_2_dt_nv_m", "solution_period_2_dt_v_nm","solution_period_2_dt_nv_nm"))
system.time(f_per1 <- parApply(cl1, total_grid_forper1, 1, solve_period_1_new_comp))
f_per1 <- do.call("cbind",f_per1)
solution_period_1 <- cbind(total_grid_forper1, t(f_per1))
solution_period_1_dt <- as.data.table(solution_period_1)


### Initial period, t=0 ####
solve_period_0 <- function()
{
  func_payoffs <- function(for_z0, for_K0)
  {
    solution_period_1_small <- solution_period_1_dt[solution_period_1_dt$z0==for_z0 & solution_period_1_dt$K0==for_K0,]
    manager_total_payoff <- (sum(as.numeric(unlist(solution_period_1_small$man_payoff))
                                 *(s_grid[2]-s_grid[1])*norm_pdf(as.numeric(unlist(solution_period_1_small$s1)),0,sqrt(sigmas^2+sigmaepsilons^2))
                                 *(d_grid[2]-d_grid[1])*norm_pdf(as.numeric(unlist(solution_period_1_small$d1)),0,sqrt(sigmad^2+sigmaepsilond^2))))
    bank_total_payoff <- (sum(as.numeric(unlist(solution_period_1_small$bank_payoff))
                              *(s_grid[2]-s_grid[1])*norm_pdf(as.numeric(unlist(solution_period_1_small$s1)),0,sqrt(sigmas^2+sigmaepsilons^2))
                              *(d_grid[2]-d_grid[1])*norm_pdf(as.numeric(unlist(solution_period_1_small$d1)),0,sqrt(sigmad^2+sigmaepsilond^2))))
    
    return(c(manager_total_payoff, bank_total_payoff))
  }
  func_payoffs_vect <- Vectorize(func_payoffs, vectorize.args = c("for_z0", "for_K0"))
  payoffs <- t(func_payoffs_vect(x_double_grid_small[,1],x_double_grid_small[,2]))
  payoffs <- data.frame(manager_total_payoff=payoffs[,1],bank_total_payoff=payoffs[,2])
  
  payoffs$to_choose <- as.numeric(payoffs$bank_total_payoff>=I)
  payoffs$to_choose <- payoffs$to_choose*payoffs$manager_total_payoff
  
  best_z0 <- x_double_grid_small[which(payoffs$to_choose==max(payoffs$to_choose)),1]
  best_K0 <- x_double_grid_small[which(payoffs$to_choose==max(payoffs$to_choose)),2]
  bank_payoff_for_bestz0K0 <- payoffs$bank_total_payoff[which(payoffs$to_choose==max(payoffs$to_choose))]
  manager_payoff_for_bestz0K0 <- payoffs$manager_total_payoff[which(payoffs$to_choose==max(payoffs$to_choose))]
  
  #bank_payoff <- payoffs[which(payoffs[,1]==max(payoffs[,1])),3]
  #bank_alt_payoff <- payoffs[which(payoffs[,1]==max(payoffs[,1])),4]
  #man_payoff <- payoffs[which(payoffs[,1]==max(payoffs[,1])),2]
  return(c(best_z0, best_K0, bank_payoff_for_bestz0K0, manager_payoff_for_bestz0K0))
}
solution_period_0 <- solve_period_0()

stopCluster(cl1)























