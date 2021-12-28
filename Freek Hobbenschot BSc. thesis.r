

library(ggplot2)
library(dplyr)
library(reshape2)
library(fitdistrplus)
library(gridExtra)
stringsAsFactors = FALSE
library(tidyverse)

################################################################################
#CALIBRATION#
################################################################################


#Here we collect data-points from the US Treasury yield curve on 1 July 2008
times <- c(1/12,3/12,6/12,1,2,3,5,7,10,20,30)
i <- c(1.92,1.87,2.13,2.38,2.63,2.90,3.33,3.62,4.01,4.60,4.55)/100

#Here we collect data-points from the US Treasury yield curve on 1 July 2021
#i <- c(0.05,0.05,0.05,0.09,0.25,0.47,0.89,1.24,1.48,2.01,2.07)/100

#See https://www.treasury.gov/resource-center/data-chart-center/interest-rates/Pages/TextView.aspx?data=yieldYear&year=2008

#Now we set up the log of the corresponding bond prices. 
p0 <- 1/(1+i)^times 
log_p0 <- log(p0)

#Given that log(p(0,t))=a+b*sqrt(t)+c*t, we estimate the corresponding coefficients 
regression <- lm(log_p0 ~ I(sqrt(times)) + times) 
pars <- unname(c(regression$coefficients)) #We remove the names such that only the coefficients are left. 

#Here we define a function for the bond prices using the obtained coefficients. 
log_p <- function(p=pars, t) {  #define function for log(p*(0,t)), using estimated coefficients
  p[1] + p[2]*sqrt(t) + p[3]*t
}

#Now we want to determine if we have a good fit. 
test_fit <- as.matrix(lapply(1:30,log_p,p=pars),ncol=1) #For times is 1 till 30 we look at the corresponding price and store this in a data-list. 
plot(1:30,test_fit,col="blue") + points(times,log_p0, col="red") #We plot both to see if our approximation comes close to the points given by US Treasury. 


#The coefficients give an approximation that is close to the observed values from US Treasury. 
#The coefficients will now be used to define the following functions: 
fitted_f <- function(p=pars, t) {  
  -(p[2]*(1/2)*(1/sqrt(t))+p[3])
}
fitted_fT <- function(p=pars, t) {
  p[2]*(1/4)*t^(-3/2)
}

fun_g <- function(sigma_r, alpha, t) {
  sigma_r^2/(2*alpha^2)*(1-exp(-alpha*t))^2
}

deriv_g <- function(sigma_r, alpha, t) {
  sigma_r^2/alpha*(1-exp(-alpha*t))*exp(-alpha*t)
}



#Now we use these functions to define theta, which is a function of time that is used to determine the short rate. 
theta <- function(sigma_r, alpha, t) {
  
  ff <- fitted_f(t=t)
  ffT <- fitted_fT(t=t)
  fg <- fun_g(sigma_r=sigma_r, alpha=alpha, t=t)
  dg <- deriv_g(sigma_r=sigma_r, alpha=alpha, t=t)
  
  outcome <- ffT + dg + alpha*(ff + fg)
  as.matrix(outcome)
}


#Now we want to look at the price of an annuity.

#We first fit a function for the bond prices and take a look at the coefficients that follow from this.
reg_p0 <- lm(p0 ~ I(sqrt(times)) + times + I(times^2))  #fit function for bond prices p(0,t)
pars_p0 <- reg_p0$coefficients

#Now we use these coefficients to define a function for these prices.
fun_p0 <- function(t) {
  pars_p0[1] + pars_p0[2]*sqrt(t) + pars_p0[3]*t + pars_p0[4]*t^2
}


#We again want to determine if these coefficients gives us a good fit. 
fx <- as.matrix(lapply(1:30,fun_p0),ncol=1)
plot(1:30,fx,col="blue") + points(times,p0, col="red")

#We argue that we have a good fit and continue with the calibration of the future bond prices. 
#We first define BtT and continue with ptT.
BtT <- function(alpha, t, m) {
  (1/alpha)*(1-exp(-alpha*(m-t)))
}


ptT <- function(alpha, sigma_r, t, m, rate) { # m = maturity t = time
  ifelse(m<=t, 1, fun_p0(t=m)/fun_p0(t=t)*exp(BtT(alpha=alpha, t=t, m=m)*fitted_f(t=t)-(sigma_r^2/(4*alpha))*(BtT(alpha=alpha, t=t, m=m))^2*(1-exp(-2*alpha*t))-BtT(alpha=alpha, t=t, m=m)*rate))
}

#The bond price is used to determine the value of an annuity. 
#This value is also depended on the short rate and therefore used within the function that we create for the GMIB. 
#For more information on calibration see Bjork (2020)

################################################################################
#FUNCTION FOR A GMIB WITHOUT RESET OPTION#
################################################################################

disc_YT = function(steps = 520, acc = 10, term = 20, g = 0.05, rg = 0.02, c = 0.02, Sf0 = 1000, initial = 1000, sigma_s = 0.2, r0 = 0.0238, alpha = 0.35, sigma_r = 0.015, dummy =0, maxdummy=0) {
  
  #Initial is the lump sum paid by the Policyholder at the start of the contract.
  #Sf0 is the value in the investment account at the start of the GMIB that is valuated (either at t=o or when valuing the reset option at t=T with T being the time of annuity)
  #dummy is the number of years between the start of the contract and the start of the GMIB. When using this function to value the reset option, this dummy is T. 
  #maxdummy is a dummy for max Sf_year before the GMIB; When valuing a reset option, we also take into account the maximum of Sf in the years upto time T. 
  #The distinction between Initial and Sf0, and the addition of dummy and maxdummy are made to use the function disc_YT also to valuate the worth of the reset option. 
  
  #(r0 is set to be equal to the Daily Treasury Par Yield Curve Rate of 1Yr for July 1, 2008. It is 0.0238.)
  #(The values of the other parameters are chosen based on what Marshall et al. 2008 have chosen. See thesis.)
  
  #'steps' is the total number of steps, which is the accumulation time times the number of intervals in each year (e.g. 10x52)
  
  #Starting Conditions
  dt <- acc/steps #accumulation period divided by steps wanted 
  t_times <- seq(0, acc, dt) #sequence creation from 0 till maturity date by dt steps
  v <- c() #empty vector to collect discounted payoff values per simulation
  at <- c() #empty vector to collect the value of the annuity related to the contract at maturity
  bb <- c() #empty vector to collect benefit base values per simulation
  yt <- c() #empty vector to collect payoff values per simulation
  Sf_year <- c() #empty vector to collect end of the year account values used to compute the benefit base and fee
  
  theta_t <- as.matrix(lapply(1:steps, theta, alpha = alpha, sigma_r = sigma_r))
  r <- c(r0, rep(NA, steps)) #vector short rate
  dSf <- c() #change in value of the investment account
  dr <-  c() #change in value of the short rate
  #starting values 
  year <- 0 #start value
  BBfee  <- 0 #start value

  #Defining Sf_year and Sf to store the investment account values. 
  Sf_year <- matrix(ncol = 1, nrow = steps + 1, data = NA) 
  Sf_year[ ,1] <- c(Sf0[1], rep(NA, steps))
  Sf <- matrix(ncol = 1, nrow = steps + 1, data = NA)
  Sf[ ,1] <- c(Sf0[1], rep(NA, steps))
  #Defining a vector to store values of the benefit base and the fees that are paid by the policyholder. 
  BBfee_vector <- rep(0, acc)
  delta_2 <- rep(0, acc)
 
  #Now we continue with a for-loop that defines the value of the investment account at each time-step. 
  
  step_year <- c()
  year_vector <- rep(0, acc)
  
  for (j in 1:steps) {
    
    
    dr[j] <- alpha*(theta_t[[j]]/alpha-r[j])*dt + sigma_r*rnorm(1, mean = 0, sd = sqrt(dt))
    r[j + 1] <- r[j] + dr[j]
    dSf[j] <- r[j]*Sf[j, 1]*dt + sigma_s*Sf[j, 1]*rnorm(1, mean = 0, sd = sqrt(dt))
    
    #Recalculating Benefit Base and subtracting annual fee
    if ((j - 1) %% (steps / acc) == 0 | j == steps) {
      step_year[j] = j
      
      for (l in 1:acc) {
        if (j - 1 > steps / acc * l) { #Fees get only yearly reduced, so we look at the values at the start of each year. 
          year_vector[l] <- year_vector[l] + 1
          BBfee_vector[l] <- max(initial * (1 + rg)^(dummy+year_vector[l]), Sf[j, 1], maxdummy, na.rm = T) 
          delta_2[l] <- c * BBfee_vector[l] #Calculating the fee before subtracting it. 
          Sf[j, 1] <- Sf[j, 1] - delta_2[l]
          Sf_year[j, 1] <- Sf[j, 1] 
        }
      }
    }
    Sf[j + 1,] <- Sf[j,] + dSf[j]
  }
  
  #This part of the codes is used to prevent negative values or NULL returns. 
  #If the value of the investment account reaches 0, it will stay 0 till the end of the contract. 
  Sf[Sf<0]<-0
  Sf_df_neg <- as.data.frame(Sf) %>% mutate(t_times = t_times)
  Sf_df <- as.data.frame(Sf_df_neg[, apply(Sf_df_neg, 2, function(x) !any(x < 0, na.rm = T))])

  #Now we price the Price the annuity.
  pTTi <- lapply((acc+dummy):((acc+dummy)+term-1),ptT,m=term,sigma_r=sigma_r,alpha=alpha, rate = tail(r,1))
  aT <- Reduce("+",pTTi)
  print(aT)
  #Now we look at the final pay-off of the contract, assuming the policyholder is rational and chooses the maximum value.
  second <- max(Sf_year[,1], na.rm=T) #'Second' is the maximum of the investment account values at the start of each year. 
  third <- max(second,maxdummy) #We make use of a dummy 'maxdummy' which is (if applicable) the maximum of the values each year in the investment account before the period of reset.
                                #'Third' maximizes over the maximum values each year in the investment account, both for the period before reset and for the period after reset.
  BBT <- max(initial*(1 + rg)^((dummy+acc)), third) 
  YT <- max(BBT*g*aT, tail(Sf_df$V1,1)) 
 
  
  disc <- exp(-acc*mean(r)) #Here we use the mean of the short rate to determine the discount factor
  Pay_Off <- YT #pay-off
  disc_YT_portf <- disc*YT #discounted pay-off
  datalist = list(aT, disc_YT_portf, tail(r,1), tail(Sf_df$V1,1), Pay_Off, disc, second) #Here we create a data-list for the output of the function
  names(datalist) = c("koek","discounted pay off", "final interest value", "final investment value", "pay off", "discount factor", "maximum value investment account")

  
  return(datalist)
} 

set.seed(8)
disc_YT(10)$"koek"

#################################################################################
#MONTE CARLO FOR A GMIB WITHOUT RESET OPTION#
#################################################################################

set.seed(8)

#number of simulated pay-offs
#The difference of 20k simulations and 5k simulations resulted in a difference of 0.2% in pay-off. 
nsim = 1000
STORE <- rep(NA,nsim)

start.time <- Sys.time()

#Performing a simple Monte Carlo
for(i in 1:nsim){
  
  #For each simulation take a new sample
  STORE[i] <- disc_YT(acc=10)$'discounted pay off'
  
}

end.time <- Sys.time()
end.time-start.time

mean(STORE)

summary(STORE)

#STORE

#################################################################################
#SEQUENTIAL MONTE CARLO FOR A GMIB WITH RESET OPTION#
#################################################################################


set.seed(8)

#number of simulated pay-offs
nsim2 = 1000

NMS <- rep(NA,nsim2)
YES <- rep(NA,nsim2)
interest <- rep(NA,nsim2)
start.time <- Sys.time()
#For each simulation we perform a Monte Carlo for the value of the reset option
for(i in 1:nsim2){

nsim1=1000

#For every simulation in nsim2, new values for r0, sf0 and maxdummy will be assigned. 
TEST <- disc_YT(acc=10)

#Define the time in years of the 'reset period' (the annuitization time of the reset period).
accY1 <- 1


#Here we have nsim samples for the reset option
Y1_valuation <- replicate(nsim1, disc_YT(steps = 52, acc=accY1, r0=TEST$'final interest value', initial = 1000, Sf0=TEST$'final investment value', dummy = 10, maxdummy = TEST$second))
 
 #We create a list with pay-poff values that followed from exercising the reset option.
Y1_payoff = rep(NA,nsim1)
 for(j in 1:nsim1){
   Y1_payoff[j]=Y1_valuation[1,j]}

 #Expected Pay Off reset option
Y1_payoff <- unlist(Y1_payoff, recursive = TRUE, use.names = TRUE)
Y1_payoff <- mean(Y1_payoff)

 #We create a list with short rate values that followed from exercising the reset option.
Y1_interest = rep(NA,nsim1)
 for(j in 1:nsim1){
   Y1_interest[j]=Y1_valuation[2,j]}

 #Expected interest reset option
Y1_interest <- unlist(Y1_interest, recursive = TRUE, use.names = TRUE)
Y1_interest <- exp(-accY1*mean(Y1_interest))

 #We take the maximum of the payoff at time 10 and the 1 year discounted payoff at time 11. Then we take the discount back to time 0. 
final_payoff <- TEST$'discount factor' * max(TEST$'pay off', Y1_interest*Y1_payoff)

 #Yes if exercising the reset option will increase the contract value, otherwise no.
yes<- if(TEST$'pay off'<Y1_interest*Y1_payoff) 1 else 0

  
 #Storing the outcomes of each sample in a vector
NMS[i] <- final_payoff
YES[i] <- yes
}
end.time <- Sys.time()

end.time-start.time


#OUTCOMES

100 * mean(YES) #The percentage of times the policyholder should use the reset option.

mean (NMS) #The expected payoff of the contract.

max(NMS) #The maximum payoff that followed from the simulations.
min(NMS) #The minimum payoff that followed from the simulations.
#NMS



##############################################################################
#PLOTS#
##############################################################################
library(ggplot2)
library(tidyr)


#We first make a plot for the GMIB without reset option using set.seed=8, nsim=10000
#From c=0 till c=0.1, stepsize 0.02. 
#I also add c=0.01 because of the relatively large change in V(c) between 0 and 0.2 and to get a better fit of the curve.  
x <-  c(0,0.01,0.02,0.04,0.06,0.08,0.1) #c-values
y1 <- c(1113.256,875.2879,789.0488,753.0723,746.1454,743.2793,741.9217) #g=4% check
y2 <- c(1167.089,955.5961,880.7809,847.1624,839.4136,836.1892,834.6619) #g=4.5% check
y3 <- c(1232.338,1042.229,974.7003,941.2831,932.6818,929.0991,927.4021) #g=5% check
y4 <- c(1309.973,1134.219,1070.127,1035.407,1025.95,1022.009,1020.142) #g=5.5% check
y5 <- c(1412.843,1233.34,1166.775,1129.535,1119.218,1114.919,1112.883) #g=6% check

df <- data.frame(x,y1,y2,y3,y4,y5)

data2 <- df %>%
  pivot_longer(y1:y5, names_to = "group", values_to = "y")

###TEST 1#####
ggplot(data2, aes(x, y, color = group, linetype=group)) +
  #geom_point(size = 3) +    # increased size for increased visibility
  geom_smooth(method = "auto", se = FALSE)+
  scale_linetype_manual("g-value", values=c(2,3,4,5,6), labels=c("4%", "4.5%", "5%", "5,5%", "6%"))+
  #scale_shape_manual("g-value", values=c(7, 8, 9, 10, 11), labels=c("4%", "4.5%", "5%", "5,5%", "6%")) +
  scale_color_manual("g-value", values=c(2, 3, 4, 5, 6), labels=c("4%", "4.5%", "5%", "5,5%", "6%"))+
  geom_line(aes(y=1000), linetype = "solid", color="darkgreen", lwd = 1) +
  ylab("V(c)") + xlab("c") + labs(title="Valuation of the GMIB without reset option")+
  theme_light() +
  theme(plot.title = element_text(color="black", size=12, face="italic", hjust = 0.5))

hist(STORE, prob=TRUE)
hist(NMS)

##########################################
#We now make a plot for the GMIB with reset option using set.seed=8, nsim1=1000, nsim2=1000
x <-  c(0.0,0.01,0.02,0.04,0.06,0.08,0.1) #c-values
z1 <- c(1119.206,882.6474,793.6509,755.8901,747.0849,744.2382,742.8573) #g=4% check
z2 <- c(1173.412,962.9609,886.099,850.0865,840.4705,837.268,835.7144) #g=4.5%
z3 <- c(1240.265,1049.009,979.9796,944.3013,933.8561,930.2978,928.5716) #g=5%
z4 <- c(1319.934,1141.204,1075.376,1038.675,1027.242,1023.328,1021.429) #g=5.5%
z5 <- c(1422.448,1240.374,1172.038,1133.1,1120.627,1116.357,1114.286) #g=6%

#set.seed-8, nsim1 = 1000 and nasim2 = 1000 resulting in a simulation time of 30 minutes per simulation
#Increasing c by 0.02, then...;

df2 <- data.frame(x,z1,z2,z3,z4,z5)

data2 <- df2 %>%
  pivot_longer(z1:z5, names_to = "group", values_to = "z")

ggplot(data2, aes(x, z, color = group, linetype=group)) +
  #geom_point(size = 3) +    # increased size for increased visibility
  geom_smooth(method = "auto", se = FALSE)+
  scale_linetype_manual("g-value", values=c(2,3,4,5,6), labels=c("4%", "4.5%", "5%", "5,5%", "6%"))+
  #scale_shape_manual("g-value", values=c(7, 8, 9, 10, 11), labels=c("4%", "4.5%", "5%", "5,5%", "6%")) +
  scale_color_manual("g-value", values=c(2, 3, 4, 5, 6), labels=c("4%", "4.5%", "5%", "5,5%", "6%"))+
  geom_line(aes(y=1000), linetype = "solid", color="darkgreen", lwd = 1) +
  ylab("V(c)") + xlab("c") + labs(title="Valuation of the GMIB with reset option")+
  theme_light() +
  theme(plot.title = element_text(color="black", size=12, face="italic", hjust = 0.5))

#COLECTING DATA MANUAL AFTER RUNNING EACH SIMULATION
# max at g = 4% goes from 7892.568 to 5021.334 to ? to 1491.886 to 1231.827 to 1082.257 to 1064.547
# min at g = 4% goes from 555.1874 540.1614 to to ? to 516.9621 to 516.9621 to 516.9621 to 516.9621
# %YES at g = 4% goes from 13.3 to 9 to 7.4 to 6.2 to 6.6 to 6.7 to 6.7 

# max at g = 4.5% goes from 7892.568 to 5021.334 to 3180.017 to 1678.372 to 1385.805 to 1217.539 to 1197.616
# min at g = 4.5% goes from 618.2458 to 607.6815 to 604.2027 to 581.5823 to 581.5823 to 581.5823 to 581.5823
# %YES at g = 4.5% goes from 13 to 8.7 to 7.2 to 6.1 to 6.6 to 6.7 to 6.7

# max at g = 5% goes from 7892.568 to 5021.334 to 3180.017 to 1864.858 to 1539.784 to 1352.822 to 1330.684
# min at g = 5% goes from 686.9398 to 675.2017 to 671.3416 to 646.2026 to 646.2026 to 646.2026 to 646.2026
# %YES at g = 5% goes from 12.7 to 7.9 to 6.8 to 6.1 to 6.6 to 6.7 to 6.7

# max at g = 5.5% goes from 7892.568 to 5021.334 to 3180.017 to 2051.344 to 1693.762 to 1488.104 to 1463.752
# min at g = 5.5% goes from 755.6337 to 742.7219 to 738.4757 to 710.8229 to 710.8229 to 710.8229 to 710.8229
# %YES at g = 5.5% goes from 16.2 to 9.3 to 6.2 to 6.1 to 6.6 to 6.7 to 6.7

# max at g = 6% goes from 8102.401 to 5147.968 to 3255,419 to 2237.83 to 1847.74 to 1623.386 to 1596.821
# min at g = 6% goes from 821.1475 to 810.242 to 805.6099 to 775.431 to 775.4431 to 775.4431 to 775.4431
# %YES at g = 6% goes from 11.3 to 7.3 to 6.1 to 6.1 to 6.6 to 6.7 to 6.7




################
#APROXIMATIONS#
################

locator()

#We have bijective data
#This approximation might slightly differ from the approximation of the line that follows from geom_smooth, 
#however this gives us a close approximation of c corresponding to different g's. 

#For the GMIB without reset option
f1 <- approxfun(df$x, df$y1) #4%
f1(0.00476)

f2 <- approxfun(df$x, df$y2) #4.5% 
f2(0.00790)

f3 <- approxfun(df$x, df$y3) #5%
f3(0.01625)

#For the GMIB with reset option (use different df)
g1 <- approxfun(df2$x, df2$z1) #4%
g1(0.00501)

#0.480% to 0.501%

g2 <- approxfun(df2$x, df2$z2) #4.5%
g2(0.00824)

# 0.790% to 0.820%

g3 <- approxfun(df2$x, df2$z3) #5% 
g3(0.0171)

#From 1.63% to 1.71%


##########################################
#SENSITIVITY ANALYSIS#
##########################################
#We now make a plot for the GMIB with reset option for yield curves 2008 and 2021, g value 5 and 5.5 using set.seed=8, nsim1=1000, nsim2=1000
#Note that when using the yield curve rate for 2021, the r0 value changes to 0.09


x <-  c(0.0,0.01,0.02,0.04,0.06,0.08,0.1) #c-values
w1 <- c(1240.27,1049.01,979.98,944.31,933.86,930.30,928.57) #g=5%, 2008
w2 <- c(1319.93,1141.20,1075.38,1038.68,1027.24,1023.33,1021.43) #g=5.5%, 2008
w3 <- c(1305.51,1115.89,1042.32,993.84,974.01,963.18,956.61) #g=5%, 2021
w4 <- c(1411.08,1221.24,1145.49,1093.22,1071.41,1059.50,1052.27) #g=5.5%, 2021


#set.seed-8, nsim1 = 1000 and nasim2 = 1000 resulting in a simulation time of 30 minutes per simulation
#Increasing c by 0.02, then...;

df3 <- data.frame(x,w1,w2,w3,w4)

data3 <- df3 %>%
  pivot_longer(w1:w4, names_to = "group", values_to = "w")

ggplot(data3, aes(x, w, color = group, linetype=group)) +
  #geom_point(size = 3) +    # increased size for increased visibility
  geom_smooth(method = "auto", se = FALSE)+
  scale_linetype_manual("g-value", values=c(2,3,4,5,6), labels=c("5%, 2008", "5.5%, 2008", "5%, 2021", "5,5%, 2021"))+
  scale_color_manual("g-value", values=c(2, 3, 4, 5, 6), labels=c("5%, 2008", "5.5%, 2008", "5%, 2021", "5,5%, 2021"))+
  geom_line(aes(y=1000), linetype = "solid", color="darkgreen", lwd = 1) +
  ylab("V(c)") + xlab("c") + labs(title="Valuation at different yield curve rates")+
  theme_light() +
  theme(plot.title = element_text(color="black", size=12, face="italic", hjust = 0.5))

