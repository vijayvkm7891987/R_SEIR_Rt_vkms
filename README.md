# R_SEIR_Rt_vkms
---
title: "covid_19_India_Rt_SEIR"
author: "Vijay_Kumar_Mishra(VKM)"
date: "01/05/2020"
output: ioslides_presentation
---

*******SB_imp_covid.india_defining_vector********
###V1<-c("2020-01-30","2020-01-31","2020-02-01","2020-02-02","2020-02-03","2020-02-04","2020-02-05","2020-02-06","2020-02-07","2020-02-08","2020-02-09","2020-02-10","2020-02-11","2020-02-12","2020-02-13","2020-02-14","2020-02-15","2020-02-16","2020-02-17","2020-02-18","2020-02-19","2020-02-20","2020-02-21","2020-02-22","2020-02-23","2020-02-24","2020-02-25","2020-02-26","2020-02-27","2020-02-28","2020-02-29","2020-03-01","2020-03-02","2020-03-03","2020-03-04","2020-03-05","2020-03-06","2020-03-07","2020-03-08","2020-03-09","2020-03-10","2020-03-11","2020-03-12","2020-03-13","2020-03-14","2020-03-15","2020-03-16","2020-03-17","2020-03-18","2020-03-19","2020-03-20","2020-03-21","2020-03-22","2020-03-23","2020-03-24","2020-03-25","2020-03-26","2020-03-27","2020-03-28","2020-03-29","2020-03-30","2020-03-31","2020-04-01","2020-04-02","2020-04-03","2020-04-04","2020-04-05","2020-04-06","2020-04-07","2020-04-08","2020-04-09","2020-04-10","2020-04-11","2020-04-12","2020-04-13","2020-04-14","2020-04-15","2020-04-16","2020-04-17","2020-04-18","2020-04-19","2020-04-20","2020-04-21","2020-04-22","2020-04-23", "2020-04-24", "2020-04-25","2020-04-25", "2020-04-26","2020-04-27","2020-04-28")

###V2<- c(1,1,1,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5,5,28,30,31,34,39,43,56,62,73,82,102,113,119,142,156,194,244,330,396,499,536,657,727,887,987,1024,1251,1397,1998,2543,2567,3082,3588,4778,5311,5916,6725,7598,8446,9205,10453,11487,12322,13430,14352,15722,17615,18539,20080,20178,23077,24530,26283,27890,29451,31324)

##covid.ind=data.frame(V1, V2)######

##View(covid.ind)

covid.ind <-c(1,1,1,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5,5,28,30,31,34,39,43,56,62,73,82,102,113,119,142,156,194,244,330,396,499,536,657,727,887,987,1024,1251,1397,1998,2543,2567,3082,3588,4778,5311,5916,6725,7598,8446,9205,10453,11487,12322,13430,14352,15722,17615,18539,20080,20178,23077,24530,26283,27890,29451,31324)

#####estimation of time dependent reproduction number using a Bayesian approach for India############
library(R0)

mGT_covid <- generation.time("gamma", c(8.4,3.8))

SB_covid<-est.R0.SB(covid.ind, mGT_covid)

SB_covid

windows()  

SB_covid[["R"]]
SB_covid[["conf.int"]]

windows()

plot(SB_covid)

plotfit(SB_covid) 

#########e@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##################################################
**************************************************
remove (list = objects() )
####SEIR_Modeling_India_VKM########################################
library(deSolve)
#####creating function to return derivatives#######################
seir_model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  S = state_values [1]        # susceptibles
  E = state_values [2]        # exposed
  I = state_values [3]        # infectious
  R = state_values [4]        # recovered
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
         {
           # compute derivatives
           dS = (-beta * S * I)
           dE = (beta * S * I) - (delta * E)
           dI = (delta * E) - (gamma * I)
           dR = (gamma * I)
           
           # combine results
           results = c (dS, dE, dI, dR)
           list (results)
         }
    )
}

####parameters#############################################################
contact_rate = 10                     # number of contacts per day
transmission_probability = 0.07       # transmission probability
infectious_period = 7                 # infectious period
latent_period = 2.3                   # latent period
####computation of beta(transmission rate), gamma(recovery rate)###########
beta_value = contact_rate * transmission_probability
gamma_value = 1 / infectious_period
delta_value = 1 / latent_period

#####computation of reproduction number####################################
Ro = beta_value / gamma_value

####parameters#############################################################
parameter_list = c (beta = beta_value, gamma = gamma_value, delta = delta_value)

#####initial value for subpopulation
W = 1330830990  # susceptible hosts
X = 1           # infectious hosts
Y = 0           # recovered hosts
Z = 9           # exposed hosts

######computation of total population##########################################
N = W + X + Y + Z

#####initial values for df
initial_values = c (S = W/N, E = X/N, I = Y/N, R = Z/N)

daypoints = seq (0, 150, by=1)

#####simulate the SEIR###################################################
output = lsoda (initial_values, daypoints, seir_model, parameter_list)

####plotting dynamics of susceptible, exposed, infected and recovered popn
windows()
# susceptible hosts over time
plot (S ~time, data = output, type='l', ylim = c(0,1), col ='orange', ylab = 'proportion', main = 'SEIR Modeling of COVID-19 outbreak,India (30Jan-27June, 2020)') 

# remain on same frame
par (new = TRUE)    

# exposed hosts over time
plot (E ~ time, data = output, type='l', ylim = c(0,1), col = 'yellow', ylab = '', axes = FALSE)

# remain on same frame
par (new = TRUE) 

# infectious hosts over time
plot (I ~ time, data = output, type='l', ylim = c(0,1), col = 'red', ylab = '', axes = FALSE) 

# remain on same frame
par (new = TRUE)  

# recovered hosts over time
plot (R ~ time, data = output, type='l', ylim = c(0,1), col = 'green', ylab = '', axes = FALSE)

legend(145,0.7,legend=c("S","E","I","R"),col=c("orange","yellow","red","green"), lty=1,cex=0.8)
##end##################################################################################
