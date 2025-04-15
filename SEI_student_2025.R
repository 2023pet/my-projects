
gc()
rm(list = ls())
gc()


# load the libraries
library(deSolve)
library(tidyverse)
library(dplyr)


# read in the climate data 'student_climate_met_df.rds'. 
# student_climate_met_df.rds data us save with columes: day, Temp, RH and country"
# Temp is temperature in °C and RH is relative humidity in %

setwd("/home/guest/Desktop/Malaria Data Science/Advanced Modelling")
#**please select the right country data to use for simulation*
met_df=readRDS("/home/guest/Desktop/Malaria Data Science/Advanced Modelling/student_climate_met_df.rds")

head(met_df)

Temp=met_df$Temp

#Temp=mean(met_df$Temp)


RH=met_df$RH
#country <- c("benin", "tanzania","nigeria")

Temp = met_df %>% filter(country=="benin") %>% select(Temp) %>% unlist() %>% unname()
#Temp <- rep(mean(Temp_seasonal) ,366)

RH= met_df %>% filter(country=="benin") %>% select(RH) %>% unlist()%>% unname()

plot(Temp)

plot(RH)
#view(Tem)

################################################################################
#                                                                              #
#              the SEI model for malaria in the mosquito population            #
#                                                                              #
################################################################################

#** the SEI model is incomplete. Complete the transmission system by adding adult mosquito mortality to the system*
SEI_mod <- function(t,y,parms) {
  
  bS<-bS[pmax(1, ceiling(t))] # mosquito birth rate
  bE<-bE[pmax(1, ceiling(t))] # mosquito contact rate
  bI<-bI[pmax(1, ceiling(t))] # mosquito EIP
  fIH<-fIH[pmax(1, ceiling(t))] # fraction of infectious humans
  muM<-muM[pmax(1, ceiling(t))] # mosquito mortality rate
  
  SM=unname(y["SM"]) ;EM=unname(y["EM"]) ;IM=unname(y["IM"]);
  
  #** build out your ODE here using the parameters and states above*
  #** use the variable name muM to save the mortality rates *
  M = SM +EM +IM
  dSM= bS*M- muM*SM - (bE *fIH)*SM#<--- insert your ode
  dEM= (bE *fIH)*SM-  muM*EM - bI *EM#<--- insert your ode
  dIM= bI *EM -   muM*IM #<--- insert your ode
    
  list(c(dSM,dEM,dIM))
  
}

################################################################################
#                                                                              #
#                       temperature-regulated parameters                       #
#                                                                              #
################################################################################

#** mosquito birth rate, based on time between oviposition and adult emergence. *
#** Since we are assessing the role of temperature, assume egg to adult survival is constant and unaffected by rainfall here*
bS=pmax(0,0.000111*Temp*(Temp-14.7)*sqrt(34-pmin(34,Temp))); # from Mordecai et al 2013  and Understanding the link between malaria risk and climate
bS
#** mosquito feed rate, based on duration of oviposition. *
bE=pmin(1,pmax(0,0.017*Temp-0.165))   # from shapiro et al 2017, see supplement

#** parasite development rate *
bI=pmin(1,pmax(0,(0.000112*Temp*(Temp-15.384)*sqrt(35-Temp)))) # from Mordecai et al 2013  and Understanding the link between malaria risk and climate

#**assume a fraction of infection from human populaton is constant= 0.2* 
#** to limit confound and effects to temperature alone*

fIH=rep(0.20,length(Temp));

################################################################################
#                                                                              #
#               Adult mosquito survival models from Lunde et al.               #
#                                                                              #
################################################################################
beta_0 = 0.00113 + RH**2 - 0.158 * RH - 6.61
beta_1 = -2.32 * 1e-4 * RH**2 + 0.0515 * RH + 1.06
beta_2 = 4.1e-6 * RH**2 - 1.09 * 1e-3 * RH - 0.0255

mart1= -0.0016 * Temp ^2 + 0.054* Temp + 0.45#martens 1, see equation 1
#mart1
mart2= exp(-1 /(-4.4+1.31*Temp-0.03*Temp^2))#martens 2, see equation 2

mart3= exp(-1/(-4.31564 + 2.19646 * Temp - 0.058276 * Temp^2))

bayoh_ermet= -2.123 * 10^-7 * Temp^5 + 1.951 * 10^-5 * Temp^4   #bayoh-ermet, see equation 3
-6.394 * 10^-4 * Temp^3 + 8.217 * 10^-3 * Temp^2 
-1.865 * 10^-2 * Temp + 7.238 * 10^-1

bayoh_Parham= -exp(1/(Temp**2 * beta_2 + Temp * beta_1 + beta_0)) # bayoh-Parham, see equation 4

bayoh_mord= -0.000828*Temp^2+0.0367*Temp+0.522

models_df=list(mart1=mart1,
            mart2=mart2,
            mart3=mart3,
            bayoh_ermet=bayoh_ermet,
            bayoh_Parham=bayoh_Parham,
            bayoh_mord=bayoh_mord)

res_df=list()
for(g in 1:length(models_df)){

#** daily mosquito mortality rate*
muM=  -log(models_df[[g]]) # Remember to convert the survival probability to mortality rates, see equation 5
plot(muM)
# set the initial state conditions
Mpop=1e4 # total mosquito population
EM0=1e3 # exposed mosquitoes at time t = 0
IM0=0 # infected & infectious mosquitoes at time t = 0
SM0=Mpop-IM0-EM0 # susceptible mosquitoes at time t = 0

# create the states variable, to store the initial state conditions
states=c(SM=SM0, EM=EM0, IM=IM0)

# set the time step of the simulation process. In this case every day
tmstep=1

# create the full times to be simulated over, in this case 366 days
numDay=365 # number of days

times=seq(1,numDay, by=tmstep)

#simulate the transmission model
sim_out=ode(y=states, times=times, func = SEI_mod, parms = NULL,method=NULL)
plot(sim_out)
#**to get the EIR number of infectious bites per human per day EIR=maZ*
 
m=rowSums(sim_out[,c("SM","EM","IM")])/(Mpop/5) # mosquito density per human. Here assume mosquitoes are 5 times the number of humans from t=0 
a=bE[1:numDay] # mosquito biting rate, i.e. the mosquito contact rate
Z=sim_out[,"IM"]/rowSums(sim_out[,c("SM","EM","IM")]) #Proportion of mosquitoes infectious
EIR=m*a*Z

#** enter and save your model results to be shared on the google sheet*
#*https://docs.google.com/spreadsheets/d/1XACnk16DjikHROdN4Rh7DEvLLxtPZkR-mPi-9Vlm7YU/edit?gid=0#gid=0

res_df[[g]]=data.frame(model= names(models_df)[g],# add the model,
                          country= "benin", # add the country,
                          IM=mean(sim_out[,"IM"]), 
                          avg_Temp=mean(Temp, na.rm = T),
                          EIR=sum(EIR, na.rm = T),
                          seasonal= 0# add wether you are using seasonal or nonseanal temperatuere
                         )

}

do.call(rbind, res_df)
################################################################################
#                                                                              #
#                              Discussion questions                            #
#                                                                              #
################################################################################

#**Q1: Plot a timeseries of the daily number of infectious mosquitoes for each country and model. How do the different temperature models compare to each other?*

# NOTE: for questions 2-4, your results output are needed in the google sheet. We will answer the questions as a class.
# Google link: https://docs.google.com/spreadsheets/d/1XACnk16DjikHROdN4Rh7DEvLLxtPZkR-mPi-9Vlm7YU/edit?gid=0#gid=0

#**Q2: Is the pattern of average number of infectious mosquitoes consistent in the direction of the authors’ finding on transmission efficiency? Remember to use non-seasonal temperatures to answer this* 
#**For Q2 enter the values for:*
#**Average number of infected mosquitoes per year for your countries in column C of the google sheet*
#**Annual average temperature in Celsius for your countries into column D of the google sheet*

#**Q3. How is regional variability affecting the spatial distribution of malaria risk (i.e. the EIR)? Remember to use non-seasonal temperatures to answer this.*
#**For Q3 enter the values for*
#**Total number of infectious bites per person (i.e., annual EIR) in column E of the google sheet for your countries (use non-seasonal temperature to derive this.)*

#**Q4. Compare the annual risk of malaria transmission (i.e., annual EIR) between seasonal and non-seasonal temperature. where is the added value of using seasonality more pronounced than others?*
#**For Q4 enter the values for*
#**Seasonal annual EIR for your countries in column F  of the google sheet *

