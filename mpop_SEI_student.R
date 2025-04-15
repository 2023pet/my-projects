gc()
rm(list=ls())
rm()

#** load needed libraries*

library(deSolve)
library(zoo)
library(tidyverse)

#** read in your climate data and select your country from the dataset*
#*
setwd("/home/guest/Desktop/Malaria Data Science/Advanced Modelling")

#*
all_met_df=readRDS("/home/guest/Desktop/Malaria Data Science/Advanced Modelling/student_met_data.rds")

getwd()

met_df=all_met_df %>% filter(country=="drc")

#** dataframe to store your scenario information *
#
scen_df=data.frame(country="drc",scenario_name=c("climatology","scen_1.9","scen_4.5","scen_8.5"),
                   Temp_scen=c("climatology_Temp","dy_Temp_1.9","dy_Temp_4.5","dy_Temp_8.5"),
                   ppt_scen=c("climatology_ppt","dy_ppt_1.9","dy_ppt_4.5","dy_ppt_8.5"),
                   start_year=c(1988,1988,1988,1988),
                   end_year=c(2023,2035,2035,2035))


res_df=list() # list to store scenario results
for(v in 1:nrow(scen_df)) try({
  country=scen_df[v,]$country
  scenario_name=scen_df[v,]$scenario_name
  ppt_scen=scen_df[v,]$ppt_scen
  Temp_scen=scen_df[v,]$Temp_scen
  start_year=scen_df[v,]$start_year
  
  end_year=scen_df[v,]$end_year

ppt=unlist(unname(met_df[met_df$year>=start_year & met_df$year<=end_year,][,ppt_scen]))  #** <- load you rainfall data, with the correct date range*
Temp=unlist(unname(met_df[met_df$year>=start_year & met_df$year<=end_year,][,Temp_scen]))  #** <- load you Temperature data, with the correct date range*

# get the year and day of the year too
days=unlist(unname(met_df[met_df$year>=start_year & met_df$year<=end_year,]$day))
years=unlist(unname(met_df[met_df$year>=start_year & met_df$year<=end_year,]$year))

#** create the mosquito population model equations from equation (1) in White et al  2011 here*
mpop <- function(t,y,parms) {
  
  beta <- beta[pmax(1,ceiling(t))] # daily egg laying rate
  d_E <- d_E[pmax(1,ceiling(t))] # duration of egg development
  d_L <- d_L[pmax(1,ceiling(t))] # duration of larvae development
  d_P <- d_P[pmax(1,ceiling(t))] # duration of pupae development
  mu0E <- mu0E[pmax(1,ceiling(t))] # baseline daily mortality rate of eggs
  K <- K[pmax(1,ceiling(t))] # carrying capacity
  mu0L <- mu0L[pmax(1,ceiling(t))] # baseline daily mortality rate of larvae
  muP <- muP[pmax(1,ceiling(t))] # daily mortality rate of pupae
  muM <- muM[pmax(1,ceiling(t))] # daily mortality rate of adult
  gamma <- gamma[pmax(1,ceiling(t))] # relative density dependent morality
  
  # transmission parameters
  bE<-bE[pmax(1, ceiling(t))] # mosquito contact rate
  bI<-bI[pmax(1, ceiling(t))] # mosquito EIP
  fIH<-fIH[pmax(1, ceiling(t))] # fraction of infectious humans
  
  # the states
  E=y["E"];L=y["L"];P=y["P"]
  SM=unname(y["SM"]) ;EM=unname(y["EM"]) ;IM=unname(y["IM"]);
  
#**build out your coupled model here*
  M = SM +EM +IM
  dE= beta*M-E/d_E-   mu0E*((1+(E+L)/K))*E#eggs
  dL=  E/d_E-L/d_L-  mu0L*((1+gamma*(E+L)/K))*L#larvae  
  dP= L/d_L-P/d_P -muP*P #pupae
  
  dSM= 1/2*(P/d_P)- muM*SM - (bE *fIH)*SM#<--- insert your ode
  
  dEM= (bE *fIH)*SM-  muM*EM - bI *EM#<--- insert your ode
  dIM= bI *EM -   muM*IM #<--- insert your ode
  
    list(c(dE,dL,dP,dSM,dEM,dIM))
}

# set the duration of simulation in days. Use the length of rain data
numDay=length(ppt)

############################# carrying capacity  #############################
#** use the carrying capacity from equation (8) in White et al 2011* 
#** you will need tau, mean_rain and lambda to derive K the carrying capacity*

tau=7 #** <-insert the length of days over which rainfall contributes to the carrying capacity*

mean_rain=rollmean(ppt, k=tau, fill=NA, align="right") #** <-insert the rolling mean rainfall during the past tau days*

lambda =1e6 #** the site scaling factor*

K=lambda*mean_rain #** compute K the site-adjusted carrying-capacity*
K[is.na(K)]=mean(K, na.rm = T) # replace missing values in K with the average

############################## development times ##############################
#** instead of the d_E,d_L and d_P in White et al, use the temperature-regulated development here*
#** See the Table 2 exercise handout for temperature formulation for each immature stage *

# temperature-regulated egg duration
a_d_E=1.3346 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out *
b_d_E=13.1687 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
c_d_E=13.8354 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
d_d_E=5.6711 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *

d_E=a_d_E+(b_d_E/(1+((Temp/c_d_E)^d_d_E))) #** <- compute temperature-regulated development duration here *

# temperature-regulated egg to larave duration
a_d_L=15.405231 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
b_d_L=-8.461096 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
c_d_L=22.166679 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
d_d_L=-16.260036 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
  
d_L=a_d_L+(b_d_L/(1+((Temp/c_d_L)^d_d_L))) #** <- compute temperature-regulated development duration here *
  
# temperature-regulated larvae to pupae duration
a_d_P=10.740867 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
b_d_P=-9.549660 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
c_d_P=12.876982 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
d_d_P=-4.231138 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
  
d_P=a_d_P+(b_d_P/(1+((Temp/c_d_P)^d_d_P))) #** <- compute temperature-regulated development duration here *
  
#################################################################################

############################ daily mortality rates ############################
#** instead of the muE0,muL and muP in White et al, use the temperature and moisture regulated mortality here*
#** See the  Table 2 exercise assignment for temperature and moisture formulations for each immature stage *

############################ moisture regulated mortality rates ############################
#** But first calculate the dessication days D and the probability immature survival due to moisture following Parham et al 2012 *

#** for dessication, i.e., drying days, how many days has it been without rain (without rain defined as rainfall less than 1 mm)*
#** run the lines of code below to calculate D*
D=numeric()
D[1]=0

for(i in 2:length(ppt)){
  # get the last rain event
  rain_event=tail(which(ppt[1:(i-1)]>=1),1)
  
  if(length(rain_event)==1){
  # how many days since it last rained
  D[i]=length(which(ppt[(rain_event):i]<1))
  }else{
    D[i]=0  
  }
  }

#** Egg *
wE=0.405 #** <- insert egg sensitivity to dessication here*
p_E_R=2*exp(-wE*D)/(1+exp(-wE*D)) #** <- insert formulation relating dessication to Egg survival probability*

#** Larvae *
wL=.855 #** <- insert larvae sensitivity to dessication here*
p_L_R=2*exp(-wL*D)/(1+exp(-wL*D)) #** <- insert formulation relating dessication to larvae survival probability*

#** Pupae *
wP=0.602 #** <- insert larvae sensitivity to dessication here*
p_P_R=2*exp(-wP*D)/(1+exp(-wP*D)) #** <- insert formulation relating dessication to larvae survival probability*
  
################################################################################
################# temperature regulated survival probability ###################

#** Egg *
a_p_E_T=-0.25223 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out *
b_p_E_T=0.0936 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
c_p_E_T=-0.002026 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *

p_E_T=(a_p_E_T+(b_p_E_T*Temp)+(c_p_E_T*Temp^2)) #** <- compute temperature-regulated Egg survival probability due to temperature here *

#** Larvae *
a_p_L_T=0.3594 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
b_p_L_T=0.04983 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
c_p_L_T=-0.0010 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
  
p_L_T=(a_p_L_T+(b_p_L_T*Temp)+(c_p_L_T*Temp^2)) #** <- compute temperature-regulated pupae survival probability due to temperature here *
  
#** Pupae *

a_p_P_T=0.19067 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
b_p_P_T=0.02643 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
c_p_P_T=-5.3e-05 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
  
p_P_T=(a_p_P_T+(b_p_P_T*Temp)+(c_p_P_T*Temp^2)) #** <- compute temperature-regulated pupae survival probability due to temperature here *
  
  
#** Adult *
a_p_M_T=-0.000828 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
b_p_M_T= 0.0367#** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
c_p_M_T=0.522 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
  
p_M_T=a_p_M_T* Temp^2 + b_p_M_T *Temp + c_p_M_T #** <- compute temperature-regulated adult survival probability due to temperature here *
  
################################### daily mortality rates #####################################

mu0E=-log(p_E_R*p_E_T) #** <- compute the daily mortality rate from probability of Egg survival due to temperature and dessication*
  
mu0L=-log(p_L_R*p_L_T) #** <- compute the daily mortality rate from probability of Larvae survival due to temperature and dessication*

#** set differential mortality due to density dependent according to Prior median in Table 1, White et al 2011*
gamma=rep(13.06, numDay) #** <- insert relative density dependent morality value here *

muP=-log(p_P_R*p_P_T)  #** <- compute the daily mortality rate from probability of Pupae survival due to temperature and dessication*

muM=-log(p_M_T) #** <- compute the daily mortality rate. See  Table 2 exercise handout for probability of adult survival due to temperature*


################################### Oviposition #####################################

#** To calculate the daily number of eggs laid per adult, you will need *
#** emax, GP (the length of the gonotrophy cycle), and muM (daily mortality rate)*
#** use the equation 3 in White el 2011 to derive a temperature-regulated daily eggs laid*
#** for this, use the temperature-regulated GP (see  Table 2 exercise handout) and temperature-dependent muM calculated above for adults*

emax=93.6 #** <- insert maximum number of eggs mosquitoes can lay during an oviposition. See table 1 prior median *

GP=1/pmin(1,pmax(0,0.017*Temp-0.165))  #** <-insert gonotrophy duration due to temperature. see  Table 2 exercise handout for formulation*

beta=emax*(muM/(exp(GP*muM)-1)) #** calculate beta, the daily number of eggs laid per adult using equation 3 in white et al 2011*
################################################################################


######################## Transmission parameters ###############################
bE=pmin(1,pmax(0,0.017*Temp-0.165))   # shapiro et al 2017, see supplement
bI=pmin(1,pmax(0,(0.000112*Temp*(Temp-15.384)*sqrt(35-Temp)))) # Mordecai et al 2013  and Understanding the link between malaria risk and climate


# assume a force of infection from human populaton that is constant= 0.2
fIH=rep(0.20,numDay)

# set the initial state conditions
#** use the following initial state conditions *
Mpop=1e5 # initial adult population
EM0=1e4 # exposed mosquitoes at time t = 0
IM0=0 # infected & infectious mosquitoes at time t = 0
SM0=Mpop-IM0-EM0 # susceptible mosquitoes at time t = 0
E0=0 # initial egg population
L0=P0=0 # initial larvae and pupae population

# the time sequence over which to simulate
times=seq(0,numDay, by=1)

# initialize the model
init_states=c(E=E0, L=L0,P=P0,SM=SM0,EM=EM0,IM=IM0)

# run the model. simulate mosquito population
sim_out=ode(y=init_states, times=times, func = mpop, parms = NULL,method=NULL)
plot(sim_out)
#**Get the EIR number of infectious bites per human per day EIR=maZ*
m=rowSums(sim_out[,c("SM","EM","IM")])/(Mpop/5) # mosquito density per human. Here assume mosquitoes are 5 times the number of humans from t=0 
a=bE[1:numDay] # mosquito biting rate, i.e. the mosquito contact rate
Z=sim_out[,"IM"]/rowSums(sim_out[,c("SM","EM","IM")]) #Proportion of mosquitoes infectious
EIR=m[-(numDay+1)]*a*Z[-(numDay+1)] 

res_df[[v]]=data.frame(country=country,scenario_name=scenario_name,day=days, year=years,
                       E=sim_out[,"E"][-(numDay+1)],L=sim_out[,"L"][-(numDay+1)],P=sim_out[,"P"][-(numDay+1)],
                       SM=sim_out[,"SM"][-(numDay+1)],EM=sim_out[,"EM"][-(numDay+1)],IM=sim_out[,"IM"][-(numDay+1)],
                       EIR=EIR,
                       oviposition=beta,sporogony=1/bI,gonotrophy_length=1/bE,
                       egg_duration=d_E,larvae_duration=d_L,pupae_duration=d_P,
                       egg_survival=exp(-mu0E),larvae_survival=exp(-mu0L),
                       pupae_survival=exp(-muP),
                       adult_survival=exp(-muM))

},silent = T)


do.call(rbind, res_df)
################################################################################
#                                                                              #
#                              presentation questions                          #
#                                                                              #
################################################################################

#**Q1: Describe the current climate conditions of your country*
#**Use the plot of the annual trend and seasonality of temperature and rainfall to describe this.*
#**Use the 1988-2023 observed_Temp and observed_ppt as your temperature and rainfall variables from the climate dataset to do this.*

#write.csv(met_df, "/home/guest/Desktop/Malaria Data Science/Advanced Modelling/met_df.csv")




## 
filtered_df <- met_df %>% filter(year >= 1988 & year <= 2023)
print(head(filtered_df))


# Check the max year
max(met_df$year)

# Annual Temperature for the year 1988 -2023 
m1 <- met_df %>% filter(year >= 1988 & year <= 2023)%>%
  group_by(year)%>% summarise(annual_trend = mean(observed_Temp))%>%
  ggplot( aes(x = year, y = annual_trend)) +
  geom_line() +  # Add line plot
  geom_point() +  # Add points
  labs(title = "Annual Temperature Over Years",
       x = "Year",
       y = "Temperature") +
  theme_minimal()

##  Daily Temperature for the year 1988 - 2023

m2 <- met_df %>% filter(year >= 1988 & year <= 2023)%>%
  group_by(day)%>% summarise(daily_trend = mean(observed_Temp))%>%
  ggplot( aes(x = day, y = daily_trend)) +
  geom_line() +  # Add line plot
  geom_point() +  # Add points
  labs(title = "Observed Temperature Over Years",
       x = "Year",
       y = "Temperature") +
  theme_minimal()

# Annual Rainfall for the year 1988 - 2023 

m3 <- met_df %>% filter(year >= 1988 & year <= 2023)%>%
  group_by(year)%>% summarise(annual_trend = mean(observed_ppt))%>%
  ggplot( aes(x = year, y = annual_trend)) +
  geom_line() +  # Add line plot
  geom_point() +  # Add points
  labs(title = "Annual Rainfall Over Years",
       x = "Year",
       y = "Rainfall") +
  theme_minimal()


#################

m4 <- met_df %>% filter(year >= 1988 & year <= 2023)%>%
  group_by(day)%>% summarise(daily_trend = mean(observed_Temp))%>%
  ggplot( aes(x = day, y = daily_trend)) +
  geom_line() +  # Add line plot
  geom_point() +  # Add points
  labs(title = "Seasonal Temperature Trend Over Years 1988 to 2023",
       x = "Year",
       y = "Temperature") +
  theme_minimal()


#######################

m5 <- met_df %>% filter(year >= 1988 & year <= 2023)%>%
   group_by(day)%>% summarise(daily_trend = mean(observed_ppt))%>%
  ggplot( aes(x = day, y = daily_trend)) +
  geom_line() +  # Add line plot
  geom_point() +  # Add points
  labs(title = "Seasonal Rainfall Trend Over Years 1988 to 2023",
       x = "Year",
       y = "Temperature") +
  theme_minimal()

grid.arrange(m1, m2, m3, m4, m5, nrow=3, ncol=3)


#**Q2: According to the climate data, what is the motivation for investigating how climate change may alter the risk of malaria transmission?*
#**To help you answer this, compare the recent annual trend of temperature and rainfall from 2018-2023 (taken as recent) to the climate average of temperature and rainfall (based on data from 1988-2023)*
#**Also, compare the seasonality of rainfall and temperature from 2018-2023 to the daily climatology of rainfall & temperature (based on data from 1988-2023)*
#**Discuss the differences between the annual trends & climate averages and also the differences between recent seasonality vs monthly climatology*
#**why would a change in the averages and seasonality be a concern to the risk of malaria transmission*


# Filter for the DRC data
met_df <- all_met_df %>% filter(country == "drc")

# Define the climate period (1988-2023) and the recent period (2018-2023)
df_climate <- met_df %>% filter(year >= 1988, year <= 2023)
df_recent  <- met_df %>% filter(year >= 2018, year <= 2023)

# ---------------------------
# Annual Means Calculation
# ---------------------------
# Calculate annual mean for temperature and precipitation for the recent period (2018-2023)
annual_recent <- df_recent %>% 
  group_by(year) %>% 
  summarise(
    observed_Temp = mean(observed_Temp, na.rm = TRUE),
    observed_ppt = mean(observed_ppt, na.rm = TRUE)
  )

# Calculate overall climate averages (1988-2023)
climate_avg_temp <- mean(df_climate$observed_Temp, na.rm = TRUE)
climate_avg_ppt  <- mean(df_climate$observed_ppt, na.rm = TRUE)

# ---------------------------
# Seasonality Calculation
# ---------------------------
# Here, we assume that the 'day' column represents either the day-of-year or a monthly index.
# Calculate seasonal averages (by day) for the recent period
seasonal_recent <- df_recent %>% 
  group_by(day) %>% 
  summarise(
    observed_Temp = mean(observed_Temp, na.rm = TRUE),
    observed_ppt = mean(observed_ppt, na.rm = TRUE)
  )

# Calculate seasonal climatology (daily averages for 1988-2023)
seasonal_climatology <- df_climate %>% 
  group_by(day) %>% 
  summarise(
    observed_Temp = mean(observed_Temp, na.rm = TRUE),
    observed_ppt = mean(observed_ppt, na.rm = TRUE)
  )

# Combine all datasets into one
rainfall_temp_trend <- met_df %>%
  filter(year >= 1988 & year <= 2023) %>%
  group_by(year) %>%
  summarise(
    annual_temp = mean(observed_Temp, na.rm = TRUE),
    annual_ppt = mean(observed_ppt, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c(annual_temp, annual_ppt),
               names_to = "variable",
               values_to = "annual_trend") %>%
  mutate(period = ifelse(year >= 2018, "2018-2023", "1988-2017"),  # Changed split point
         variable = recode(variable,
                           "annual_temp" = "Temperature",
                           "annual_ppt" = "Rainfall"))

# Plot both Temperature and Rainfall on one panel
p1 <- ggplot(rainfall_temp_trend, aes(x = year, y = annual_trend, color = period, group = 1)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ variable, scales = "free_y") +  # Separate panels for Temperature & Rainfall
  labs(title = "Annual Temperature and Rainfall Trends (1988-2017 & 2018-2023)",
       x = "Year",
       y = "Annual Trend",
       color = "Period") +
  theme_minimal() +
  scale_color_manual(values = c("1988-2017" = "blue", "2018-2023" = "red")) + 
  theme(strip.text = element_text(size = 14, face = "bold"))# Assign colors manually



# Seasonality Plot - Temperature with Legend
p3 <- ggplot() +
  geom_line(data = seasonal_climatology, 
            aes(x = day, y = observed_Temp, color = "1988-2023"), 
            linetype = "dashed", size = 1) +
  geom_line(data = seasonal_recent, 
            aes(x = day, y = observed_Temp, color = "2018-2023"), 
            size = 1) +
  scale_color_manual(name = "Data Source", 
                     values = c("1988-2023" = "red", "2018-2023" = "blue")) +
  labs(title = "Seasonal Cycle of Temperature",
       x = "Day (or Month Index)",
       y = "Temperature (°C)") +
  theme_minimal()

# Seasonality Plot - Precipitation with Legend
p4 <- ggplot() +
  geom_line(data = seasonal_climatology, 
            aes(x = day, y = observed_ppt, color = "1988-2023"), 
            linetype = "dashed", size = 1) +
  geom_line(data = seasonal_recent, 
            aes(x = day, y = observed_ppt, color = "2018-2023"), 
            size = 1) +
  scale_color_manual(name = "Data Source", 
                     values = c("1988-2023" = "red", "2018-2023" = "darkgreen")) +
  labs(title = "Seasonal Cycle of Precipitation",
       x = "Day (or Month Index)",
       y = "Precipitation") +
  theme_minimal()

# Display the plots
#print(p1)
#print(p2)
print(p3)
print(p4)

grid.arrange(p1,p3,p4, ncol=2, nrow=2)





#**Q3: Build out your coupled model of mosquito population and malaria transmission. And describe the system*


$$
M = S_M + E_M + I_M
$$
  
$$
\begin{aligned}
\frac{dE}{dt} &= \beta M - \frac{E}{d_E} - \mu_{0E}\left(1+\frac{E+L}{K}\right)E, \\\\[8pt]
\frac{dL}{dt} &= \frac{E}{d_E} - \frac{L}{d_L} - \mu_{0L}\left(1+\gamma \frac{E+L}{K}\right)L, \\\\[8pt]
\frac{dP}{dt} &= \frac{L}{d_L} - \frac{P}{d_P} - \mu_P P, \\\\[8pt]
\frac{dS_M}{dt} &= \frac{1}{2}\frac{P}{d_P} - \mu_M S_M - \left(b_E \, f_{IH}\right) S_M, \\\\[8pt]
\frac{dE_M}{dt} &= \left(b_E \, f_{IH}\right) S_M - \mu_M E_M - b_I E_M, \\\\[8pt]
\frac{dI_M}{dt} &= b_I E_M - \mu_M I_M.
\end{aligned}
$$
  


  \noindent
\textbf{Model Description:} \\
This model describes the dynamics of mosquito populations and malaria transmission using a system of coupled ordinary differential equations (ODEs). The system is divided into two main parts:
  
  $$
  \textbf{Immature Stages:}
$$
  
  $$
  \begin{array}{rcl}
\textbf{Eggs (E):} & \quad & \text{Eggs are produced by the adult mosquito population at a rate } \beta, \\
&      & \text{and then develop into larvae over a period } d_E, \\
&      & \text{while experiencing mortality } \mu_{0E} \text{ that also accounts for density-dependent effects via the carrying capacity } K. \\
\end{array}
$$
  
  $$
  \begin{array}{rcl}
\textbf{Larvae (L):} & \quad & \text{Eggs hatch into larvae, which then develop into pupae over a period } d_L. \\
&      & \text{Their loss is governed by the developmental process and a density-dependent mortality } \mu_{0L}, \\
&      & \text{influenced by a factor } \gamma \text{ related to the overall immature population.} \\
\end{array}
$$
  
  $$
  \begin{array}{rcl}
\textbf{Pupae (P):} & \quad & \text{Larvae transition into pupae, which then develop into adults over a period } d_P, \\
&      & \text{while incurring mortality at a rate } \mu_P. \\
\end{array}
$$
  
  $$
  \textbf{Adult Mosquito Dynamics and Malaria Transmission:}
$$
  
  $$
  \begin{array}{rcl}
\textbf{Susceptible Adults (S\_M):} & \quad & \text{A fraction of the emerging adults from the pupal stage enter the susceptible compartment.} \\
&      & \text{They are subject to natural mortality at a rate } \mu_M \text{ and may become infected upon biting infectious humans.} \\
\end{array}
$$
  
  $$
  \begin{array}{rcl}
\textbf{Exposed Adults (E\_M):} & \quad & \text{When susceptible mosquitoes bite infectious humans (at a rate determined by } b_E \text{ and a force of infection } f_{IH}), \\
  &      & \text{they move to the exposed class. They remain here until the malaria parasite develops within them.} \\
  \end{array}
  $$
    
    $$
    \begin{array}{rcl}
  \textbf{Infectious Adults (I\_M):} & \quad & \text{Exposed mosquitoes become infectious after the extrinsic incubation period,} \\
  &      & \text{governed by the rate } b_I. \text{ These infectious mosquitoes contribute to malaria transmission.} \\
  \end{array}
  $$





#**Q4: Describe the mosquito lifecycle and malaria transmission parameters based on temperature and moisture regulation entomology in your country*
#**Be sure to plot & discuss how your local climate conditions of temperature and rainfall dictate the oviposition rate, the sporogonic cycle, the gonotrophic cycle, the development time taken for egg, larvae and pupae to develop, and the survival of egg, larvae and pupae and adults*

    library(gridExtra)
  clim_df <- res_df[[which(scen_df$scenario_name == "climatology")]] #* <- select the climatology scenario from the results
  
  a1 <- ggplot(clim_df %>% group_by(day)%>% summarise(oviposition = mean(oviposition)), aes(x=day, y=oviposition)) +
    geom_line() +
    labs(title="Daily Oviposition Rate", y="Eggs per mosquito", x="Days of year") +
    theme_minimal()
  
  a2 <- ggplot(clim_df %>% group_by(day)%>% summarise(sporogony = mean(sporogony)), aes(x=day, y=sporogony)) +
    geom_line() +
    labs(title="Sporogonic Cycle Duration", y="Sporogonic Cycle", x="Days of year") +
    theme_minimal()
  
  a3 <- ggplot(clim_df %>% group_by(day)%>% summarise(gonotrophy_length = mean(gonotrophy_length)), aes(x=day, y=gonotrophy_length)) +
    geom_line() +
    labs(title="Gonotrophic Cycle Duration", y="Gonotrophic Cycle", x="Days of year") +
    theme_minimal()
  
  a4 <- ggplot(clim_df %>% group_by(day)%>% summarise(egg_duration = mean(egg_duration)), aes(x=day, y=egg_duration)) +
    geom_line() +
    labs(title="Egg Duration", y="Egg development duration ", x="Days of year") +
    theme_minimal()
  
  grid.arrange(a1,a2,a3,a4, ncol=2, nrow=2)
  
  a5 <- ggplot(clim_df %>% group_by(day)%>% summarise(larve_duration = mean(larvae_duration)), aes(x=day, y=larve_duration)) +
    geom_line() +
    labs(title="Larve Development DDuration", y="Larve development duration", x="Days of year") +
    theme_minimal()
  
  
  a6 <- ggplot(clim_df %>% group_by(day)%>% summarise(pupae_duration = mean(pupae_duration)), aes(x=day, y=pupae_duration)) +
    geom_line() +
    labs(title="Pupae Development Duration", y="Pupae development duration", x="Days of year") +
    theme_minimal()
  
  a7 <- ggplot(clim_df %>% group_by(day)%>% summarise(egg_survival = mean(egg_survival)), aes(x=day, y=egg_survival)) +
    geom_line() +
    labs(title="Egg Survival Rate", y="Survival Probability of the eggs", x="Days of year") +
    theme_minimal()
  
  a8 <- ggplot(clim_df %>% group_by(day)%>% summarise(larvae_survival = mean(larvae_survival)), aes(x=day, y=larvae_survival)) +
    geom_line() +
    labs(title="Larvae Survival Rate", y="Survival Probability of the larve", x="Days of year") +
    theme_minimal()
  
  grid.arrange(a5,a6,a7,a8, ncol=2, nrow=2)
  a9 <- ggplot(clim_df %>% group_by(day)%>% summarise(pupae_survival = mean(pupae_survival)), aes(x=day, y=pupae_survival)) +
    geom_line() +
    labs(title="Pupae Survival Rate", y="Survival Probability of the pupae", x="Days of year") +
    theme_minimal()
  
  a10 <- ggplot(clim_df %>% group_by(day)%>% summarise(adult_survival = mean(adult_survival)), aes(x=day, y=adult_survival)) +
    geom_line() +
    labs(title="Adult Mosquiotes Survival Rate  ", y="Survival Probability of the adults", x="Days of year") +
    theme_minimal()
  
  grid.arrange(a9,a10, ncol=2, nrow=1)
  
# use model simulations to answer the questions below
#**Q5:According to the daily climatological conditions (calculated based on 1988-2023 data), what is the expected annual EIR?*
#**To do answer this, run the model using the climatology_Temp and climatology_ppt as inputs for your Temp and ppt variables.*
#**Now Let the model run from 1988 to 2023, then calculate the annual EIR of the last 10 years (2014 to 2023) and then take the average of the annual EIRs over these 10 years.*
#**Use this as the expected annual EIR according to the climatology*


  # Filter climatology data for 1988-2023:
 # df_clim <- met_df %>% filter(year >= 1988, year <= 2023)
  
  # Use climatological temperature and rainfall:
 # Temp <- unlist(df_clim$climatology_Temp)
  #ppt  <- unlist(df_clim$climatology_ppt)
  #days <- unlist(df_clim$day)
  #years <- unlist(df_clim$year)
  #numDay <- length(ppt)
  
  # Set initial state conditions
  Mpop <- 1e5       # initial adult mosquito population
  EM0  <- 1e4       # initial exposed mosquitoes
  IM0  <- 0         # initial infectious mosquitoes
  SM0  <- Mpop - EM0 - IM0  # initial susceptible adults
  E0   <- 0        # initial eggs
  L0   <- 0        # initial larvae
  P0   <- 0        # initial pupae
  
  init_states <- c(E = E0, L = L0, P = P0, SM = SM0, EM = EM0, IM = IM0)
  times <- seq(0, numDay, by = 1)
  
  # Run the model simulation using your ODE system (mpop)
  sim_out <- ode(y = init_states, times = times, func = mpop, parms = NULL, method = "lsoda")
  
  # Calculate total adult mosquito population from simulation:
  total_adults <- rowSums(sim_out[, c("SM", "EM", "IM")])
  # Mosquito density per human (assume 5 mosquitoes per human initially)
  m <- total_adults / (Mpop / 5)
  # Biting rate (a) from precomputed beta_vec (assumed to be computed using climatology data)
  a_rate <- bE[1:numDay]
  # Proportion of infectious mosquitoes:
  Z <- sim_out[,"IM"] / total_adults
  
  #m=rowSums(sim_out[,c("SM","EM","IM")])/(Mpop/5) # mosquito density per human. Here assume mosquitoes are 5 times the number of humans from t=0 
  a=bE[1:numDay] # mosquito biting rate, i.e. the mosquito contact rate
  Z=sim_out[,"IM"]/rowSums(sim_out[,c("SM","EM","IM")]) #Proportion of mosquitoes infectious
  EIR=m[-(numDay+1)]*a*Z[-(numDay+1)] 
  
  
  
  
  
  # Compute daily EIR (adjust dimensions by dropping last time point):
  EIR_daily <- m[-(numDay + 1)] * a_rate * Z[-(numDay + 1)]
  
  # Create a data frame with day, year, and daily EIR:
  sim_df <- data.frame(
    day = days[-(numDay + 1)],
    year = years[-(numDay + 1)],
    EIR = EIR_daily
  )
  
  # Aggregate daily EIR to obtain annual EIR (using the mean daily EIR per year)
  annual_EIR <- sim_df %>%
    group_by(year) %>%
    summarise(annual_EIR = mean(EIR, na.rm = TRUE))
  
  
  sim_df <- data.frame(
    day = days[-(numDay + 1)],
    year = years[-(numDay + 1)],
    EIR = EIR_daily
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Filter climatology data for 1988-2023:
  #df_clim <- met_df %>% filter(year >= 1988, year <= 2023)
  df_clim <- res_df[[which(scen_df$scenario_name == "climatology")]]
  
  # Aggregate daily EIR to obtain annual EIR (using the mean daily EIR per year)
  annual_EIR <- df_clim %>%
    group_by(year) %>%
    summarise(annual_EIR = sum(EIR, na.rm = TRUE))
  
   #Focus on the last 10 years (2014-2023)
  annual_EIR_recent <- annual_EIR %>% filter(year >= 2014, year <= 2023)
  expected_annual_EIR <- mean(annual_EIR_recent$annual_EIR, na.rm = TRUE)
  
  print(annual_EIR_recent)
  print(paste("Expected Annual EIR (2014-2023):", round(expected_annual_EIR, 2)))
  

  ##############################  
  
  eir_df <- df_clim %>% filter(year >= 2014 & year <= 2023)
  
  annual_EIR <- eir_df %>% group_by(year) %>%
    summarise(total_EIR = sum(EIR, na.rm=TRUE))
  
  mean_EIR_2014_2023 <- mean(annual_EIR$total_EIR)
  print(mean_EIR_2014_2023)
  
  
#**Q6:Plot & examine the average daily EIR (i.e., daily seasonality) calculated from these last 10 years.*
#**What is the average and the range of number of infectious bites per day (the EIR) one would expect to receive in a season?*
#**Also according to the climatological EIR, when is transmission most/least likely to happen*

  
  sim_df <- data.frame(
    day = days[-(numDay + 1)],
    year = years[-(numDay + 1)],
    EIR = EIR_daily
  )
  
  
  sim_df_recent <- sim_df %>% filter(year >= 2014, year <= 2023)
  
  #df_clim <- met_df %>% filter(year >= 1988, year <= 2023)
  
  # Calculate average daily EIR (i.e., seasonal cycle) over the last 10 years
  daily_seasonality <- sim_df_recent %>%
    group_by(day) %>%
    summarise(
      mean_EIR = mean(EIR, na.rm = TRUE),
      min_EIR = min(EIR, na.rm = TRUE),
      max_EIR = max(EIR, na.rm = TRUE)
    )
  
  # Plot the seasonal cycle: Mean daily EIR with range (min-max)
  p <- ggplot(daily_seasonality, aes(x = day, y = mean_EIR)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = min_EIR, ymax = max_EIR), fill = "blue", alpha = 0.3) +
    labs(title = "Average Daily EIR (Seasonality) for 2014-2023",
         x = "Day of Year",
         y = "EIR (infectious bites per person per day)") +
    theme_minimal()
  
  print(p)
  
  # Calculate overall average and range (min and max) of daily EIR across the season
  overall_avg_EIR <- mean(daily_seasonality$mean_EIR, na.rm = TRUE)
  overall_min_EIR <- min(daily_seasonality$min_EIR, na.rm = TRUE)
  overall_max_EIR <- max(daily_seasonality$max_EIR, na.rm = TRUE)
  
  cat("Overall Average Daily EIR:", round(overall_avg_EIR, 2), "\n")
  cat("Range of Daily EIR: [", round(overall_min_EIR, 2), ",", round(overall_max_EIR, 2), "]\n")
  
  
  ####################### Another plot
  
  eir_df$day_of_year <-df_clim$day[df_clim$year >= 2014 & df_clim$year <= 2023]
  
  daily_EIR <- eir_df %>% group_by(day) %>%
    summarise(mean_EIR = mean(EIR, na.rm = TRUE))
  
  ggplot(daily_EIR, aes(x=day, y=mean_EIR)) +
    geom_line() +
    labs(title="Seasonal Daily EIR (Climatology 2014–2023)",
         x="Day of Year", y="EIR (Infectious bites/day)") +
    theme_minimal()
  
  
  
  
# use future climate simulations to answer the questions below
#**Q7:Now run the model using the future climate scenario dy_Temp_1.9(for temperature), and dy_ppt_1.9(for rainfall). Let the model run from 1988 to 2035*
#**Compute and plot the average annual EIR over the years 2030-2035, compare it to the expected annual EIR from Q5.Discuss any differences in risk*
#**Plot & examine the average daily EIR(i.e., daily seasonality) under this future climate (2030-2035). Plot & discuss the seasonality in the EIR compared to the climatology EIR from Q6.*
#**What is the average and the range of infectious bites one would expect to receive in a season?*
#**Also under this future scenario, when is transmission most/least likely to happen*

  ## Compute and plot the average annual EIR over the years 2030-2035 under Scenario 1.9
  #Extract the data for each scenario
  
  
  df_climatology <- res_df[[1]]
  df_scen_1.9 <- res_df[[2]]
  df_scen_4.5 <- res_df[[3]]
  df_scen_8.5 <- res_df[[4]]
  
  # Filter for years 2030-2035
  filtered_data <- df_scen_1.9 %>%
    filter(year >= 2030 & year <= 2035)
  
  # Compute the average annual EIR
  annual_EIR <- filtered_data %>%
    group_by(year) %>%
    summarise(mean_EIR = sum(EIR, na.rm = TRUE))
  
  annual_EIR$difference <- annual_EIR$mean_EIR - 3256.65
  # Plot the results
  q1 <- ggplot(annual_EIR, aes(x = year, y = difference)) +
    geom_line(color = "blue") +
    geom_point(size = 3, color = "red") +
    labs(title = "Average Annual EIR (2030-2035)",
         x = "Year",
         y = "Difference in EIRfrom the baseline") +
    theme_minimal()
  
  
  ## Averag daily plot
  # Extract the data for the future scenario (2030-2035)
  df_future <-df_scen_1.9 %>%
    filter(year >= 2030 & year <= 2035)
  
  # Compute the average daily EIR over 2030-2035 for future scenario
  daily_EIR_1.9 <- df_future %>%
    group_by(day) %>%
    summarise(mean_EIR = mean(EIR, na.rm = TRUE))
  
  ## Comparison
  q2 <- ggplot() +
    geom_line(data = daily_EIR, aes(x = day, y = mean_EIR, color = "Climatology"), size = 1, linetype = "dashed") +
    geom_line(data = daily_EIR_1.9, aes(x = day, y = mean_EIR, color = "Scenario 1.9"), size = 1) +
    labs(title = "Comparison of Daily EIR Seasonality (Climatology vs. Scenario 1.9)",
         x = "Day of Year",
         y = "Mean Daily EIR",
         color = "Scenario") +  # Legend title
    scale_color_manual(values = c("Climatology" = "blue", "Scenario 1.9" = "red")) +
    theme_minimal()
  
  
  # Compute the range of infectious bites (EIR) in a season (min & max)
  eir_range_1.9 <- tibble(
    min = min(daily_EIR_1.9$mean_EIR, na.rm = TRUE),
    mean = mean(daily_EIR_1.9$mean_EIR, na.rm = TRUE),
    max = max(daily_EIR_1.9$mean_EIR, na.rm = TRUE)
  )
  
  print(eir_range_1.9)
  
  
  # Identify peak & lowest transmission periods
  peak_day <- daily_EIR_1.9 %>% filter(mean_EIR == max(daily_EIR_1.9$mean_EIR)) %>% select(day)
  low_day <- daily_EIR_1.9 %>% filter(mean_EIR == min(daily_EIR_1.9$mean_EIR)) %>% select(day)
  
  print(paste("Peak Transmission Days:", peak_day$day))
  print(paste("Lowest Transmission Days:", low_day$day))
  
  grid.arrange(q1, q2, nrow=1, ncol=2)
  
#**Q8:Run the model now, but using the future climate scenario dy_Temp_4.5, and dy_ppt_4.5, let the model run from 1988 to 2035*
#**Compute and plot the average annual EIR over the years 2030-2035, compare it to the expected annual EIR from Q5.Discuss any differences in risk*
#**Plot & examine the average daily EIR(i.e., daily seasonality) under this future climate (2030-2035). Plot & discuss the seasonality in the EIR compared to the climatology EIR from Q6.*
#**What is the average and the range of infectious bites one would expect to receive in a season?*
#**Also under this future scenario, when is transmission most/least likely to happen*

  # ---------------------------
  # Data Setup: Future Scenario (dy_Temp_4.5, dy_ppt_4.5, 1988-2035)
  # ---------------------------
  df_future45 <- met_df %>% filter(year >= 1988, year <= 2035)
  
  # Extract future scenario variables using dy_Temp_4.5 and dy_ppt_4.5
  Temp_future45 <- unlist(df_future45$dy_Temp_4.5)
  ppt_future45  <- unlist(df_future45$dy_ppt_4.5)
  days_future45 <- unlist(df_future45$day)
  years_future45 <- unlist(df_future45$year)
  numDay_future45 <- length(ppt_future45)
  
  # ---------------------------
  # Compute Future Parameter: beta_vec_future45
  # (Using the same coefficients as before)
  emax <- 93.6
  a_p_M_T <- -0.000828; b_p_M_T <- 0.0367; c_p_M_T <- 0.522
  
  GP_vec_future45 <- 1 / pmin(1, pmax(0, 0.017 * Temp_future45 - 0.165))
  mu_M_vec_future45 <- -log(a_p_M_T * Temp_future45^2 + b_p_M_T * Temp_future45 + c_p_M_T)
  beta_vec_future45 <- emax * (mu_M_vec_future45 / (exp(GP_vec_future45 * mu_M_vec_future45) - 1))
  
  # ---------------------------
  # Set initial conditions (same as before)
  # ---------------------------
  Mpop <- 1e5       # initial adult mosquito population
  EM0  <- 1e4       # initial exposed mosquitoes
  IM0  <- 0         # initial infectious mosquitoes
  SM0  <- Mpop - EM0 - IM0  # initial susceptible adults
  E0   <- 0        # initial eggs
  L0   <- 0        # initial larvae
  P0   <- 0        # initial pupae
  
  init_states <- c(E = E0, L = L0, P = P0, SM = SM0, EM = EM0, IM = IM0)
  times_future45 <- seq(0, numDay_future45, by = 1)
  
  # ---------------------------
  # Run the model with future climate scenario (dy_Temp_4.5, dy_ppt_4.5)
  # ---------------------------
  sim_out_future45 <- ode(y = init_states, times = times_future45, func = mpop, parms = NULL, method = "lsoda")
  
  # ---------------------------
  # Calculate Daily EIR for Future Scenario (dy_Temp_4.5)
  # ---------------------------
  total_adults_future45 <- rowSums(sim_out_future45[, c("SM", "EM", "IM")])
  m_future45 <- total_adults_future45 / (Mpop / 5)
  a_rate_future45 <- beta_vec_future45[1:numDay_future45]  # future biting rate vector
  Z_future45 <- sim_out_future45[,"IM"] / total_adults_future45
  
  EIR_daily_future45 <- m_future45[-(numDay_future45 + 1)] * a_rate_future45 * Z_future45[-(numDay_future45 + 1)]
  
  sim_df_future45 <- data.frame(
    day = days_future45[-(numDay_future45 + 1)],
    year = years_future45[-(numDay_future45 + 1)],
    EIR = EIR_daily_future45
  )
  
  # ---------------------------
  # Annual EIR Calculation for Future Scenario (2030-2035)
  # ---------------------------
  annual_EIR_future45 <- sim_df_future45 %>%
    group_by(year) %>%
    summarise(annual_EIR = mean(EIR, na.rm = TRUE))
  
  annual_EIR_future45_recent <- annual_EIR_future45 %>% filter(year >= 2030, year <= 2035)
  expected_annual_EIR_future45 <- mean(annual_EIR_future45_recent$annual_EIR, na.rm = TRUE)
  
  print(annual_EIR_future45_recent)
  print(paste("Expected Annual EIR (2030-2035) [Future Scenario dy_Temp_4.5]:", round(expected_annual_EIR_future45, 2)))
  
  # ---------------------------
  # Daily Seasonality for Future Scenario (2030-2035)
  # ---------------------------
  sim_df_future45_recent <- sim_df_future45 %>% filter(year >= 2030, year <= 2035)
  
  daily_seasonality_future45 <- sim_df_future45_recent %>%
    group_by(day) %>%
    summarise(
      mean_EIR = mean(EIR, na.rm = TRUE),
      min_EIR = min(EIR, na.rm = TRUE),
      max_EIR = max(EIR, na.rm = TRUE)
    )
  
  p_future45 <- ggplot(daily_seasonality_future45, aes(x = day, y = mean_EIR)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = min_EIR, ymax = max_EIR), fill = "blue", alpha = 0.3) +
    labs(title = "Average Daily EIR (2030-2035) - Future Scenario dy_Temp_4.5",
         x = "Day of Year",
         y = "EIR (infectious bites per person per day)") +
    theme_minimal()
  
  print(p_future45)
  
  overall_avg_EIR_future45 <- mean(daily_seasonality_future45$mean_EIR, na.rm = TRUE)
  overall_min_EIR_future45 <- min(daily_seasonality_future45$min_EIR, na.rm = TRUE)
  overall_max_EIR_future45 <- max(daily_seasonality_future45$max_EIR, na.rm = TRUE)
  
  cat("Future Scenario dy_Temp_4.5 - Overall Average Daily EIR:", round(overall_avg_EIR_future45, 2), "\n")
  cat("Future Scenario dy_Temp_4.5 - Range of Daily EIR: [", round(overall_min_EIR_future45, 2), ",", round(overall_max_EIR_future45, 2), "]\n")
  
  
  
  # Annual EIR Plot with X-Axis as Years
  annual_EIR_future45_recent <- annual_EIR_future45 %>% filter(year >= 2030, year <= 2035)
  
  p_annual_future45 <- ggplot(annual_EIR_future45_recent, aes(x = year, y = annual_EIR)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "blue", size = 3) +
    labs(title = "Average Annual EIR (2030-2035) - Future Scenario dy_Temp_4.5",
         x = "Year",
         y = "Annual EIR (infectious bites per person per year)") +
    theme_minimal()
  
  print(p_annual_future45)
  
  
  
  
  
  
  
  
  
  
  
  
  #####################################
  
  ## Compute and plot the average annual EIR over the years 2030-2035 under Scenario 1.9
  #Extract the data for each scenario
  df_climatology <- res_df[[1]]
  df_scen_1.9 <- res_df[[2]]
  df_scen_4.5 <- res_df[[3]]
  df_scen_8.5 <- res_df[[4]]
  
  # Filter for years 2030-2035
  filtered_data <- df_scen_4.5 %>%
    filter(year >= 2030 & year <= 2035)
  
  
  
  # Compute the average annual EIR
  annual_EIR <- filtered_data %>%
    group_by(year) %>%
    summarise(mean_EIR = sum(EIR, na.rm = TRUE))
  
  annual_EIR$difference <- annual_EIR$mean_EIR - 3256.65
  
  # Plot the results
  l1 <- ggplot(annual_EIR, aes(x = year, y = difference)) +
    geom_line(color = "blue") +
    geom_point(size = 3, color = "red") +
    labs(title = "Difference Annual EIR (2030-2035)",
         x = "Year",
         y = "Difference EIR from the baseline") +
    theme_minimal()
  
  
  ## Averag daily plot
  # Extract the data for the future scenario (2030-2035)
  df_future <-df_scen_4.5 %>%
    filter(year >= 2030 & year <= 2035)
  
  # Compute the average daily EIR over 2030-2035 for future scenario
  daily_EIR_4.5 <- df_future %>%
    group_by(day) %>%
    summarise(mean_EIR = mean(EIR, na.rm = TRUE))
  
  ## Comparison
  l2 <- ggplot() +
    geom_line(data = daily_EIR, aes(x = day, y = mean_EIR, color = "Climatology"), size = 1, linetype = "dashed") +
    geom_line(data = daily_EIR_4.5, aes(x = day, y = mean_EIR, color = "Scenario 4.5"), size = 1) +
    labs(title = "Comparison of Daily EIR Seasonality",
         x = "Day of Year",
         y = "Mean Daily EIR",
         color = "Scenario") +  # Legend title
    scale_color_manual(values = c("Climatology" = "blue", "Scenario 4.5" = "red")) +
    theme_minimal()
  
  
  # Compute the range of infectious bites (EIR) in a season (min & max)
  eir_range_4.5 <- tibble(
    min = min(daily_EIR_4.5$mean_EIR, na.rm = TRUE),
    mean = mean(daily_EIR_4.5$mean_EIR, na.rm = TRUE),
    max = max(daily_EIR_4.5$mean_EIR, na.rm = TRUE)
  )
  
  print(eir_range_4.5)
  
  
  # Identify peak & lowest transmission periods
  peak_day <- daily_EIR_4.5 %>% filter(mean_EIR == max(daily_EIR_4.5$mean_EIR)) %>% select(day)
  low_day <- daily_EIR_4.5 %>% filter(mean_EIR == min(daily_EIR_4.5$mean_EIR)) %>% select(day)
  
  print(paste("Peak Transmission Days:", peak_day$day))
  print(paste("Lowest Transmission Days:", low_day$day))  
  
  
  grid.arrange(l1, l2, nrow=1, ncol=2)
  
  
  
  
  
  
#**Q9:Again, run the model but using the future climate scenario dy_Temp_8.5, and dy_ppt_8.5. Let the model run from 1988 to 2035*
#**Compute & plot the average annual EIR over the years 2030-2035, compare it to the expected annual EIR from Q5.Discuss any differences in risk*
#**Plot & examine the average daily EIR(i.e., daily seasonality) under this future climate (2030-2035).Plot & discuss the seasonality in the EIR compared to the climatology EIR from Q6.*
#**What is the average and the range of infectious bites one would expect to receive in a season?*
#**Also under this future scenario, when is transmission most/least likely to happen*

  library(deSolve)
  library(dplyr)
  library(ggplot2)
  
  # ---------------------------
  # Data Setup: Future Scenario (dy_Temp_8.5, dy_ppt_8.5; 1988-2035)
  # ---------------------------
  df_future85 <- met_df %>% filter(year >= 1988, year <= 2035)
  
  Temp_future85 <- unlist(df_future85$dy_Temp_8.5)
  ppt_future85  <- unlist(df_future85$dy_ppt_8.5)
  days_future85 <- unlist(df_future85$day)
  years_future85 <- unlist(df_future85$year)
  numDay_future85 <- length(ppt_future85)
  
  times_future85 <- seq(0, numDay_future85, by = 1)
  
  # ---------------------------
  # Compute beta_vec_future85 using future scenario variables
  # (Using existing coefficients: emax, a_p_M_T, b_p_M_T, c_p_M_T)
  GP_vec_future85 <- 1 / pmin(1, pmax(0, 0.017 * Temp_future85 - 0.165))
  mu_M_vec_future85 <- -log(a_p_M_T * Temp_future85^2 + b_p_M_T * Temp_future85 + c_p_M_T)
  beta_vec_future85 <- emax * (mu_M_vec_future85 / (exp(GP_vec_future85 * mu_M_vec_future85) - 1))
  
  # ---------------------------
  # Set initial conditions (assumed same as before)
  # ---------------------------
  init_states <- c(E = 0, L = 0, P = 0, SM = (1e5 - 1e4), EM = 1e4, IM = 0)
  
  # ---------------------------
  # Run the Model Simulation with Future Scenario 8.5
  # ---------------------------
  sim_out_future85 <- ode(y = init_states, times = times_future85, func = mpop, parms = NULL, method = "lsoda")
  
  # ---------------------------
  # Calculate Daily EIR for Future Scenario 8.5
  # ---------------------------
  total_adults_future85 <- rowSums(sim_out_future85[, c("SM", "EM", "IM")])
  m_future85 <- total_adults_future85 / (1e5 / 5)
  a_rate_future85 <- beta_vec_future85[1:numDay_future85]
  Z_future85 <- sim_out_future85[,"IM"] / total_adults_future85
  EIR_daily_future85 <- m_future85[-(numDay_future85 + 1)] * a_rate_future85 * Z_future85[-(numDay_future85 + 1)]
  
  sim_df_future85 <- data.frame(
    day = days_future85[-(numDay_future85 + 1)],
    year = years_future85[-(numDay_future85 + 1)],
    EIR = EIR_daily_future85
  )
  
  # ---------------------------
  # Annual EIR Calculation for 2030-2035 (Future Scenario 8.5)
  # ---------------------------
  annual_EIR_future85 <- sim_df_future85 %>%
    group_by(year) %>%
    summarise(annual_EIR = mean(EIR, na.rm = TRUE))
  
  annual_EIR_future85_recent <- annual_EIR_future85 %>% filter(year >= 2030, year <= 2035)
  expected_annual_EIR_future85 <- mean(annual_EIR_future85_recent$annual_EIR, na.rm = TRUE)
  
  print(annual_EIR_future85_recent)
  print(paste("Expected Annual EIR (2030-2035) [Future Scenario 8.5]:", round(expected_annual_EIR_future85, 2)))
  
  # ---------------------------
  # Daily Seasonality for 2030-2035 (Future Scenario 8.5)
  # ---------------------------
  sim_df_future85_recent <- sim_df_future85 %>% filter(year >= 2030, year <= 2035)
  daily_seasonality_future85 <- sim_df_future85_recent %>%
    group_by(day) %>%
    summarise(
      mean_EIR = mean(EIR, na.rm = TRUE),
      min_EIR = min(EIR, na.rm = TRUE),
      max_EIR = max(EIR, na.rm = TRUE)
    )
  
  p_future85 <- ggplot(daily_seasonality_future85, aes(x = day, y = mean_EIR)) +
    geom_line(color = "red", size = 1) +
    geom_ribbon(aes(ymin = min_EIR, ymax = max_EIR), fill = "red", alpha = 0.3) +
    labs(title = "Average Daily EIR (2030-2035) - Future Scenario 8.5",
         x = "Day of Year",
         y = "EIR (infectious bites per person per day)") +
    theme_minimal()
  
  print(p_future85)
  
  overall_avg_EIR_future85 <- mean(daily_seasonality_future85$mean_EIR, na.rm = TRUE)
  overall_min_EIR_future85 <- min(daily_seasonality_future85$min_EIR, na.rm = TRUE)
  overall_max_EIR_future85 <- max(daily_seasonality_future85$max_EIR, na.rm = TRUE)
  
  cat("Future Scenario 8.5 - Overall Average Daily EIR:", round(overall_avg_EIR_future85, 2), "\n")
  cat("Future Scenario 8.5 - Range of Daily EIR: [", round(overall_min_EIR_future85, 2), ",", round(overall_max_EIR_future85, 2), "]\n")
  
  
  
  
  
  
  
  
  # ---------------------------
  # Annual EIR Plot with X-Axis in Years
  # ---------------------------
  annual_EIR_future85_recent <- annual_EIR_future85 %>% filter(year >= 2030, year <= 2035)
  
  p_annual <- ggplot(annual_EIR_future85_recent, aes(x = year, y = annual_EIR)) +
    geom_line(color = "red", size = 1) +
    geom_point(color = "red", size = 3) +
    labs(title = "Average Annual EIR (2030-2035) - Future Scenario 8.5",
         x = "Year",
         y = "Annual EIR (infectious bites per person per year)") +
    theme_minimal()
  
  print(p_annual)
  
  
  
  
  
  
  
  
  #####################################
  
  ## Compute and plot the average annual EIR over the years 2030-2035 under Scenario 1.9
  #Extract the data for each scenario
  df_climatology <- res_df[[1]]
  df_scen_1.9 <- res_df[[2]]
  df_scen_4.5 <- res_df[[3]]
  df_scen_8.5 <- res_df[[4]]
  
  # Filter for years 2030-2035
  filtered_data <- df_scen_8.5 %>%
    filter(year >= 2030 & year <= 2035)
  
  # Compute the average annual EIR
  annual_EIR <- filtered_data %>%
    group_by(year) %>%
    summarise(mean_EIR = sum(EIR, na.rm = TRUE))
  
  annual_EIR$difference <- annual_EIR$mean_EIR - 3256.65
  
  # Plot the results
  N1 <- ggplot(annual_EIR, aes(x = year, y = difference)) +
    geom_line(color = "blue") +
    geom_point(size = 3, color = "red") +
    labs(title = "Difference in Annual EIR (2030-2035) under df_scen_8.5 vs baseline",
         x = "Year",
         y = "Difference in EIR from the baseline") +
    theme_minimal()
  
  
  ## Averag daily plot
  # Extract the data for the future scenario (2030-2035)
  df_future <-df_scen_8.5 %>%
    filter(year >= 2030 & year <= 2035)
  
  # Compute the average daily EIR over 2030-2035 for future scenario
  daily_EIR_8.5 <- df_future %>%
    group_by(day) %>%
    summarise(mean_EIR = mean(EIR, na.rm = TRUE))
  
  ## Comparison
 N2 <-  ggplot() +
    geom_line(data = daily_EIR, aes(x = day, y = mean_EIR, color = "Climatology"), size = 1, linetype = "dashed") +
    geom_line(data = daily_EIR_8.5, aes(x = day, y = mean_EIR, color = "Scenario 8.5"), size = 1) +
    labs(title = "Comparison of Daily EIR Seasonality (Climatology vs. Scenario 8.5)",
         x = "Day of Year",
         y = "Mean Daily EIR",
         color = "Scenario") +  # Legend title
    scale_color_manual(values = c("Climatology" = "blue", "Scenario 4.5" = "red")) +
    theme_minimal()
  
  
  # Compute the range of infectious bites (EIR) in a season (min & max)
  eir_range_4.5 <- tibble(
    min = min(daily_EIR_8.5$mean_EIR, na.rm = TRUE),
    mean = mean(daily_EIR_8.5$mean_EIR, na.rm = TRUE),
    max = max(daily_EIR_8.5$mean_EIR, na.rm = TRUE)
  )
  
  print(eir_range_4.5)
  
  
  # Identify peak & lowest transmission periods
  peak_day <- daily_EIR_8.5 %>% filter(mean_EIR == max(daily_EIR_8.5$mean_EIR)) %>% select(day)
  low_day <- daily_EIR_8.5 %>% filter(mean_EIR == min(daily_EIR_8.5$mean_EIR)) %>% select(day)
  
  print(paste("Peak Transmission Days:", peak_day$day))
  print(paste("Lowest Transmission Days:", low_day$day))  
  
  grid.arrange(N1, N2, nrow=1, ncol=2)
  
  
  
  
  
  
  
  
  
  
  
  
  
# Conclusions
  
#**Q10: Reflect over the future risk predicted under the different the three scenarios. Are the risk of transmission, as you would expect?.*
#**How similar are they across future conditions? Which scenario is more concerning for malaria transmission?*

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  