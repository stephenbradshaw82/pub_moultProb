#' Author: Stephen Bradshaw
#' University of Tasmania
#' Date: Jan 2025
#' Title: SRL modelling - Moult prediction
#' Details:
#'    - Create simulation data
#'    - Establish STAN model inputs
#'    - Run STAN model to predict moult probability
#'    - [20250123] Data includes 13 simulation datasets and STAN model to save outputs
#'    
#### Package / Library Requirements ####
rm(list=ls())
options(stringsAsFactors = FALSE)
options(scipen = 999)

#' Packages (not in environment) -------------------------------------------
list.of.packages <- c("magrittr", "tidyr", "dplyr", "stringr", "purrr"
                      , "gam", "rebus", "lubridate", "modeest"
                      , "rstan", "shinystan", "stringr", "StanHeaders"
                      , "overlapping", "stringi", "cowplot"
                      , "ggmcmc")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#' Libraries ---------------------------------------------------------------
req_packages <- list.of.packages 
sapply(req_packages, require, character.only = TRUE, quietly=TRUE)
#####

#### FUNCTIONS ####
source("src/functions.R")
#####

#### Create Empty Data list ####
ldata <- list()
#####

##############################################
############ SCENARIO 1 GENERAL ##############
##############################################

#### [NEW GEN SB: "gen30lib18grow52"] As per above with Generated ERROR (and some other distribution mods) #####
if (100==100){
  #--> Initial Features ####
  set.seed(1)
  endYear<-2 #capture and recapture between time 0 to endYear
  nLobsters<-3000 #Number of lobsters
  annualSteps<-12 #How many time steps per year (12 for roughly months)
  propAnimalsExtraInfo<-0.3 #Proportion of animals with extra moult evidence
  upperLimAtLiberty <- 18 #Maximum number of months at Liberty
  
  startDates<-runif(nLobsters,0,endYear-1/annualSteps) #Uniform across period
  libertyLength<-runif(nLobsters, 1,upperLimAtLiberty)/annualSteps
  
  libertyLength[libertyLength<0]<-1/annualSteps #No negative times
  
  endDates<-startDates+libertyLength
  
  #Keep only lobsters recaptured before endYear
  select<-endDates<endYear
  
  select1k <- which(select)
  select1k <- sample(select1k, size=1000, replace=FALSE)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  #Moulting probability by time period
  moultP <- c(0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.15,0.05,0.02,0.02,0.01)
  moultP<-moultP[1:annualSteps]
  
  ##DAMAGE
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveDamageRegen <- c(0.1, 0.8)
  
  ##PLEOPOD
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObservePleopodRegen <- c(0.4, 0.99)
  
  
  ##MATURITY
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveMaturity <- c(0.2, 0.5)
  
  
  ##For model that relies on growth increments to estimate moulting probability
  #Growth increment distribution parameters
  growth_mu<-5
  growth_sigma<-2 #Includes variability in growth and measurement
  
  #Variation in size if they didn't moult (just observation error)
  nogrowth_sigma<-1
  
  
  #--> Generate Data ####
  #Contains the number of each time period that each lobster was at liberty
  atLiberty<-matrix(0,nLobsters,annualSteps)
  
  #Calculating proportion/number of steps that each lobster was at liberty
  for (iLobster in 1:nLobsters){
    started <- FALSE;
    ended<-FALSE;
    
    for (iYear in 1:endYear)
    {
      for (iStep in 1:annualSteps)
      {
        if (!started & startDates[iLobster]<((iYear-1)+(iStep)/annualSteps))
        {
          started <- TRUE
          atLiberty[iLobster,iStep]<-1-(startDates[iLobster]%%(1/annualSteps)*annualSteps)
        } else if (started & (endDates[iLobster]>((iYear-1)+(iStep)/annualSteps)))
        {
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+1 
        } else if (started & !ended)
        {
          ended<-TRUE;
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+(endDates[iLobster]%%(1/annualSteps)*annualSteps)
        }
      }
    }
  } ##end iLobster
  
  #--> Calculate moulting observations ####
  moultedP<-c() #Probability it moulted
  moultedOutcome<-c() #Did it moult?
  damagedecreaseObservation<-c()
  pleopoddecreaseObservation<-c()
  maturitychangeObservation<-c()
  growth<-c() #The growth increment
  
  for (iLobster in 1:nLobsters){
    
    moultedP[iLobster] <- 1-prod((1-moultP)^atLiberty[iLobster,])
    moultedOutcome[iLobster] <- rbinom(1,1,moultedP[iLobster])

    damagedecreaseObservation[iLobster] <- rbinom(1,1,chanceObserveDamageRegen[moultedOutcome[iLobster]+1])
    pleopoddecreaseObservation[iLobster] <- rbinom(1,1,chanceObservePleopodRegen[moultedOutcome[iLobster]+1])
    maturitychangeObservation[iLobster] <- rbinom(1,1,chanceObserveMaturity[moultedOutcome[iLobster]+1])
    
    ##calculate the growth
    if (moultedOutcome[iLobster]){
      growth[iLobster]<-rnorm(1,growth_mu,growth_sigma)
    }else{
      growth[iLobster]<-rnorm(1,0,nogrowth_sigma)
    }
  }
  
  ##Assign only some valid records (makes some data invalid as no information from those animals --> proportion of animals not-informative)
  damagedecreaseObservation[sample(1:length(damagedecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(damagedecreaseObservation)), replace=FALSE)] <- -1
  pleopoddecreaseObservation[sample(1:length(pleopoddecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(pleopoddecreaseObservation)), replace=FALSE)] <- -1
  maturitychangeObservation[sample(1:length(maturitychangeObservation), size=floor((1-propAnimalsExtraInfo)*length(maturitychangeObservation)), replace=FALSE)] <- -1 ##NEW
  
  #--> Set data list ####
  sdata <-list(
    name = "gen30lib18grow52"
    ,nLobsters=nLobsters
    ,nPeriods=annualSteps
    ,obsGI=growth
    ,obsDam = damagedecreaseObservation
    ,obsPle = pleopoddecreaseObservation
    ,obsMat = maturitychangeObservation
    ,atLiberty=atLiberty
    ,isFem = 1
  )
  
  ldata[[length(ldata)+1]] <- sdata
  
  # sdata$nLobsters
  # sdata$obsGI %>% mean()
  # sdata$obsGI %>% sd()
  # sdata$obsDam %>% table()
  # sdata$obsPle %>% table()
  # sdata$obsMat %>% table()
  
  #####
} #end 100==100
#####

##############################################
############ SCENARIO 2 (% DATA) #############
##############################################

#### [NEW GEN SB: "gen05lib18grow52"] AS ABOVE --> CHANGED to (5% info) 0.95 UNINFORMATIVE MOULTING [SCENARIO 1] #####
if (200==200){
  #--> Initial Features ####
  set.seed(1)
  endYear<-2 #capture and recapture between time 0 to endYear
  nLobsters<-3000 #Number of lobsters
  annualSteps<-12 #How many time steps per year (12 for roughly months)
  propAnimalsExtraInfo<-0.05 #Proportion of animals with extra moult evidence
  upperLimAtLiberty <- 18 #Maximum number of months at Liberty
  
  startDates<-runif(nLobsters,0,endYear-1/annualSteps) #Uniform across period
  libertyLength<-runif(nLobsters, 1,upperLimAtLiberty)/annualSteps
  
  libertyLength[libertyLength<0]<-1/annualSteps #No negative times
  
  endDates<-startDates+libertyLength
  
  #Keep only lobsters recaptured before endYear
  select<-endDates<endYear
  
  select1k <- which(select)
  select1k <- sample(select1k, size=1000, replace=FALSE)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  #Moulting probability by time period
  moultP <- c(0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.15,0.05,0.02,0.02,0.01)
  moultP<-moultP[1:annualSteps]
  
  ##DAMAGE
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveDamageRegen <- c(0.1, 0.8)
  
  ##PLEOPOD
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObservePleopodRegen <- c(0.4, 0.99)
  
  
  ##MATURITY
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveMaturity <- c(0.2, 0.5)
  
  
  ##For model that relies on growth increments to estimate moulting probability
  #Growth increment distribution parameters
  growth_mu<-5
  growth_sigma<-2 #Includes variability in growth and measurement
  
  #Variation in size if they didn't moult (just observation error)
  nogrowth_sigma<-1
  
  #--> Generate Data ####
  #Contains the number of each time period that each lobster was at liberty
  atLiberty<-matrix(0,nLobsters,annualSteps)
  
  #Calculating proportion/number of steps that each lobster was at liberty
  for (iLobster in 1:nLobsters){
    started <- FALSE;
    ended<-FALSE;
    
    for (iYear in 1:endYear)
    {
      for (iStep in 1:annualSteps)
      {
        if (!started & startDates[iLobster]<((iYear-1)+(iStep)/annualSteps))
        {
          started <- TRUE
          atLiberty[iLobster,iStep]<-1-(startDates[iLobster]%%(1/annualSteps)*annualSteps)
        } else if (started & (endDates[iLobster]>((iYear-1)+(iStep)/annualSteps)))
        {
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+1 
        } else if (started & !ended)
        {
          ended<-TRUE;
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+(endDates[iLobster]%%(1/annualSteps)*annualSteps)
        }
      }
    }
  } ##end iLobster
  
  #--> Calculate moulting observations ####
  moultedP<-c() #Probability it moulted
  moultedOutcome<-c() #Did it moult?
  damagedecreaseObservation<-c()
  pleopoddecreaseObservation<-c()
  maturitychangeObservation<-c() ##NEW
  growth<-c() #The growth increment
  
  for (iLobster in 1:nLobsters){
    
    moultedP[iLobster] <- 1-prod((1-moultP)^atLiberty[iLobster,])
    moultedOutcome[iLobster] <- rbinom(1,1,moultedP[iLobster])
    
    damagedecreaseObservation[iLobster] <- rbinom(1,1,chanceObserveDamageRegen[moultedOutcome[iLobster]+1])
    pleopoddecreaseObservation[iLobster] <- rbinom(1,1,chanceObservePleopodRegen[moultedOutcome[iLobster]+1])
    maturitychangeObservation[iLobster] <- rbinom(1,1,chanceObserveMaturity[moultedOutcome[iLobster]+1]) ##NEW
    
    ##calculate the growth
    if (moultedOutcome[iLobster]){
      growth[iLobster]<-rnorm(1,growth_mu,growth_sigma)
    }else{
      growth[iLobster]<-rnorm(1,0,nogrowth_sigma)
    }
  }
  
  ##Assign only some valid records (makes some data invalid as no information from those animals --> proportion of animals not-informative)
  damagedecreaseObservation[sample(1:length(damagedecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(damagedecreaseObservation)), replace=FALSE)] <- -1
  pleopoddecreaseObservation[sample(1:length(pleopoddecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(pleopoddecreaseObservation)), replace=FALSE)] <- -1
  maturitychangeObservation[sample(1:length(maturitychangeObservation), size=floor((1-propAnimalsExtraInfo)*length(maturitychangeObservation)), replace=FALSE)] <- -1
  
  moultedP[atLiberty[,12]>0]
  
  #--> Set data list ####
  sdata <-list(
    name = "gen05lib18grow52"
    ,nLobsters=nLobsters
    ,nPeriods=annualSteps
    ,obsGI=growth
    ,obsDam = damagedecreaseObservation
    ,obsPle = pleopoddecreaseObservation
    ,obsMat = maturitychangeObservation
    ,atLiberty=atLiberty
    ,isFem = 1
  )
  
  
  ldata[[length(ldata)+1]] <- sdata
  
  # sdata$nLobsters
  # sdata$obsGI %>% mean()
  # sdata$obsGI %>% sd()
  # sdata$obsDam %>% table()
  # sdata$obsPle %>% table()
  # sdata$obsMat %>% table()
  #####
} #end 200==200
#####

#### [NEW GEN SB: "gen50lib18grow52"] AS ABOVE--> CHANGED to (50% info) 0.5 UNINFORMATIVE MOULTING [SCENARIO 1] #####
if (210==210){
  #--> Initial Features ####
  set.seed(1)
  endYear<-2 #capture and recapture between time 0 to endYear
  nLobsters<-3000 #Number of lobsters
  annualSteps<-12 #How many time steps per year (12 for roughly months)
  propAnimalsExtraInfo<-0.50 #Proportion of animals with extra moult evidence
  upperLimAtLiberty <- 18 #Maximum number of months at Liberty
  
  startDates<-runif(nLobsters,0,endYear-1/annualSteps) #Uniform across period
  libertyLength<-runif(nLobsters, 1,upperLimAtLiberty)/annualSteps
  
  libertyLength[libertyLength<0]<-1/annualSteps #No negative times
  
  endDates<-startDates+libertyLength
  
  #Keep only lobsters recaptured before endYear
  select<-endDates<endYear
  
  select1k <- which(select)
  select1k <- sample(select1k, size=1000, replace=FALSE)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  #Moulting probability by time period
  moultP <- c(0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.15,0.05,0.02,0.02,0.01)
  moultP<-moultP[1:annualSteps]
  
  ##DAMAGE
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveDamageRegen <- c(0.1, 0.8)
  
  ##PLEOPOD
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObservePleopodRegen <- c(0.4, 0.99)
  
  ##MATURITY
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveMaturity <- c(0.2, 0.5)
  
  ##For model that relies on growth increments to estimate moulting probability
  #Growth increment distribution parameters
  growth_mu<-5
  growth_sigma<-2 #Includes variability in growth and measurement
  
  #Variation in size if they didn't moult (just observation error)
  nogrowth_sigma<-1
  
  #--> Generate Data ####
  #Contains the number of each time period that each lobster was at liberty
  atLiberty<-matrix(0,nLobsters,annualSteps)
  
  #Calculating proportion/number of steps that each lobster was at liberty
  for (iLobster in 1:nLobsters){
    started <- FALSE;
    ended<-FALSE;
    
    for (iYear in 1:endYear)
    {
      for (iStep in 1:annualSteps)
      {
        if (!started & startDates[iLobster]<((iYear-1)+(iStep)/annualSteps))
        {
          started <- TRUE
          atLiberty[iLobster,iStep]<-1-(startDates[iLobster]%%(1/annualSteps)*annualSteps)
        } else if (started & (endDates[iLobster]>((iYear-1)+(iStep)/annualSteps)))
        {
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+1 
        } else if (started & !ended)
        {
          ended<-TRUE;
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+(endDates[iLobster]%%(1/annualSteps)*annualSteps)
        }
      }
    }
  } ##end iLobster
  
  #--> Calculate moulting observations ####
  moultedP<-c() #Probability it moulted
  moultedOutcome<-c() #Did it moult?
  damagedecreaseObservation<-c()
  pleopoddecreaseObservation<-c()
  maturitychangeObservation<-c() ##NEW
  growth<-c() #The growth increment
  
  for (iLobster in 1:nLobsters){
    
    moultedP[iLobster] <- 1-prod((1-moultP)^atLiberty[iLobster,])
    moultedOutcome[iLobster] <- rbinom(1,1,moultedP[iLobster])

    damagedecreaseObservation[iLobster] <- rbinom(1,1,chanceObserveDamageRegen[moultedOutcome[iLobster]+1])
    pleopoddecreaseObservation[iLobster] <- rbinom(1,1,chanceObservePleopodRegen[moultedOutcome[iLobster]+1])
    maturitychangeObservation[iLobster] <- rbinom(1,1,chanceObserveMaturity[moultedOutcome[iLobster]+1]) ##NEW
    
    ##calculate the growth
    if (moultedOutcome[iLobster]){
      growth[iLobster]<-rnorm(1,growth_mu,growth_sigma)
    }else{
      growth[iLobster]<-rnorm(1,0,nogrowth_sigma)
    }
  }
  
  ##Assign only some valid records (makes some data invalid as no information from those animals --> proportion of animals not-informative)
  damagedecreaseObservation[sample(1:length(damagedecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(damagedecreaseObservation)), replace=FALSE)] <- -1
  pleopoddecreaseObservation[sample(1:length(pleopoddecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(pleopoddecreaseObservation)), replace=FALSE)] <- -1
  maturitychangeObservation[sample(1:length(maturitychangeObservation), size=floor((1-propAnimalsExtraInfo)*length(maturitychangeObservation)), replace=FALSE)] <- -1
  
  #--> Set data list ####
  sdata <-list(
    name = "gen50lib18grow52"
    ,nLobsters=nLobsters
    ,nPeriods=annualSteps
    ,obsGI=growth
    ,obsDam = damagedecreaseObservation
    ,obsPle = pleopoddecreaseObservation
    ,obsMat = maturitychangeObservation
    ,atLiberty=atLiberty
    ,isFem = 1
  )
  
  ldata[[length(ldata)+1]] <- sdata
  
  # sdata$nLobsters
  # sdata$obsGI %>% mean()
  # sdata$obsGI %>% sd()
  # sdata$obsDam %>% table()
  # sdata$obsPle %>% table()
  # sdata$obsMat %>% table()
  #####
} #end 210==210
#####

#### [NEW GEN SB: "gen80lib18grow52"] AS ABOVE--> CHANGED to (80% info) 0.2 UNINFORMATIVE MOULTING [SCENARIO 1] #####
if (220==220){
  #--> Initial Features ####
  set.seed(1)
  endYear<-2 #capture and recapture between time 0 to endYear
  nLobsters<-3000 #Number of lobsters
  annualSteps<-12 #How many time steps per year (12 for roughly months)
  propAnimalsExtraInfo<-0.80 #Proportion of animals with extra moult evidence
  upperLimAtLiberty <- 18 #Maximum number of months at Liberty
  
  startDates<-runif(nLobsters,0,endYear-1/annualSteps) #Uniform across period
  libertyLength<-runif(nLobsters, 1,upperLimAtLiberty)/annualSteps
  
  libertyLength[libertyLength<0]<-1/annualSteps #No negative times
  
  endDates<-startDates+libertyLength
  
  #Keep only lobsters recaptured before endYear
  select<-endDates<endYear
  
  select1k <- which(select)
  select1k <- sample(select1k, size=1000, replace=FALSE)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  #Moulting probability by time period
  moultP <- c(0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.15,0.05,0.02,0.02,0.01)
  moultP<-moultP[1:annualSteps]
  
  ##DAMAGE
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveDamageRegen <- c(0.1, 0.8)
  
  ##PLEOPOD
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObservePleopodRegen <- c(0.4, 0.99)
  
  ##MATURITY
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveMaturity <- c(0.2, 0.5)
  
  
  ##For model that relies on growth increments to estimate moulting probability
  #Growth increment distribution parameters
  growth_mu<-5
  growth_sigma<-2 #Includes variability in growth and measurement
  
  #Variation in size if they didn't moult (just observation error)
  nogrowth_sigma<-1
  
  #--> Generate Data ####
  #Contains the number of each time period that each lobster was at liberty
  atLiberty<-matrix(0,nLobsters,annualSteps)
  
  #Calculating proportion/number of steps that each lobster was at liberty
  for (iLobster in 1:nLobsters){
    started <- FALSE;
    ended<-FALSE;
    
    for (iYear in 1:endYear)
    {
      for (iStep in 1:annualSteps)
      {
        if (!started & startDates[iLobster]<((iYear-1)+(iStep)/annualSteps))
        {
          started <- TRUE
          atLiberty[iLobster,iStep]<-1-(startDates[iLobster]%%(1/annualSteps)*annualSteps)
        } else if (started & (endDates[iLobster]>((iYear-1)+(iStep)/annualSteps)))
        {
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+1 
        } else if (started & !ended)
        {
          ended<-TRUE;
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+(endDates[iLobster]%%(1/annualSteps)*annualSteps)
        }
      }
    }
  } ##end iLobster
  
  #--> Calculate moulting observations ####
  moultedP<-c() #Probability it moulted
  moultedOutcome<-c() #Did it moult?
  damagedecreaseObservation<-c()
  pleopoddecreaseObservation<-c()
  maturitychangeObservation<-c() ##NEW
  growth<-c() #The growth increment
  
  for (iLobster in 1:nLobsters){
    
    moultedP[iLobster] <- 1-prod((1-moultP)^atLiberty[iLobster,])
    moultedOutcome[iLobster] <- rbinom(1,1,moultedP[iLobster])
    
    damagedecreaseObservation[iLobster] <- rbinom(1,1,chanceObserveDamageRegen[moultedOutcome[iLobster]+1])
    pleopoddecreaseObservation[iLobster] <- rbinom(1,1,chanceObservePleopodRegen[moultedOutcome[iLobster]+1])
    maturitychangeObservation[iLobster] <- rbinom(1,1,chanceObserveMaturity[moultedOutcome[iLobster]+1]) ##NEW
    
    ##calculate the growth
    if (moultedOutcome[iLobster]){
      growth[iLobster]<-rnorm(1,growth_mu,growth_sigma)
    }else{
      growth[iLobster]<-rnorm(1,0,nogrowth_sigma)
    }
  }
  
  ##Assign only some valid records (makes some data invalid as no information from those animals --> proportion of animals not-informative)
  damagedecreaseObservation[sample(1:length(damagedecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(damagedecreaseObservation)), replace=FALSE)] <- -1
  pleopoddecreaseObservation[sample(1:length(pleopoddecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(pleopoddecreaseObservation)), replace=FALSE)] <- -1
  maturitychangeObservation[sample(1:length(maturitychangeObservation), size=floor((1-propAnimalsExtraInfo)*length(maturitychangeObservation)), replace=FALSE)] <- -1

  #--> Set data list ####
  sdata <-list(
    name = "gen80lib18grow52"
    ,nLobsters=nLobsters
    ,nPeriods=annualSteps
    ,obsGI=growth
    ,obsDam = damagedecreaseObservation
    ,obsPle = pleopoddecreaseObservation
    ,obsMat = maturitychangeObservation
    ,atLiberty=atLiberty
    ,isFem = 1
  )
  
  ldata[[length(ldata)+1]] <- sdata

  # sdata$nLobsters
  # sdata$obsGI %>% mean()
  # sdata$obsGI %>% sd()
  # sdata$obsDam %>% table()
  # sdata$obsPle %>% table()
  # sdata$obsMat %>% table()
  #####
} #end 220==220
#####

##############################################
############ SCENARIO 3 (TIME) ###############
##############################################

#### [NEW GEN SB: "gen30lib24grow52"] CHANGES TO TIME AT LIBERTY (18months --> 24months) As per above with Generated ERROR #####
if (300==300){
  #--> Initial Features ####
  set.seed(1)
  endYear<-2 #capture and recapture between time 0 to endYear
  nLobsters<-3000 #Number of lobsters
  annualSteps<-12 #How many time steps per year (12 for roughly months)
  propAnimalsExtraInfo<-0.3 #Proportion of animals with extra moult evidence
  upperLimAtLiberty <- 24 #Maximum number of months at Liberty
  
  startDates<-runif(nLobsters,0,endYear-1/annualSteps) #Uniform across period
  libertyLength<-runif(nLobsters, 1,upperLimAtLiberty)/annualSteps
  
  libertyLength[libertyLength<0]<-1/annualSteps #No negative times
  
  endDates<-startDates+libertyLength
  
  #Keep only lobsters recaptured before endYear
  select<-endDates<endYear
  select1k <- which(select)

  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  #Moulting probability by time period
  moultP <- c(0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.15,0.05,0.02,0.02,0.01)
  moultP<-moultP[1:annualSteps]
  
  ##DAMAGE
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveDamageRegen <- c(0.1, 0.8)
  
  ##PLEOPOD
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObservePleopodRegen <- c(0.4, 0.99)
  
  ##MATURITY
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveMaturity <- c(0.2, 0.5)
  
  ##For model that relies on growth increments to estimate moulting probability
  #Growth increment distribution parameters
  growth_mu<-5
  growth_sigma<-2 #Includes variability in growth and measurement
  
  #Variation in size if they didn't moult (just observation error)
  nogrowth_sigma<-1
  
  #--> Generate Data ####
  #Contains the number of each time period that each lobster was at liberty
  atLiberty<-matrix(0,nLobsters,annualSteps)
  
  #Calculating proportion/number of steps that each lobster was at liberty
  for (iLobster in 1:nLobsters){
    started <- FALSE;
    ended<-FALSE;
    
    for (iYear in 1:endYear)
    {
      for (iStep in 1:annualSteps)
      {
        if (!started & startDates[iLobster]<((iYear-1)+(iStep)/annualSteps))
        {
          started <- TRUE
          atLiberty[iLobster,iStep]<-1-(startDates[iLobster]%%(1/annualSteps)*annualSteps)
        } else if (started & (endDates[iLobster]>((iYear-1)+(iStep)/annualSteps)))
        {
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+1 
        } else if (started & !ended)
        {
          ended<-TRUE;
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+(endDates[iLobster]%%(1/annualSteps)*annualSteps)
        }
      }
    }
  } ##end iLobster
  
  #--> Modification of data (removal of all zeros) ####
  
  select1k <- atLiberty %>% rowSums() %>% "!="(0) %>% which()
  select1k <- sample(select1k, size=1000, replace=FALSE)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  atLiberty<-atLiberty[select1k,]
  
  #--> Calculate moulting observations ####
  moultedP<-c() #Probability it moulted
  moultedOutcome<-c() #Did it moult?
  damagedecreaseObservation<-c()
  pleopoddecreaseObservation<-c()
  maturitychangeObservation<-c() ##NEW
  growth<-c() #The growth increment
  
  for (iLobster in 1:nLobsters){
    
    moultedP[iLobster] <- 1-prod((1-moultP)^atLiberty[iLobster,])
    moultedOutcome[iLobster] <- rbinom(1,1,moultedP[iLobster])

    damagedecreaseObservation[iLobster] <- rbinom(1,1,chanceObserveDamageRegen[moultedOutcome[iLobster]+1])
    pleopoddecreaseObservation[iLobster] <- rbinom(1,1,chanceObservePleopodRegen[moultedOutcome[iLobster]+1])
    maturitychangeObservation[iLobster] <- rbinom(1,1,chanceObserveMaturity[moultedOutcome[iLobster]+1])
    
    ##calculate the growth
    if (moultedOutcome[iLobster]){
      growth[iLobster]<-rnorm(1,growth_mu,growth_sigma)
    }else{
      growth[iLobster]<-rnorm(1,0,nogrowth_sigma)
    }
  }
  
  ##Assign only some valid records (makes some data invalid as no information from those animals --> proportion of animals not-informative)
  damagedecreaseObservation[sample(1:length(damagedecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(damagedecreaseObservation)), replace=FALSE)] <- -1
  pleopoddecreaseObservation[sample(1:length(pleopoddecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(pleopoddecreaseObservation)), replace=FALSE)] <- -1
  maturitychangeObservation[sample(1:length(maturitychangeObservation), size=floor((1-propAnimalsExtraInfo)*length(maturitychangeObservation)), replace=FALSE)] <- -1
  
  #--> Set data list ####
  sdata <-list(
    name = "gen30lib24grow52"
    ,nLobsters=nLobsters
    ,nPeriods=annualSteps
    ,obsGI=growth
    ,obsDam = damagedecreaseObservation
    ,obsPle = pleopoddecreaseObservation
    ,obsMat = maturitychangeObservation
    ,atLiberty=atLiberty
    ,isFem = 1
  )
  
  ldata[[length(ldata)+1]] <- sdata

  # sdata$nLobsters
  # sdata$obsGI %>% mean()
  # sdata$obsGI %>% sd()
  # sdata$obsDam %>% table()
  # sdata$obsPle %>% table()
  # sdata$obsMat %>% table()
  #####
} #end 300==300
#####

#### [NEW GEN SB: "gen30lib12grow52"] CHANGES TO TIME AT LIBERTY (18months --> 12months) As per above with Generated ERROR #####
if (310==310){
  #--> Initial Features ####
  set.seed(1)
  endYear<-2 #capture and recapture between time 0 to endYear
  nLobsters<-3000 #Number of lobsters
  annualSteps<-12 #How many time steps per year (12 for roughly months)
  propAnimalsExtraInfo<-0.3 #Proportion of animals with extra moult evidence
  upperLimAtLiberty <- 12 #Maximum number of months at Liberty
  
  startDates<-runif(nLobsters,0,endYear-1/annualSteps) #Uniform across period
  libertyLength<-runif(nLobsters, 1,upperLimAtLiberty)/annualSteps
  
  libertyLength[libertyLength<0]<-1/annualSteps #No negative times
  
  endDates<-startDates+libertyLength
  
  #Keep only lobsters recaptured before endYear
  select<-endDates<endYear
  
  select1k <- which(select)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  #Moulting probability by time period
  moultP <- c(0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.15,0.05,0.02,0.02,0.01)
  moultP<-moultP[1:annualSteps]
  
  ##DAMAGE
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveDamageRegen <- c(0.1, 0.8)
  
  ##PLEOPOD
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObservePleopodRegen <- c(0.4, 0.99)
  
  ##MATURITY
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveMaturity <- c(0.2, 0.5)
  
  ##For model that relies on growth increments to estimate moulting probability
  #Growth increment distribution parameters
  growth_mu<-5
  growth_sigma<-2 #Includes variability in growth and measurement
  
  #Variation in size if they didn't moult (just observation error)
  nogrowth_sigma<-1

  #--> Generate Data ####
  #Contains the number of each time period that each lobster was at liberty
  atLiberty<-matrix(0,nLobsters,annualSteps)
  
  #Calculating proportion/number of steps that each lobster was at liberty
  for (iLobster in 1:nLobsters){
    started <- FALSE;
    ended<-FALSE;
    
    for (iYear in 1:endYear)
    {
      for (iStep in 1:annualSteps)
      {
        if (!started & startDates[iLobster]<((iYear-1)+(iStep)/annualSteps))
        {
          started <- TRUE
          atLiberty[iLobster,iStep]<-1-(startDates[iLobster]%%(1/annualSteps)*annualSteps)
        } else if (started & (endDates[iLobster]>((iYear-1)+(iStep)/annualSteps)))
        {
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+1 
        } else if (started & !ended)
        {
          ended<-TRUE;
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+(endDates[iLobster]%%(1/annualSteps)*annualSteps)
        }
      }
    }
  } ##end iLobster
  
  #--> Modification of data (removal of all zeros) ####
  select1k <- atLiberty %>% rowSums() %>% "!="(0) %>% which()
  select1k <- sample(select1k, size=1000, replace=FALSE)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  atLiberty<-atLiberty[select1k,]
  
  #--> Calculate moulting observations ####
  moultedP<-c() #Probability it moulted
  moultedOutcome<-c() #Did it moult?
  damagedecreaseObservation<-c()
  pleopoddecreaseObservation<-c()
  maturitychangeObservation<-c() ##NEW
  growth<-c() #The growth increment
  
  for (iLobster in 1:nLobsters){
    
    moultedP[iLobster] <- 1-prod((1-moultP)^atLiberty[iLobster,])
    moultedOutcome[iLobster] <- rbinom(1,1,moultedP[iLobster])

    damagedecreaseObservation[iLobster] <- rbinom(1,1,chanceObserveDamageRegen[moultedOutcome[iLobster]+1])
    pleopoddecreaseObservation[iLobster] <- rbinom(1,1,chanceObservePleopodRegen[moultedOutcome[iLobster]+1])
    maturitychangeObservation[iLobster] <- rbinom(1,1,chanceObserveMaturity[moultedOutcome[iLobster]+1])
    
    ##calculate the growth
    if (moultedOutcome[iLobster]){
      growth[iLobster]<-rnorm(1,growth_mu,growth_sigma)
    }else{
      growth[iLobster]<-rnorm(1,0,nogrowth_sigma)
    }
  }
  
  ##Assign only some valid records (makes some data invalid as no information from those animals --> proportion of animals not-informative)
  damagedecreaseObservation[sample(1:length(damagedecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(damagedecreaseObservation)), replace=FALSE)] <- -1
  pleopoddecreaseObservation[sample(1:length(pleopoddecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(pleopoddecreaseObservation)), replace=FALSE)] <- -1
  maturitychangeObservation[sample(1:length(maturitychangeObservation), size=floor((1-propAnimalsExtraInfo)*length(maturitychangeObservation)), replace=FALSE)] <- -1
  
  #--> Set data list ####
  sdata <-list(
    name = "gen30lib12grow52"
    ,nLobsters=nLobsters
    ,nPeriods=annualSteps
    ,obsGI=growth
    ,obsDam = damagedecreaseObservation
    ,obsPle = pleopoddecreaseObservation
    ,obsMat = maturitychangeObservation#rep(-1, length(pleopoddecreaseObservation))
    ,atLiberty=atLiberty
    ,isFem = 1
  )
  
  ldata[[length(ldata)+1]] <- sdata

  # sdata$nLobsters
  # sdata$obsGI %>% mean()
  # sdata$obsGI %>% sd()
  # sdata$obsDam %>% table()
  # sdata$obsPle %>% table()
  # sdata$obsMat %>% table()
  #####
} #end 310==310
#####

#### [NEW GEN SB: "gen30lib06grow52"] CHANGES TO TIME AT LIBERTY (18months --> 6months) As per above with Generated ERROR #####
if (320==320){
  #--> Initial Features ####
  set.seed(1)
  endYear<-2 #capture and recapture between time 0 to endYear
  nLobsters<-3000 #Number of lobsters
  annualSteps<-12 #How many time steps per year (12 for roughly months)
  propAnimalsExtraInfo<-0.3 #Proportion of animals with extra moult evidence
  upperLimAtLiberty <- 6 #Maximum number of months at Liberty
  
  startDates<-runif(nLobsters,0,endYear-1/annualSteps) #Uniform across period
  libertyLength<-runif(nLobsters, 1,upperLimAtLiberty)/annualSteps
  
  libertyLength[libertyLength<0]<-1/annualSteps #No negative times
  
  endDates<-startDates+libertyLength
  
  #Keep only lobsters recaptured before endYear
  select<-endDates<endYear
  
  select1k <- which(select)

  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  #Moulting probability by time period
  moultP <- c(0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.15,0.05,0.02,0.02,0.01)
  moultP<-moultP[1:annualSteps]
  
  ##DAMAGE
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveDamageRegen <- c(0.1, 0.8)
  
  ##PLEOPOD
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObservePleopodRegen <- c(0.4, 0.99)
  
  ##MATURITY
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveMaturity <- c(0.2, 0.5)
  
  ##For model that relies on growth increments to estimate moulting probability
  #Growth increment distribution parameters
  growth_mu<-5
  growth_sigma<-2 #Includes variability in growth and measurement
  
  #Variation in size if they didn't moult (just observation error)
  nogrowth_sigma<-1
  
  #--> Generate Data ####
  #Contains the number of each time period that each lobster was at liberty
  atLiberty<-matrix(0,nLobsters,annualSteps)
  
  #Calculating proportion/number of steps that each lobster was at liberty
  for (iLobster in 1:nLobsters){
    started <- FALSE;
    ended<-FALSE;
    
    for (iYear in 1:endYear)
    {
      for (iStep in 1:annualSteps)
      {
        if (!started & startDates[iLobster]<((iYear-1)+(iStep)/annualSteps))
        {
          started <- TRUE
          atLiberty[iLobster,iStep]<-1-(startDates[iLobster]%%(1/annualSteps)*annualSteps)
        } else if (started & (endDates[iLobster]>((iYear-1)+(iStep)/annualSteps)))
        {
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+1 
        } else if (started & !ended)
        {
          ended<-TRUE;
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+(endDates[iLobster]%%(1/annualSteps)*annualSteps)
        }
      }
    }
  } ##end iLobster
  
  
  #--> Modification of data (removal of all zeros) ####
  select1k <- atLiberty %>% rowSums() %>% "!="(0) %>% which()
  select1k <- sample(select1k, size=1000, replace=FALSE)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  atLiberty<-atLiberty[select1k,]
  
  #--> Calculate moulting observations ####
  moultedP<-c() #Probability it moulted
  moultedOutcome<-c() #Did it moult?
  damagedecreaseObservation<-c()
  pleopoddecreaseObservation<-c()
  maturitychangeObservation<-c()
  growth<-c() #The growth increment
  
  for (iLobster in 1:nLobsters){
    
    moultedP[iLobster] <- 1-prod((1-moultP)^atLiberty[iLobster,])
    moultedOutcome[iLobster] <- rbinom(1,1,moultedP[iLobster])
    
    damagedecreaseObservation[iLobster] <- rbinom(1,1,chanceObserveDamageRegen[moultedOutcome[iLobster]+1])
    pleopoddecreaseObservation[iLobster] <- rbinom(1,1,chanceObservePleopodRegen[moultedOutcome[iLobster]+1])
    maturitychangeObservation[iLobster] <- rbinom(1,1,chanceObserveMaturity[moultedOutcome[iLobster]+1])
    
    ##calculate the growth
    if (moultedOutcome[iLobster]){
      growth[iLobster]<-rnorm(1,growth_mu,growth_sigma)
    }else{
      growth[iLobster]<-rnorm(1,0,nogrowth_sigma)
    }
  }
  
  ##Assign only some valid records (makes some data invalid as no information from those animals --> proportion of animals not-informative)
  damagedecreaseObservation[sample(1:length(damagedecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(damagedecreaseObservation)), replace=FALSE)] <- -1
  pleopoddecreaseObservation[sample(1:length(pleopoddecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(pleopoddecreaseObservation)), replace=FALSE)] <- -1
  maturitychangeObservation[sample(1:length(maturitychangeObservation), size=floor((1-propAnimalsExtraInfo)*length(maturitychangeObservation)), replace=FALSE)] <- -1 ##NEW
  
  #--> Set data list ####
  sdata <-list(
    name = "gen30lib06grow52"
    ,nLobsters=nLobsters
    ,nPeriods=annualSteps
    ,obsGI=growth
    ,obsDam = damagedecreaseObservation
    ,obsPle = pleopoddecreaseObservation
    ,obsMat = maturitychangeObservation
    ,atLiberty=atLiberty
    ,isFem = 1
  )
  
  ldata[[length(ldata)+1]] <- sdata

  # sdata$nLobsters
  # sdata$obsGI %>% mean()
  # sdata$obsGI %>% sd()
  # sdata$obsDam %>% table()
  # sdata$obsPle %>% table()
  # sdata$obsMat %>% table()
  #####
} #end 320==320
#####

##############################################
############ SCENARIO 4 (GROWTH) #############
##############################################

#### [NEW GEN SB: "gen30lib18grow72"] CHANGES Growth N(7,2) #####
if (400==400){
  #--> Initial Features ####
  set.seed(1)
  endYear<-2 #capture and recapture between time 0 to endYear
  nLobsters<-3000 #Number of lobsters
  annualSteps<-12 #How many time steps per year (12 for roughly months)
  propAnimalsExtraInfo<-0.3 #Proportion of animals with extra moult evidence
  upperLimAtLiberty <- 18
  grow_mu <- 7
  grow_sigma <- 2
  
  startDates<-runif(nLobsters,0,endYear-1/annualSteps) #Uniform across period
  libertyLength<-runif(nLobsters, 1,upperLimAtLiberty)/annualSteps
  
  libertyLength[libertyLength<0]<-1/annualSteps #No negative times
  
  endDates<-startDates+libertyLength
  
  #Keep only lobsters recaptured before endYear
  select<-endDates<endYear
  
  select1k <- which(select)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  #Moulting probability by time period
  moultP <- c(0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.15,0.05,0.02,0.02,0.01)
  moultP<-moultP[1:annualSteps]
  
  ##DAMAGE
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveDamageRegen <- c(0.1, 0.8)
  
  ##PLEOPOD
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObservePleopodRegen <- c(0.4, 0.99)
  
  ##MATURITY
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveMaturity <- c(0.2, 0.5)
  
  ##For model that relies on growth increments to estimate moulting probability
  #Growth increment distribution parameters
  growth_mu<-grow_mu#5
  growth_sigma<-grow_sigma#2 #Includes variability in growth and measurement
  
  #Variation in size if they didn't moult (just observation error)
  nogrowth_sigma<-1

  #--> Generate Data ####
  #Contains the number of each time period that each lobster was at liberty
  atLiberty<-matrix(0,nLobsters,annualSteps)
  
  #Calculating proportion/number of steps that each lobster was at liberty
  for (iLobster in 1:nLobsters){
    started <- FALSE;
    ended<-FALSE;
    
    for (iYear in 1:endYear)
    {
      for (iStep in 1:annualSteps)
      {
        if (!started & startDates[iLobster]<((iYear-1)+(iStep)/annualSteps))
        {
          started <- TRUE
          atLiberty[iLobster,iStep]<-1-(startDates[iLobster]%%(1/annualSteps)*annualSteps)
        } else if (started & (endDates[iLobster]>((iYear-1)+(iStep)/annualSteps)))
        {
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+1 
        } else if (started & !ended)
        {
          ended<-TRUE;
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+(endDates[iLobster]%%(1/annualSteps)*annualSteps)
        }
      }
    }
  } ##end iLobster
  
  #--> Modification of data (removal of all zeros) ####
  select1k <- atLiberty %>% rowSums() %>% "!="(0) %>% which()
  select1k <- sample(select1k, size=1000, replace=FALSE)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  atLiberty<-atLiberty[select1k,]
  
  #--> Calculate moulting observations ####
  moultedP<-c() #Probability it moulted
  moultedOutcome<-c() #Did it moult?
  damagedecreaseObservation<-c()
  pleopoddecreaseObservation<-c()
  maturitychangeObservation<-c()
  growth<-c() #The growth increment
  
  for (iLobster in 1:nLobsters){
    
    moultedP[iLobster] <- 1-prod((1-moultP)^atLiberty[iLobster,])
    moultedOutcome[iLobster] <- rbinom(1,1,moultedP[iLobster])

    damagedecreaseObservation[iLobster] <- rbinom(1,1,chanceObserveDamageRegen[moultedOutcome[iLobster]+1])
    pleopoddecreaseObservation[iLobster] <- rbinom(1,1,chanceObservePleopodRegen[moultedOutcome[iLobster]+1])
    maturitychangeObservation[iLobster] <- rbinom(1,1,chanceObserveMaturity[moultedOutcome[iLobster]+1])
    
    ##calculate the growth
    if (moultedOutcome[iLobster]){
      growth[iLobster]<-rnorm(1,growth_mu,growth_sigma)
    }else{
      growth[iLobster]<-rnorm(1,0,nogrowth_sigma)
    }
  }
  
  ##Assign only some valid records (makes some data invalid as no information from those animals --> proportion of animals not-informative)
  damagedecreaseObservation[sample(1:length(damagedecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(damagedecreaseObservation)), replace=FALSE)] <- -1
  pleopoddecreaseObservation[sample(1:length(pleopoddecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(pleopoddecreaseObservation)), replace=FALSE)] <- -1
  maturitychangeObservation[sample(1:length(maturitychangeObservation), size=floor((1-propAnimalsExtraInfo)*length(maturitychangeObservation)), replace=FALSE)] <- -1 ##NEW

  #--> Set data list ####
  sdata <-list(
    name = "gen30lib18grow72"
    ,nLobsters=nLobsters
    ,nPeriods=annualSteps
    ,obsGI=growth
    ,obsDam = damagedecreaseObservation
    ,obsPle = pleopoddecreaseObservation
    ,obsMat = maturitychangeObservation#rep(-1, length(pleopoddecreaseObservation))
    ,atLiberty=atLiberty
    ,isFem = 1
  )
  
  ldata[[length(ldata)+1]] <- sdata

  # sdata$nLobsters
  # sdata$obsGI %>% mean()
  # sdata$obsGI %>% sd()
  # sdata$obsDam %>% table()
  # sdata$obsPle %>% table()
  # sdata$obsMat %>% table()
  #####
} #end 400==400
#####

#### [NEW GEN SB: "gen30lib18grow32"] CHANGES Growth N(3,2) #####
if (410==410){
  #--> Initial Features ####
  set.seed(1)
  endYear<-2 #capture and recapture between time 0 to endYear
  nLobsters<-3000 #Number of lobsters
  annualSteps<-12 #How many time steps per year (12 for roughly months)
  propAnimalsExtraInfo<-0.3 #Proportion of animals with extra moult evidence
  upperLimAtLiberty <- 18
  grow_mu <- 3
  grow_sigma <- 2
  
  startDates<-runif(nLobsters,0,endYear-1/annualSteps) #Uniform across period
  libertyLength<-runif(nLobsters, 1,upperLimAtLiberty)/annualSteps
  
  libertyLength[libertyLength<0]<-1/annualSteps #No negative times
  
  endDates<-startDates+libertyLength
  
  #Keep only lobsters recaptured before endYear
  select<-endDates<endYear
  
  select1k <- which(select)

  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  #Moulting probability by time period
  moultP <- c(0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.15,0.05,0.02,0.02,0.01)
  moultP<-moultP[1:annualSteps]
  
  ##DAMAGE
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveDamageRegen <- c(0.1, 0.8)
  
  ##PLEOPOD
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObservePleopodRegen <- c(0.4, 0.99)
  
  ##MATURITY
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveMaturity <- c(0.2, 0.5)
  
  ##For model that relies on growth increments to estimate moulting probability
  #Growth increment distribution parameters
  growth_mu<-grow_mu#5
  growth_sigma<-grow_sigma#2 #Includes variability in growth and measurement
  
  #Variation in size if they didn't moult (just observation error)
  nogrowth_sigma<-1

  #--> Generate Data ####
  #Contains the number of each time period that each lobster was at liberty
  atLiberty<-matrix(0,nLobsters,annualSteps)
  
  #Calculating proportion/number of steps that each lobster was at liberty
  for (iLobster in 1:nLobsters){
    started <- FALSE;
    ended<-FALSE;
    
    for (iYear in 1:endYear)
    {
      for (iStep in 1:annualSteps)
      {
        if (!started & startDates[iLobster]<((iYear-1)+(iStep)/annualSteps))
        {
          started <- TRUE
          atLiberty[iLobster,iStep]<-1-(startDates[iLobster]%%(1/annualSteps)*annualSteps)
        } else if (started & (endDates[iLobster]>((iYear-1)+(iStep)/annualSteps)))
        {
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+1 
        } else if (started & !ended)
        {
          ended<-TRUE;
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+(endDates[iLobster]%%(1/annualSteps)*annualSteps)
        }
      }
    }
  } ##end iLobster
  
  
  #--> Modification of data (removal of all zeros) ####
  select1k <- atLiberty %>% rowSums() %>% "!="(0) %>% which()
  select1k <- sample(select1k, size=1000, replace=FALSE)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  atLiberty<-atLiberty[select1k,]
  
  #--> Calculate moulting observations ####
  moultedP<-c() #Probability it moulted
  moultedOutcome<-c() #Did it moult?
  damagedecreaseObservation<-c()
  pleopoddecreaseObservation<-c()
  maturitychangeObservation<-c()
  growth<-c() #The growth increment
  
  for (iLobster in 1:nLobsters){
    
    moultedP[iLobster] <- 1-prod((1-moultP)^atLiberty[iLobster,])
    moultedOutcome[iLobster] <- rbinom(1,1,moultedP[iLobster])

    damagedecreaseObservation[iLobster] <- rbinom(1,1,chanceObserveDamageRegen[moultedOutcome[iLobster]+1])
    pleopoddecreaseObservation[iLobster] <- rbinom(1,1,chanceObservePleopodRegen[moultedOutcome[iLobster]+1])
    maturitychangeObservation[iLobster] <- rbinom(1,1,chanceObserveMaturity[moultedOutcome[iLobster]+1]) ##NEW
    
    ##calculate the growth
    if (moultedOutcome[iLobster]){
      growth[iLobster]<-rnorm(1,growth_mu,growth_sigma)
    }else{
      growth[iLobster]<-rnorm(1,0,nogrowth_sigma)
    }
  }
  
  ##Assign only some valid records (makes some data invalid as no information from those animals --> proportion of animals not-informative)
  damagedecreaseObservation[sample(1:length(damagedecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(damagedecreaseObservation)), replace=FALSE)] <- -1
  pleopoddecreaseObservation[sample(1:length(pleopoddecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(pleopoddecreaseObservation)), replace=FALSE)] <- -1
  maturitychangeObservation[sample(1:length(maturitychangeObservation), size=floor((1-propAnimalsExtraInfo)*length(maturitychangeObservation)), replace=FALSE)] <- -1
  
  #--> Set data list ####
  sdata <-list(
    name = "gen30lib18grow32"
    ,nLobsters=nLobsters
    ,nPeriods=annualSteps
    ,obsGI=growth
    ,obsDam = damagedecreaseObservation
    ,obsPle = pleopoddecreaseObservation
    ,obsMat = maturitychangeObservation
    ,atLiberty=atLiberty
    ,isFem = 1
  )
  
  ldata[[length(ldata)+1]] <- sdata

  # sdata$nLobsters
  # sdata$obsGI %>% mean()
  # sdata$obsGI %>% sd()
  # sdata$obsDam %>% table()
  # sdata$obsPle %>% table()
  # sdata$obsMat %>% table()
  #####
} #end 410==410
#####

#### [NEW GEN SB: "gen30lib18grow11"] CHANGES Growth N(1,1) #####
if (420==420){
  #--> Initial Features ####
  set.seed(1)
  endYear<-2 #capture and recapture between time 0 to endYear
  nLobsters<-3000 #Number of lobsters
  annualSteps<-12 #How many time steps per year (12 for roughly months)
  propAnimalsExtraInfo<-0.3 #Proportion of animals with extra moult evidence
  upperLimAtLiberty <- 18
  grow_mu <- 1
  grow_sigma <- 1
  
  startDates<-runif(nLobsters,0,endYear-1/annualSteps) #Uniform across period
  libertyLength<-runif(nLobsters, 1,upperLimAtLiberty)/annualSteps
  
  libertyLength[libertyLength<0]<-1/annualSteps #No negative times
  
  endDates<-startDates+libertyLength
  
  #Keep only lobsters recaptured before endYear
  select<-endDates<endYear
  
  select1k <- which(select)

  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  #Moulting probability by time period
  moultP <- c(0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.15,0.05,0.02,0.02,0.01)
  moultP<-moultP[1:annualSteps]
  
  ##DAMAGE
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveDamageRegen <- c(0.1, 0.8)
  
  ##PLEOPOD
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObservePleopodRegen <- c(0.4, 0.99)
  
  ##MATURITY
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveMaturity <- c(0.2, 0.5)
  
  ##For model that relies on growth increments to estimate moulting probability
  #Growth increment distribution parameters
  growth_mu<-grow_mu#5
  growth_sigma<-grow_sigma#2 #Includes variability in growth and measurement
  
  #Variation in size if they didn't moult (just observation error)
  nogrowth_sigma<-1

  #--> Generate Data ####
  #Contains the number of each time period that each lobster was at liberty
  atLiberty<-matrix(0,nLobsters,annualSteps)
  
  #Calculating proportion/number of steps that each lobster was at liberty
  for (iLobster in 1:nLobsters){
    started <- FALSE;
    ended<-FALSE;
    
    for (iYear in 1:endYear)
    {
      for (iStep in 1:annualSteps)
      {
        if (!started & startDates[iLobster]<((iYear-1)+(iStep)/annualSteps))
        {
          started <- TRUE
          atLiberty[iLobster,iStep]<-1-(startDates[iLobster]%%(1/annualSteps)*annualSteps)
        } else if (started & (endDates[iLobster]>((iYear-1)+(iStep)/annualSteps)))
        {
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+1 
        } else if (started & !ended)
        {
          ended<-TRUE;
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+(endDates[iLobster]%%(1/annualSteps)*annualSteps)
        }
      }
    }
  } ##end iLobster
  
  #--> Modification of data (removal of all zeros) ####
  select1k <- atLiberty %>% rowSums() %>% "!="(0) %>% which()
  select1k <- sample(select1k, size=1000, replace=FALSE)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  atLiberty<-atLiberty[select1k,]
  
  #--> Calculate moulting observations ####
  moultedP<-c() #Probability it moulted
  moultedOutcome<-c() #Did it moult?
  damagedecreaseObservation<-c()
  pleopoddecreaseObservation<-c()
  maturitychangeObservation<-c()
  growth<-c() #The growth increment
  
  for (iLobster in 1:nLobsters){
    
    moultedP[iLobster] <- 1-prod((1-moultP)^atLiberty[iLobster,])
    moultedOutcome[iLobster] <- rbinom(1,1,moultedP[iLobster])
    
    damagedecreaseObservation[iLobster] <- rbinom(1,1,chanceObserveDamageRegen[moultedOutcome[iLobster]+1])
    pleopoddecreaseObservation[iLobster] <- rbinom(1,1,chanceObservePleopodRegen[moultedOutcome[iLobster]+1])
    maturitychangeObservation[iLobster] <- rbinom(1,1,chanceObserveMaturity[moultedOutcome[iLobster]+1]) ##NEW
    
    ##calculate the growth
    if (moultedOutcome[iLobster]){
      growth[iLobster]<-rnorm(1,growth_mu,growth_sigma)
    }else{
      growth[iLobster]<-rnorm(1,0,nogrowth_sigma)
    }
  }
  
  ##Assign only some valid records (makes some data invalid as no information from those animals --> proportion of animals not-informative)
  damagedecreaseObservation[sample(1:length(damagedecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(damagedecreaseObservation)), replace=FALSE)] <- -1
  pleopoddecreaseObservation[sample(1:length(pleopoddecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(pleopoddecreaseObservation)), replace=FALSE)] <- -1
  maturitychangeObservation[sample(1:length(maturitychangeObservation), size=floor((1-propAnimalsExtraInfo)*length(maturitychangeObservation)), replace=FALSE)] <- -1

  #--> Set data list ####
  sdata <-list(
    name = "gen30lib18grow11"
    ,nLobsters=nLobsters
    ,nPeriods=annualSteps
    ,obsGI=growth
    ,obsDam = damagedecreaseObservation
    ,obsPle = pleopoddecreaseObservation
    ,obsMat = maturitychangeObservation
    ,atLiberty=atLiberty
    ,isFem = 1
  )
  
  ldata[[length(ldata)+1]] <- sdata

  # sdata$nLobsters
  # sdata$obsGI %>% mean()
  # sdata$obsGI %>% sd()
  # sdata$obsDam %>% table()
  # sdata$obsPle %>% table()
  # sdata$obsMat %>% table()
  #####
} #end 420==420
#####

####################################################
############ SCENARIO 5 (RM EVIDENCE) ##############
####################################################

#### [NEW GEN SB: "gen30lib18grow52_grDamage"] CHANGES Growth + Damage #####
if (500==500){
  #--> Initial Features ####
  set.seed(1)
  endYear<-2 #capture and recapture between time 0 to endYear
  nLobsters<-3000 #Number of lobsters
  annualSteps<-12 #How many time steps per year (12 for roughly months)
  propAnimalsExtraInfo<-0.3 #Proportion of animals with extra moult evidence
  upperLimAtLiberty <- 18
  grow_mu <- 5
  grow_sigma <- 2
  
  startDates<-runif(nLobsters,0,endYear-1/annualSteps) #Uniform across period
  libertyLength<-runif(nLobsters, 1,upperLimAtLiberty)/annualSteps
  
  libertyLength[libertyLength<0]<-1/annualSteps #No negative times
  
  endDates<-startDates+libertyLength
  
  #Keep only lobsters recaptured before endYear
  select<-endDates<endYear
  
  select1k <- which(select)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  #Moulting probability by time period
  moultP <- c(0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.15,0.05,0.02,0.02,0.01)
  moultP<-moultP[1:annualSteps]
  
  ##DAMAGE
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveDamageRegen <- c(0.1, 0.8)
  
  ##PLEOPOD
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObservePleopodRegen <- c(0.4, 0.99)
  
  ##MATURITY
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveMaturity <- c(0.2, 0.5)
  
  ##For model that relies on growth increments to estimate moulting probability
  #Growth increment distribution parameters
  growth_mu<-grow_mu#5
  growth_sigma<-grow_sigma#2 #Includes variability in growth and measurement
  
  #Variation in size if they didn't moult (just observation error)
  nogrowth_sigma<-1

  #--> Generate Data ####
  #Contains the number of each time period that each lobster was at liberty
  atLiberty<-matrix(0,nLobsters,annualSteps)
  
  #Calculating proportion/number of steps that each lobster was at liberty
  for (iLobster in 1:nLobsters){
    started <- FALSE;
    ended<-FALSE;
    
    for (iYear in 1:endYear)
    {
      for (iStep in 1:annualSteps)
      {
        if (!started & startDates[iLobster]<((iYear-1)+(iStep)/annualSteps))
        {
          started <- TRUE
          atLiberty[iLobster,iStep]<-1-(startDates[iLobster]%%(1/annualSteps)*annualSteps)
        } else if (started & (endDates[iLobster]>((iYear-1)+(iStep)/annualSteps)))
        {
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+1 
        } else if (started & !ended)
        {
          ended<-TRUE;
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+(endDates[iLobster]%%(1/annualSteps)*annualSteps)
        }
      }
    }
  } ##end iLobster
  
  #--> Modification of data (removal of all zeros) ####
  select1k <- atLiberty %>% rowSums() %>% "!="(0) %>% which()
  select1k <- sample(select1k, size=1000, replace=FALSE)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  atLiberty<-atLiberty[select1k,]
  
  #--> Calculate moulting observations ####
  moultedP<-c() #Probability it moulted
  moultedOutcome<-c() #Did it moult?
  damagedecreaseObservation<-c()
  pleopoddecreaseObservation<-c()
  maturitychangeObservation<-c()
  growth<-c() #The growth increment
  
  for (iLobster in 1:nLobsters){
    
    moultedP[iLobster] <- 1-prod((1-moultP)^atLiberty[iLobster,])
    moultedOutcome[iLobster] <- rbinom(1,1,moultedP[iLobster])

    damagedecreaseObservation[iLobster] <- rbinom(1,1,chanceObserveDamageRegen[moultedOutcome[iLobster]+1])
    pleopoddecreaseObservation[iLobster] <- rbinom(1,1,chanceObservePleopodRegen[moultedOutcome[iLobster]+1])
    maturitychangeObservation[iLobster] <- rbinom(1,1,chanceObserveMaturity[moultedOutcome[iLobster]+1])
    
    ##calculate the growth
    if (moultedOutcome[iLobster]){
      growth[iLobster]<-rnorm(1,growth_mu,growth_sigma)
    }else{
      growth[iLobster]<-rnorm(1,0,nogrowth_sigma)
    }
  }
  
  ##Assign only some valid records (makes some data invalid as no information from those animals --> proportion of animals not-informative)
  damagedecreaseObservation[sample(1:length(damagedecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(damagedecreaseObservation)), replace=FALSE)] <- -1
  pleopoddecreaseObservation[sample(1:length(pleopoddecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(pleopoddecreaseObservation)), replace=FALSE)] <- -1
  maturitychangeObservation[sample(1:length(maturitychangeObservation), size=floor((1-propAnimalsExtraInfo)*length(maturitychangeObservation)), replace=FALSE)] <- -1
  
  #--> Set data list ####
  sdata <-list(
    name = "gen30lib18grow52_grDamage"
    ,nLobsters=nLobsters
    ,nPeriods=annualSteps
    ,obsGI=growth
    ,obsDam = damagedecreaseObservation
    # ,obsDam = rep(-1, length(damagedecreaseObservation))
    # ,obsPle = pleopoddecreaseObservation
    ,obsPle = rep(-1, length(pleopoddecreaseObservation))
    # ,obsMat = maturitychangeObservation
    ,obsMat = rep(-1, length(maturitychangeObservation))
    ,atLiberty=atLiberty
    ,isFem = 1
  )
  
  ldata[[length(ldata)+1]] <- sdata

  sdata$nLobsters
  sdata$obsGI %>% mean()
  sdata$obsGI %>% sd()
  sdata$obsDam %>% table()
  sdata$obsPle %>% table()
  sdata$obsMat %>% table()
  #####
} #end 510==510
#####

#### [NEW GEN SB: "gen30lib18grow52_grPleopod"] CHANGES Growth + Pleopod #####
if (510==510){
  #--> Initial Features ####
  set.seed(1)
  endYear<-2 #capture and recapture between time 0 to endYear
  nLobsters<-3000 #Number of lobsters
  annualSteps<-12 #How many time steps per year (12 for roughly months)
  propAnimalsExtraInfo<-0.3 #Proportion of animals with extra moult evidence
  upperLimAtLiberty <- 18
  grow_mu <- 5
  grow_sigma <- 2
  
  startDates<-runif(nLobsters,0,endYear-1/annualSteps) #Uniform across period
  libertyLength<-runif(nLobsters, 1,upperLimAtLiberty)/annualSteps
  
  libertyLength[libertyLength<0]<-1/annualSteps #No negative times
  
  endDates<-startDates+libertyLength
  
  #Keep only lobsters recaptured before endYear
  select<-endDates<endYear
  
  select1k <- which(select)

  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  #Moulting probability by time period
  moultP <- c(0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.15,0.05,0.02,0.02,0.01)
  moultP<-moultP[1:annualSteps]
  
  ##DAMAGE
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveDamageRegen <- c(0.1, 0.8)
  
  ##PLEOPOD
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObservePleopodRegen <- c(0.4, 0.99)
  
  ##MATURITY
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveMaturity <- c(0.2, 0.5)
  
  ##For model that relies on growth increments to estimate moulting probability
  #Growth increment distribution parameters
  growth_mu<-grow_mu#5
  growth_sigma<-grow_sigma#2 #Includes variability in growth and measurement
  
  #Variation in size if they didn't moult (just observation error)
  nogrowth_sigma<-1

  #--> Generate Data ####
  #Contains the number of each time period that each lobster was at liberty
  atLiberty<-matrix(0,nLobsters,annualSteps)
  
  #Calculating proportion/number of steps that each lobster was at liberty
  for (iLobster in 1:nLobsters){
    started <- FALSE;
    ended<-FALSE;
    
    for (iYear in 1:endYear)
    {
      for (iStep in 1:annualSteps)
      {
        if (!started & startDates[iLobster]<((iYear-1)+(iStep)/annualSteps))
        {
          started <- TRUE
          atLiberty[iLobster,iStep]<-1-(startDates[iLobster]%%(1/annualSteps)*annualSteps)
        } else if (started & (endDates[iLobster]>((iYear-1)+(iStep)/annualSteps)))
        {
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+1 
        } else if (started & !ended)
        {
          ended<-TRUE;
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+(endDates[iLobster]%%(1/annualSteps)*annualSteps)
        }
      }
    }
  } ##end iLobster
  
  #--> Modification of data (removal of all zeros) ####
  select1k <- atLiberty %>% rowSums() %>% "!="(0) %>% which()
  select1k <- sample(select1k, size=1000, replace=FALSE)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  atLiberty<-atLiberty[select1k,]
  
  #--> Calculate moulting observations ####
  moultedP<-c() #Probability it moulted
  moultedOutcome<-c() #Did it moult?
  damagedecreaseObservation<-c()
  pleopoddecreaseObservation<-c()
  maturitychangeObservation<-c()
  growth<-c() #The growth increment
  
  for (iLobster in 1:nLobsters){
    
    moultedP[iLobster] <- 1-prod((1-moultP)^atLiberty[iLobster,])
    moultedOutcome[iLobster] <- rbinom(1,1,moultedP[iLobster])

    damagedecreaseObservation[iLobster] <- rbinom(1,1,chanceObserveDamageRegen[moultedOutcome[iLobster]+1])
    pleopoddecreaseObservation[iLobster] <- rbinom(1,1,chanceObservePleopodRegen[moultedOutcome[iLobster]+1])
    maturitychangeObservation[iLobster] <- rbinom(1,1,chanceObserveMaturity[moultedOutcome[iLobster]+1])
    
    ##calculate the growth
    if (moultedOutcome[iLobster]){
      growth[iLobster]<-rnorm(1,growth_mu,growth_sigma)
    }else{
      growth[iLobster]<-rnorm(1,0,nogrowth_sigma)
    }
  }
  
  ##Assign only some valid records (makes some data invalid as no information from those animals --> proportion of animals not-informative)
  damagedecreaseObservation[sample(1:length(damagedecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(damagedecreaseObservation)), replace=FALSE)] <- -1
  pleopoddecreaseObservation[sample(1:length(pleopoddecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(pleopoddecreaseObservation)), replace=FALSE)] <- -1
  maturitychangeObservation[sample(1:length(maturitychangeObservation), size=floor((1-propAnimalsExtraInfo)*length(maturitychangeObservation)), replace=FALSE)] <- -1
  
  #--> Set data list ####
  sdata <-list(
    name = "gen30lib18grow52_grPleopod"
    ,nLobsters=nLobsters
    ,nPeriods=annualSteps
    ,obsGI=growth
    # ,obsDam = damagedecreaseObservation
    ,obsDam = rep(-1, length(damagedecreaseObservation))
    ,obsPle = pleopoddecreaseObservation
    # ,obsPle = rep(-1, length(pleopoddecreaseObservation))
    # ,obsMat = maturitychangeObservation
    ,obsMat = rep(-1, length(maturitychangeObservation))
    ,atLiberty=atLiberty
    ,isFem = 1
  )
  
  ldata[[length(ldata)+1]] <- sdata
  
  # sdata$nLobsters
  # sdata$obsGI %>% mean()
  # sdata$obsGI %>% sd()
  # sdata$obsDam %>% table()
  # sdata$obsPle %>% table()
  # sdata$obsMat %>% table()
  #####
} #end 510==510
#####

#### [NEW GEN SB: "gen30lib18grow52_grMaturity"] CHANGES Growth + Maturity #####
if (520==520){
  #--> Initial Features ####
  set.seed(1)
  endYear<-2 #capture and recapture between time 0 to endYear
  nLobsters<-3000 #Number of lobsters
  annualSteps<-12 #How many time steps per year (12 for roughly months)
  propAnimalsExtraInfo<-0.3 #Proportion of animals with extra moult evidence
  upperLimAtLiberty <- 18
  grow_mu <- 5
  grow_sigma <- 2
  
  startDates<-runif(nLobsters,0,endYear-1/annualSteps) #Uniform across period
  libertyLength<-runif(nLobsters, 1,upperLimAtLiberty)/annualSteps
  
  libertyLength[libertyLength<0]<-1/annualSteps #No negative times
  
  endDates<-startDates+libertyLength
  
  #Keep only lobsters recaptured before endYear
  select<-endDates<endYear
  
  select1k <- which(select)

  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  #Moulting probability by time period
  moultP <- c(0.05,0.05,0.1,0.4,0.7,0.4,0.2,0.15,0.05,0.02,0.02,0.01)
  moultP<-moultP[1:annualSteps]
  
  ##DAMAGE
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveDamageRegen <- c(0.1, 0.8)
  
  ##PLEOPOD
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObservePleopodRegen <- c(0.4, 0.99)
  
  ##MATURITY
  #First element is it didn't moult but data indicates it did
  #Second element is it did moult and data indicates it did
  chanceObserveMaturity <- c(0.2, 0.5)
  
  ##For model that relies on growth increments to estimate moulting probability
  #Growth increment distribution parameters
  growth_mu<-grow_mu#5
  growth_sigma<-grow_sigma#2 #Includes variability in growth and measurement
  
  #Variation in size if they didn't moult (just observation error)
  nogrowth_sigma<-1
  
  #--> Generate Data ####
  #Contains the number of each time period that each lobster was at liberty
  atLiberty<-matrix(0,nLobsters,annualSteps)
  
  #Calculating proportion/number of steps that each lobster was at liberty
  for (iLobster in 1:nLobsters){
    started <- FALSE;
    ended<-FALSE;
    
    for (iYear in 1:endYear)
    {
      for (iStep in 1:annualSteps)
      {
        if (!started & startDates[iLobster]<((iYear-1)+(iStep)/annualSteps))
        {
          started <- TRUE
          atLiberty[iLobster,iStep]<-1-(startDates[iLobster]%%(1/annualSteps)*annualSteps)
        } else if (started & (endDates[iLobster]>((iYear-1)+(iStep)/annualSteps)))
        {
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+1 
        } else if (started & !ended)
        {
          ended<-TRUE;
          atLiberty[iLobster,iStep]<-atLiberty[iLobster,iStep]+(endDates[iLobster]%%(1/annualSteps)*annualSteps)
        }
      }
    }
  } ##end iLobster
  
  
  #--> Modification of data (removal of all zeros) ####
  select1k <- atLiberty %>% rowSums() %>% "!="(0) %>% which()
  select1k <- sample(select1k, size=1000, replace=FALSE)
  
  startDates<-startDates[select1k]
  endDates<-endDates[select1k]
  libertyLength<-libertyLength[select1k]
  nLobsters<-length(select1k)
  
  atLiberty<-atLiberty[select1k,]
  
  #--> Calculate moulting observations ####
  moultedP<-c() #Probability it moulted
  moultedOutcome<-c() #Did it moult?
  damagedecreaseObservation<-c()
  pleopoddecreaseObservation<-c()
  maturitychangeObservation<-c()
  growth<-c() #The growth increment
  
  for (iLobster in 1:nLobsters){
    
    moultedP[iLobster] <- 1-prod((1-moultP)^atLiberty[iLobster,])
    moultedOutcome[iLobster] <- rbinom(1,1,moultedP[iLobster])

    damagedecreaseObservation[iLobster] <- rbinom(1,1,chanceObserveDamageRegen[moultedOutcome[iLobster]+1])
    pleopoddecreaseObservation[iLobster] <- rbinom(1,1,chanceObservePleopodRegen[moultedOutcome[iLobster]+1])
    maturitychangeObservation[iLobster] <- rbinom(1,1,chanceObserveMaturity[moultedOutcome[iLobster]+1])
    
    ##calculate the growth
    if (moultedOutcome[iLobster]){
      growth[iLobster]<-rnorm(1,growth_mu,growth_sigma)
    }else{
      growth[iLobster]<-rnorm(1,0,nogrowth_sigma)
    }
  }
  
  ##Assign only some valid records (makes some data invalid as no information from those animals --> proportion of animals not-informative)
  damagedecreaseObservation[sample(1:length(damagedecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(damagedecreaseObservation)), replace=FALSE)] <- -1
  pleopoddecreaseObservation[sample(1:length(pleopoddecreaseObservation), size=floor((1-propAnimalsExtraInfo)*length(pleopoddecreaseObservation)), replace=FALSE)] <- -1
  maturitychangeObservation[sample(1:length(maturitychangeObservation), size=floor((1-propAnimalsExtraInfo)*length(maturitychangeObservation)), replace=FALSE)] <- -1

  #--> Set data list ####
  sdata <-list(
    name = "gen30lib18grow52_grMaturity"
    ,nLobsters=nLobsters
    ,nPeriods=annualSteps
    ,obsGI=growth
    # ,obsDam = damagedecreaseObservation
    ,obsDam = rep(-1, length(damagedecreaseObservation))
    # ,obsPle = pleopoddecreaseObservation
    ,obsPle = rep(-1, length(pleopoddecreaseObservation))
    ,obsMat = maturitychangeObservation
    # ,obsMat = rep(-1, length(maturitychangeObservation))
    ,atLiberty=atLiberty
    ,isFem = 1
  )
  
  ldata[[length(ldata)+1]] <- sdata
  
  # sdata$nLobsters
  # sdata$obsGI %>% mean()
  # sdata$obsGI %>% sd()
  # sdata$obsDam %>% table()
  # sdata$obsPle %>% table()
  # sdata$obsMat %>% table()
  #####
} #end 520==520
#####

################################################
############ EXECUTE STAN MODELS ###############
################################################

#--> Test Stan works ####
# func_testStanWorks()

#--> R Stan ####
##Simulation data
for (mods in 1:length(ldata)){

  sdata <- ldata[[mods]]
  
  sdata$isFem
  sdata$obsGI %>% hist(main="obsGI", breaks=100)
  sdata$obsDam %>% table()
  sdata$obsPle %>% table()
  sdata$obsMat %>% table()
  sdata$atLiberty %>% range() #(2 times occurred)
  sdata$atLiberty %>% hist(main="atlib")
  sdata$atLiberty %>% rowSums() %>% hist(main="rowSums atLiberty")
  sdata$atLiberty %>% colMeans()
  
  mod_model <- "models/rSTAN_moultProb_np.stan" #No informed priors used for simulation data
  # mod_model <- "models/rSTAN_moultProb_priors.stan"  #Structured model demonstrating prior implementation
  
  ##Specify specs for model
  iter = 1e5
  model <- rstan::stan_model(mod_model, auto_write=FALSE)
  
  set.seed(1)
  
  #--> Running ####
  fit_mcmc = rstan::sampling(
    object=model
    , data = sdata
    , control = list(adapt_delta=0.8) #default
    , chains = 4
    , cores = 4
    , iter = iter
    , warmup = iter/2
    , thin = 4
  )
  
  ##--> Checking results ####
  fit_mcmc
  capture.output(print(fit_mcmc))

  ##checks
  if(1==100){
    sum(sdata$obsPle) / sdata$nLobsters
    sum(sdata$obsDam) / sdata$nLobsters
    
    fit <- rstan::extract(fit_mcmc)
    
    par(mfrow=c(2,1))
    
    fit$growth_mu %>% hist()
    fit$growth_sigma %>% hist()
    fit$nogrowth_sigma %>% hist()
    
    fit$damageIndicates_moulted %>% hist()
    fit$damageIndicates_nonmoulted %>% hist()
    
    fit$pleopodIndicates_moulted %>% hist()
    fit$pleopodIndicates_nonmoulted %>% hist()
    
    fit$maturityIndicates_moulted %>% hist()
    fit$maturityIndicates_nonmoulted %>% hist()
    
  }
  
  if(1==100){
    fit <- rstan::extract(fit_mcmc)
    
    #--> Basic plot ####
    tmp <- data.frame(
      p = fit$p %>% colMeans()
      ,sd = fit$p %>% apply(2, sd)
    )
    
    # Create a Line plot with ggplot, connecting the points
    ggplot(tmp, aes(x = factor(seq_along(p)), y = p)) +
      # geom_violin(fill = "skyblue", color = "black", alpha = 0.7) +
      geom_point(position = position_jitter(width = 0), color = "black", alpha = 0.5) +
      geom_line(aes(x = factor(seq_along(p)), y = p, group = 1), color = "black", alpha = 0.5) +
      geom_errorbar(aes(ymin = p - 2*sd, ymax = p + 2*sd), width = 0.2, position = position_dodge(0.75)) +
      geom_line(aes(x = factor(seq_along(p)), y = c(0.05, 0.05, 0.1, 0.4, 0.7, 0.4, 0.2, 0.15, 0.05, 0.02, 0.02, 0.01), group = 2),
                # moultP <- c(0.05, 0.05, 0.1, 0.4, 0.7, 0.4, 0.2, 0.15, 0.05, 0.02, 0.02, 0.01)
                color = "red", linetype = "dashed") +
      labs(x = "Month", y = "p", title = "Violin Plot of moult prob (red shows expected)") +
      theme_minimal()
    
    #--> Other params ####
    fit$growth_mu %>% hist(main="growth_mu", breaks=50)
    # abline(v=5, col="red", lwd=3)
    fit$growth_sigma %>% hist(main="growth_sigma", breaks=50)
    # abline(v=2, col="red", lwd=3)
    fit$nogrowth_sigma %>% hist(main="nogrowth_sigma", breaks=50)
    # abline(v=1, col="red", lwd=3)
  }
  #####
  
  #--> SAVED OUTPUTS ####
  ## Specify the folder name
  folder_name <- "outputs"
  
  ## Check if the folder exists, and create it if it doesn't
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
    cat("Folder 'outputs' created.\n")
  } else {
    cat("Folder 'outputs' already exists.\n")
  }
  
  samples <- ggmcmc::ggs(fit_mcmc)
  
  ##--> Write out ####
  gc()
  
  ##Save FIT
  fit <- rstan::extract(fit_mcmc)
  
  ##Condition for this file to save as size can be large
  if(mods <= length(ldata)){
    saveRDS(fit_mcmc, file = paste("outputs/out_moult", sdata$name, "FITMCMC.rds", sep="_"))
  }
  
  saveRDS(fit, file = paste("outputs/out_moult", sdata$name, "FIT.rds", sep="_"))
  
  ##Save DATA
  saveRDS(sdata, file = paste("outputs/out_moult", sdata$name, "DATA.rds", sep="_"))
  
  ##Save Summary
  out_summ <- capture.output(print(fit_mcmc))
  cat(out_summ, file=paste("outputs/out_moult", sdata$name, "SUMMARY.txt", sep="_"),sep="\n",append=TRUE)
  

  ##Save Histograms
  samples <- ggmcmc::ggs(fit_mcmc)
  samples_filtered <- samples
  
  ##Save PLOTS
  png(paste("outputs/out_moult", sdata$name, "HIST.png", sep="_"))
  plot(
    ggplot(samples_filtered, aes(x=value, color=Parameter, fill=Parameter)) +
      geom_histogram(alpha=0.6) +
      theme(
        legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)
      ) +
      xlab("Moult Probability") +
      ylab("Months") +
      facet_wrap(~Parameter)
  )
  dev.off()
  
  png(paste("outputs/out_moult", sdata$name, "TRACEPLOT.png", sep="_"))
  plot(
    ggplot(samples_filtered, aes(x=Iteration, y=value, color=as.factor(Chain), fill=Parameter)) +
      geom_line(alpha=0.6) +
      theme(
        legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)
      ) +
      xlab("Mixing") +
      ylab("Months") +
      facet_wrap(~Parameter)
  )
  dev.off()
  #####
}

#####

################################################
################################################
################################################

