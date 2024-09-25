#############FLMB Stocking Model##############
library(tidyverse)


#input parameters
burn_in<-20#let the model settle in before stocking starts for this many years
nfish<-500#number of starting superindividuals
nstock<-7#change this to reflect stocking rate
years_stocked<-8#number of years stocked
years_post_stock<-5#number of years after stocking to continue running model
years<-burn_in+years_stocked+years_post_stock#add together all years modeled
maxfish<-1000#max starting number fish within a superindividual
imms<-2#number of immigrants
perc_min<-0#min starting %FLMB
perc_max<-12#max starting #FLMB
n_sims<-10#number of times to replicate the simulation

#stock_regime<-burn_in+c(1:9,14:22)#stock every other year- change based on years stocked if needed
stock_regime<-burn_in+(1:years_stocked)#stock every year- use this instead of above if stocking every year


#results array for plotting and visualizing model results
results<-data.frame(matrix(ncol=10, nrow=years*n_sims));colnames(results)<-c("Year","nFish","nadult_inds","perc_FLMB","nnlmb","nhybrids","nF1","nFxF","nFLMB","nstock")

######Start the model######
for(t in 1:n_sims){  #repeat simulation multiple times in order to produce an average

  #create data
  Fish<-data.frame(matrix(ncol=5,nrow=nfish))
  colnames(Fish)<-c("Age", "Per_FLMB","n_inds","sex","fertile")
  Fish[,1]<-round(runif(nfish, 0, 10))#age between 0-10
  Fish[,2]<-round(runif(nfish,perc_min,perc_max))#Random uniform distribution of %FLMB genes 
  Fish[,3]<-round(runif(nfish, maxfish/2, maxfish))#start with random number of individuals between max and max/2
  Fish[,4]<-sample(c(0,1),replace=TRUE, size=nfish)#0=male, 1=female
  Fish[,5]<-0#0 is not fertile, 1 is fertile- set during "reproduction" below
  
  #immigrant and stocked fish matrices
  Stockers<-data.frame(matrix(ncol=5,nrow=nstock));colnames(Stockers)<-c("Age","Per_FLMB","n_inds","sex", "fertile")
  immigrant<-data.frame(matrix(ncol=5,nrow=imms));colnames(immigrant)<-c("Age","Per_FLMB","n_inds","sex", "fertile")
  
  
  for (j in 1:years){ #Number of years for the simulation
    
    
    #Stocked fish data-recreate each year
    Stockers[,1]<-0; Stockers[,2]<-100; Stockers[,3]<-maxfish; Stockers[,4]<-sample(c(0,1),replace=TRUE, size=nstock);Stockers[,5]<-0
    #immigrant fish data
    immigrant[,1]<-0; immigrant[,2]<-0; immigrant[,3]<-maxfish; immigrant[,4]<-sample(c(0,1),replace=TRUE, size=imms);immigrant[,5]<-0
    
    #set up reproduction
    newoffspring<-data.frame(matrix(ncol=5,nrow=nrow(Fish)));colnames(newoffspring)<-c("Age","Per_FLMB","n_inds","sex","fertile") 
    #Reproduction
    for (i in 1:nrow(Fish)){ 
      if(Fish[i,1]>=3){Fish[i,5]<-1}
      if(Fish[i,1]==2){Fish[i,5]<-round(runif(1,0,1))}#all age3+ set as fertile, 50% of age2 are fertile for both male and female
    }
    fspawners<-filter(Fish, fertile==1 & sex==0)#select fertile females
    
    for (i in 1:nrow(Fish)){ 
      #each is spawned with a random other fish to get the %FLMB, a single offspring results)
      if(Fish[i,5]==1 & Fish[i,4]==0){#proceed if they are a fertile male
        randfish<-sample_n(fspawners, 1)#pick a random adult female to mate with
        newoffspring[i,2]<-(Fish[i,2]+randfish[1,2])/2#take mean of the %FLMB  
        newoffspring[i,1]<-0#age=0
        newoffspring[i,3]<-round((Fish[i,3]*(Fish[i,1]-1))+(100*(Fish[i,1]-1))); if(newoffspring[i,3]>maxfish){newoffspring[i,3]<-maxfish}#older ages produce more per individual fish
        newoffspring[i,4]<-sample(0:1,1)#randomly make male or female
        newoffspring[i,5]<-0
      }}
    newoffspring<-na.omit(newoffspring)#remove na's from offspring matrix
    if(j %in% stock_regime){if(nrow(newoffspring)>148-nstock){newoffspring<-sample_n(newoffspring,148-nstock)}}#during stocking years, limit the number of fish lower so that total with stocking = 150
    if(nrow(newoffspring)>148){newoffspring<-sample_n(newoffspring,148)}#Randomly select 148 "successful" spawns
    if(j<=burn_in){newoffspring[,2]<-runif(nrow(newoffspring),perc_min, perc_max)}#during burn in, offspring get random %FLB alleles
    
    
    #mortality - is even across age groups
    surv<-rnorm(nrow(Fish),0.6,0.1)#survival random normal distribution
    Fish[,3]<-round(Fish[,3]*surv)#take out some n_inds based on surv
    Fish<-Fish[!Fish$n_inds<=5,]#if n_inds <=5 the fish is removed from the array
    Fish[,1]<-Fish[,1]+1#add a year to age  
    Fish<-Fish[!Fish$Age==13,]#if age=13, take out that row
    
    
    #Add offspring to the array
    Fish<-bind_rows(Fish,newoffspring)
    
    #add stocked fish to the array (after the burn in period and before the post_stocked years, or just in stock regime years)
    if(j %in% stock_regime){  #use this to match to stock regime
      #if(j > burn_in & j < years_stocked+burn_in){  #use this to base off of stocking years
      Fish<-bind_rows(Fish,Stockers)  
    }
    
    #add immigrants to the array
    Fish<-bind_rows(Fish, immigrant)
    
    #save results
    results[((years*t)-(years)+j),1]<-j-burn_in
    results[((years*t)-(years)+j),2]<-nrow(Fish)#save n Fish
    results[((years*t)-(years)+j),3]<-sum(Fish$n_inds)#n individuals
    results[((years*t)-(years)+j),4]<-mean(Fish[,2])#save mean %FLMB for the year
    results[((years*t)-(years)+j),5]<-nrow(filter(Fish,Per_FLMB<=5))#NLMB
    results[((years*t)-(years)+j),6]<-nrow(filter(Fish,Per_FLMB>5 & Per_FLMB<50))#FxNLMB
    results[((years*t)-(years)+j),7]<-0#nrow(filter(Fish,Per_FLMB>=45 & Per_FLMB<55))#F1
    results[((years*t)-(years)+j),8]<-nrow(filter(Fish,Per_FLMB>=50 & Per_FLMB<95))#FxFLMB
    results[((years*t)-(years)+j),9]<-nrow(filter(Fish,Per_FLMB>=95))#save number of FLMB
    results[((years*t)-(years)+j),10]<-nstock#keep track of simulations with multiple stocking levels
    results[((years*t)-(years)+j),11]<-t#Simulation # for replications
  }  }#####End of simulation
