#Packages
library(data.table)
library(dplyr)
library(tibble)
library(tidyverse)
library(rtracklayer)
library(biomaRt)

#Function to get the metapopulation a population belongs to
PopToMetapop <- function(pop){
  #pops by metapop
  AFRpops<-c('ACB','ASW','ESN','GWD','LWK','MSL','YRI', 'MKK')
  EASpops<-c('CDX','CHB','CHS','JPT','KHV', 'CHD')
  EURpops<-c('CEU','FIN','GBR','IBS','TSI')
  SASpops<-c('BEB','GIH','ITU','PJL','STU')
  AMRpops<-c('CLM','MXL','PEL','PUR', 'MEX')
  #Getting the metapop
  metapop <- case_when(pop %in% AFRpops ~ 'AFR',
                       pop %in% EURpops ~ 'EUR',
                       pop %in% EASpops ~ 'EAS',
                       pop %in% SASpops ~ 'SAS',
                       pop %in% AMRpops ~ 'AMR')
  if (!is.na(metapop)){
    return(metapop)
  }
}

#Function to get the evidence type (LD, SFS or PC) based on the stat
StatToEvidence <- function(stat){
  LD <- c('iHS', 'XPEHH', 'Rsb', 'XPEHH/Rsb', 'LRH', 'Hudson_HapTest', 'DepaulisVeuille_K_H', 'LD', 'HS', 'nSL')
  DFS <- c('FuLi_D','FuLi_F','Tajima_D', 'MAFthreshold', 'EwensWatterson_F', 'Fu_Fs', 'HKA', 'Heterozygosity', 'SF2')
  HFDA <- c('FayWu_H')
  PD <- c('Fst', 'GeographicCorrelation', 'Pexcess', 'Snn')
  PC <- c('Ka/Ks', 'dN/dS','PAML', 'MK', 'Bn/Bs', 'Pa/Ps' , 'PRFModel', 'pNC/pNR')
  #Getting the evidence type
  ev <- case_when(stat %in% LD ~ 'LD',
                  stat %in% DFS ~ 'DFS',
                  stat %in% HFDA ~ 'HFDA',
                  stat %in% PD ~ 'PD',
                  stat %in% PC ~ 'PC')
  if (!is.na(ev)){
    return(ev)
  }
}


#Obtaining metapopulations from list of populations (vectorized function)
GetMetapops <- function(popsvector){
  resultvector <- c()
  for (i in 1:length(popsvector)){ #For each row
    popsi <- str_split_1(popsvector[i], ',') #they are separated by commas
    result <- c() #Result for each record
    for (pop in popsi){ #For each population in a record, I check what metapopulation it belongs to
      metapop <- PopToMetapop(pop) #I get the metapop 
      if (!metapop %in% result){ #I check if that metapop is already in results for that record. 
        result <- append(result, metapop) #If not, I add it
      }
    }
    #I add the result for the iteration (separating each metapop by a comma) to the final result vector
    resultvector <- append(resultvector, paste(result, collapse=','))
  }
  return(resultvector)
}

#Obtaining the list of evidence
GetEvidence <- function(stats){
  resultvector <- c() 
  for (i in 1:length(stats)){ #For each row 
    statsi <- str_split_1(stats[i], ',') #Stats are separated by commas
    result <- c() #Result for each record
    for (stat in statsi){ #Checking each stat
      evidence <- StatToEvidence(stat)
      if (!evidence %in% result){
        result <- append(result, evidence)
      }
    }
    #I add the result for the iteration (separating each evidence by a comma) to the final result vector
    resultvector <- append(resultvector, paste(result, collapse=','))
  }
  return(resultvector)
}