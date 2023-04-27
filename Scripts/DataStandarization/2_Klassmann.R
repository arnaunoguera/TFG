source('C:/Users/Usuari/OneDrive - Universitat Autònoma de Barcelona/TFG/Bases de dades/data/GlobalFunctions.R')

#There's several txt files with the detected regions, which contain the tables copied from the original pdf.
#Columns are separated by a space
#The first step is to import all the different tables 

#The pops object will be used to acces the different directories. The first 6 have to be combines with iHS (i.e. '2_iHS_CEU.txt')
#The others have to be combined with XPEHHRsb: XPEHHRsb_JPTvsYRI
pops <- c('CEU', 'CHB', 'JPT', 'YRI', 'CEUGBR', 'CHBCHS',
          'CEUvsCHB', 'CEUvsJPT', 'CEUvsYRI', 'CHBvsJPT', 'CHBvsYRI', 'JPTvsYRI', 'CEUGBRvsCHBCHS')

regions <- list() #List where I will store all the tables

#Loop to access all the different 
stat <- 'iHS' #Initially, the stat is iHS (to be changed during the iterations)
statstr <- 'iHS'
for (i in 1:length(pops)){
  #Just to annotate correctly the populations in the table: 
  pop <- c(pops[i])
  if (pop == 'CEUGBR' | pop == 'CHBCHS'){
    pop <- c(substr(pop,1,3), substr(pop,4,6)) #Now pop will be CEU,GBR or CHB, CHS
  } #I have to differentiate the first 6 cases (iHS) from the last (XPEHH/Rsb)
  else if (i > 6){ #This will be true for the records pertaining to XPEHH/Rsb
    stat <- 'XPEHHRsb' #This is for the path
    statstr <- 'XPEHH/Rsb' #This is for the resulting table
    if (pop != 'CEUGBRvsCHBCHS'){
      pop <- str_split_1(pop, 'vs')
    } else {
      pop <- c('CEU','GBR','CHB','CHS')
    }
  }
  path <- paste0('2_', stat, '_', pops[i], '.txt') #the path to the file
  regions[[i]] <- read.table(path, header = FALSE, sep = ' ', skip = 1, #I skip 1 line to skip the heather (it pastes weird from the pdf).
                             col.names = c('id', 'chr', 'start', 'end', 'length', 'mrk', 'mean', 'max',
                                           'extr_mrk', 'percent')) %>% #And then I add column names manually)
    mutate(Pop = paste0(pop, collapse = ','), Stats = statstr)
}#The result is a list with 13 elements, corresponding to the 13 tables 

#I obtain a table joining all the separate tables
BigTable <- data.table()
for (regiontable in regions){
  BigTable <- rbind(BigTable, regiontable)
}

#Obtaining metapopulations from list of populations (vectorized function)
GetMetapops <- function(popsvector){
  resultvector <- c()
  for (i in 1:length(popsvector)){ #For each row
    popsi <- str_split_1(popsvector[i], ',') #they are separated by commas
    result <- c() #Result for each record
    for (pop in popsi){ #For each population in a record, I check what metapopulation it belongs to
      if (grepl('2', pop, fixed = TRUE)){ #If pop contains a '2', that's a population differentation test (i.e. FIN2YRI).
        #I want this to return 'EUR2AFR
        pops <- str_split_1(pop, '2') #Getting the 2 pops involved: 'FIN', 'YRI'
        metapops <- c()
        metapops[1] <- PopToMetapop(pops[1]) #1st metapop
        metapops[2] <- PopToMetapop(pops[2]) #2nd metapop
        metapop <- paste(metapops, collapse='2') #Saving both metapops in str metapop, separated by a '2'
      }
      else {
        metapop <- PopToMetapop(pop) #If it contains no '2', we get the metapop directly
      }
      if (!metapop %in% result){ #I check if that metapop is already in results for that record. 
        result <- append(result, metapop) #If not, I add it
      }
    }
    #I add the result for the iteration (separating each metapop by a comma) to the final result vector
    resultvector <- append(resultvector, paste(result, collapse=','))
  }
  return(resultvector)
}
  

#For the final table, I need to fix the IDs, convert the starts and ends to bases (it's in Mb), get the metapops, and add the evidence and source
Out <- BigTable %>% rename(Pops = Pop, Chr = chr) %>% 
  mutate(ID = rownames(BigTable), Start = start*1000000, End = end*1000000, Metapops = GetMetapops(Pops), 
         Evidence = 'LD', Source = 'Klassmann et al.', GenesAnnotation = NA) %>% 
  select(ID, Chr, Start, End, Pops, Metapops, Stats, Evidence, Source, GenesAnnotation)

#Saving the resulting table to a file
write.csv(Out, '../../ProcessedData/2_KlassmannEtAl.csv',row.names = FALSE)