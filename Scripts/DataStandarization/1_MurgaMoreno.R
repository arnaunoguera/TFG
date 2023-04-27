#Running file with necessary functions
source(file = 'C:/Users/Usuari/OneDrive - Universitat Autònoma de Barcelona/TFG/Bases de dades/data/GlobalFunctions.R')

#Importing PopHumanScan dataset
PHS <- read.csv("1_rawPophumanscanTable.tab", header = TRUE, sep = "\t")
#Coordinates are in hg19

#Obtaining metapopulations from list of populations (vectorized function)
GetMetapopsPHS <- function(popsvector){
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

#Obtaining types of evidence from the stats column (and label column, to get significant alpha)
#Alpha isn't written in the stat column when it is significant (there is evidence of Protein Changes)
#Instead, it's in the label column only
GetEvidencePHS <- function(stats, label){
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
    #Now we have to see if there is evidence of protein changes (significant alpha), checking the label column
    if (grepl('Alpha', label[i], fixed = TRUE)){
      result <- append(result, 'PC') #IF that's the case, I add the 5th evidence type (PC)
    }
    #I add the result for the iteration (separating each evidence by a comma) to the final result vector
    resultvector <- append(resultvector, paste(result, collapse=','))
  }
  return(resultvector)
}

#I have to convert the coordinates from hg19 to hg38. I will do it using Liftover online.
#I need to export the coordinates to a BED file
granges19 <- GRanges(seqnames = PHS$chr, ranges = IRanges(start = PHS$start, end = PHS$end)) #Genomic ranges object
#Note: chromosomes have to be as 'chr7', not '7'
export.bed(object = granges19, con = 'coordinatesHg19.bed')

#I do the conversion at http://genome.ucsc.edu/cgi-bin/hgLiftOver. Importing the results: 
granges38 <- import.bed('coordinatesHg38.bed')

#One record couldn't be converted, so I will eliminate it manually before replacing the coordinates. It's chr2	172742418	172880598
PHS <- PHS %>% filter(chr != 'chr2' | start != 172742419) #Coordinate has to be +1 because it UCSC analyzes it as 0-based
#I'll convert the GRanges object into a table
granges38table <- as.data.table(granges38)
#Getting the start and end vectors. These 2 will be as long as the data table, once the record that wasn't transformed is removed
startvector <- granges38table$start
endvector <- granges38table$end

#I get the standarized table I use for all datasets
Out <- PHS %>% dplyr::rename(Pops = pops, Stats = stat) %>%
  separate(chr, sep = "chr", into = c(NA, "Chr")) %>%
  mutate(ID = rownames(PHS), Metapops = GetMetapopsPHS(Pops), Evidence = GetEvidencePHS(Stats, label),
         Source = 'PopHumanScan', GenesAnnotation = NA, Start = startvector, End = endvector) %>%
  select(ID, Chr, Start, End, Pops, Metapops, Stats, Evidence, Source, GenesAnnotation)

#I add 'Alpha' to the Stats of those records with Protein Changes (because it wasn't in the PHS table)
Out[grepl('PC',Out$Evidence),'Stats'] <-  paste0(Out[grepl('PC',Out$Evidence),'Stats'], ',Alpha')

write.csv(Out, '../../ProcessedData/1_PHS.csv',row.names = FALSE)

