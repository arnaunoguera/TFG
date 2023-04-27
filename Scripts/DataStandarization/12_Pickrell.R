#Running file with necessary functions
source(file = 'C:/Users/Usuari/OneDrive - Universitat Autònoma de Barcelona/TFG/Bases de dades/data/GlobalFunctions.R')

#The data is in 16 files extracted from the supplementary materials pdf, separated by test (iHS or XP-EHH) and population
tests <- c('iHS', 'XPEHH')
metapops <- c('AMR', 'BPY', 'BSP', 'EAS', 'EUR', 'MDE', 'OCN', 'SAS')

#I use a loop to access all the files
Data <- data.table()
for (test in tests){
  for (metapop in metapops){ #Iteration for every metapop. BPY and BSP are both pops from AFR
    #I get the path
    path <- paste0(test, '_', metapop, '.txt')
    if (metapop %in% c('BPY', 'BSP')){
      pop <- metapop #I put BPY or BSP as the population, and the metapopulation is AFR
      metapop <- 'AFR'
    } else{
      pop <- NA #There is no population for these
    }
    #I add the new table to Data
    Data <- rbind(Data, read.table(path, header = FALSE, sep = ' ', col.names = c('chr', 'range', NA), fill = TRUE) %>%
                    mutate(Metapops = metapop, Pops = pop, Stats = test))
  }
}
#I convert the coordinates to bases (they are in Mb) and fix the chromosomes (chr7 -> 7)
Data <- Data %>% separate(range, into = c('startMb', 'endMb'), sep = ':', remove = TRUE, convert = TRUE) %>%
  mutate(start = startMb * 1000000, end = endMb * 1000000, Chr = str_replace(chr, 'chr', ''))
  
#I get the gRanges object to convert the coordinates hg17->hg19->hg38
granges17 <- GRanges(seqnames = Data$chr, ranges = IRanges(start = Data$start, end = Data$end))
export.bed(granges17, 'coordinateshg17.bed')

#24 records failed conversion to hg19
#3 failed conversion to hg38:
  #Split in new: chr19	7248999	7694000	.	0	. This corresponds to 7200000-7600000
  #Split in new: chr19	7248999	7694000	.	0	.
  #Split in new: chr11	87322351	87760352	.	0	. This corresponds to 87000000-87400000
#I saved the 27 failed records in a file
failed <- read.table(file = 'Failedhg17.txt', header = FALSE, sep = '\t', comment.char = '#',
                     col.names = c('chr', 'start_1', 'end', NA, NA, NA))
#I add a column to Data to filter those in failed
Data <- Data %>% mutate(failed = FALSE)
for (i in 1:nrow(Data)){
  pos_start <- round(Data[i, 'start']) #For some reason, 3 failed records aren't recognized without round
  pos_end <- round(Data[i, 'end'])
  chrom <- Data[i,'chr']
  for (n in 1:nrow(failed)){
    if (pos_start == failed[n,'start_1']+1 & pos_end == failed[n,'end'] & chrom == failed[n,'chr']){
      Data[i, 'failed'] <- TRUE
      break
    }
  }
}
#I filter these records
Data <- Data %>% filter(failed == FALSE)

#I get the converted coordinates
granges38 <- import.bed('coordinateshg38.bed')
startvector <- as.data.table(granges38)$start
endvector <- as.data.table(granges38)$end

Out <- Data %>% mutate(ID = rownames(Data), Start = startvector, End = endvector, Evidence = 'LD', 
                       Source = 'Pickrell et al.', GenesAnnotation = NA) %>% 
  dplyr::select(ID, Chr, Start, End, Pops, Metapops, Stats, Evidence, Source, GenesAnnotation)

write.csv(Out, '../../ProcessedData/12_Pickrell.csv',row.names = FALSE)
