source('C:/Users/Usuari/OneDrive - Universitat Autònoma de Barcelona/TFG/Bases de dades/data/GlobalFunctions.R')

#Importing the table
Data <- read.csv(file = 'JohnsonData.csv', header = TRUE, colClasses = c(NA, 'character', rep(NA, 7)))
#I use colclasses to define the class of the CHR column to be str, so I can change chr 23 to X

#I change chr 23 to X for the coordinates change to work
Data[Data$CHR == '23', 'CHR'] <- 'X'

#I have to convert the coordinates from hg19 to hg38. I will do it using LiftOver online.
#I need to export the coordinates to a BED file
granges19 <- GRanges(seqnames = paste0('chr',Data$CHR), 
                     ranges = IRanges(start = Data$START_POS, end = Data$STOP_POS)) #Genomic ranges object
#Note: chromosomes have to be as 'chr7', not '7'
export.bed(object = granges19, con = 'coordinatesHg19.bed')

#I do the conversion at http://genome.ucsc.edu/cgi-bin/hgLiftOver. Importing the results: 
granges38 <- import.bed('coordinatesHg38.bed')
#I'll convert the GRanges object into a table
granges38table <- as.data.table(granges38)
#Getting the start and end vectors. These 2 will be as long as the data table, once the record that wasn't transformed is removed
startvector <- granges38table$start
endvector <- granges38table$end

#98 records failed conversion. As they are so many, I downloaded it as a file
FailedConversion <- read.table('FailedConversion.txt', header = FALSE, sep = '\t', col.names = c('chr', 'start1', 'end', rep(NA,3)),
                               comment.char = '#') %>%
  separate(col = chr, into = c(NA, "chr"), sep = 'chr', remove = TRUE) %>%
  mutate(start = start1 + 1) #I add 1 to account for 0-based indexing

#To delete the failed regions, I delete those where the start and end positions are found in the FailedConversion object
#It theoretically could delete records accidentally, but exactly 98 records are deleted, so I keep it that way
PrunedData <- Data %>% filter(!START_POS %in% FailedConversion$start | !STOP_POS %in% FailedConversion$end)

#I prepare the table to export
Out <- PrunedData %>% dplyr::rename(Pops = POP, Chr = CHR) %>%
  mutate(ID = rownames(PrunedData), Start = startvector, End = endvector, Metapops = GetMetapops(Pops), 
         Stats = 'iHS', Evidence = 'LD', Source = 'Johnson et al.', GenesAnnotation = NA) %>%
  select(ID, Chr, Start, End, Pops, Metapops, Stats, Evidence, Source, GenesAnnotation)

write.csv(Out, '../../ProcessedData/5_JohnsonEtAl.csv',row.names = FALSE)

