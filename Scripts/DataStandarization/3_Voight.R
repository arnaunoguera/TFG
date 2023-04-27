source('C:/Users/Usuari/OneDrive - Universitat Autònoma de Barcelona/TFG/Bases de dades/data/GlobalFunctions.R')

#There are 3 tables, 1 for CEU, 1 for YRI and 1 for CHB and JPT together
metapops <- c('EAS','EUR','AFR')
pops <- list(c('JPT','CHB'), c('CEU'), c('YRI'))

#Table to add all tables
JointTable <- data.table()
#Loop to access all the files
for (n in 1:3){
  metapop <- metapops[n]
  pop <- pops[[n]]
  path <- paste0("3_", metapop, "pbio.0040072.sd00", as.character(n), '.txt') #Name of the files
  JointTable <- rbind(JointTable, read.table(path, header=FALSE, sep='\t', comment.char = 'r', #Comment char = r to ignore all inconsistent columns after rs....
                                             col.names = c('Chr', 'start', 'end', 'iHS', 'V1', 'V2', 'V3', 'V4', 'V5')) %>% #I don't know the meaning of all columns, but I have the necessary ones
                        mutate(Pops = paste0(pop, collapse=','), Metapops = metapop))
}

#The coordinates are in hg16 and need to be converted to hg38. Hg16 -> hg38 can't be done directly. It's hg16 -> hg19 -> hg38
#I need to export the coordinates to a BED file. Chromosomes have to be as 'chr7', not '7'
granges16 <- GRanges(seqnames = paste0('chr',JointTable$Chr),
                     ranges = IRanges(start = JointTable$start, end = JointTable$end)) #Genomic ranges object
export.bed(object = granges16, con = 'coordinatesHg16.bed')

#I do the conversions at http://genome.ucsc.edu/cgi-bin/hgLiftOver.
#I convert from hg16 to hg19. It fails on 4 records: 
  #Split in new chr2	126999999	127100000	.	0	.
  #Partially deleted in new chr8	48099999	48200000	.	0	.
  #Partially deleted in new chr15	26899999	27000000	.	0	.
  #Split in new chr1	2299999	2400000	.	0	
#I then do the conversion hg19 -> hg38. Failed on 1 record:
  #Partially deleted in new chr5	138723783	138871499	.	0	.
  #I checked manually and this corresponds to the range 138799999	138900000 in hg16

#Importing the results: 
granges38 <- import.bed('coordinatesHg38.bed')

#I filter the 5 failed regions from the JointTable object. +1 at start because it uses 0-based indexing
JointTable <- JointTable %>% filter((Chr != '2' | start != 126999999+1) & (Chr != '8' | start != 48099999+1) 
                                    & (Chr != '15' | start != 26899999+1) & (Chr != '1' | start != 2299999+1)
                                    & (Chr != '5' | start != 138799999+1))
#I'll convert the GRanges object into a table
granges38table <- as.data.table(granges38)
#Getting the start and end vectors. These 2 will be as long as the data table, once the record that wasn't transformed is removed
startvector <- granges38table$start
endvector <- granges38table$end

#Finish up the table:
Out <- JointTable %>% mutate(ID = rownames(JointTable), Stats = 'iHS', Evidence = 'LD',
                             Source = 'Voight et al.', GenesAnnotation = NA, 
                             Start = startvector, End = endvector) %>%
  select(ID, Chr, Start, End, Pops, Metapops, Stats, Evidence, Source, GenesAnnotation)

#Saving the resulting table to a file
write.csv(Out, '../../ProcessedData/3_VoightEtAl.csv',row.names = FALSE)
