#Running file with necessary functions
source(file = 'C:/Users/Usuari/OneDrive - Universitat Autònoma de Barcelona/TFG/Bases de dades/data/GlobalFunctions.R')

#Getting the Data:
AFR <- read.table('ASW.txt', header = FALSE, sep = '\t', col.names = c('Chr', 'pos', 'DFS')) %>%
  mutate(Metapops = 'AFR', Pops = 'ASW')
EUR <- read.table('CEU.txt', header = FALSE, sep = '\t', col.names = c('Chr', 'pos', 'DFS')) %>%
  mutate(Metapops = 'EUR', Pops = 'CEU')
EAS <- read.table('CHB.txt', header = FALSE, sep = '\t', col.names = c('Chr', 'pos', 'DFS')) %>%
  mutate(Metapops = 'EAS', Pops = 'CHB')
Combined <- read.table('Combined.txt', header = FALSE, sep = '\t', col.names = c('Chr', 'pos', 'DFS')) %>%
  mutate(Metapops = 'AFR,EUR,EAS', Pops = 'ASW,CEU,CHB')

#Joining the data and making the position into a range of 10kb
Data <- rbind(rbind(AFR,EUR), rbind(EAS,Combined)) %>% mutate(starthg16 = pos -5000, endhg16 = pos+5000)

#I create the gRanges object to transform the coordinates
granges16 <- GRanges(seqnames = paste0('chr',Data$Chr), ranges = IRanges(start = Data$starthg16, end = Data$endhg16))
export.bed(granges16, con = 'hg16coordinated.bed')

#Conversion was successful in all 181 records. I import them
granges38 <- import.bed(con = 'hg38coordinated.bed')
startvector <- as.data.table(granges38)$start
endvector <- as.data.table(granges38)$end

Out <- Data %>% mutate(ID = rownames(Data), Start = startvector, End = endvector, Stats = 'CLR', Evidence = 'DFS',
                       Source = 'Williamson et al.', GenesAnnotation = NA) %>% 
  dplyr::select(ID, Chr, Start, End, Pops, Metapops, Stats, Evidence, Source, GenesAnnotation)

write.csv(Out, '../../ProcessedData/11_Williamson.csv',row.names = FALSE)
