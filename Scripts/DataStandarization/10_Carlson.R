#Running file with necessary functions
source(file = 'C:/Users/Usuari/OneDrive - Universitat Autònoma de Barcelona/TFG/Bases de dades/data/GlobalFunctions.R')

#I get the data. I say '/' is the comment character to ignore all that follows (since it's heterogeneous)
AFR <- read.table('AFR.txt', header = FALSE, sep = ' ', comment.char = '/', col.names = c('chrom', 'range', NA)) %>%
  mutate(Pops = 'ASW', Metapops = 'AFR')
EUR <- read.table('EUR.txt', header = FALSE, sep = ' ', comment.char = '/', col.names = c('chrom', 'range', NA)) %>%
  mutate(Pops = 'CEU', Metapops = 'EUR')
EAS <- read.table('EAS.txt', header = FALSE, sep = ' ', comment.char = '/', col.names = c('chrom', 'range', NA)) %>%
  mutate(Pops = 'CHD', Metapops = 'EAS')

#I join the tables, fix the chr and separate the range
Data <- rbind(AFR, rbind(EUR, EAS)) %>% mutate(chr = str_replace(chrom, ':', '')) %>%
  mutate(Chr = str_replace(chr, 'chr', '')) %>%
  separate(range, into = c('starthg17', 'endhg17'), sep = '–', remove = TRUE, convert = TRUE)

#I save the ranges in a gRanges object to transform the coordinates
granges17 <- GRanges(seqnames = Data$chr, ranges = IRanges(start = Data$starthg17, end = Data$endhg17))
export.bed(granges17, 'hg17coordinates.bed')

#Conversion to hg19 failed in 1 record:
  #Partially deleted in new: chr15	82779999	83180000	.	0	.
#Conversion to hg38 successful on the remaining 58 records
#I filter out the failed record
Data <- Data %>% filter(Chr != '15' | endhg17 != 83180000)

#I import the converted coordidnates
granges38 <- import.bed('hg38coordinates.bed')
startvector <- as.data.table(granges38)$start
endvector <- as.data.table(granges38)$end

#I create the final table
Out <- Data %>% mutate(ID = rownames(Data), Start = startvector, End = endvector, Stats = 'Tajima_D', 
                       Evidence = 'DFS', Source = 'Carlson et al.', GenesAnnotation = NA) %>%
  dplyr::select(ID, Chr, Start, End, Pops, Metapops, Stats, Evidence, Source, GenesAnnotation)

write.csv(Out, '../../ProcessedData/10_Carlson.csv',row.names = FALSE)
