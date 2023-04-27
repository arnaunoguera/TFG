#Running file with necessary functions
source(file = 'C:/Users/Usuari/OneDrive - Universitat Autònoma de Barcelona/TFG/Bases de dades/data/GlobalFunctions.R')

Data <- read.table('data.txt', header = TRUE, sep = '\t')

#I get the granges object to convert the coordinates (hg19->hg38)
granges19 <- GRanges(seqnames = paste0('chr', Data$Chr), ranges = IRanges(start= Data$Sw.Start, end=Data$Sw.End))
export.bed(granges19, 'hg19coordinates.bed')

#I get the transformed coordinates, All were successful
granges38 <- import.bed('hg38coordinates.bed')
startvector <- as.data.table(granges38)$start
endvector <- as.data.table(granges38)$end

Out <- Data %>% mutate(ID = rownames(Data), Start = startvector, End = endvector, Pops = 'earlyEUR', Metapops = 'EUR', 
                       Stats = 'SF2', Evidence = 'DFS', Source = 'Souilmi et al.', GenesAnnotation = NA) %>%
  select(ID, Chr, Start, End, Pops, Metapops, Stats, Evidence, Source, GenesAnnotation)

write.csv(Out, '../../ProcessedData/8_Souilmi.csv',row.names = FALSE)
