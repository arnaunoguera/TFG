#Running file with necessary functions
source(file = 'C:/Users/Usuari/OneDrive - Universitat Autònoma de Barcelona/TFG/Bases de dades/data/GlobalFunctions.R')

#I read the table 
Data <- read.table('Data.csv', header = TRUE, sep = ';', na.strings = '-') %>% separate(chr, into = c(NA, 'Chr'),sep = 'chr',convert = TRUE)

#I save the coordinates to convert them to hg38 (they are in hg18)
granges18 <- GRanges(seqnames=paste0('chr',Data$Chr), ranges = Data$Possition.bp.)
export.bed(granges18, 'hg18coordinates.bed')

#I import the converted coordinates. All of them worked
granges38 <- import.bed('hg38coordinates.bed')
coordinates <- as.data.table(granges38)$start

#I get a mutated table with the new coordinates, and add a column for start and end (+5000 and -5000)
DataHg38 <- Data %>% mutate(Possition.bp. = coordinates, Start = coordinates - 5000, End = coordinates + 5000, overlap = FALSE)
#Now I check if any positions are already included in the 10kb intervals I created on another position (per population)
for (n in 1:nrow(DataHg38)){
  pop <- DataHg38$Population[n]
  pos <- DataHg38$Possition.bp.[n]
  comparison <- DataHg38[1:n-1,] %>% filter(Population == pop) #I compare it to all previous SNPs for that population
  if (nrow(comparison) == 0){
    next
  }
  IntervalsStart <- comparison$Start
  IntervalsEnd <- comparison$End
  for (i in 1:length(IntervalsStart)){ #I check if the position is in any of the intervals to compare with
    if (pos %in% IntervalsStart[i]:IntervalsEnd[i]){
      DataHg38[n,'overlap'] <- TRUE
      break
    }
  }
}
DataHg38 <- DataHg38 %>% filter(overlap == FALSE)

#I create the final table
Out <- DataHg38 %>% dplyr::rename(Pops = Population) %>% 
  mutate(ID = rownames(DataHg38), Chr = as.character(Chr), Metapops = GetMetapops(Pops), Stats = 'nSL', 
         Evidence = 'LD', Source = 'Ferrer-Admetlla et al.', GenesAnnotation = NA) %>% 
  select(ID, Chr, Start, End, Pops, Metapops, Stats, Evidence, Source, GenesAnnotation)

write.csv(Out, '../../ProcessedData/6_nSL.csv',row.names = FALSE)
