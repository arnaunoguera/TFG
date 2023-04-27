source('C:/Users/Usuari/OneDrive - Universitat Autònoma de Barcelona/TFG/Bases de dades/data/GlobalFunctions.R')

#I import the data and ignore those regions where balancing selection has been detected (I'm studying positive selection)
InitialData <- read.table('data.csv', header = TRUE, sep = ';') 
#I also filter a record of mitochondrial DNA and 5 that contain different locations together
InitialData <- InitialData %>% filter(!grepl('multiple', InitialData$Chromosomal.position..HG16.),
                                      Type.of.selection != 'Balancing', 
                                      Chromosomal.position..HG16. != 'mtDNA')

#The regions are in hg16. I save them in a file to convert them
write(InitialData$Chromosomal.position..HG16., file = 'hg16coordinates.txt')

#Conversion to hg19 fails on 2 records: 
  #Split in new: chr1:25359819-25419545
  #Partially deleted in new: chr14:105161971-106270487
#Further conversion to hg38 fails on 1 record: 
  #Split in new: chr19:54799856-54804238. This corresponds to chr19:59491668-59496050 in hg16
#I filter the records that couldn't be transformed to hg38
InitialData <- InitialData %>% filter(Chromosomal.position..HG16. != 'chr1:25359819-25419545', 
                                      Chromosomal.position..HG16. != 'chr14:105161971-106270487',
                                      Chromosomal.position..HG16. != 'chr19:59491668-59496050')
#I import the new coordinates. Now they correspond to my remaining records
coordinates38 <- read.table('hg38coordinates.bed', header=FALSE, sep = ':', col.names = c('chr', 'region')) %>% 
  separate(region, into = c('Start', 'End'), sep = '-', convert = TRUE) %>% separate(chr, into = c(NA, 'Chr'), sep = 'chr')

#I get the final table
Out <- InitialData %>% dplyr::rename(Stats = Tests) %>% 
  mutate(ID = rownames(InitialData), Chr = coordinates38$Chr, Start = coordinates38$Start, End = coordinates38$End,
         Pops = NA, Metapops=NA, Evidence = GetEvidence(Stats), Source = 'Sabeti el al. (Multiple)',
         GenesAnnotation = NA) %>%
  dplyr::select(ID, Chr, Start, End, Pops, Metapops, Stats, Evidence, Source, GenesAnnotation)

write.csv(Out, '../../ProcessedData/4_Sabeti.csv',row.names = FALSE)
