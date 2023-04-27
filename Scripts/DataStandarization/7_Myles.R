#Running file with necessary functions
source(file = 'C:/Users/Usuari/OneDrive - Universitat Autònoma de Barcelona/TFG/Bases de dades/data/GlobalFunctions.R')

#There's 3 imput tables, and I add a column with their metapop:
AFR_CHB <- read.table(file = 'AFR_CHB.txt', header = TRUE, sep = ' ') %>% mutate(Metapops = 'AFR2EAS', Pops = 'ASW2CHB')
AFR_EUR <- read.table(file = 'AFR_EUR.txt', header = TRUE, sep = ' ') %>% mutate(Metapops = 'AFR2EUR', Pops = 'ASW2CEU')
EUR_CHB <- read.table(file = 'EUR_CHB.txt', header = TRUE, sep = ' ') %>% mutate(Metapops = 'EUR2EAS', Pops = 'CEU2CHB')

#I join the tables:
JointTable <- rbind(AFR_CHB, rbind(AFR_EUR, EUR_CHB))

#I create the gRAnges object to convert the coordinates (hg16->hg19->hg38) in LiftOver
granges16 <- GRanges(seqnames = paste0('chr', JointTable$chromosome), ranges = IRanges(start = JointTable$Start_position, 
                                                                                       end = JointTable$End_position))
export.bed(granges16, 'hg16coordinates.bed')

#Conversion to hg19 failed on 6 regions. 
  #Split in new: chr17	36379440	36868166	.	0	.
  #Split in new: chr2	96827767	97911239	.	0	.
  #Partially deleted in new: chr5	91624002	91967865	.	0	.
  #Partially deleted in new: chr8	42580305	47884542	.	0	.
  #Split in new: chr15	25885799	27121431	.	0	.
  #Split in new: chr17	62923123	63599553	.	0	.
#Conversion to hg38 successful in the remaining 1068 records
#I filter out the failed conversions. 
JointTable <- JointTable %>% filter(chromosome != 17 | End_position != 36868166,
                                    chromosome!= 2 | End_position != 97911239, 
                                    chromosome!= 5 | End_position != 91967865,
                                    chromosome!= 8 | End_position != 47884542,
                                    chromosome!= 15 | End_position != 27121431,
                                    chromosome!= 17 | End_position != 63599553)

granges38 <- import.bed('hg38coordinates.bed')
startvector <- as.data.table(granges38)$start
endvector <- as.data.table(granges38)$end

#The evidence in this study is differentation among populations
Out <- JointTable %>% dplyr::rename(Chr = chromosome) %>%
  mutate(ID = rownames(JointTable), Start = startvector, End = endvector, Stats = 'Fst', Evidence = 'PD', 
         Source = 'Myles et al.', GenesAnnotation = NA) %>%
  select(ID, Chr, Start, End, Pops, Metapops, Stats, Evidence, Source, GenesAnnotation)

write.csv(Out, '../../ProcessedData/7_Myles.csv',row.names = FALSE)
    
    

