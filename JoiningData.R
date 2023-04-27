source('../data/GlobalFunctions.R')

#Getting the table
Data <- rbind(read.csv('1_PHS.csv', colClasses = c('character', 'factor', 'numeric', 'numeric', rep('character', 4), 'factor')) %>% mutate(ID = paste0('01_', ID)), 
              rbind(read.csv('2_KlassmannEtAl.csv', colClasses = c('character', 'factor', 'numeric', 'numeric', rep('character', 4), 'factor')) %>% mutate(ID = paste0('02_', ID)), 
                    rbind(read.csv('3_VoightEtAl.csv', colClasses = c('character', 'factor', 'numeric', 'numeric', rep('character', 4), 'factor')) %>% mutate(ID = paste0('03_', ID)), 
                          rbind(read.csv('4_Sabeti.csv', colClasses = c('character', 'factor', 'numeric', 'numeric', rep('character', 4), 'factor')) %>% mutate(ID = paste0('04_', ID)), 
                                rbind(read.csv('5_JohnsonEtAl.csv', colClasses = c('character', 'factor', 'numeric', 'numeric', rep('character', 4), 'factor')) %>% mutate(ID = paste0('05_', ID)), 
                                      rbind(read.csv('6_nSL.csv', colClasses = c('character', 'factor', 'numeric', 'numeric', rep('character', 4), 'factor')) %>% mutate(ID = paste0('06_', ID)), 
                                            rbind(read.csv('7_Myles.csv', colClasses = c('character', 'factor', 'numeric', 'numeric', rep('character', 4), 'factor')) %>% mutate(ID = paste0('07_', ID)), 
                                                  rbind(read.csv('8_Souilmi.csv', colClasses = c('character', 'factor', 'numeric', 'numeric', rep('character', 4), 'factor')) %>% mutate(ID = paste0('08_', ID)), 
                                                        rbind(Data <- read.csv('9_Wang.csv', colClasses = c('character', 'factor', 'numeric', 'numeric', rep('character', 4), 'factor')) %>% mutate(ID = paste0('09_', ID)), 
                                                              rbind(read.csv('10_Carlson.csv', colClasses = c('character', 'factor', 'numeric', 'numeric', rep('character', 4), 'factor')) %>% mutate(ID = paste0('10_', ID)), 
                                                                    rbind(read.csv('11_Williamson.csv', colClasses = c('character', 'factor', 'numeric', 'numeric', rep('character', 4), 'factor')) %>% mutate(ID = paste0('11_', ID)), 
                                                                          read.csv('12_Pickrell.csv', colClasses = c('character', 'factor', 'numeric', 'numeric', rep('character', 4), 'factor')) %>% mutate(ID = paste0('12_', ID)))))))))))))
  
#Gene annotation
ensembl <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl') #I access genome build hg38
print('Connected to Ensembl')
for (n in 1:nrow(Data)) {
  chr <- Data[n, 'Chr']
  inici <- Data[n, 'Start']
  final <- Data[n, 'End']
  genes <- getBM(attributes = 'ensembl_gene_id', filters =  c('chromosome_name', 'start', 'end'), 
                 values = list(chr, inici, final), mart = ensembl) %>% pull(ensembl_gene_id)
  Data[n,"Genes"] <- paste(genes, collapse=',') 
  #With each iteration, I get the genes from 1 region and save them to the table, as a str separated by commas
  print(n)
}

write.csv(Data %>% dplyr::select(ID, Chr, Start, End, Pops, Metapops, Stats, Evidence, Source, Genes), 'AllData.csv',row.names = FALSE)
