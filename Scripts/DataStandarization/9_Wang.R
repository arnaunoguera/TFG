#Running file with necessary functions
source(file = 'C:/Users/Usuari/OneDrive - Universitat Autònoma de Barcelona/TFG/Bases de dades/data/GlobalFunctions.R')

#I get the data. The input, in this case, is just a list of genes where selection was detected. 
Genes <- read.table('Data.txt', header = FALSE, col.names = c('Gene'))

#I connect to ensembl and get the data
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
Annotation <- getBM(attributes=c('hgnc_symbol','chromosome_name', 'start_position', 'end_position', 'strand'),
      filters=c('hgnc_symbol'),
      values=list(Genes$Gene),
      mart=ensembl)

#2 records are duplicated, because ensembl has mapped them onto a non-canon chromosome too. I will remove them
#Checking for duplicates:
Annotation[duplicated(Annotation$hgnc_symbol),]
#removing the 2 particular records:
Annotation <- Annotation %>% filter(chromosome_name != 'CHR_HG2114_PATCH', chromosome_name != 'CHR_HG142_HG150_NOVEL_TEST') 

#I create the final table
Out <- Annotation %>% dplyr::rename(Chr = chromosome_name, Start = start_position, End = end_position, 
                                    GenesAnnotation = hgnc_symbol) %>%
  mutate(ID = rownames(Annotation), Pops = 'CEU,CHB,JPT,YRI', Metapops = 'EUR,EAS,AFR', Stats = 'LDD',
         Evidence = 'LD', Source = 'Wang et al.') %>%
  dplyr::select(ID, Chr, Start, End, Pops, Metapops, Stats, Evidence, Source, GenesAnnotation)

write.csv(Out, '../../ProcessedData/9_Wang.csv',row.names = FALSE)
