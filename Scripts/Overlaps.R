source('../data/GlobalFunctions.R')
#Importing the data
Data <- read.csv(file= 'AllData.csv', header = TRUE)
#I will study the overlap by genes. So, first, I want a vector with all the genes in the Data table
GenesVector <- unique(unlist(strsplit(Data$Genes, split = ',')))
#There are 27162 unique genes in the table
GenesData <- data.table(GeneID = GenesVector, Count = as.numeric(NA), Chr = as.character(NA), Pops = as.character(NA), 
                        Metapops = as.character(NA),Stats = as.character(NA), Evidence = as.character(NA), 
                        SourceIDs = as.character(NA))
#Function to summarize information for a gene based on a vector (ex: c('AFR', 'EUR,AFR', 'EAS') -> 'AFR,EUR,EAS')
summarizeInfo <- function(inputvector){
  elements <- unique(unlist(strsplit(inputvector, ',')))
  return(paste0(elements, collapse = ','))
}
#I summarize each gene
for (n in 1:nrow(GenesData)){
  geneID <- GenesData$GeneID[n]
  info <- Data[grep(geneID, Data$Genes, fixed = TRUE),]
  GenesData[n,'Count'] <- nrow(info)
  chr <- info$Chr[1] #All chr should be the same because it's the same gene
  if (any(info$Chr!=chr)){
    print(paste0('Chromosome discrepancies in ', geneID))
  }
  GenesData[n,'Chr'] <- chr
  GenesData[n,'Pops'] <- summarizeInfo(info$Pops)
  GenesData[n,'Metapops'] <- summarizeInfo(info$Metapops)
  GenesData[n,'Stats'] <- summarizeInfo(info$Stats)
  GenesData[n,'Evidence'] <- summarizeInfo(info$Evidence)
  GenesData[n,'SourceIDs'] <- summarizeInfo(info$ID)
  if (n%%100 == 0){
    print(n)
  }
}
#I connect to Ensembl and add the gene type for each gene
ensembl <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl') #I access genome build hg38
print('Connected to Ensembl')
genesType <- getBM(attributes = c('gene_biotype','ensembl_gene_id'), filters =  'ensembl_gene_id', 
                                  values = GenesData$GeneID, mart = ensembl)
GenesData_2 <- merge(x = GenesData, y = genesType, by.x = 'GeneID', by.y = 'ensembl_gene_id')
#I save the file, ordered by the ratio Count/Number of metapopulatios. Number of pops is  str_count(Metapops, ',') + 1
write.csv(GenesData_2 %>% mutate(ratio = Count/(str_count(Metapops, ',') + 1)) %>% arrange(desc(ratio)), 
          file = 'GenesData.csv', row.names = FALSE)
