source('../data/GlobalFunctions.R')
Data <- read.csv('../ProcessedData/AllData.csv')
GenesData <- read.csv('../ProcessedData/GenesData.csv')
HumanGenome <- read.table(file = 'HumanGenome.txt', header = TRUE, sep = '\t')

################# Gene type pie chart ##############
GeneTypesCount <- GenesData %>% group_by(gene_biotype) %>% summarise(count = n()) %>% arrange(desc(count))
GeneTypesKeep <- GeneTypesCount %>% filter(count > 500) %>%
  add_row(gene_biotype = 'other', count = sum(GeneTypesCount %>% filter(count < 500) %>% pull(count)))
ggplot(data= GeneTypesKeep, mapping = aes(x = '', y = count, fill = gene_biotype)) +
  geom_bar(stat = 'identity', width = 1, color = 'white') + 
  coord_polar('y', start = 0) + theme_void()

pie(GeneTypesKeep$count, labels = GeneTypesKeep$gene_biotype, radius = -1, 
    col = c("#FDE725", "#B4DE2C", "#7FBC41", "#4D9221", "#35873D", "#276419", "#1D3409", "#0C170D"))

GeneTypesKeep <- GeneTypesKeep %>% mutate(biotype = case_when(gene_biotype == 'protein_coding' ~ 'Protein coding',
                                                              gene_biotype == 'lncRNA' ~ 'lncRNA',
                                                              gene_biotype == 'processed_pseudogene' ~ 'Processed pseudogene',
                                                              gene_biotype == 'misc_RNA' ~ 'miscRNA',
                                                              gene_biotype == 'unprocessed_pseudogene' ~ 'Unprocessed pseudogene',
                                                              gene_biotype == 'snRNA' ~ 'snRNA',
                                                              gene_biotype == 'miRNA' ~ 'miRNA',
                                                              gene_biotype == 'other' ~ 'Other'))
pie(GeneTypesKeep$count, labels = GeneTypesKeep$biotype, radius = -1, 
    col = c("#75E6DA", "#A3EBB1", "#349B89","#189AB4", "#006A8A", "#2F5496", "#05445E", '#0C2D48'))
pie(GeneTypesKeep$count, labels = GeneTypesKeep$biotype, radius = -1, 
    col = c("#FDE725", "#B4DE2C", "#7FBC41", "#4D9221", "#35873D", "#276419", "#1D3409", "#0C170D"))


################## Chromosome counts ###########################
ggplot(GenesData %>% mutate(coding = case_when(gene_biotype == 'protein_coding' ~ TRUE,
                                               TRUE ~ FALSE)), 
       mapping = aes(x= factor(Chr, levels = c(1:22, 'X', 'Y')), fill = coding)) + geom_bar(stat ='count', position = 'stack') +
  scale_fill_manual(values = c("#276419","#B4DE2C"))

################## Chromosome counts amb els gens de cada chr (gens = codificant + gens RNA) ###########################
GenesDataGrouped <- GenesData %>% mutate(coding = case_when(gene_biotype == 'protein_coding' ~ 'Protein coding genes',
                                                            gene_biotype %in% c('translated_unprocessed_pseudogene', 'TR_J_pseudogene', 'IG_J_pseudogene', 'IG_C_pseudogene', 'pseudogene', 'TR_V_pseudogene', 'unitary_pseudogene', 'IG_V_pseudogene', 'transcribed_unitary_pseudogene', 'rRNA_pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unprocessed_pseudogene', 'unprocessed_pseudogene', 'processed_pseudogene') ~ 'Pseudogenes',
                                                            TRUE ~ 'RNA genes'))

ggplot(GenesDataGrouped) + 
  geom_bar(data = GenesDataGrouped, mapping = aes(x= factor(Chr, levels = c(1:22, 'X', 'Y')), fill = factor(coding, levels=c('Pseudogenes','RNA genes', 'Protein coding genes'))), stat ='count', position = 'stack') +
  scale_fill_manual(values = c("#276419","#B4DE2C", "#FDE725"), name = 'Candidates of positive selecion') + 
  labs(x = 'Chromosome', y = 'Number of genes', colour = 'Total genes per chromosome') +
  geom_line(data = HumanGenome, mapping = aes(x = Name, y = Gene, group = 'drop', color = 'Number of genes'), linewidth = 1.1) + 
  geom_point(data = HumanGenome, mapping = aes(x = Name, y = Gene, color = 'Number of genes'), size = 3)+ 
  scale_color_manual(values = '#35873D') + theme_light() + theme(axis.title = element_text(size = 15), 
                                                                 legend.text = element_text(size= 13), 
                                                                 legend.title = element_text(size= 13))
######Correlation
GeneCorr <- left_join(x = HumanGenome %>% dplyr::select(Name, Gene), y = GenesData %>% group_by(Chr) %>% summarise(GenesRegion = n()), 
                      by = c('Name' = 'Chr')) %>% 
  dplyr::rename('Chr' = 'Name', 'GenesChr' = 'Gene')
cor.test(GeneCorr$GenesRegion, GeneCorr$GenesChr, method = 'pearson')
cor.test(GeneCorr$GenesRegion, GeneCorr$GenesChr, method = 'spearman')


###################### Regions detected in each chromosome
ggplot(data = Data, mapping = aes(x = factor(Chr, levels = c(1:22, 'X', 'Y')))) + 
  geom_bar(stat = 'count', position = 'dodge', fill= '#2F5496') + labs(y = 'Number of candidate regions' , x = 'Chromosome', 
                                                                      tag = 'B') + theme_light() +
  theme(axis.title = element_text(size = 15), plot.tag = element_text(size= 16))

################### Plot to determine the amount of bp detected as suffering natural selection in each chromosome
#Function to input a vector of starts and one of ends: 07329887-107339887, etc.and get the amount of bp represented, considering overlaps
MeasureBP <- function(startvector, endvector){ #The input is the start and end vectors for ONE chromosome only
  regions <- data.table(Start = startvector, End = endvector, Checked = FALSE)
  totalregions <- data.table(Start = numeric(), End = numeric())
  for (n0 in 1:nrow(regions)){
    if (regions[n0,'Checked'] == TRUE){
      next #If that row has already been checked previously, skip to the next
    }
    print(paste0(n0, ' of ', nrow(regions)))
    startbp <- as.numeric(regions[n0,'Start']) #This will be replaced by the new start and end of the interval if it overlaps any regions
    endbp <- as.numeric(regions[n0,'End'])
    if (n0 == nrow(regions)){
      break #If it's the last region, there's no need to compare against others
    }
    n1 <- n0
    while (TRUE){ #I check the region against all others
      n1 <- n1 + 1
      if (regions[n1,'Checked'] == TRUE){
        if (n1 == nrow(regions)){
          break #If this is the last row, we end the loop
        }
        next #If that row has already been checked previously, skip to the next
      }
      startbp1 <- as.numeric(regions[n1,'Start'])
      endbp1 <- as.numeric(regions[n1,'End'])
      if ((startbp1 >= startbp & startbp1 <= endbp)|(endbp1 >= startbp & endbp1 <= endbp)){ #This means they overlap
        #We have to keep the lowest start and the highest end to get the whole region
        startbp <- min(c(startbp, startbp1))
        endbp <- max(c(endbp, endbp1))
        regions[n1, 'Checked'] <- TRUE #This region has already been compared
        #We will repeat this untill all regions have been included in the overlap or do not overlap (we have to repeat from the beginning because some regions might not overlap before but may do now that it is bigger)
        n1 <- n0
      }
      if (n1 == nrow(regions)){
        break #If this is the last row, we end the loop
      }
    }#When the for loop ends, startbp and endbp have the start and end coordinates of he whole overlap. I save these.
    totalregions <- totalregions %>% add_row(Start = startbp, End = endbp)
  }#When the for loop ends, the whole table has been analyzed
  #The total amount of bp is the sum of all end - start coordinates
  return(totalregions %>% mutate(bp = End - Start) %>% summarise(BP =sum(bp)))
}
#I summarise every chr like that
ChromosomeBP <- Data %>% group_by(Chr) %>% summarise(BP = MeasureBP(Start, End))
write.csv(ChromosomeBP, file = 'BPCountsByChr.csv')
#Summarised counts: 
ChromosomeBP <- data.table(Chr = c(1,10:19,2,20:22,3:9,'X','Y'), BP = c(95339371, 58662407, 62107280, 63719479, 43275984, 37717935,
                                                                    37478299, 37273279, 34744306, 31892375, 22432739, 110437724,
                                                                    29202991, 17859042, 16558885, 91080865, 77957934, 77102804,
                                                                    79623363, 73592863, 68551167, 51002302, 54917080, 144569))
#Summarising for each chr and plot
ggplot(as.data.table(rbind(ChromosomeBP %>% mutate(Type = 'Candidate regions') %>% dplyr::rename(Name = Chr),
             HumanGenome %>% mutate(Type = 'Human Genome', 'BP' = Size..Mb. * 1000000) %>% dplyr::select(Type, Name, BP))),
       mapping = aes(x = factor(Name, levels = c(1:22, 'X', 'Y')), y = as.numeric(BP)/1000000, fill = Type)) +
  geom_col(position = 'dodge') + scale_fill_manual(values = c("#276419","#B4DE2C"), name = '') +
  labs(x = 'Chromosome', y = 'Size (Mb)') + theme_light() + theme(axis.title = element_text(size = 15))
######Correlation test between chromosome size and the size of the candidate region in that chromosome
SizeTable <- left_join(x = HumanGenome %>% dplyr::select(Name, Size..Mb.), y = ChromosomeBP, by = c('Name' = 'Chr')) %>% 
  dplyr::rename('Chr' = 'Name', 'ChrMb' = 'Size..Mb.', 'CandidateBp' = 'BP') %>% mutate('ChrBp' = ChrMb * 1000000)
cor.test(SizeTable$ChrMb, SizeTable$CandidateBp, method = 'pearson')
cor.test(SizeTable$ChrMb, SizeTable$CandidateBp, method = 'spearman')


###############Histogram of the detected region size
Data %>% mutate(Size = End - Start) %>% ggplot(mapping= aes(x = log10(Size))) + geom_histogram(fill = '#B4DE2C', color = 'black')

##############Histogram of sources
library(ggbreak)
Data %>%
  ggplot(mapping = aes(x= Source)) + geom_bar(fill = '#2F5496') + scale_y_break(c(3000,14500)) + 
  scale_y_continuous(limits = c(0, 15500), breaks = seq(0, 15500, 500)) + labs(y = 'Number of candidate regions', x = 'Source', tag = 'A') +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust =1), axis.title.x = element_text(margin = margin(t = -20)), 
                        axis.title = element_text(size = 15), plot.tag = element_text(size= 16))

#############Region count by metapopulation
metapops <- unlist(strsplit(unlist(strsplit(Data$Metapops, split = ',')), split = '2'))
ggplot(data.table(Metapop = metapops) %>% filter(!is.na(Metapop)), 
       mapping = aes(x = factor(Metapop, levels = c('AFR','EAS','EUR','SAS','AMR','MDE','OCN')), fill = Metapop)) + 
  geom_bar() + scale_fill_manual(values = c('AFR' = "#F7F14A", 'EAS' = "#33B033", 'EUR'="#5691C4", 'SAS' = "#A965BA", 
                                            'AMR' = '#AB3B35', 'OCN' = '#F28500', 'MDE' = '#C0C2C9')) +
  labs(x = 'Metapopulation', y = 'Number of candidate regions', tag = 'C') + theme_light() +
  theme(axis.title = element_text(size = 15), plot.tag = element_text(size= 16), legend.position = 'none') 


###################Genes analysis
nSource <- function(SourceIDsVector){
  nVector <- c()
  for (register in 1:length(SourceIDsVector)){
    SourceIDs <- strsplit(unlist(strsplit(SourceIDsVector[register], ',')), '_')
    IDsfinal <- c()
    for (ID in SourceIDs){
      if (!ID[1] %in% IDsfinal){
        IDsfinal <- append(IDsfinal, ID[1])
      }
    }
    nVector <- append(nVector, length(IDsfinal))
  }
  return(nVector)
}

GenesData <- GenesData %>% mutate(nSource = nSource(SourceIDs))

TopSourceGenes <- GenesData %>% arrange(desc(nSource)) %>% head(31) 

ensembl <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl') #I access genome build hg38
print('Connected to Ensembl')
genesname <- getBM(attributes = c('external_gene_name','ensembl_gene_id'), filters =  'ensembl_gene_id', 
                   values = TopSourceGenes$GeneID, mart = ensembl)
TopSourceGenes <- merge(x = TopSourceGenes, y = genesname, by.x = 'GeneID', by.y = 'ensembl_gene_id')
TopSourceGenes[21, "external_gene_name"] <- 'ENSG00000233005'
TopSourceGenes[22, "external_gene_name"] <- 'COX7B Pseudogene'
TopSourceGenes[30, "external_gene_name"] <- 'ENSG00000286810' 
TopSourceGenes[23, "gene_biotype"] <- 'unprocessed_pseudogene'
TopSourceGenes <- TopSourceGenes %>% filter(external_gene_name != 'SULT1C2P1' | gene_biotype != 'lncRNA')

ggplot(TopSourceGenes, mapping = aes(x = nSource, y = factor(external_gene_name, levels = TopSourceGenes %>% arrange(nSource) %>% pull(external_gene_name)), 
                                     fill = gene_biotype)) + geom_col() +
  scale_fill_manual(values = c('protein_coding' = '#FDE725','lncRNA' = '#B4DE2C', 'misc_RNA' = '#7FBC41', 'processed_pseudogene' = '#4D9221',  
                               'unprocessed_pseudogene' = '#276419'), 
                    breaks = c('protein_coding','lncRNA','misc_RNA', 'processed_pseudogene','unprocessed_pseudogene'),
                    labels = c('protein_coding' = 'Protein-coding','lncRNA' = 'lncRNA', 'misc_RNA' = 'miscRNA', 'processed_pseudogene' = 'Processed pseudogene',  
                               'unprocessed_pseudogene' = 'Unprocessed pseudogene'), name = 'Gene biotype') +
  labs(y = 'Gene name', x = 'Number of papers', tag = '') + theme_light() +
  theme(axis.title = element_text(size = 15), plot.tag = element_text(size= 16)) +
  scale_x_continuous(breaks = seq(0,10,2))


###################Genes analysis 2
TopMetapopGenes <- GenesData %>% arrange(desc(ratio)) %>% head(36) 
ensembl <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl') #I access genome build hg38
print('Connected to Ensembl')
genesname <- getBM(attributes = c('external_gene_name','ensembl_gene_id'), filters =  'ensembl_gene_id', 
                   values = TopMetapopGenes$GeneID, mart = ensembl)
TopMetapopGenes <- merge(x = TopMetapopGenes, y = genesname, by.x = 'GeneID', by.y = 'ensembl_gene_id') %>% filter(gene_biotype != 'TEC') %>% 
  mutate(GeneName = case_when(external_gene_name == ''~GeneID, TRUE~external_gene_name))
TopMetapopGenes[22, "external_gene_name"] <- 'RPL21P125'
TopMetapopGenes[23, "external_gene_name"] <- 'NDUFB8'


ggplot(TopMetapopGenes, mapping = aes(x = Count, y = factor(GeneName, levels = TopMetapopGenes %>% arrange(ratio) %>% pull(GeneName)), 
                                     color = Metapops, fill = gene_biotype)) + geom_col(size = 1.5) +
  labs(y = 'Gene name', x = 'Number of records', tag = 'B') + 
  scale_color_manual(name = 'Metapopulations', values = c('AFR' = "#F7F14A", 'EAS' = "#33B033", 'SAS,AFR,EAS,EUR,AMR' = 'orange'),
                    labels = c('SAS,AFR,EAS,EUR,AMR' = 'AFR, EUR, EAS, SAS, AMR')) +
  scale_fill_manual(values = c('protein_coding' = '#75E6DA','lncRNA' = '#349B89','processed_pseudogene' = '#006A8A'), 
                    breaks = c('protein_coding','lncRNA','processed_pseudogene'),
                    labels = c('protein_coding' = 'Protein-coding','lncRNA' = 'lncRNA', 'processed_pseudogene' = 'Processed pseudogene'), 
                    name = 'Gene biotype') +
  scale_x_break(c(13,55), scales = c(0.2,1)) +  scale_x_continuous(breaks = c(0,5,10,55,60), limits = c(0,60)) +
  theme_light() +
  theme(axis.title = element_text(size = 15), plot.tag = element_text(size= 16)) 
