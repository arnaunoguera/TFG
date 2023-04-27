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
  scale_fill_manual(values = c("#276419","#7FBC41", "#FDE725"), name = 'Candidates of natural selecion') + 
  labs(x = 'Chromosome', y = 'Number of genes', colour = 'Total genes per chromosome') +
  geom_line(data = HumanGenome, mapping = aes(x = Name, y = Gene, group = 'drop', color = '# of genes'), linewidth = 1.1) + 
  geom_point(data = HumanGenome, mapping = aes(x = Name, y = Gene, color = '# of genes'), size = 3)+ 
  scale_color_manual(values = '#64B87B') + theme_light()

###################### Regions detected in each chromosome
ggplot(data = Data, mapping = aes(x = factor(Chr, levels = c(1:22, 'X', 'Y')))) + geom_bar(stat = 'count', position = 'dodge', fill= '#64C23A')

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
ggplot(as.data.table(rbind(ChromosomeBP %>% mutate(Type = 'Detected regions') %>% dplyr::rename(Name = Chr),
             HumanGenome %>% mutate(Type = 'Human Genome', 'BP' = Size..Mb. * 1000000) %>% dplyr::select(Type, Name, BP))),
       mapping = aes(x = factor(Name, levels = c(1:22, 'X', 'Y')), y = as.numeric(BP)/1000000, fill = Type)) +
  geom_col(position = 'dodge') + scale_fill_manual(values = c("#276419", "#B4DE2C"), name = '') +
  labs(x = 'Chromosome', y = 'Size (Mb)') + theme_light()

######################## Plot representing each detected gene count of regions and metapopulations:
ggplot(GenesData %>% mutate(MetapopCount = str_count(Metapops, ',')+1), mapping = aes(x = MetapopCount, y = Count, color = Evidence)) +
  geom_point()

###############Histogram of the detected region size
Data %>% mutate(Size = End - Start) %>% ggplot(mapping= aes(x = log10(Size))) + geom_histogram(fill = '#B4DE2C', color = 'black')
