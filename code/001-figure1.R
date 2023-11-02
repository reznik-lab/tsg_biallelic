library(ggplot2)
library(patchwork)
library(stringr)
library(ggnewscale)
library(RColorBrewer)
library(ggpubr)
library(maftools)

###Plotting functions used in making the panels:
#Driver/VUS zygosity barplot (f1a):
zygosityBarplot <- function(zygosityData, colorvec, mainsize=12, xsize=6){
  tab <- as.data.frame(table(zygosityData[which(zygosityData$mechanism %in% names(colorvec) & zygosityData$driver_oncogenic),c('mechanism','type')]))
  tab <- tab[which(tab$type != 'Other'),]
  tab$N <- sapply(1:nrow(tab),function(i){sum(tab$Freq[which(tab$type == tab$type[i])])})
  tab$Fraction <- tab$Freq/tab$N
  colnames(tab)[1] <- c('Category')
  tab$Category <- factor(as.character(tab$Category), levels = names(colorvec))
  tab$type <- as.character(tab$type)
  nTSGs <- length(unique(zygosityData$Hugo_Symbol[which(zygosityData$type == 'TSG')]))
  nOGs <- length(unique(zygosityData$Hugo_Symbol[which(zygosityData$type == 'Oncogene')]))
  tab$type[which(tab$type == 'TSG')] <- paste('TSG (',nTSGs,')',sep='')
  tab$type[which(tab$type == 'Oncogene')] <- paste('Oncogene (',nOGs,')',sep='')
  tab$type <- factor(tab$type, levels = unique(tab$type))
  p <- ggplot(aes(x=type,y=Fraction,fill=Category), data=tab) + 
    geom_bar(stat = 'identity', position = 'stack') + theme_bw() + scale_fill_manual(values=colorvec[which(!names(colorvec) %in% 
                                                                                                             c('LOH','Mutation + Gain of Mutant'))]) + xlab('') + 
    ggtitle('Alteration zygosity across drivers') + xlab('') + theme(legend.position="left") + 
    theme(axis.text.x= element_text(size = xsize)) + theme(plot.title = element_text(size = mainsize)) + 
    theme(legend.position="right")
  p
}

#Mechanism barplot (f1c):
mechanismBarplot <- function(table_s5, mainsize=12, subtypeTable, xsize=6){
  #Reformatting table S5 for plot
  subtype_input <- table_s5[which(table_s5$pvalue.corrected < .05 & table_s5$Total_Biallelic_Events >= 10 & 
                                    table_s5$Dominant_Mechanism_Pancan != table_s5$preferred_mechanism),]
  subtype_input <- Reduce(rbind, lapply(c(4:7), function(i){
    cbind(subtype_input[,c(1,2,9)],'Fraction'=subtype_input[,i],'Mechanism'=strsplit(colnames(subtype_input)[i],'Fraction_')[[1]][2])}))
  #Adding pancancer rows:
  pancan <- unique(table_s5[which(table_s5$pvalue.corrected < .05 & table_s5$Total_Biallelic_Events >= 10 & 
                      table_s5$Dominant_Mechanism_Pancan != table_s5$preferred_mechanism),c(2,9,14:17)])
  pancan_input <- Reduce(rbind, lapply(c(3:6), function(i){cbind('Disease'='Pan-Cancer', pancan[,c(1,2)],'Fraction'=pancan[,i],
                                                                 'Mechanism'=strsplit(colnames(pancan)[i],"Fraction_")[[1]][2])}))
  pancan_input$pvalue.corrected <- NA
  pancan_input <- unique(pancan_input)
  pancan_input$Mechanism <- str_split_fixed(pancan_input$Mechanism, " ", 2)[,1]
  df <- rbind(subtype_input, pancan_input)
  #Only plotting significant:
  #x axis
  df$abbrev <- subtypeTable$ONCOTREE_CODE[match(df$Disease, subtypeTable$CANCER_TYPE_DETAILED)]
  df$n <- subtypeTable$n[match(df$Disease, subtypeTable$CANCER_TYPE_DETAILED)]
  df$abbrev_n <- str_c(df$abbrev, ' (', df$n, ')')
  df$abbrev_n[which(df$Disease == 'Pan-Cancer')] <- 'Pan-Cancer'
  #Ordering, setting factor levels (pan-cancer n is set to -1 to make sure it comes first):
  df$n[which(df$Disease == 'Pan-Cancer')] <- -1
  df <- df[order(df$Gene, as.numeric(df$n), decreasing = TRUE),]
  df$abbrev_n <- factor(df$abbrev_n, levels = c('Pan-Cancer', unique(df$abbrev_n)[-which(unique(df$abbrev_n)=='Pan-Cancer')]))
  df$Gene <- factor(df$Gene, levels = unique(df$Gene))
  df$Fraction <- as.numeric(df$Fraction)

  p <- ggplot() + geom_bar(aes(x=abbrev_n,y=Fraction, fill=Mechanism),data=df, position = 'stack',stat = 'identity') + 
    scale_fill_manual(values = c('Mut_LOH'="#58B9F4",'Mut_Fusion'="#9E65AA",'Homdel'="#242860",'Compound'="#8CA1D3")) + theme(axis.title.x = element_text(size=4)) + 
    theme(axis.text.x = element_text(size=xsize,angle=90, hjust=1), axis.title.x = element_blank()) +
    facet_wrap(~ Gene, scales = "free_x",ncol = length(unique(df$Gene))) + theme(legend.position="none")
  p
}

#Generating figure 1 A, B, C:
generateFigure1abc <- function(drivers, drivers.with_fusions, xsize=9, mainsize=12, tsgdf, subtypeTable){
  #f1a and b:
  colorvec <- c('Amplification'="#D25047",'Mutation'="#CED1CE",'Mutation + Gain of Mutant'="#EDBE59",'Fusion'='plum1','Compound'= "#8CA1D3",'Homdel'="#242860", 
                'Mut + Fusion'="#9E65AA",'Mut + LOH'="#58B9F4")
  data <- drivers
  data.with_fusions <- drivers.with_fusions
  if('Mut + fusion' %in% data.with_fusions$mechanism){data.with_fusions$mechanism[which(data.with_fusions$mechanism == 'Mut + fusion')] <- 'Mut + Fusion'}
  f1b <- zygosityBarplot(zygosityData=data.with_fusions, colorvec=colorvec, xsize = xsize)
  
  #f1c:
  f1c <- mechanismBarplot(table_s5, mainsize=mainsize, subtypeTable=subtypeTable, xsize = xsize, filterMutations = 0)
  layout <- 
    "
  ###ACCC
  "
  
  f1abc <- f1b + f1c + 
    plot_layout(design = layout, heights=c(2.7,9)) + plot_annotation(tag_levels = list(c('B','C')))
  
  #Save figure:
  pdf(paste(figurePath,'/f1abc.pdf',sep=''),width = 23.4, height = 7.6)
  print(f1abc)
  dev.off()
}
generateFigure1abc(drivers.impact, drivers.impact.with_fusions, xsize=10, mainsize=12, tsgdf, subtypeTable=select_ctds)

