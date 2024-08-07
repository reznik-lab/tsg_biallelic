library(stringr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(patchwork)

biallelicRateTileplot <- function(data, subtypeTable, colorvec, ncolors=4, xsize=14, mainsize=14){
  dat <- data
  #% rather than fraction:
  dat$`Biallelic %` <- round(100*(dat$`Biallelic Fraction`))
  dat$`Alteration %` <- round(100*(dat$Alteration_Rate))
  
  ###Adding number of samples and disease abbrev. + n to be x-axis values:
  dat$n <- subtypeTable$n[match(dat$Disease, subtypeTable$Disease)]
  dat$abbrev <- subtypeTable$CTD_ABBR[match(dat$Disease, subtypeTable$CANCER_TYPE_DETAILED)]
  dat$abbrev_n <- str_c(dat$abbrev, ' (', dat$n, ')')
  
  #Gene-level object to be used to make barplot:
  genedf <- merge(unique(mechdf[,c('Hugo_Symbol','Pathway','gene_alts')]), dat[which(dat$Disease == 'all'),], by.x='Hugo_Symbol', by.y='Gene')
  genedf$`Alteration %` <- 100*genedf$Alteration_Rate
  
  #Adding gene_alt category to 'dat' for plotting purposes:
  dat$gene_alts <- mechdf$gene_alts[match(dat$Gene, mechdf$Hugo_Symbol)]

  #Adding pathway info to inputs:
  dat$Pathway <- factor(mechdf$Pathway[match(dat$Gene, mechdf$Hugo_Symbol)], levels = levels(mechdf$Pathway))
  genedf$Pathway <- factor(mechdf$Pathway[match(genedf$Hugo_Symbol, mechdf$Hugo_Symbol)], levels = levels(mechdf$Pathway))
  
  #Ordering organs:
  organOrder <- unique(subtypeTable$organ[order(subtypeTable$organ_order)])
  dat$organ <- factor(subtypeTable$organ[match(dat$Disease, subtypeTable$CANCER_TYPE_DETAILED)], levels = organOrder)
  dat$organ_color <- as.character(subtypeTable$organ_color[match(dat$Disease, subtypeTable$CANCER_TYPE_DETAILED)])
  
  #Making sure abbrev_n levels are in right order:
  dat$abbrev_n <- factor(dat$abbrev_n, levels = unique(dat$abbrev_n))
  
  #Making sure other factor levels are in right order:
  dat$Gene <- factor(dat$Gene, levels = levels(mechdf$Hugo_Symbol))
  genedf$Hugo_Symbol <- factor(genedf$Hugo_Symbol, levels = levels(mechdf$Hugo_Symbol))
  mechdf$Hugo_Symbol <- factor(mechdf$Hugo_Symbol, levels = levels(mechdf$Hugo_Symbol))
  mechdf$Mechanism <- factor(mechdf$Mechanism, levels = c("Mutation","Mutation + Gain of Mutant","Amplification","Compound",
                                                          "Homdel","Mut + Fusion","Mut + LOH"))
  
  #Getting rid of gene/diseases with < 4 mutations:
  dat$`Alteration %`[which(dat$`Alteration %` == 0)] <- ''
  dat$`Biallelic %`[which(dat$`Alteration %` == 0)] <- 0
  dat$`Biallelic %`[which(dat$Alterations < 4)] <- 0
  dat$`Biallelic %`[which(dat$`Alteration %` == '')] <- 0
  dat$`Alteration %`[which(dat$Num_Alterations < 4 & dat$Num_Alterations > 0 | dat$`Alteration %` == '')] <- ''
  
  #Color gradient:
  colfunc <- colorRampPalette(c("white", "yellow", "red"))
  grad <- colfunc(ncolors)
  
  grad = c("#FFFFFF", "#6CC4AD", "#CBDC52", "#F8D846", "#EF273E")
  grad = c("#FFFFFF", "#6CC4AD", "yellow", "#EF273E")
  
  colfunc <- colorRampPalette(c("white", "#FFFB66", "#EF273E"))
  grad <- colfunc(ncolors)
  
  #Max alt rate plotted: 0.15
  genedf$Alteration_Rate[which(genedf$Alteration_Rate >= 0.15)] <- 0.15
  p_alt <- ggplot(aes(x=gene_alts, y=Alteration_Rate), data=genedf %>% select(gene_alts, Alteration_Rate, Pathway) %>% unique) + 
    geom_bar(stat='identity', fill='black') + 
    theme_bw() + coord_flip() + xlab('') + ylab('Alteration Fraction') + theme(axis.text.y = element_text(size=9)) + ggtitle('') + 
    theme(plot.title  = element_text(size=mainsize), axis.text.x = element_text(size=xsize)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    facet_grid(Pathway ~ ., scales = "free", space = "free") + 
    theme(strip.text = element_text(
      size = 6),
      legend.position = 'none')
  
  #X axis text has to be aded manually:
  lastGene <- head(genedf[order(genedf$Pathway, decreasing = TRUE),]$Hugo_Symbol, 1)
  dat$xlabel <- ''
  dat$xlabel[which(dat$Gene == lastGene)] <- dat$abbrev_n[which(dat$Gene == lastGene)]
  
  #Getting rid of any genes we don't have pathway assignement for:
  dat <- dat[which(!is.na(dat$Pathway)),]
  
  p_heatmap <- ggplot(data=dat,aes(x=abbrev_n,y=gene_alts,fill=`Biallelic %`)) + geom_tile() + theme_bw() + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),axis.title.y = element_blank(), plot.title = element_text(size=mainsize),
          legend.position = 'none') + 
    theme(plot.margin = unit(c(1,1,9.7,1.5), "lines")) +
    geom_text(aes(label = `Alteration %`), color = "black", size = 3) + scale_fill_gradientn(colours = grad,values = seq(from=0, to=1, length=ncolors)) + 
    labs(x='',y='Gene') + ylab('') + ggtitle('') + 
    geom_point(data=dat[which(dat$Num_Alterations < 4 & dat$Num_Alterations > 0 | dat$`Alteration %` == ''),], shape=5, aes(x=abbrev_n,y=gene_alts,fill=`Biallelic %`)) + 
    facet_grid(Pathway ~ organ, scales = "free", space = "free") + 
    theme(strip.text = element_text(size = 6)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_blank(),
          strip.background = element_blank())
  
  tmb <- ggplot(data=fga_and_tmb_dat) + geom_violin(aes(x=abbrev_n, y=TMB), fill='#a0bcbb') + theme(axis.text.x = element_blank(),axis.text.y = element_text(size = xsize),axis.title.y = element_text(size = xsize),
                                                                                                    axis.title.x = element_blank()) + xlab('Disease') + 
    geom_text(data=fga_and_tmb_dat,aes(x=abbrev_n,y=-1.5,colour=organ_color,label=abbrev_n),angle=90,size=4,vjust=.5,hjust=1) + theme(legend.position='none') +
    coord_cartesian(ylim=c(0,20), clip = 'off') + theme(plot.margin = unit(c(1,1,9.7,1.5), "lines")) + 
    facet_grid(. ~ organ, scales = "free", space = "free") + 
    theme(strip.text = element_text(size = 6)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_blank(),
          strip.background = element_blank()) 
  
  fga <- ggplot() + geom_violin(data=fga_and_tmb_dat, aes(x=abbrev_n, y=FGA),fill='#a492c5') + theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
                                                                                                     axis.text.y = element_text(size = xsize),axis.title.y = element_text(size = xsize)) + 
    facet_grid(. ~ organ, scales = "free", space = "free")+ 
    theme(strip.text = element_text(size = 6)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_blank(),
          strip.background = element_blank()) 
  
  p_bi <- ggplot(aes(x=gene_alts,y=Fraction_of_Alterations, fill=Mechanism), data=mechdf) + geom_bar(stat = 'identity',position = 'stack') + 
    coord_flip() + theme_bw() + xlab('') + ylab('Biallelic Fraction') + 
    ggtitle('') + theme(axis.text.y = element_blank(), plot.title = element_text(size=mainsize), axis.title.x = element_text(size = xsize)) + 
    scale_fill_manual(values=colorvec) + theme(legend.position = "none") + 
    facet_grid(Pathway ~ ., scales = "free", space = "free") + 
    theme(strip.text = element_text(size = 6)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_blank(),
          strip.background = element_blank()) 
  
  ###Sideways stacked barplot, showing 1) how many/fraction (?) subtypes >x% mutated; 2) how many preferentially biallelic/monoallelic.
  dat$Status <- NA
  dat$Status[which(dat$Num_Alterations >= 10 & dat$`Alteration %` >= 1)] <- 'Not Biallelic'
  dat$Status[which(dat$Num_Alterations >= 10 & dat$`Alteration %` >= 1 & dat$`Biallelic %` > 80)] <- 'Biallelic'
  dtab <- as.data.frame(table(dat[,c('Gene','Status')]))
  dtab$Status <- factor(dtab$Status, levels=c('Not Biallelic','Biallelic'))
  genedf <- merge(genedf, data.frame('Hugo_Symbol'=dtab$Gene, 'Status'=dtab$Status, 'Subtypes'=dtab$Freq), by='Hugo_Symbol')
  p_subtypes <- ggplot(data=genedf, aes(x=gene_alts,y=Subtypes, fill=Status)) + geom_bar(stat = 'identity', position='stack') + theme_bw() + 
    xlab('') + ylab('Disease Subtypes') + ggtitle('')  + coord_flip() + scale_fill_manual(values = c('Biallelic'='#1b8843', 'Not Biallelic'='#cbd8c0')) + 
    theme(legend.position = 'bottom', axis.text.x = element_blank()) + facet_grid(Pathway ~ ., scales = "free", space = "free")  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_blank(),
          strip.background = element_blank())
  
  layout <- 
    "
  ABCCCCD
  ##EEEE#
  ##FFFF#
  ##GGGG#
  "
  p <- p_alt + p_bi + p_heatmap + p_subtypes + plot_spacer() + fga + tmb + plot_layout(design = layout, widths = c(.7,.7,.7,9), heights = c(8,-1.63,.5,.5))
  p
}

#Generating figure 1:
#'data' should be supp. table S6
generateFigure2 <- function(data, xsize=9, mainsize=12, subtypeTable){
  colorvec <- c('Amplification'="#D25047",'Mutation'="#CED1CE",'Mutation + Gain of Mutant'="#EDBE59",'Fusion'='plum1','Compound'= "#8CA1D3",'Homdel'="#242860", 
                'Mut + Fusion'="#9E65AA",'Mut + LOH'="#58B9F4")
  f2 <- biallelicRateTileplot(data=data, colorvec=colorvec, subtypeTable=subtypeTable)
  f2
}

