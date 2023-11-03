#Gene expression boxplots for fig. 4:
expressionBoxplots <- function(counts, annotation=NULL, genes, counts.prot, proteins, hugo, textsize=12, mainsize=14, levels, levelNames, remove=NULL, colorvec=NULL, indices=NULL){
  if(is.null(indices)){
    if(length(levels) == 2){indices <- 1:2}else{
      indices <- c(which(levels == 'a_wt_wt'), which(levels == 'd_biallelic_wt'))
      if(length(indices) == 0){indices <- c(which(levels == 'a_wt'), which(levels == 'c_biallelic'))}
    }
  }
  if(is.null(counts$group)){
    counts$group <- annotation
    counts.prot$group <- annotation
  }
  if(!is.null(remove)){
    counts <- counts[which(!rownames(counts) %in% remove),]
    counts.prot <- counts.prot[which(!rownames(counts.prot) %in% remove),]
  }
  
  #RNA Box plots:
  if(!is.null(counts)){
    ens_boxplot <- names(hugo)[match(genes, hugo)]
    df.barplot <- counts[,c(ncol(counts),match(ens_boxplot, colnames(counts)))]
    df.barplot <- as.data.frame(Reduce(rbind, lapply(2:(length(genes)+1), function(i){as.matrix(df.barplot[,c(1,i)])})))
    df.barplot$Gene <- as.vector(sapply(genes, function(g){rep(g, nrow(counts))}))
    colnames(df.barplot)[1:2] <- c('Zygosity','Normalized_Counts')
    df.barplot$Normalized_Counts <- as.numeric(df.barplot$Normalized_Counts)
    if(is.null(colorvec)){
      colorvec <- c('#ed2225','#3a53a4','#69bd45','#8151a1','darkorange','cyan1')[1:length(levels)]
      names(colorvec) <- levelNames
    }
    df.barplot <- df.barplot[which(df.barplot$Zygosity %in% levels),]
    df.barplot$Zygosity <- levelNames[match(df.barplot$Zygosity, levels)]
    df.barplot$Zygosity <- factor(df.barplot$Zygosity, levels = levelNames)
    pvals.rna <- try(sapply(genes,function(g){signif(wilcox.test(x=df.barplot$Normalized_Counts[which(df.barplot$Gene == g & df.barplot$Zygosity == levelNames[indices[1]])],
                                                                 y=df.barplot$Normalized_Counts[which(df.barplot$Gene == g & df.barplot$Zygosity == levelNames[indices[2]])])$p.value,digits = 3)}), silent = TRUE)
    #ti.rna <- paste('p=',pvals.rna[g], ' (',levelNames[indices[1]],' vs. ',levelNames[indices[2]],')',sep='')
    ti.rna <- paste('Gene Expression (p=',pvals.rna,')',sep='')
    boxplots.rna <- lapply(genes, function(g){ggplot() + geom_boxplot(data=df.barplot[which(df.barplot$Gene == g),], 
                                                                      aes(x=Zygosity,y=Normalized_Counts,fill=Zygosity), position = 'dodge', show.legend = g == tail(genes,1)) + 
        theme_classic() + ggtitle(ti.rna) + scale_fill_manual(values=colorvec) + theme(axis.text.x = element_text(size = textsize, angle = 90), 
                                                                                       plot.title = element_text(hjust=.5,size=mainsize)) + 
        ylab(paste(g, ' Normalized Counts', sep='')) + ylim(0,75000)
    })
    if(length(genes) == 1){boxplots.rna <- boxplots.rna[[1]]}else if(length(genes) == 2){
      boxplots.rna <- ggpubr::ggarrange(boxplots.rna[[1]], boxplots.rna[[2]], ncol=2, nrow=1,labels = c('',''))
    }else if(length(genes) == 3){
      boxplots.rna <- ggpubr::ggarrange(boxplots.rna[[1]], boxplots.rna[[2]], boxplots.rna[[3]], ncol=3, nrow=1,labels = c('','',''))
    }
  }
  
  #Protein boxplots:
  if(length(strsplit(as.character(counts.prot$group[1]), split='_')[[1]]) < 2){counts.prot$group <- sapply(counts.prot$group, function(x){
    if(x=='wt'){'a_wt'}else if(x=='het'){'b_het'}else if(x=='biallelic'){'c_biallelic'}})}
  if(!is.null(counts.prot)){
    df.barplot.prot <- counts.prot[,c(ncol(counts.prot),match(proteins, colnames(counts.prot)))]
    df.barplot.prot <- as.data.frame(Reduce(rbind, lapply(2:ncol(df.barplot.prot),function(i){as.matrix(df.barplot.prot[,c(1,i)])})))
    df.barplot.prot$Gene <- unlist(lapply(proteins, function(g){rep(g, nrow(counts.prot))}))
    colnames(df.barplot.prot)[1:2] <- c('Zygosity','Protein_Level')
    df.barplot.prot$Protein_Level <- as.numeric(df.barplot.prot$Protein_Level)
    df.barplot.prot <- df.barplot.prot[which(!is.na(df.barplot.prot$Zygosity)),]
    df.barplot.prot <- df.barplot.prot[which(df.barplot.prot$Zygosity %in% levels),]
    df.barplot.prot$Zygosity <- levelNames[match(df.barplot.prot$Zygosity, levels)]
    df.barplot.prot$Zygosity <- factor(df.barplot.prot$Zygosity, levels = levelNames)
    
    pvals.prot <- try(sapply(proteins,function(p){signif(wilcox.test(x=df.barplot.prot$Protein_Level[which(df.barplot.prot$Gene == p & df.barplot.prot$Zygosity == levelNames[indices[1]])],
                                                                     y=df.barplot.prot$Protein_Level[which(df.barplot.prot$Gene == p & df.barplot.prot$Zygosity == levelNames[indices[2]])])$p.value,digits = 3)}), silent = TRUE)
    ti.prot <- paste('Protein Expression (p=',pvals.prot,')', sep='')
    boxplots.prot <- lapply(proteins, function(p){ggplot() + geom_boxplot(data=df.barplot.prot[which(df.barplot.prot$Gene == p),], 
                                                                          aes(x=Zygosity,y=Protein_Level,fill=Zygosity), position = 'dodge', show.legend = p == tail(proteins,1)) + 
        theme_classic() + ggtitle(ti.prot) + scale_fill_manual(values=colorvec) + theme(axis.text.x = element_text(size = textsize, angle = 90), 
                                                                                        plot.title = element_text(hjust=.5,size=mainsize)) + 
        ylab('Beta-catenin Protein Level') + ylim(-2,3)
    })
    if(length(genes) == 1){boxplots.prot <- boxplots.prot[[1]]}else if(length(genes) == 2){
      boxplots.prot <- ggpubr::ggarrange(boxplots.prot[[1]], boxplots.prot[[2]], ncol=2, nrow=1,labels = c('',''))
    }else if(length(genes) == 3){
      boxplots.prot <- ggpubr::ggarrange(boxplots.prot[[1]], boxplots.prot[[2]], boxplots.prot[[3]], ncol=3, nrow=1,labels = c('','',''))
    }
  }
  
  #Combining boxplots:
  p <- ggarrange(boxplots.rna, boxplots.prot, ncol=2, common.legend = TRUE, legend = 'bottom')
  p
}

#VUS vs. Driver expression plots for fig. 5:
expressionBoxplots2 <- function(res.vus_v_wt, res.vus_v_driver, counts.prot, counts, gene, disease, path,
                                genes_boxplot, genes_boxplot.prot, geneInTitle=FALSE, mainsize=15){
  res.vus_v_wt$lab <- res.vus_v_driver$lab <- ''
  labelled <- gene
  res.vus_v_wt$lab[which(res.vus_v_wt$hugo %in% labelled)] <- res.vus_v_wt$hugo[which(res.vus_v_wt$hugo %in% labelled)]
  res.vus_v_driver$lab[which(res.vus_v_driver$hugo %in% labelled)] <- res.vus_v_driver$hugo[which(res.vus_v_driver$hugo %in% labelled)]
  
  ###RNA expression:
  if(geneInTitle){ti <- ggtitle(paste(gene,' - ',disease,' DGE (VUS v \n WT)',sep=''))}else{ti <- ggtitle(paste(disease,' DGE (VUS v \n WT)',sep=''))}
  res.vus_v_wt$Significant <- as.character(res.vus_v_wt$padj < .05)
  res.vus_v_wt$Significant[which(res.vus_v_wt$hugo == gene)] <- gene
  scale <- c('blue','gray','red')
  names(scale) <- c('TRUE','FALSE',gene)
  
  if(geneInTitle){ti <- ggtitle(paste(gene,' - ',disease,' DGE (Drivers v \n VUS)',sep=''))}else{ti <- ggtitle(paste(disease,' DGE (Drivers v \n VUS)',sep=''))}
  res.vus_v_driver$Significant <- as.character(res.vus_v_driver$padj < .05)
  res.vus_v_driver$Significant[which(res.vus_v_driver$hugo == gene)] <- gene
  
  #RNA Box plots:
  ens_boxplot <- names(hugo)[match(genes_boxplot, hugo)]
  df.barplot <- counts[,c(ncol(counts),match(ens_boxplot, colnames(counts)))]
  df.barplot <- as.data.frame(rbind(as.matrix(df.barplot[,c(1,2)]), as.matrix(df.barplot[,c(1,3)])))
  df.barplot$Gene <- as.vector(sapply(genes_boxplot, function(g){rep(g, nrow(counts))}))
  colnames(df.barplot)[1:2] <- c('Mutation_Class','Normalized_Counts')
  df.barplot$Normalized_Counts <- as.numeric(df.barplot$Normalized_Counts)
  pvals.rna <- sapply(c(genes_boxplot),function(g){signif(wilcox.test(x=df.barplot$Normalized_Counts[which(df.barplot$Gene == g & df.barplot$Mutation_Class == 'wt')],
                                                                      y=df.barplot$Normalized_Counts[which(df.barplot$Gene == g & df.barplot$Mutation_Class == 'VUS')])$p.value,digits = 3)})
  pvals.rna.drivers <- sapply(c(genes_boxplot),function(g){signif(wilcox.test(x=df.barplot$Normalized_Counts[which(df.barplot$Gene == g & df.barplot$Mutation_Class == 'wt')],
                                                                              y=df.barplot$Normalized_Counts[which(df.barplot$Gene == g & df.barplot$Mutation_Class == 'Driver')])$p.value,digits = 3)})
  
  minim <- min(na.omit(df.barplot[which(df.barplot$Gene == genes_boxplot[1]),]$Normalized_Counts))
  maxim <- 1.2*max(na.omit(df.barplot[which(df.barplot$Gene == genes_boxplot[1]),]$Normalized_Counts))
  boxplot1 <- ggplot() + geom_boxplot(data=df.barplot[which(df.barplot$Gene == genes_boxplot[1]),], 
                                      aes(x=Gene,y=Normalized_Counts,fill=Mutation_Class), position = 'dodge', show.legend = FALSE) + ggtitle(genes_boxplot[1]) + theme_classic() + 
    ylab('Gene expression \n (Normalized Counts)') + ylim(0,20000) + xlab('') + theme(axis.text.x = element_blank()) + 
    guides(fill=guide_legend(title=paste(gene,' Mut. Class',sep=''))) + theme(plot.title = element_text(hjust=.5, size=mainsize)) + 
    annotate("text", x = 1, y = 20000, label = paste('p=',signif(pvals.rna[1]),' (WT v. VUS)', sep=''), size=4)
  
  minim <- min(na.omit(df.barplot[which(df.barplot$Gene == genes_boxplot[2]),]$Normalized_Counts))
  maxim <- 1.1*max(na.omit(df.barplot[which(df.barplot$Gene == genes_boxplot[2]),]$Normalized_Counts))
  boxplot2 <- ggplot() + geom_boxplot(data=df.barplot[which(df.barplot$Gene == genes_boxplot[2]),], 
                                      aes(x=Gene,y=Normalized_Counts,fill=Mutation_Class), position = 'dodge', show.legend = FALSE) + ggtitle(genes_boxplot[2]) + theme_classic() +
    guides(fill=guide_legend(title=paste(gene,' Mut. Class',sep=''))) + theme(plot.title = element_text(hjust=.5, size=mainsize)) + ylim(0,100000) + 
    annotate("text", x = 1, y = 100000, label = paste('p=',signif(pvals.rna[2]),' (WT v. VUS)', sep=''), size=4)
  
  #Protein boxplots:
  df.barplot.prot <- counts.prot[,c(ncol(counts.prot),match(genes_boxplot.prot, colnames(counts.prot)))]
  df.barplot.prot <- as.data.frame(Reduce(rbind, lapply(2:ncol(df.barplot.prot),function(i){as.matrix(df.barplot.prot[,c(1,i)])})))
  df.barplot.prot$Gene <- unlist(lapply(genes_boxplot.prot, function(g){rep(g, nrow(counts.prot))}))
  colnames(df.barplot.prot)[1:2] <- c('Mutation_Class','Protein_Level')
  df.barplot.prot$Protein_Level <- as.numeric(df.barplot.prot$Protein_Level)
  pvals.prot <- signif(sapply(genes_boxplot.prot,function(g){wilcox.test(x=df.barplot.prot$Protein_Level[which(df.barplot.prot$Gene == g & df.barplot.prot$Mutation_Class == 'wt')],
                                                                         y=df.barplot.prot$Protein_Level[which(df.barplot.prot$Gene == g & df.barplot.prot$Mutation_Class == 'VUS')])$p.value}),digits = 3)
  pvals.prot.driver <- signif(sapply(genes_boxplot.prot,function(g){wilcox.test(x=df.barplot.prot$Protein_Level[which(df.barplot.prot$Gene == g & df.barplot.prot$Mutation_Class == 'wt')],
                                                                                y=df.barplot.prot$Protein_Level[which(df.barplot.prot$Gene == g & df.barplot.prot$Mutation_Class == 'Driver')])$p.value}),digits = 3)
  df.barplot.prot <- df.barplot.prot[which(!is.na(df.barplot.prot$Mutation_Class)),]
  minim <- min(na.omit(df.barplot.prot[which(df.barplot.prot$Gene == genes_boxplot.prot[1]),]$Protein_Level))
  maxim <- 1.1*max(na.omit(df.barplot.prot[which(df.barplot.prot$Gene == genes_boxplot.prot[1]),]$Protein_Level))
  boxplot3 <- ggplot() + geom_boxplot(data=df.barplot.prot[which(df.barplot.prot$Gene == genes_boxplot.prot[1]),], 
                                      aes(x=Gene,y=Protein_Level,fill=Mutation_Class), position = 'dodge', show.legend = TRUE) + theme(plot.title = element_blank()) + theme_classic() + 
    xlab('') + ylab('Protein Level') + guides(fill=guide_legend(title=paste(gene,' Mut. Class',sep=''))) + ylim(minim,maxim) + 
    annotate("text", x = 1, y = maxim, label = paste('p=',signif(pvals.prot[1]),' (WT v. VUS)', sep=''), size=4)
  
  minim <- min(na.omit(df.barplot.prot[which(df.barplot.prot$Gene == genes_boxplot.prot[1]),]$Protein_Level))
  maxim <- 1.1*max(na.omit(df.barplot.prot[which(df.barplot.prot$Gene == genes_boxplot.prot[2]),]$Protein_Level))
  boxplot4 <- ggplot() + geom_boxplot(data=df.barplot.prot[which(df.barplot.prot$Gene == genes_boxplot.prot[2]),], 
                                      aes(x=Gene,y=Protein_Level,fill=Mutation_Class), position = 'dodge', show.legend = TRUE) + theme(plot.title = element_blank()) + theme_classic() + 
    xlab('') + ylab('Protein Level') + guides(fill=guide_legend(title=paste(gene,' Mut. Class',sep=''))) + ylim(minim,maxim) + 
    annotate("text", x = 1, y = maxim, label = paste('p=',signif(pvals.prot[2]),' (WT v. VUS)', sep=''), size=4)
  
  #Making p-value etc. summary table:
  names(pvals.prot) <- names(pvals.prot.driver) <- genes_boxplot.prot
  names(pvals.rna) <- names(pvals.rna.drivers) <- genes_boxplot
  stats <- list('protein_vus'=pvals.prot, 'rna_vus'=pvals.rna, 'protein_drivers'=pvals.prot.driver, 'rna_drivers'=pvals.rna.drivers)
  plots <- ggarrange(boxplot1, boxplot2, boxplot3, boxplot4,ncol=2,nrow=2,common.legend = TRUE,legend = 'right',
                     labels = c('','',''))
  list('plots'=plots, 'stats'=stats)
}

#Co-mutation plots for figure 5.
mutationBarplot.class <- function(data, gene, disease, excludeHet=TRUE, cutoff=.05, textSize=12, geneset=NULL, dataset=NULL, returnTable=FALSE){
  mutcats <- c('Biallelic - Mut + LOH','Biallelic - Mut + fusion','Biallelic - compound','Heterozygous (Mutation)','Gain-of-mutant - Mut + copy gain')
  data$mutated <- data$zygosity_call %in% mutcats
  data$zygosity_call <- as.character(data$zygosity_call)
  #Het LOH cases will be grouped with WT for this analysis:
  data$zygosity_call[which(data$zygosity_call %in% c('Amplification - CNA','Heterozygous - LOH (gene-level)') | data$zygosity_call == 'wt' | is.na(data$zygosity_call))] <- 'WT'
  data$zygosity_call <- factor(data$zygosity_call, levels = unique(rev(c(mutcats,'WT'))))
  data <- data[which(data$disease_subtype == disease),]
  xgeneData <- unique(data[which(data$Hugo_Symbol == gene),c('Tumor_Sample_Barcode','mutation_class','zygosity_call')])
  xgeneData$mutation_class[which(xgeneData$zygosity_call == 'WT')] <- 'WT'
  if(excludeHet){xgeneData$mutation_class[which(xgeneData$zygosity == 'het')] <- NA}
  xgeneData <- xgeneData[,c('Tumor_Sample_Barcode','mutation_class')]
  colnames(xgeneData)[2] <- 'mutation_class_x'
  xgeneData$mutation_class_x <- factor(xgeneData$mutation_class_x, levels = c('WT','VUS','Driver'))
  data <- data[which(data$zygosity_call %in% c(mutcats,'WT')),]
  data <- merge(data, xgeneData, by='Tumor_Sample_Barcode')
  tab <- table(data[which(data$Hugo_Symbol != gene),c('Hugo_Symbol','mutation_class_x','mutated')])
  tab <- data.frame('gene'=dimnames(tab)[1], 'n_driver'=rowSums(tab[,'Driver',]),'mutated_driver'=tab[,'Driver','TRUE'], 
                    'n_vus'=rowSums(tab[,'VUS',]),'mutated_vus'=tab[,'VUS','TRUE'],'n_wt'=rowSums(tab[,'WT',]),'mutated_wt'=tab[,'WT','TRUE'])
  geneOrder <- tab$Hugo_Symbol[order(tab$mutated_driver + tab$mutated_vus + tab$mutated_wt + tab$mutated_silent, decreasing = TRUE)]
  
  #Stratified by mutation class:
  tab2 <- Reduce(rbind, lapply(unique(data$Hugo_Symbol), function(g){
    cbind('Hugo_Symbol'=g,as.data.frame(table(data[which(data$Hugo_Symbol == g),c('mutation_class_x','zygosity_call')])))}))
  tab2 <- merge(tab2, tab, by=c('Hugo_Symbol'))
  tab2$Fraction <- NA
  tab2$Fraction[which(tab2$mutation_class_x == 'WT')] <- tab2$Freq[which(tab2$mutation_class_x == 'WT')]/tab2$n_wt[which(tab2$mutation_class_x == 'WT')]
  tab2$Fraction[which(tab2$mutation_class_x == 'Silent')] <- tab2$Freq[which(tab2$mutation_class_x == 'Silent')]/tab2$n_silent[which(tab2$mutation_class_x == 'Silent')]
  tab2$Fraction[which(tab2$mutation_class_x == 'VUS')] <- tab2$Freq[which(tab2$mutation_class_x == 'VUS')]/tab2$n_vus[which(tab2$mutation_class_x == 'VUS')]
  tab2$Fraction[which(tab2$mutation_class_x == 'Driver')] <- tab2$Freq[which(tab2$mutation_class_x == 'Driver')]/tab2$n_driver[which(tab2$mutation_class_x == 'Driver')]
  
  #Getting confidence intervals for overall mutation frac.
  tabDriver <- tab[,c(c('Hugo_Symbol','n_driver','mutated_driver'))]
  tabVUS <- tab[,c(c('Hugo_Symbol','n_vus','mutated_vus'))]
  tabWT <- tab[,c(c('Hugo_Symbol','n_wt','mutated_wt'))]
  colnames(tabDriver)[2:3] <- colnames(tabVUS)[2:3] <- colnames(tabWT)[2:3] <- c('n','mutated')
  tabDriver$Fraction <- tabDriver$mutated/tabDriver$n
  tabVUS$Fraction <- tabVUS$mutated/tabVUS$n
  tabWT$Fraction <- tabWT$mutated/tabWT$n
  tabDriver <- cbind(tabDriver, t(sapply(1:nrow(tabDriver),function(i){
    ci <- try(prop.test(n=tabDriver$n[i], x=tabDriver$mutated[i])$conf.int, silent = TRUE)
    if(class(ci) == 'try-error'){c(NA,NA)}else{ci}})))
  tabVUS <- cbind(tabVUS, t(sapply(1:nrow(tabVUS),function(i){
    ci <- try(prop.test(n=tabVUS$n[i], x=tabVUS$mutated[i])$conf.int, silent = TRUE)
    if(class(ci) == 'try-error'){c(NA,NA)}else{ci}})))
  tabWT <- cbind(tabWT, t(sapply(1:nrow(tabWT),function(i){
    ci <- try(prop.test(n=tabWT$n[i], x=tabWT$mutated[i])$conf.int, silent = TRUE)
    if(class(ci) == 'try-error'){c(NA,NA)}else{ci}})))
  
  #Since we only care about alteration fraction, not specific zygosity calls:
  tabDriver$mutation_class_x <- 'Driver'
  tabVUS$mutation_class_x <- 'VUS'
  tabWT$mutation_class_x <- 'WT'
  merged <- rbind(tabDriver, tabVUS, tabWT)
  colnames(merged)[5:6] <- c('Lower','Upper')
  merged$pvalue_driver_vus <- sapply(unique(merged$Hugo_Symbol),function(g){
    pval <- try(prop.test(n=merged$n[which(merged$Hugo_Symbol == g & merged$mutation_class_x != 'WT')], 
                          x=merged$mutated[which(merged$Hugo_Symbol == g & merged$mutation_class_x != 'WT')])$p.value,silent=TRUE)
    if(class(pval) == 'try-error'){NA}else{pval}})
  merged$pvalue_vus_wt <- sapply(unique(merged$Hugo_Symbol),function(g){
    pval <- try(prop.test(n=merged$n[which(merged$Hugo_Symbol == g & merged$mutation_class_x != 'Driver')], 
                          x=merged$mutated[which(merged$Hugo_Symbol == g & merged$mutation_class_x != 'Driver')])$p.value,silent=TRUE)
    if(class(pval) == 'try-error'){NA}else{pval}})
  merged$total_mutated <- sapply(merged$Hugo_Symbol, function(g){sum(merged$mutated[which(merged$Hugo_Symbol == g)])})
  merged$total_samples <- nrow(xgeneData)
  if(cutoff < 1 & cutoff > 0){
    merged <- merged[which(merged$total_mutated/merged$total_samples >= cutoff),]
  }else{merged <- merged[which(merged$total_mutated >= cutoff),]}
  fdr <- p.adjust(c(merged$pvalue_driver_vus, merged$pvalue_vus_wt), method = 'fdr')
  merged$qvalue_driver_vus <- fdr[1:(length(fdr)/2)]
  merged$qvalue_vus_wt <- fdr[(length(fdr)/2 + 1):length(fdr)]
  if(is.null(geneset)){
    merged <- merged[which(merged$qvalue_driver_vus < .05 | merged$qvalue_vus_wt < .05),]
    merged <- merged[order(merged$total_mutated, decreasing = TRUE),]
    merged$Hugo_Symbol <- factor(merged$Hugo_Symbol, levels = unique(merged$Hugo_Symbol))
  }else{
    merged <- merged[which(merged$Hugo_Symbol %in% geneset),]
    merged$Hugo_Symbol <- factor(merged$Hugo_Symbol, levels = geneset)
  }
  merged$Percent <- merged$Fraction*100
  merged$Upper <- merged$Upper*100
  merged$Lower <- merged$Lower*100
  colorvec <- c('Driver'='#6795ad','VUS'='#a8c956','WT'='#b4b4b4')
  p <- ggplot(data=merged, aes(x=Hugo_Symbol, y=Percent, fill=mutation_class_x)) + geom_bar(stat = 'identity',position = position_dodge()) + 
    geom_errorbar(mapping=aes(ymin=Lower, ymax=Upper), width=0.2, size=1, position = position_dodge(.9)) + theme_classic() + scale_fill_manual(values = colorvec) +
    ggtitle(paste(gene , ' - ',disease,sep='')) + xlab(gene) + ylab('Alteration rate') + guides(fill=guide_legend(title="Mutation status")) + 
    theme(plot.title = element_blank(), axis.text.x = element_text(angle = 90, hjust=1,vjust=.5))
  if(returnTable){list('plot'=p,'table'=merged)}else{p}
}

