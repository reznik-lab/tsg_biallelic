### Loading necessary packages and sourcing necessary functions:
library(maftools)
library(readxl)
library(stringr)
library(ggpubr)
library(patchwork)
library(ggrepel)
library(dplyr)

#Loading a few auxiliary functions needed:
source('code/auxiliary_functions.R')

#Figure 5: VUSs:
#table_s8 should be used for enrichment_driver_vus
generateFigure5 <- function(enrichment_driver_vus, keap1_driver_maf, keap1_vus_maf, keap1_mut_tab, expr_res_luad_keap1, 
                            drivers_and_vus.luad, mainsize=15){
  #Merging drivers and vuses:
  mutcats <- c('Heterozygous (Mutation)','Biallelic - Mut + LOH','Gain-of-mutant - Mut + copy gain','Biallelic - compound','Biallelic - Mut + fusion')
  drivers_and_vus <- as.data.frame(rbind(drivers, 
                           vus %>% 
                             filter(grepl('Mutation|Biallelic|Gain', zygosity_call)) %>%
                             filter(!Hugo_Symbol %in% blacklist$Hugo_Symbol)))
  
  ###Plotting driver enrichment against vus enrichment
  enrdat <- enrichment_driver_vus
  enrdat$Mutations <- enrdat$drivers_num_MutOnly + enrdat$driver_num_MutLOH + enrdat$vus_num_MutOnly + enrdat$vus_num_MutLOH
  colnames(enrdat)[c(6,10)] <- c('P_Value_Driver','P_Value_VUS')
  enrdat$label <- sapply(enrdat$Disease, function(d){
    str <- strsplit(d, split='[ ]')[[1]]
    l <- length(str)
    if(l == 1){out <- d}else{
      a <- floor(l/2)
      paste(paste(str[1:a], collapse=' '),'\n', paste(str[(a+1):length(str)],collapse=' '), collapse='')
    }
  })
  
  enrdat <-
    enrdat %>%
    mutate(log10p_driver = -log10(P_Value_Driver),
           log10p_vus = -log10(P_Value_VUS)) %>%
    mutate(log10p_driver_sig = P_Value_Driver < 0.05,
           log10p_vus_sig = P_Value_VUS < 0.05) %>% 
    filter(Gene %in% ((enrdat %>% filter(log10p_vus_sig, vus_num_MutLOH + vus_num_MutOnly >= 5))$Gene)) %>%
    mutate(to_show = (log10p_vus_sig | (log10p_driver > 18))) 
  
  enrdat <-
    enrdat %>% 
    left_join(select_ctds %>% select(Disease=CANCER_TYPE_DETAILED, ONCOTREE_CODE)) %>%
    filter(!is.na(ONCOTREE_CODE))
  
  enrdat$vus_OR <- as.numeric(enrdat$vus_OR)
  enrdat <-
    enrdat %>%
    mutate(display_text = paste0(ONCOTREE_CODE, ' - ', Gene)) %>%
    mutate(display_text = paste0(ONCOTREE_CODE)) %>%
    mutate(logOR_vus = log(vus_OR)) %>%
    mutate(logOR_vus = ifelse(logOR_vus > 5, 6, logOR_vus))

  ### Plotting vus enrichment p-values against driver enrichment p-values
  #select_genes <- c('KEAP1','STK11',names(sort(table(enrdat[which(!enrdat$Gene %in% c('KEAP1','STK11')),c('Gene','log10p_vus_sig')])[,'TRUE'], 
  #                                             decreasing = TRUE)[1:10]))
  select_genes <- c('ATM','BAP1','CASP8','CBFB','DAXX','KEAP1','MEN1','NF2','PTEN','STK11','TSC2')
  enrdat <- enrdat[which(enrdat$Gene %in% select_genes),]
  f5a <- ggplot(enrdat, aes(x = log10p_vus, y = log10p_driver, fill=Gene, color=Gene, size=logOR_vus)) + 
    geom_hline(yintercept = 1.31, linewidth=1, color='darkred', linetype='dotted') +
    geom_vline(xintercept = 1.31, linewidth=1, color='darkred', linetype='dotted') +
    geom_point(pch=21) +
    geom_text_repel(data=enrdat %>% filter(to_show), 
                    aes(x=log10p_vus, y=log10p_driver, label=display_text), 
                    size=6, show.legend = F, box.padding = unit(1, "lines"),
                    max.overlaps = 17) +
    scale_fill_manual(values=c(ggthemes::tableau_color_pal(palette= 'Tableau 20')(20) [c(1,3,5,7,9,11,13,15,17,19,8,12)])) +
    scale_color_manual(values=c(ggthemes::tableau_color_pal(palette= 'Tableau 20')(20) [c(1,3,5,7,9,11,13,15,17,19,8,12)])) +
    scale_y_continuous(breaks=c(0:10)* 20) +
    scale_x_continuous(breaks=c(0:10)* 15) +
    theme(axis.text=element_text(size=14),
          axis.title.x  = element_text(size=14, face='plain'),
          axis.title.y  = element_text(size=14, face='plain'),
          axis.text.x = element_text(hjust=0.5, vjust=0.5),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(face = 'bold',hjust = 0.5),
          panel.grid.major  = element_line(colour=NA),
          panel.grid.minor = element_line(colour=NA),
          panel.spacing = unit(0.2, "lines"),
          panel.background = element_rect(fill=NA, color=NA),
          legend.text = element_text(size = 14, face="bold"),
          legend.position="right",
          legend.justification = c(0, 1),
          legend.direction="vertical",
          legend.key = element_rect(fill = "transparent", colour = "transparent")
    ) +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    xlab("-log10(p-value VUS)") +
    ylab("-log10(p-value driver)")
  
  #Lollipop plot of KEAP1 driver and VUS mutations (stored as mafs)
  #f5b, generated separately using maftools:
  lollipopPlot2(m1 = keap1_driver_maf, m2 = keap1_vus_maf,gene = 'KEAP1', m1_name = 'Drivers',m2_name = 'VUSs')
  
  ###LUAD KEAP1 barplot:
  f5ca <- ggplot(aes(x=Class,y=LOH_Fraction, fill=Class), data=keap1_mut_tab[keap1_mut_tab$Dataset == 'MSK-IMPACT',]) + 
    geom_bar(stat = 'identity', position = 'stack') + theme_bw() + scale_fill_manual(values=c('Driver'='#ED6A50','VUS'='#F9CBA8','WT'='#7284A1')) + 
    ylab('% with LOH') + ggtitle('') + theme(legend.position="none") + xlab('MSK-IMPACT') + 
    theme(axis.text.x= element_blank()) + 
    geom_signif(comparisons = list(c("Driver", "VUS"),c("Driver","WT")), map_signif_level=TRUE, annotations = c('***','***')) + 
    geom_linerange(stat = 'identity',data=keap1_mut_tab[keap1_mut_tab$Dataset == 'MSK-IMPACT',], 
                   mapping=aes(ymin=Upper, ymax=Lower), size=1)
  
  #KEAP1 TCGA:
  f5cb <- ggplot(aes(x=Class,y=LOH_Fraction, fill=Class), data=keap1_mut_tab[keap1_mut_tab$Dataset == 'TCGA',]) + 
    geom_bar(stat = 'identity', position = 'stack') + theme_bw() + scale_fill_manual(values=c('Driver'='#ED6A50','VUS'='#F9CBA8','WT'='#7284A1')) + 
    ylab('% with LOH') + xlab('TCGA') + ggtitle('TCGA') + theme(legend.position="none") + 
    theme(axis.text.x= element_blank()) + 
    geom_signif(comparisons = list(c("Driver", "VUS"),c("Driver","WT")), map_signif_level=TRUE, annotations = c('***','***')) + 
    geom_linerange(stat = 'identity',data=keap1_mut_tab[keap1_mut_tab$Dataset == 'TCGA',],
                   mapping=aes(ymin=Upper, ymax=Lower), size=1)
  
  f5c <- ggarrange(f5ca, f5cb, nrow=2, common.legend = TRUE)
  
  #Expression (RNA and protein) plots:
  f5d <- expressionBoxplots2(res.vus_v_wt=expr_res_luad_keap1$res.vus_v_wt, res.vus_v_driver=expr_res_luad_keap1$res.vus_v_driver,gene='KEAP1', 
              disease='LUAD', path=NULL, counts = expr_res_luad_keap1$counts, counts.prot=expr_res_luad_keap1$counts.prot, 
              genes_boxplot=c('NFE2L2','NQO1'), genes_boxplot.prot=c('NRF2','NQO1'))$plots
  
  #LUAD comutation plot:
  f5e <- mutationBarplot.class(data=drivers_and_vus.luad, gene='KEAP1', disease='Lung Adenocarcinoma',cutoff = .05)
  
  #Combining panels:
  layout <- 
    "
  AAAA#####
  AAAABBCCC
  DDDD#####
  "
  f5 <- f5a + f5c + f5d + f5e + 
    plot_layout(design = layout, heights = c(.75,.75,1.1,1)) + plot_annotation(tag_levels = 'A')
  print(f5)
}


