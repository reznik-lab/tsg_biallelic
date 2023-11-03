### Loading necessary packages and sourcing necessary functions:
library(maftools)
library(readxl)
library(stringr)
library(RColorBrewer)
library(ggpubr)
library(DESeq2)
library(patchwork)
library(ggrepel)

generateFigure4a <- function(enrichmentRes, subtypeTable, mainsize=14, textsize=12){
  select_genes_4a = c("APC", "BAP1", "CDH1", "MEN1", "NF1", "NF2", "SMAD4", "STK11")
  if(is.null(enrichmentRes$ONCOTREE_CODE)){enrichmentRes$ONCOTREE_CODE <- subtypeTable$ONCOTREE_CODE[match(enrichmentRes$Disease, subtypeTable$CANCER_TYPE_DETAILED)]}
  enrichmentRes$alt_rate <- (enrichmentRes$Num_MutLOH + enrichmentRes$Num_Mut_only)/enrichmentRes$N
  colnames(enrichmentRes)[which(colnames(enrichmentRes) == 'mut_fraction')] <- 'alt_rate'
  colnames(enrichmentRes)[which(colnames(enrichmentRes) == 'gene')] <- 'Gene'
  colnames(enrichmentRes)[which(colnames(enrichmentRes) == 'disease')] <- 'Disease'
  enrichmentRes$neg_log10_p_value <- -log10(enrichmentRes$p_value_corrected)

  d4a <- 
    enrichmentRes %>% 
    mutate(gene = Gene, disease = Disease) %>% 
    filter(gene %in% select_genes_4a) %>% 
    mutate(display_text = paste0(Gene, ' - ', ONCOTREE_CODE),
           is_significant = !is.na(p_value_corrected)  & p_value_corrected < 0.05) %>%
    filter(Disease != 'all') %>% 
    mutate(log10p = neg_log10_p_value,
           log_or = log_or) %>%
    mutate(log_or = ifelse(log_or > 4.8, 5, 
                          ifelse(log_or < -2, -2, log_or))) %>%
    mutate(to_display = (log10p > 5 | (alt_rate > 0.1 & p_value_corrected < 0.05))) %>%
    mutate(Gene = ifelse(is.na(is_significant) | !is_significant, 'zzzDummy', Gene))
  
  d4a$Gene = factor(d4a$Gene, levels=sort(d4a$Gene %>% unique))
  
  # `palette` must be one of Tableau 10, Tableau 20, Color Blind, Seattle Grays, Traffic, Miller Stone, 
  # Superfishel Stone, Nuriel Stone, Jewel Bright, Summer, Winter, Green-Orange-Teal, Red-Blue-Brown, Purple-Pink-Gray, 
  # Hue Circle, Classic 10, Classic 10 Medium, Classic 10 Light, Classic 20, Classic Gray 5, Classic Color Blind, 
  # Classic Traffic Light, Classic Purple-Gray 6, Classic Purple-Gray 12, Classic Green-Orange 6, Classic Green-Orange 12, 
  # Classic Blue-Red 6, Classic Blue-Red 12, Classic Cyclic.

  d4a_to_plot <- d4a[which(d4a$Gene != 'zzzDummy'),]
  f4a <-
    ggplot(d4a_to_plot %>% filter(!is_significant),
         aes(x = alt_rate, y = log10p, fill=Gene, color=Gene, size=log_or)) +
    geom_hline(yintercept = 1.31, linewidth=1, linetype='dotted') +
    geom_point(pch=21, alpha=0.5) +
    geom_point(data=d4a_to_plot %>% filter(is_significant), pch=21, alpha=1, color='black') +
    geom_text_repel(data=d4a_to_plot, 
                    aes(x=alt_rate, y=log10p, label=display_text), 
                    size=6, show.legend = F, box.padding = unit(1, "lines"),
                    max.overlaps = 17) +
    scale_fill_manual(values=c(ggthemes::tableau_color_pal(palette= 'Tableau 20')(14))) +
    scale_color_manual(values=c(ggthemes::tableau_color_pal(palette= 'Tableau 20')(14))) +
    scale_y_continuous(breaks=c(0:20)* 10) +
    scale_x_continuous(breaks=c(0:10)* 0.1) +
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
    xlab("Alteration rate") +
    ylab("-log10(adj. p-value)")
    
  d4a_to_plot <- d4a %>% filter(alt_rate < 0.1 & log10p >= 1 & log10p < 5) %>% mutate(to_display=T) %>% 
    filter(Gene != 'zzzDummy')
  f4a_inset <-
    ggplot(d4a_to_plot %>% filter(!is_significant),
           aes(x = alt_rate, y = log10p, fill=Gene, color=Gene, size=log_or)) +
      geom_point(pch=21, alpha=0.5) +
      geom_point(data=d4a_to_plot %>% filter(is_significant), pch=21, alpha=1, color='black') +
      geom_text_repel(data=d4a_to_plot %>% filter(to_display), 
                      aes(x=alt_rate, y=log10p, label=display_text), 
                      size=6, show.legend = F, box.padding = unit(1, "lines"),
                      max.overlaps = 17) +
      scale_fill_manual(values=c(ggthemes::tableau_color_pal(palette= 'Tableau 20')(14))) +
      scale_color_manual(values=c(ggthemes::tableau_color_pal(palette= 'Tableau 20')(14))) +
      scale_y_continuous(breaks=c(0:10), limits=c(1,5)) +
      scale_x_continuous(breaks=c(0:10)* 0.01) +
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
      xlab("mutation rate") +
      ylab("-log10(adj. p-value)")
  
  f4a
}

#Figure 3: Rare TSGs:
generateFigure4 <- function(apc_mut_tab, enrichmentRes, expr.res.luad, subtypeTable, mainsize=15, textsize=14){
  #f4a:
  ###Selection vs. mut. rate:
  f4a <- generateFigure4a(enrichmentRes=enrichmentRes, subtypeTable = subtypeTable, mainsize=14, textsize=12)
  
  ### Plotting biallelic frequency in TCGA vs. impact (LUAD and PRAD):
  p1 <- ggplot(data = apc_mut_tab[apc_mut_tab$Disease == 'LUAD' & apc_mut_tab$Dataset == 'IMPACT',], aes(x=Class, y=Fraction)) + 
    geom_bar(stat = 'identity', position = 'dodge', fill='black') + theme_classic() + ggtitle('') + 
    geom_signif(comparisons = list(c("Mutated", "WT")), map_signif_level=TRUE, annotations = '***') + 
    theme(legend.position='none', axis.text = element_text(size = textsize, angle = 90), axis.title = element_text(size = textsize)) + 
    ylab('LOH Fraction') + xlab('IMPACT') + ylim(0,1) + geom_linerange(stat = 'identity',
        data=apc_mut_tab[apc_mut_tab$Disease == 'LUAD' & apc_mut_tab$Dataset == 'IMPACT',], mapping=aes(ymin=Upper, ymax=Lower), size=1,col='gray47')
  p2 <- ggplot(data = apc_mut_tab[apc_mut_tab$Disease == 'LUAD' & apc_mut_tab$Dataset == 'TCGA',], aes(x=Class, y=Fraction)) + 
    geom_bar(stat = 'identity', position = 'dodge', fill='#6F8093') + theme_classic() + ggtitle('') + 
    geom_signif(comparisons = list(c("Mutated", "WT")), map_signif_level=TRUE, annotations = 'NS.') + 
    theme(legend.position='none', axis.text = element_text(size = textsize, angle = 90), axis.title = element_text(size = textsize)) + 
    ylab('') + xlab('TCGA') + ylim(0,1) + geom_linerange(stat = 'identity',
         data=apc_mut_tab[apc_mut_tab$Disease == 'LUAD' & apc_mut_tab$Dataset == 'TCGA',], mapping=aes(ymin=Upper, ymax=Lower), size=1,col='black')
  p3 <- ggplot(data = apc_mut_tab[apc_mut_tab$Disease == 'PRAD' & apc_mut_tab$Dataset == 'IMPACT',], aes(x=Class, y=Fraction)) + 
    geom_bar(stat = 'identity', position = 'dodge', fill='black') + theme_classic() + ggtitle('') + 
    geom_signif(comparisons = list(c("Mutated", "WT")), map_signif_level=TRUE, annotations = '**') +
    theme(legend.position='none', axis.text = element_text(size = textsize, angle = 90), axis.title = element_text(size = textsize)) + 
    ylab('') + xlab('IMPACT') + ylim(0,1) + geom_linerange(stat = 'identity',
          data=apc_mut_tab[apc_mut_tab$Disease == 'PRAD' & apc_mut_tab$Dataset == 'IMPACT',], mapping=aes(ymin=Upper, ymax=Lower), size=1,col='gray47')
  p4 <- ggplot(data = apc_mut_tab[apc_mut_tab$Disease == 'PRAD' & apc_mut_tab$Dataset == 'TCGA',], aes(x=Class, y=Fraction)) + 
    geom_bar(stat = 'identity', position = 'dodge', fill='#6F8093') + theme_classic() + ggtitle('') + 
    geom_signif(comparisons = list(c("Mutated", "WT")), map_signif_level=TRUE, annotations = '**') + 
    theme(legend.position='none', axis.text = element_text(size = textsize, angle = 90), axis.title = element_text(size = textsize)) + 
    ylab('') + xlab('TCGA') + ylim(0,1) + geom_linerange(stat = 'identity',
           data=apc_mut_tab[apc_mut_tab$Disease == 'PRAD' & apc_mut_tab$Dataset == 'TCGA',], mapping=aes(ymin=Upper, ymax=Lower), size=1,col='black')
  
  f4b <- ggarrange(annotate_figure(ggarrange(p1, p2, ncol=2), top = text_grob("LUAD", size = 20)),
            annotate_figure(ggarrange(p3, p4, ncol=2), top = text_grob("PRAD", size = 20)), ncol=2)
  
  #Expression (RNA and protein) boxplot for APC-associated genes:
  f4c <- expressionBoxplots(expr.res.luad$counts, annotation=expr.res.luad$columnData$category, genes=c('CTNNB1'), counts.prot=expr.res.luad$counts.prot, 
        proteins=c('BETACATENIN'), hugo,textsize=12, levels=c('a_wt_wt','b_het_wt','c_wt_mut','d_biallelic_wt'), 
        levelNames = c('APC-WT/CTNNB1-WT','APC-Het/CTNNB1-WT','APC-WT/CTNNB1-Mut','APC-Biallelic/CTNNB1-WT'), remove=NULL)

  
  #Combining panels:
  layout <- 
    "
  ABC
  "
  f4 <- f4a + f4b + f4c  + 
    plot_layout(design = layout, widths = c(4,3,3)) + plot_annotation(tag_levels = 'A')
  print(f4)
}

