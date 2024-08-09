
load("f2_input.Rdata")

################################
#### Gene alteration rates panel
################################
{
  f2_alt_rates_to_plot <-
    f2_alt_rates %>%
    filter(Hugo_Symbol %in% f2_select_genes) %>%
    mutate(Hugo_Symbol=factor(Hugo_Symbol, f2_select_genes))
  
  
  f2_alt_rates_to_plot$gene_label = factor(f2_alt_rates_to_plot$gene_label, 
                                           levels=(f2_alt_rates_to_plot %>% 
                                                     ungroup %>% 
                                                     select(Hugo_Symbol, Pathway, gene_label) %>%
                                                     unique %>%
                                                     arrange(Pathway, Hugo_Symbol))$gene_label)
  
  f2_alteration_rate_plot <-
    ggplot(f2_alt_rates_to_plot %>% 
             mutate(gene_alt_rate = pmin(gene_alt_rate, 0.15)),  
           aes(x=gene_label, y=gene_alt_rate)) + 
    geom_bar(stat='identity', fill='black') + 
    facet_grid(. ~ Pathway, scales = "free", space = "free") +
    theme(axis.text=element_text(size=10, color='black'),
          #axis.text.x =  element_text(size=14,face='plain', angle=90), #element_blank(),
          axis.text.x =  element_blank(),
          axis.title.x  = element_text(size=14,face='plain'),
          axis.title.y  = element_text(size=14,face='plain'),
          plot.title = element_text(face = 'bold',hjust = 0.5),
          panel.grid.major = element_line(colour=NA),
          panel.grid.minor = element_line(colour=NA),
          panel.background = element_rect(fill=NA, color='black'),
          panel.spacing = unit(0.2, "lines"),
          legend.text = element_text(size = 14, face="bold"),
          legend.position="right",
          legend.justification = c(0, 1),
          legend.direction="vertical",
          strip.background = element_rect(colour=NA, fill=NA),
          strip.text = element_text(size=11, face='bold', 
                                    vjust= 0.5, hjust=0, angle=90)) +
    xlab("") +
    ylab("")
}

################################
#### Heatmap panel
################################
{
  biallelic_heatmap_colors <- colorRampPalette(c("white", "#FFFB66", "#EF273E"))(4)
  
  f2_biallelic_rates_to_plot <-
    f2_biallelic_rates %>% 
    filter(CANCER_TYPE_DETAILED %in% select_ctds_to_show,
           Hugo_Symbol %in% f2_select_genes) %>% 
    mutate(Pathway = factor(Pathway, levels=pathway_order)) 
  
  f2_biallelic_rates_to_plot$gene_label = 
    factor(f2_biallelic_rates_to_plot$gene_label, levels=(f2_alt_rates_to_plot$gene_label %>% sort))
  
  f2_biallelic_rates_matrix <-
    cross_join(f2_biallelic_rates_to_plot %>% select(CANCER_TYPE_DETAILED, organ, ctd_abbrev_n) %>% unique, 
               f2_biallelic_rates_to_plot %>% select(Pathway, gene_label) %>% unique) %>%
    mutate(f_biallelic_rate = 0,
           f_alteration_rate_txt = "") 
  
  ctd_order = (f2_biallelic_rates_to_plot %>% arrange((N_ctd)))$ctd_abbrev_n %>% unique
  
  f2_biallelic_rates_matrix$organ = factor(f2_biallelic_rates_matrix$organ, levels=organ_order)
  f2_biallelic_rates_matrix$ctd_abbrev_n = factor(f2_biallelic_rates_matrix$ctd_abbrev_n, levels=ctd_order)
  
  f2_biallelic_rates_to_plot$organ = factor(f2_biallelic_rates_to_plot$organ, levels=organ_order)
  f2_biallelic_rates_to_plot$ctd_abbrev_n = factor(f2_biallelic_rates_to_plot$ctd_abbrev_n, levels=ctd_order)
  
  f2_heatmap_plot <-
    ggplot(f2_biallelic_rates_matrix,
           aes(x=gene_label, y=ctd_abbrev_n, fill=f_biallelic_rate)) + 
    geom_tile(fill='white', color=NA)  + 
    geom_tile(data=f2_biallelic_rates_to_plot %>% filter(show_alteration_rate)) +
    geom_text(data=f2_biallelic_rates_to_plot %>% filter(show_alteration_rate), 
              aes(label = f_alteration_rate_txt), 
              color = "black", size = 4) + 
    geom_point(data=f2_biallelic_rates_to_plot %>% filter(show_diamond), 
               aes(x=gene_label, y=ctd_abbrev_n), 
               color='gray70', size=0.7, shape=5) +
    scale_fill_gradientn(colours = biallelic_heatmap_colors,
                         values = seq(from=0, to=1, length=4)) + 
    facet_grid(organ ~ Pathway, scales = "free", space = "free") +
    theme(axis.text=element_text(size=10, color='black'),
          #axis.text.x = element_text(size=14,face='plain', angle=90), #element_blank(),
          axis.text.x =  element_blank(),
          axis.title.x  = element_text(size=14,face='plain'),
          axis.title.y  = element_text(size=14,face='plain'),
          #axis.line = element_line(colour = "black"),
          plot.title = element_text(face = 'bold',hjust = 0.5),
          panel.grid.major = element_line(colour=NA),
          panel.grid.minor = element_line(colour=NA),
          panel.background = element_rect(fill=NA, color='black'),
          panel.spacing = unit(0.2, "lines"),
          legend.text = element_text(size = 14, face="bold"),
          legend.position="right",
          legend.justification = c(0, 1),
          legend.direction="vertical",
          strip.background = element_rect(colour=NA, fill=NA),
          strip.text = element_blank()) +
    xlab("") +
    ylab("")
}


################################
#### biallelic mechanisms plot
################################
{
  colorvec <- c('Biallelic - Mut + LOH' = "#58B9F4",
                'Biallelic - Homdel'="#242860", 
                'Biallelic - compound'= "#8CA1D3",
                'Heterozygous'="gray80")
  
  f2_mechanism_rates_to_plot <-
    f2_mechanism_rates %>%
    filter(Hugo_Symbol %in% f2_select_genes) %>% 
    mutate(Pathway = factor(Pathway, levels=pathway_order)) %>%
    mutate(mechanism = factor(mechanism, levels= rev(colorvec %>% names)))
  
  f2_mechanism_rates_to_plot$gene_label = 
    factor(f2_mechanism_rates_to_plot$gene_label, levels=(f2_alt_rates_to_plot$gene_label %>% sort))
  
  f2_mechanisms_plot <-
    ggplot(f2_mechanism_rates_to_plot, aes(x=gene_label, y=f_mutated_mechansim, fill=as.factor(mechanism))) + 
    geom_bar(stat='identity') + 
    scale_fill_manual(values=colorvec) +
    facet_grid(. ~ Pathway, scales = "free", space = "free") +
    theme(axis.text=element_text(size=10, color='black'),
          axis.text.x =  element_blank(),
          axis.title.x  = element_text(size=14,face='plain'),
          axis.title.y  = element_text(size=14,face='plain'),
          plot.title = element_text(face = 'bold',hjust = 0.5),
          panel.grid.major = element_line(colour=NA),
          panel.grid.minor = element_line(colour=NA),
          panel.background = element_rect(fill=NA, color='black'),
          panel.spacing = unit(0.2, "lines"),
          legend.text = element_text(size = 14, face="bold"),
          legend.position="right",
          legend.justification = c(0, 1),
          legend.direction="vertical",
          strip.background = element_rect(colour=NA, fill=NA),
          strip.text = element_blank()) +
    xlab("") +
    ylab("")
  
}           

################################
#### summary of high biallelic rates across subtypes
################################
{
  f2_biallelic_ctds_summary <-
    rbind(f2_alt_rates %>%
            select(Pathway, Hugo_Symbol, gene_label) %>% 
            mutate(variable = 'n_high_biallelic'),
          f2_alt_rates %>%
            select(Pathway, Hugo_Symbol, gene_label) %>% 
            mutate(variable = 'n_not_high_biallelic')) %>%  
    filter(Hugo_Symbol %in% f2_select_genes) %>% 
    left_join(
      f2_biallelic_rates %>%
        filter(n_altered >= 10, f_alteration_rate >= 0.01) %>%
        group_by(Pathway, gene_label) %>%
        summarise(n_high_biallelic = length(which(f_biallelic_rate >= 0.8)),
                  n_not_high_biallelic = length(which(f_biallelic_rate < 0.8)))  %>%
        melt(id.vars=c("Pathway", "gene_label"))) %>%
    mutate(value = ifelse(is.na(value), 0, value)) %>% 
    mutate(variable = factor(variable, levels=c("n_not_high_biallelic", 
                                                "n_high_biallelic")))
  
  f2_biallelic_ctds_summary$gene_label = 
    factor(f2_biallelic_ctds_summary$gene_label, levels=(f2_alt_rates_to_plot$gene_label %>% sort))
  
  
  f2_biallelic_ctds_summary_plot <-
    ggplot(f2_biallelic_ctds_summary, aes(x=gene_label, y=value, fill=variable)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values = c('n_high_biallelic'='#1b8843', 
                                 'n_not_high_biallelic'='#cbd8c0')) + 
    facet_grid(. ~ Pathway, scales = "free", space = "free") +
    theme(axis.text=element_text(size=10, color='black'),
          axis.text.x = element_text(hjust=1, size=10, angle=90, vjust=0.5, face='bold'),
          axis.title.x  = element_text(size=14,face='plain'),
          axis.title.y  = element_text(size=14,face='plain'),
          plot.title = element_text(face = 'bold',hjust = 0.5),
          panel.grid.major = element_line(colour=NA),
          panel.grid.minor = element_line(colour=NA),
          panel.background = element_rect(fill=NA, color='black'),
          panel.spacing = unit(0.2, "lines"),
          legend.text = element_text(size = 14, face="bold"),
          legend.position="right",
          legend.justification = c(0, 1),
          legend.direction="vertical",
          strip.background = element_rect(colour=NA, fill=NA),
          strip.text = element_blank()) +
    xlab("") +
    ylab("")
}  

library(patchwork)
f2_combined_plot <-
  f2_alteration_rate_plot +
  f2_mechanisms_plot +
  f2_heatmap_plot +
  f2_biallelic_ctds_summary_plot +
  plot_layout(ncol=1, heights=c(0.35, 0.35, 8, 0.35))

pdf(paste0('figure_2.pdf'), width=18, height=19)
f2_combined_plot
dev.off()

