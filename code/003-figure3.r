library(ggplot2)
library(patchwork)
library(stringr)
library(ggnewscale)
library(RColorBrewer)
library(ggtext)
library(binom)
library(tidyr)

generateFigure3 <- function(enrichmentRes, subtypeTable, clinicalData, mainsize=12, borders=TRUE){-

  ### gene selection criteria for figure 3: if it is significant in at least one cancertype, or if it is
  select_genes_to_plot_fig3 <-
    (enrichmentRes %>% 
       mutate(n_mut = mut + mutloh) %>% 
       mutate(is_sig = p_value_corrected < 0.05) %>%
       group_by(gene) %>% 
       dplyr::summarise(n_sig_all = length(which(is_sig)),
                        n_sig_without_pancan = length(which(is_sig & disease !='all')),
                        has_pancan_pos_sel = length(which(disease=='all' & is_sig & or > 1)),
                        has_pancan_neg_sel = length(which(disease=='all' & is_sig & or < 1)),
                        n_mut_all = sum(n_mut),
                        n_mut_ctd = length(which(disease !='all' & n_mut >= 20))) %>% 
       filter(n_mut_ctd > 0) %>% arrange(desc(n_sig_without_pancan)))$gene
  
  fig3_data <- 
    enrichmentRes %>%
    filter(gene %in% select_genes_to_plot_fig3) %>%
    left_join(classes_by_heuristic_criteria) %>%
    mutate(class = factor(class, levels=c("Class 1", "Class 2", "Class 3", "Class 4"))) %>%
    mutate(disease = ifelse(disease == 'all', "Pancan.", disease)) %>%
    left_join(subtypeTable %>% select(disease = Disease, organ, organ_color, organ_order, CTD_ABBR)) %>%
    mutate(organ = ifelse(disease == 'Pancan.', "Pancan.", organ),
           organ_color = ifelse(disease == 'Pancan.', "black", organ_color),
           organ_order = ifelse(disease == 'Pancan.', 0, organ_order),
           CTD_ABBR = ifelse(disease == 'Pancan.', 'Pancan.', CTD_ABBR)) %>%
    mutate(CANCER_TYPE_DETAILED = disease) %>%
    left_join(clinicalData %>% 
                group_by(ONCOTREE_CODE) %>% 
                summarise(N_otc = n())) %>%
    left_join(f1d_table %>%
                data.table %>%
                janitor::clean_names() %>%
                select(gene, disease, alterations, biallelic_fraction)) %>%
    replace_na(list(biallelic_fraction = 0, alterations = 0)) %>%
    mutate(ctd_abbreviation = ifelse(disease == 'Pancan.', 'Pancan.', paste0(CTD_ABBR, ' (n=', N_otc, ')'))) %>%
    mutate(mutation_rate = (mut + mutloh)/(mut + mutloh + loh + wt)) %>%
    mutate(log_or = ifelse(log_or > 4, 4, log_or)) %>% 
    mutate(significance_label = ifelse(p_value_corrected < 0.05 & log_or > 0, "Sig-pos", 
                                       ifelse(p_value_corrected < 0.05 & log_or < 0, "Sig-neg",
                                              ifelse(biallelic_fraction >= 0.8 & alterations >=20, 'NotSig-but-High', "N.S")))) %>% 
    #ifelse(mutloh_rate >= 0.8 & n_mutations >=20, 'NotSig-but-High', "N.S")))) %>% 
    mutate(significant = p_value_corrected < 0.05) %>%
    mutate(gene = factor(gene, levels = rev(select_genes_to_plot_fig3))) 
  
  fig3_data$ctd_abbreviation = factor(fig3_data$ctd_abbreviation, levels=c("Pancan.", 
                                                                           (fig3_data %>% 
                                                                              filter(ctd_abbreviation != 'Pancan.') %>%
                                                                              select(ctd_abbreviation, N_otc) %>% 
                                                                              unique %>% 
                                                                              arrange(desc(N_otc)))$ctd_abbreviation))
  
  fig3_data$organ = factor(fig3_data$organ, levels=c("Pancan.", (fig3_data %>% 
                                                                   filter(organ != 'Pancan.') %>%
                                                                   select(organ, organ_order) %>% 
                                                                   unique %>% 
                                                                   arrange((organ_order)))$organ))
  
  #####
  fig3_gene_x_ctd_plot <-
    ggplot(fig3_data,aes(x=ctd_abbreviation,y=gene)) + 
    geom_point(aes(size = mutation_rate, fill = log_or, colour = significance_label),pch=21, stroke=1) +
    facet_grid(class ~ organ, scales='free', space='free') +
    theme_bw() + 
    scale_fill_gradientn(breaks=c(-2, 0, 2, 4), na.value='gray90', colors=c("#2541a9", "white", "#ff9e83", "#FF0000" ), limits=c(-2,4)) +
    scale_color_manual(values=c("NotSig-but-High" = "brown", "Sig-pos" = "#010001", "Sig-neg" = "#1767dd", "N.S" = "gray76")) +
    scale_size_area(breaks=c(0.2, 0.4, 0.6, 0.8)) +
    theme(
      axis.text.x = element_text(hjust=1, vjust=0.5, angle=90),
      axis.text.y = element_text(size = 13),
      plot.title = element_text(hjust = 0.5, size=mainsize),
      panel.grid.major = element_line(color='gray90', size=0.2),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 0)) + 
    theme(plot.margin = unit(c(1,1,9.7,1.5), "lines")) + 
    theme(strip.text.x = element_blank(),
          panel.spacing = unit(0, "lines"))
  
  
  fig3_data %>%
    filter(disease != 'Pancan.') %>% 
    left_join(f1d_table %>% 
                filter(Alterations >= 20, 
                       Disease !='all') %>% 
                select(disease = Disease, 
                       gene = Gene, 
                       overall_biallelic_rate = `Biallelic Fraction`,
                       overall_alterations = Alterations)) %>% 
    replace_na(list(overall_biallelic_rate = 0, overall_alterations = 0))  %>% 
    #filter(gene %in% c("SMARCA4", "SMAD3", "RASA1", "KMT2C")) %>% 
    filter(gene %in% c("PBRM1")) %>% 
    mutate(fig3_mutations_N = mut + mutloh,
           fig2_alterations_N = overall_alterations,
           fig3_biallelic_rate = (mutloh/(mut+mutloh))) %>% 
    select(gene, disease, fig2_alterations_N, fig2_biallelic_rate = overall_biallelic_rate, fig3_mutations_N, fig3_biallelic_rate, p_value_corrected) %>% 
    arrange(desc(fig2_alterations_N), 
            desc(fig3_mutations_N)) %>% 
    data.frame 
  filter(fig2_alterations_N >= 20 | fig3_mutations_N >=20) %>% 
    arrange(gene)
  
  
  
  fig3_data_summary <-
    fig3_data %>%
    filter(disease != 'Pancan.') %>% 
    left_join(f1d_table %>% 
                filter(Alterations >= 20, 
                       `Biallelic Fraction` >= 0.8, 
                       Disease !='all') %>% 
                select(disease = Disease, 
                       gene = Gene, 
                       overall_biallelic_rate = `Biallelic Fraction`,
                       overall_alterations = Alterations)) %>% 
    replace_na(list(overall_biallelic_rate = 0, overall_alterations = 0)) %>% 
    filter(mut+mutloh >= 20 | overall_alterations >= 20) %>% 
    mutate(sig_pos = p_value_corrected < 0.05 & or > 1) %>%
    mutate(sig_neg = p_value_corrected < 0.05 & or < 1) %>%
    mutate(high_biallelic = !sig_pos & !sig_neg & overall_biallelic_rate >= 0.8) %>%
    mutate(not_sig_gt20 = !sig_pos & !sig_neg & !high_biallelic & (overall_alterations >=20 | mut + mutloh >= 20)) %>% 
    group_by(gene, class) %>%
    summarise(total = n(),
              n_sig_high_biallelic = length(which(sig_pos | high_biallelic)),
              #n_sig_pos = length(which(sig_pos)),
              n_sig_neg = length(which(sig_neg)),
              #high_biallelic = length(which(high_biallelic)),
              not_sig_gt20 = length(which(not_sig_gt20))) %>%
    mutate(n_sig_high_biallelic = ifelse(n_sig_high_biallelic > 25, 25, n_sig_high_biallelic))
  
  
  fig3_data_summary_melt <-
    fig3_data_summary %>%
    melt(id.vars=c('gene', 'class', 'total'))
  
  fig3_data_summary_melt$variable = factor(fig3_data_summary_melt$variable, levels=rev(c("n_sig_high_biallelic","n_sig_neg", "not_sig_gt20")))
  
  fig3_data_summary_melt$gene = factor(fig3_data_summary_melt$gene, levels=rev(select_genes_to_plot_fig3))
  
  fig3_barplot <-
    ggplot(fig3_data_summary_melt %>% mutate(total = ifelse(total > 25, 25, total)), aes(x=gene, y=value, fill=variable)) +
    geom_bar(stat='identity') +
    geom_text(aes(x=gene, y=(total + 2), label=total), hjust=0, size=2) +
    facet_grid(class ~ ., space='free', scales='free') +
    scale_fill_manual(values=c("n_sig_high_biallelic" = "#0a9396", "n_sig_neg" = "#e63946", "not_sig_gt20" = "#d1d1d1")) +
    theme(axis.text=element_text(size=8, color='black'),
          axis.text.x = element_text(hjust=1, angle=90, vjust=0.5),
          plot.title = element_text(face = 'bold',hjust = 0.5),
          panel.grid.major  = element_line(colour=NA),
          panel.grid.minor = element_line(colour=NA),
          #axis.line = element_line(),
          plot.margin = margin(0, 0, 0, 0, "pt"),
          panel.spacing = unit(0.1, "lines"),
          panel.background = element_rect(fill='white', color='black'),
          legend.text = element_text(size = 14, face="bold"),
          legend.position="none",
          legend.justification = c(0, 1),
          legend.direction="vertical",
          strip.background = element_blank(),
          strip.text = element_blank()
    ) +
    xlab("") +
    ylab("")  + coord_flip()
  
  
  fig3_data_summary_melt_prop_chart <-
    fig3_data_summary %>%
    mutate(total_gt20 = n_sig_high_biallelic + n_sig_neg + not_sig_gt20) %>%
    mutate(f_sig_high_biallelic = n_sig_high_biallelic/total_gt20,
           f_sig_neg = n_sig_neg/total_gt20,
           f_not_sig_gt20 = not_sig_gt20/total_gt20) %>%
    select(gene, class, total, starts_with("f_")) %>%
    melt(id.vars=c("gene", "class", "total")) %>%
    mutate(total_label = paste0(gene, " (n=", total, ")")) %>% 
    mutate(total_label = ifelse(variable=='f_sig_neg', total_label, ''))
  
  fig3_data_summary_melt_prop_chart$variable = factor(fig3_data_summary_melt_prop_chart$variable, levels=rev(c("f_sig_high_biallelic","f_sig_neg", "f_not_sig_gt20")))
  fig3_data_summary_melt_prop_chart$gene = factor(fig3_data_summary_melt_prop_chart$gene, levels=rev(select_genes_to_plot_fig3))
  
  fig3_barplot2 <-
    ggplot(fig3_data_summary_melt_prop_chart, aes(x=gene, y=value, fill=variable)) +
    geom_bar(stat='identity') +
    geom_text(aes(x=gene, y=1.05, label = total_label), hjust = 0, size=2) +
    facet_grid(class ~ ., space='free', scales='free') +
    scale_fill_manual(values=c("f_sig_high_biallelic" = "#0a9396", "f_sig_neg" = "#e63946", "f_not_sig_gt20" = "#d1d1d1")) +
    theme(axis.text=element_text(size=8, color='black'),
          axis.text.x = element_text(hjust=1, angle=90, vjust=0.5),
          plot.title = element_text(face = 'bold',hjust = 0.5),
          panel.grid.major  = element_line(colour=NA),
          panel.grid.minor = element_line(colour=NA),
          #axis.line = element_line(),
          plot.margin = margin(0, 0, 0, 0, "pt"),
          panel.spacing = unit(0.1, "lines"),
          panel.background = element_rect(fill='white', color='black'),
          legend.text = element_text(size = 14, face="bold"),
          legend.position="right",
          legend.justification = c(0, 1),
          legend.direction="vertical",
          strip.background = element_blank(),
          strip.text = element_blank()
    ) +
    xlab("") +
    ylab("")  + coord_flip()
  
  pdf("./fig3.pdf", width=15, height=16)
  fig3_gene_x_ctd_plot + fig3_barplot + fig3_barplot2 + 
    plot_layout(widths = c(8,0.6, 0.6), guides = 'collect') + theme(legend.position = 'right')
  dev.off()
  
}
#generateFigure3(enrichmentRes, subtypeTable, clinicalData=clinData.filtered, mainsize=15)



