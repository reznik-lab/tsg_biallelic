

fig1b_colorevec <- c('Biallelic - Mut + LOH' = "#58B9F4",
                     'Biallelic - Mut + fusion'="#9E65AA",
                     'Biallelic - Homdel'="#242860", 
                     'Biallelic - compound'= "#8CA1D3",
                     'Het - Fusion'='plum1',
                     'Heterozygous'="gray80",
                     'Het - Focal'='#749367',
                     'Amplification'="#D25047")


f1b_gxb_drivers_zygosity_stats <-
  gxb_cohort_drivers %>% 
  filter(Tumor_Sample_Barcode %in% table_s1$SAMPLE_ID) %>%
  mutate(gene_type = ifelse(Hugo_Symbol %in% tsgs, "TSG",
                            ifelse((Hugo_Symbol %in% oncogenes), "Oncogene", "Other"))) %>% 
  filter(gene_type != 'Other') %>% 
  filter(CANCER_TYPE_DETAILED != 'all') %>% 
  
  # exclude fusions or mutations for which zygosity is not determinable
  filter(!grepl("Indeterminate", zygosity_call)) %>% 
  
  # exclude instances where a tumor with only a heterozygous driver mutation also has a VUS mutation
  left_join(gene_ctd_pairs_to_exclude %>%
              filter(ignore_this_gene_ctd_pair_for_driver_analysis) %>%
              select(Tumor_Sample_Barcode, ONCOTREE_CODE, Hugo_Symbol, to_ignore=ignore_this_gene_ctd_pair_for_driver_analysis)) %>%
  filter(is.na(to_ignore)) %>%
  data.table  %>%  
  mutate(zygosity = ifelse(grepl("Biallelic", zygosity_call), zygosity_call,
                           ifelse(zygosity_call %in% c("Heterozygous (Mutation)", "Gain-of-mutant - Mut + copy gain"), "Heterozygous", 
                                  ifelse(zygosity_call == 'Heterozygous (Fusion)', 'Het - Fusion', 
                                         ifelse(zygosity_call == 'Amplification - CNA', 'Amplification',
                                                ifelse(zygosity_call == "Heterozygous - focal LOH (gene-level)", "Het - Focal", "wt")))))) %>%
  filter(zygosity != 'wt') %>% 
  group_by(gene_type, zygosity) %>% 
  summarise(total = n()) %>%
  mutate(zygosity = factor(zygosity, levels=rev(fig1b_colorevec %>% names))) %>% 
  data.frame

f1b_gxb_drivers_zygosity_summary <- 
  f1b_gxb_drivers_zygosity_stats %>%
  left_join(
    f1b_gxb_drivers_zygosity_stats %>% 
      filter(zygosity != 'Het - Focal') %>% 
      group_by(gene_type) %>%
      summarise(N_altered_without_focal = sum(total)) %>%
      left_join(f1b_gxb_drivers_zygosity_stats %>% 
                  group_by(gene_type) %>%
                  summarise(N_altered_with_focal = sum(total)))) %>%
  mutate(f_without_focal = round(total / N_altered_without_focal, 6)) %>% 
  mutate(f_with_focal = round(total / N_altered_with_focal, 6))
  


f1b_plot <-
  ggplot(f1b_gxb_drivers_zygosity_summary %>%
           filter(zygosity != 'Het - Focal'), aes(x=gene_type, y=f_without_focal, fill=zygosity)) + 
  geom_bar(stat='identity', color='black', size=0.1) + 
  scale_fill_manual(values=fig1b_colorevec) +
  scale_y_continuous(label=percent) +
  theme(axis.text=element_text(size=10, color='black'),
        axis.text.x = element_blank(),
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

f1b_supplementary_plot <-
  ggplot(f1b_gxb_drivers_zygosity_summary, aes(x=gene_type, y=f_with_focal, fill=zygosity)) + 
  geom_bar(stat='identity', color='black', size=0.1) + 
  scale_fill_manual(values=fig1b_colorevec) +
  scale_y_continuous(label=percent) +
  theme(axis.text=element_text(size=10, color='black'),
        axis.text.x = element_blank(),
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
  ylab("") + 
  ggtitle("with focal deletions")

pdf("figure_1b_supplement_with_focal.pdf", width=5, height=6)
f1b_supplementary_plot
dev.off()
