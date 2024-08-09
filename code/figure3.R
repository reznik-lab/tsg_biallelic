
load('f3_input.Rdata')

f3_select_genes <-
  gxb_drivers_logistf %>%
  mutate(is_sig = logistf_qval < 0.05) %>%
  group_by(Hugo_Symbol) %>%
  summarise(n_sig_all = length(which(is_sig)),
            n_sig_without_pancan = length(which(is_sig & ONCOTREE_CODE != 'all')),
            n_mut_ctd = length(which(ONCOTREE_CODE !='all' & n_mutated >= 20))) %>% 
  filter(n_mut_ctd > 0) %>% 
  arrange(desc(n_sig_without_pancan))


f3_diseases =
  rbind(table_s1 %>% 
          filter(N_ctd >= 50) %>% 
          select(CANCER_TYPE_DETAILED, organ, organ_color, N_ctd, CTD_ABBR, ctd_abbrev_n),
        table_s1 %>% 
          mutate(CANCER_TYPE_DETAILED='all', organ='all', organ_color='black', N_ctd=nrow(table_s1), CTD_ABBR = 'all', ctd_abbrev_n='all (n=23713)') %>%
          select(CANCER_TYPE_DETAILED, organ, organ_color, N_ctd, CTD_ABBR, ctd_abbrev_n)) %>%
  unique

f3_diseases$organ = factor(f3_diseases$organ, levels=c('all', organ_order))
f3_diseases$ctd_abbrev_n = factor(f3_diseases$ctd_abbrev_n, levels=(f3_diseases %>% arrange(desc(N_ctd)))$ctd_abbrev_n)

f3_to_plot <- 
  gxb_drivers_logistf %>%
  left_join(f3_diseases) %>% 
  filter(Hugo_Symbol %in% f3_select_genes$Hugo_Symbol) %>%
  left_join(f2_biallelic_rates %>% select(Hugo_Symbol, ONCOTREE_CODE, n_altered, f_alteration_rate, f_biallelic_rate)) %>%
  mutate(logistf_coef = has_lohTRUE_coeff) %>%
  mutate(logistf_coef = ifelse(logistf_coef > 4, 4, ifelse(logistf_coef < -4, -4, logistf_coef))) %>%
  mutate(significant = logistf_qval < 0.05) %>% 
  mutate(significance_label = ifelse(significant & logistf_coef > 0, "Sig-pos",
                                     ifelse(significant & logistf_coef < 0, "Sig-neg",
                                            ifelse(f_biallelic_rate >= 0.8 & n_altered >= 20, "NotSig-but-High", "N.S")))) %>%
  left_join(tsg_classes %>% select(Hugo_Symbol, tsg_class)) %>% 
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = rev(f3_select_genes$Hugo_Symbol))) 

f3_circle_plot <-
  ggplot(f3_to_plot,aes(x=ctd_abbrev_n,y=Hugo_Symbol)) + 
  geom_point(aes(size = f_alteration_rate, fill = logistf_coef, colour = significance_label),pch=21, stroke=1) +
  facet_grid(tsg_class ~ organ, scales='free', space='free') +
  theme_bw() + 
  scale_fill_gradientn(breaks=c(-2, 0, 2, 4), na.value='gray90', colors=c("#2541a9", "white", "#ff9e83", "#FF0000" ), limits=c(-2,4)) +
  scale_color_manual(values=c("NotSig-but-High" = "#a31e22", "Sig-pos" = "#010001", "Sig-neg" = "#3c68b2", "N.S" = "#a9a8a8")) +
  scale_size_area(breaks=c(0.2, 0.4, 0.6, 0.8)) +
  theme(
    axis.text.x = element_text(hjust=1, vjust=0.5, angle=90),
    axis.text.y = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, size=12),
    panel.grid.major = element_line(color='gray90', size=0.2),
    axis.title.x = element_text(size = 0),
    axis.title.y = element_text(size = 0)) + 
  #theme(plot.margin = unit(c(1,1,9.7,1.5), "lines")) + 
  theme(strip.text.x = element_blank(),
        panel.spacing = unit(0, "lines"))


f3_data_summary <-
  f3_to_plot %>%
  filter(ONCOTREE_CODE != 'all') %>% 
  #filter(n_altered >= 20) %>% 
  filter(n_mutated >= 20) %>% 
  mutate(sig_pos = significance_label == 'Sig-pos',
         sig_neg = significance_label == 'Sig-neg') %>% 
  mutate(high_biallelic = !sig_pos & !sig_neg & f_biallelic_rate >= 0.8) %>%
  mutate(not_sig_gt20 = !sig_pos & !sig_neg & !high_biallelic) %>% 
  group_by(Hugo_Symbol, tsg_class) %>%
  summarise(total = n(),
            n_sig_high_biallelic = length(which(sig_pos | high_biallelic)),
            n_sig_neg = length(which(sig_neg)),
            not_sig_gt20 = length(which(not_sig_gt20))) %>%
  arrange(desc(total))

f3_data_summary_melt <- f3_data_summary %>% melt(id.vars=c('Hugo_Symbol', 'tsg_class', 'total'))
f3_data_summary_melt$variable = factor(f3_data_summary_melt$variable, levels=rev(c("n_sig_high_biallelic","n_sig_neg", "not_sig_gt20")))
f3_data_summary_melt$Hugo_Symbol = factor(f3_data_summary_melt$Hugo_Symbol, levels=rev(f3_select_genes$Hugo_Symbol))

f3_barplot <-
  ggplot(f3_data_summary_melt %>% mutate(total = ifelse(total > 25, 25, total)), 
         aes(x=Hugo_Symbol, y=value, fill=variable)) +
  geom_bar(stat='identity') +
  geom_text(data=f3_data_summary_melt %>% 
              mutate(total_txt = total) %>%
              mutate(total = ifelse(total > 25, 25, total)) %>%
              filter(variable=='n_sig_high_biallelic') ,
            aes(x=Hugo_Symbol, y=(total + 2), label=total_txt), hjust=0, size=2) +
  facet_grid(tsg_class ~ ., space='free', scales='free') +
  scale_fill_manual(values=c("n_sig_high_biallelic" = "#c16f70", "n_sig_neg" = "#5c92cd", "not_sig_gt20" = "#e0e1e0")) +
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



f3_data_summary_melt_prop_chart <-
  f3_data_summary %>%
  mutate(total_gt20 = n_sig_high_biallelic + n_sig_neg + not_sig_gt20) %>%
  mutate(f_sig_high_biallelic = n_sig_high_biallelic/total_gt20,
         f_sig_neg = n_sig_neg/total_gt20,
         f_not_sig_gt20 = not_sig_gt20/total_gt20) %>%
  select(Hugo_Symbol, tsg_class, total, starts_with("f_")) %>%
  melt(id.vars=c("Hugo_Symbol", "tsg_class", "total")) %>%
  mutate(total_label = paste0(Hugo_Symbol, " (n=", total, ")")) %>% 
  mutate(total_label = ifelse(variable=='f_sig_neg', total_label, ''))

f3_data_summary_melt_prop_chart$variable = factor(f3_data_summary_melt_prop_chart$variable, levels=rev(c("f_sig_high_biallelic","f_sig_neg", "f_not_sig_gt20")))
f3_data_summary_melt_prop_chart$Hugo_Symbol = factor(f3_data_summary_melt_prop_chart$Hugo_Symbol, levels=rev(f3_select_genes$Hugo_Symbol))

f3_barplot2 <-
  ggplot(f3_data_summary_melt_prop_chart, aes(x=Hugo_Symbol, y=value, fill=variable)) +
  geom_bar(stat='identity') +
  geom_text(aes(x=Hugo_Symbol, y=1.05, label = total_label), hjust = 0, size=2) +
  facet_grid(tsg_class ~ ., space='free', scales='free') +
  scale_fill_manual(values=c("f_sig_high_biallelic" = "#c16f70", "f_sig_neg" = "#5c92cd", "f_not_sig_gt20" = "#e0e1e0")) +
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


pdf("figure_3.pdf", width=17, height=16)
f3_circle_plot + f3_barplot + f3_barplot2 + 
  plot_layout(widths = c(8,0.6, 0.6), guides = 'collect') + theme(legend.position = 'right')
dev.off()

