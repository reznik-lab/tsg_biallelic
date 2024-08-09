 
# figure 5A

f5a_data_to_plot <-
  rbind(gxb_drivers_logistf %>% select(Hugo_Symbol, ONCOTREE_CODE, CANCER_TYPE_DETAILED),
      gxb_vus_logistf %>% select(Hugo_Symbol, ONCOTREE_CODE, CANCER_TYPE_DETAILED)) %>%
  unique %>%
  filter(ONCOTREE_CODE !='all') %>% 
  left_join(gxb_drivers_logistf %>% select(Hugo_Symbol, ONCOTREE_CODE, drivers_n_mutated = n_mutated, 
                                           drivers_n_mutloh = n_biallelic, drivers_n_mut_only = n_het, drivers_n_loh = n_loh, drivers_n_wt = n_wt,
                                           logistf_drivers_qval = logistf_qval, logistf_drivers_coeff = has_lohTRUE_coeff)) %>%
  left_join(gxb_vus_logistf %>% select(Hugo_Symbol, ONCOTREE_CODE, vus_n_mutated = n_mutated,
                                       vus_n_mutloh = n_biallelic, vus_n_mut_only = n_het, vus_n_loh = n_loh, vus_n_wt = n_wt,
                                       logistf_vus_qval = logistf_qval, logistf_vus_coeff = has_lohTRUE_coeff)) %>% 
  filter(drivers_n_mutated >=10 & vus_n_mutated >=10) %>% 
  filter(logistf_drivers_coeff > 0, logistf_vus_coeff > 0) %>% 
  mutate(both_sig = logistf_drivers_qval < 0.05 & logistf_vus_qval < 0.05) %>%
  mutate(neg_log10_drivers_qval = -log10(logistf_drivers_qval),
         neg_log10_vus_qval = -log10(logistf_vus_qval)) %>%
  mutate(neg_log10_drivers_qval = ifelse(neg_log10_drivers_qval > 20, 20, neg_log10_drivers_qval)) %>%
  mutate(neg_log10_vus_qval = ifelse(neg_log10_vus_qval > 10, 10, neg_log10_vus_qval))


f5a_select_genes = (f5a_data_to_plot %>% filter(both_sig))$Hugo_Symbol %>% unique

f5a_plot <-
  ggplot(f5a_data_to_plot %>% 
         filter(Hugo_Symbol %in% f5a_select_genes), aes(x=neg_log10_vus_qval, y=neg_log10_drivers_qval, fill=Hugo_Symbol, color=Hugo_Symbol)) +
  geom_hline(yintercept = 1.31, linewidth=1, color='darkred', linetype='dotted') +
  geom_vline(xintercept = 1.31, linewidth=1, color='darkred', linetype='dotted') +
  geom_point(pch=21, size=5, color='black') +
  scale_y_continuous(breaks=c(0:10)*5) +
  scale_x_continuous(breaks=c(0:10)*2) +
  scale_fill_manual(values=c(ggthemes::tableau_color_pal(palette= 'Tableau 20')(20) [c(1,3,5,7,9,11,13,15,17,19,8,12,2)])) +
  scale_color_manual(values=c(ggthemes::tableau_color_pal(palette= 'Tableau 20')(20) [c(1,3,5,7,9,11,13,15,17,19,8,12,2)])) +
  geom_text_repel(data=f5a_data_to_plot %>% 
                    filter(Hugo_Symbol %in% f5a_select_genes), 
                  aes(x=neg_log10_vus_qval, y=neg_log10_drivers_qval, label=ONCOTREE_CODE), 
                  size=6, show.legend = F, box.padding = unit(1, "lines"),
                  max.overlaps = 17) +
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


pdf('figure_5a.pdf',width=10,height=7)
print(f5a_plot)
dev.off() 

