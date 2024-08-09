
# Figure 4A
### Genes to include in Figure 4A: Class I TSGs that are significant in at least 2 cancer types in which they are mutated in >10% of tumors
### Those would be: APC, BAP1, CDH1, NF1, RB1, SMAD4, STK11 (pan-cancer TSGs that were found significant in 15 or more cancer types were excluded from this plot.
### TP53, CDKN2A, PTEN and RB1 were excluded)
f4a_select_genes <-
  gxb_drivers_logistf %>%
  filter(logistf_qval < 0.05) %>% 
  mutate(f_mutated = n_mutated/n) %>% 
  filter(ONCOTREE_CODE != 'all') %>%
  left_join(tsg_classes %>% select(Hugo_Symbol, tsg_class) %>% unique) %>% 
  #exclude unclassified TSGs & restrict to class 1 tsgs
  filter(!is.na(tsg_class), tsg_class == "Class 1") %>% 
  #exclude pancancer TSGs
  filter(!(Hugo_Symbol %in% c("TP53", "CDKN2A", "PTEN", "RB1"))) %>%
  group_by(Hugo_Symbol, tsg_class) %>% 
  summarise(total = length(which(f_mutated >= 0.1)),
            n_lt_10 = length(which(f_mutated < 0.1))) %>%
  filter(total > 1)


f4a_data_to_plot <-
  gxb_drivers_logistf %>%
  mutate(f_mutated = n_mutated/n) %>% 
  filter(ONCOTREE_CODE != 'all') %>%
  filter(Hugo_Symbol %in% f4a_select_genes$Hugo_Symbol) %>% 
  mutate(neg_log10_qval = -log10(logistf_qval)) %>% 
  mutate(logistf_coef = has_lohTRUE_coeff) %>%
  mutate(logistf_coef = ifelse(logistf_coef > 6, 6, ifelse(logistf_coef < -6, -6, logistf_coef)))  %>%
  mutate(is_sig = logistf_qval < 0.05) %>% 
  mutate(to_show = neg_log10_qval > 5 | (f_mutated > 0.1 & logistf_qval < 0.05) )

f4a_gene_colors <-c("")

f4a_plot <-
  ggplot(f4a_data_to_plot %>% filter(!is_sig),
       aes(x = f_mutated, y = neg_log10_qval, fill=Hugo_Symbol, color=Hugo_Symbol)) +
  geom_hline(yintercept = 1.31, linewidth=0.15, linetype='dashed') +
  geom_vline(xintercept = 0.1, linewidth=0.15, linetype='dashed') +
  geom_point(pch=19, alpha=0.4, size=3) +
  geom_point(data=f4a_data_to_plot %>% filter(is_sig), pch=21, alpha=0.75,  color='black', size=4, stroke=0.4) +
  geom_text_repel(data=f4a_data_to_plot %>% filter(is_sig & (f_mutated > 0.1 | neg_log10_qval >= 9)),
                  aes(x=f_mutated, y=neg_log10_qval, label=ONCOTREE_CODE),
                  size=6, show.legend = F, box.padding = unit(0.45, "lines"),
                  max.overlaps = 17) +
  scale_fill_manual(values=c(ggthemes::tableau_color_pal(palette= 'Classic 10')(6))) +
  scale_color_manual(values=c(ggthemes::tableau_color_pal(palette= 'Classic 10')(6))) +
  scale_y_continuous(breaks=c(0:20)* 10) +
  scale_x_continuous(breaks=c(0:10)* 0.1, labels=percent) +
  theme(axis.text=element_text(size=14, color='black'),
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

pdf('figure_4a.pdf',width=10,height=7)
print(f4a_plot)
dev.off() 


f4a_inset_data_to_plot <- f4a_data_to_plot %>% filter(f_mutated < 0.1 & neg_log10_qval >= 1.5 & neg_log10_qval < 9) %>% mutate(to_display=T)

f4a_inset <-
  ggplot(f4a_inset_data_to_plot,
         aes(x = f_mutated, y = neg_log10_qval, fill=Hugo_Symbol, color=Hugo_Symbol)) +
  geom_point(pch=21, alpha=0.75, size=4) +
  geom_text_repel(data=f4a_inset_data_to_plot, 
                  aes(x = f_mutated, y = neg_log10_qval, label=ONCOTREE_CODE), 
                  size=6, show.legend = F, box.padding = unit(0.5, "lines"),
                  max.overlaps = 17) +
  scale_fill_manual(values=c(ggthemes::tableau_color_pal(palette= 'Classic 10')(6))) +
  scale_color_manual(values=c(ggthemes::tableau_color_pal(palette= 'Classic 10')(6))) +
  scale_y_continuous(breaks=c(0:10), limits=c(1,5)) +
  scale_x_continuous(breaks=c(0:10)* 0.01, labels=percent) +
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
