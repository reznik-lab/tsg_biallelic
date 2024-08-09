
### Evaluate selection for mut-loh among drivers and VUSs.

####
#### When evaluating signal for selection for mut-loh, we want to exclude samples with VUS
#### mutations. However, there are some instances where VUSs arise in patients that also have  
#### oncogenic mutaitons. Therefore, we need to identify cases for exclusion from "driver" analysis
#### where there are: (1) only VUS mutations or, (2) if there are VUS mutations as well as another 
#### heterozygous mutation. The latter might be the case of two hits where the VUS may be oncogenic. We
#### just exclude these VUS cases to avoid any confounding effects. Similarly, for "vus" analysis,
#### we exclude any tumor that has any oncogenic mutation in that sample.
####
{
  gene_ctd_pairs_to_exclude <-
    gxb_cohort_combined %>% 
    data.table %>% 
    filter(mutation_class != 'silent') %>% 
    filter(grepl("Biallelic|Mut", zygosity_call)) %>% 
    mutate(zyg = ifelse(grepl("Biallelic", zygosity_call), "biallelic", "het")) %>%
    mutate(tag = paste0(mutation_class, "_", zyg)) %>% 
    group_by(Tumor_Sample_Barcode, ONCOTREE_CODE, Hugo_Symbol) %>% 
    summarise(vus_only = !any(mutation_class == 'driver'),
              driver_only = !any(mutation_class == 'vus'),
              has_driver_het = any(tag == 'driver_het'),
              has_any_vus = any(mutation_class=="vus"),
              has_any_driver = any(mutation_class=="driver")) %>%
    ungroup %>% 
    rowwise %>% 
    mutate(ignore_this_gene_ctd_pair_for_driver_analysis = (vus_only | (has_driver_het & has_any_vus))) %>%
    mutate(ignore_this_gene_ctd_pair_for_vus_analysis = has_any_driver) %>%
    filter(ignore_this_gene_ctd_pair_for_driver_analysis | ignore_this_gene_ctd_pair_for_vus_analysis)
}

### function to do regression analysis and gather output
{
  get_logistf_output <- function(gxb_, gene_, otc_, analysis_type = "default"){
    
    if(analysis_type == 'default') {
      mm2 = logistf(has_mut ~ has_loh + CVR_TMB_SCORE + fga + PRIM_MET,
                   data=(gxb_ %>% filter(ONCOTREE_CODE==otc_, Hugo_Symbol==gene_)))
    } else if (analysis_type == 'prim_met') {
      mm2 = logistf(has_mut ~ has_loh + CVR_TMB_SCORE + fga,
                   data=(gxb_ %>% filter(ONCOTREE_CODE==otc_, Hugo_Symbol==gene_)))
    } else {
      stop("")
    }
    
    logistf_coeffs = mm2$coefficients[-1]
    names(logistf_coeffs) = paste0(names(logistf_coeffs), "_coeff")
    
    logistf_pvals = mm2$prob[-1]
    # get exact pvalues
    {
      # Extract the Wald statistics
      wald_stats <- mm2$coefficients / sqrt(diag(mm2$var))
      
      # Compute the two-tailed p-values from the Wald statistics
      exact_p_values <- (2 * pnorm(-abs(wald_stats)))[-1]
    }
    logistf_pvals[which(logistf_pvals == 0)] = exact_p_values[which(logistf_pvals==0)]
    names(logistf_pvals) = paste0(names(logistf_pvals), "_pvalue")
    
    logistf_out = c(logistf_coeffs, logistf_pvals)
    logistf_out = logistf_out[order(names(logistf_out))]
  
    return(logistf_out)
  }
  
  ### just some hack to run the regression model with adply function
  dummy_out_logistf = data.frame(as.list(get_logistf_output(gxb_drivers, "KEAP1", "LUAD")))[] 
  dummy_out_logistf[] <- -1
  
  #get_logistf_output(gxb_drivers, "KEAP1", "LUAD")
  
  summary(logistf(has_mut ~ has_loh + CVR_TMB_SCORE + fga + PRIM_MET,
          data=(gxb_drivers %>% filter(ONCOTREE_CODE=='LUAD', Hugo_Symbol=='STK11'))))
  
  library(doParallel)
  registerDoParallel(cores=6)
}

fisher_func <- function(vec, full=F) {
  f = fisher.test(matrix(vec, byrow=T, nrow=2))
  if(full) {
    return(f)
  } else{
    return(f$p.value)
  }
}

#### Construct data table for selection analysis. 
#### Limit to specific zygosity states (exclude fusions, indeterminate, compound, gains, etc)
{
  gxb_drivers <-
    gxb_cohort_drivers %>% 
    ### Restrict to zygosity-state pairs to these:
    filter(zygosity_call %in% c('WT',
                                'Biallelic - Mut + LOH',
                                'Heterozygous (Mutation)',
                                'Gain-of-mutant - Mut + copy gain',
                                'Heterozygous - LOH (gene-level)',
                                'Heterozygous - focal LOH (gene-level)') ) %>%
    filter(CANCER_TYPE_DETAILED %in% select_diseases) %>%
    ### exclude VUS only cases
    left_join(gene_ctd_pairs_to_exclude %>% 
                filter(ignore_this_gene_ctd_pair_for_driver_analysis) %>% 
                select(Tumor_Sample_Barcode, ONCOTREE_CODE, Hugo_Symbol, to_ignore=ignore_this_gene_ctd_pair_for_driver_analysis)) %>%
    filter(is.na(to_ignore)) %>% 
    mutate(zygosity = ifelse(grepl("Biallelic", zygosity_call), "biallelic",
                             ifelse(zygosity_call %in% c("Heterozygous (Mutation)", 
                                                         "Gain-of-mutant - Mut + copy gain"), "het", 
                                    ifelse(zygosity_call %in% c("Heterozygous - LOH (gene-level)", 
                                                                "Heterozygous - focal LOH (gene-level)"), "loh", "wt")))) %>%
    filter(Hugo_Symbol %in% (oncokb %>% filter(is_tsg))$hugo_symbol) %>%
    mutate(has_loh = grepl("Biallelic|LOH", zygosity_call),
           has_mut = grepl("Biallelic|Mutation", zygosity_call)) %>%
    left_join(table_s1 %>% select(Tumor_Sample_Barcode = SAMPLE_ID, CTD = CANCER_TYPE_DETAILED, PRIM_MET:fga)) 
  
  gxb_drivers_zyg_summary <-
    gxb_drivers %>%
    group_by(ONCOTREE_CODE, CANCER_TYPE_DETAILED, Hugo_Symbol) %>%
    summarise(n=n(),
              n_mutated = length(which(grepl("Biallelic|Heterozygous \\(Mutation|Gain-of-mutant", zygosity_call))),
              n_biallelic = length(which(zygosity == 'biallelic')),
              n_het = length(which(zygosity == 'het')),
              n_loh = length(which(zygosity == 'loh')),
              n_wt = length(which(zygosity == 'wt'))) %>%
    ungroup %>% 
    filter(n_mutated >=3) %>% 
    rowwise %>% 
    mutate(pval = fisher_func(c(n_biallelic, n_het, n_loh, n_wt))) %>%
    arrange(desc(n)) 
  
  gxb_drivers_logistf <-
    plyr::adply(gxb_drivers_zyg_summary %>% filter(n_mutated>2), 1, function(x) {
      logistf_output = get_logistf_output(gxb_drivers, x$Hugo_Symbol, x$ONCOTREE_CODE)
      bind_rows(dummy_out_logistf, data.frame(as.list(logistf_output))) %>% tail(n=1)
    }, .parallel=T)
  
  gxb_drivers_logistf <-
    gxb_drivers_logistf %>%
    mutate(logistf_pval = has_lohTRUE_pvalue, 
           logistf_selection_type = ifelse(has_lohTRUE_coeff > 0, "Positive-Sel", "Negative-Sel"))
  gxb_drivers_logistf$logistf_qval = p.adjust(gxb_drivers_logistf$logistf_pval, method='BH')
  gxb_drivers_logistf$fisher_qval = p.adjust(gxb_drivers_logistf$pval, method='BH')
}

{
  gxb_vus <-
    gxb_cohort_vus %>% 
    ### Restrict to zygosity-state pairs to these:
    filter(zygosity_call %in% c('WT',
                                'Biallelic - Mut + LOH',
                                'Heterozygous (Mutation)',
                                'Gain-of-mutant - Mut + copy gain',
                                'Heterozygous - LOH (gene-level)',
                                'Heterozygous - focal LOH (gene-level)') ) %>%
    filter(CANCER_TYPE_DETAILED %in% select_diseases) %>%
    ### exclude VUS only cases
    left_join(gene_ctd_pairs_to_exclude %>% 
                filter(ignore_this_gene_ctd_pair_for_vus_analysis) %>% 
                select(Tumor_Sample_Barcode, ONCOTREE_CODE, Hugo_Symbol, to_ignore=ignore_this_gene_ctd_pair_for_driver_analysis)) %>%
    filter(is.na(to_ignore)) %>% 
    mutate(zygosity = ifelse(grepl("Biallelic", zygosity_call), "biallelic",
                             ifelse(zygosity_call %in% c("Heterozygous (Mutation)", 
                                                         "Gain-of-mutant - Mut + copy gain"), "het", 
                                    ifelse(zygosity_call %in% c("Heterozygous - LOH (gene-level)", 
                                                                "Heterozygous - focal LOH (gene-level)"), "loh", "wt")))) %>%
    filter(Hugo_Symbol %in% (oncokb %>% filter(is_tsg))$hugo_symbol) %>%
    mutate(has_loh = grepl("Biallelic|LOH", zygosity_call),
           has_mut = grepl("Biallelic|Mutation", zygosity_call)) %>%
    left_join(table_s1 %>% select(Tumor_Sample_Barcode = SAMPLE_ID, CTD = CANCER_TYPE_DETAILED, PRIM_MET:fga)) 
  
  
  gxb_vus_zyg_summary <-
    gxb_vus %>%
    group_by(ONCOTREE_CODE, CANCER_TYPE_DETAILED, Hugo_Symbol) %>%
    summarise(n=n(),
              n_mutated = length(which(grepl("Biallelic|Heterozygous \\(Mutation|Gain-of-mutant", zygosity_call))),
              n_biallelic = length(which(zygosity == 'biallelic')),
              n_het = length(which(zygosity == 'het')),
              n_loh = length(which(zygosity == 'loh')),
              n_wt = length(which(zygosity == 'wt'))) %>%
    ungroup %>% 
    filter(n_mutated >=3) %>% 
    rowwise %>% 
    mutate(pval = fisher_func(c(n_biallelic, n_het, n_loh, n_wt))) %>%
    arrange(desc(n)) 
  
  gxb_vus_logistf <-
    plyr::adply(gxb_vus_zyg_summary %>% 
            filter(n_mutated>2) %>%
            filter(n_biallelic > 0 | n_loh > 0), 1, function(x) {
      logistf_output = get_logistf_output(gxb_vus, x$Hugo_Symbol, x$ONCOTREE_CODE)
      bind_rows(dummy_out_logistf, data.frame(as.list(logistf_output))) %>% tail(n=1)
      
    }, .parallel=T)

  gxb_vus_logistf <-
    gxb_vus_logistf %>%
    mutate(logistf_pval = has_lohTRUE_pvalue, 
           logistf_selection_type = ifelse(has_lohTRUE_coeff > 0, "Positive-Sel", "Negative-Sel"))
  gxb_vus_logistf$logistf_qval = p.adjust(gxb_vus_logistf$logistf_pval, method='BH')
  gxb_vus_logistf$fisher_qval = p.adjust(gxb_vus_logistf$pval, method='BH')

}


### QC: compare between fishers test (Table S7) and regression (for drivers)
{
  table_s7 = fread('Table_S7_mutloh_fishers.txt') 
  
    rbind(table_s7 %>% select(CANCER_TYPE_DETAILED, Hugo_Symbol), 
          gxb_drivers_logistf %>% select(CANCER_TYPE_DETAILED, Hugo_Symbol)) %>%
      unique %>%
      left_join(table_s7) %>%
      left_join(gxb_drivers_logistf) %>%
      mutate(fisher_qval = p_value_corrected,
             fisher_sig = p_value_corrected < 0.05,
             fisher_sig = ifelse(is.na(fisher_qval), "not-evaluated", ifelse(fisher_qval < 0.05, "Yes", "No")),
             fisher_selection_type = ifelse(OR > 1, "Positive-Sel", "Negative-Sel"),
             logistf_sig = ifelse(is.na(logistf_qval), "not-evaluated", ifelse(logistf_qval < 0.05, "Yes", "No"))) %>%
      mutate(fisher_sig = factor(fisher_sig, levels=c("Yes", "No", "not-evaluated")),
             logistf_sig = factor(logistf_sig, levels=c("Yes", "No", "not-evaluated"))) %>%
      mutate(discordant_sig = ifelse((fisher_sig == "Yes" & logistf_sig == "No") | (fisher_sig == 'No' & logistf_sig == 'Yes'),
                                     "Yes", "")) %>% 
      #filter(fisher_sig == 'No', logistf_sig == 'Yes') %>% View
      #filter(fisher_sig == 'Yes', logistf_sig == 'No') %>% View
      #mutate(f_biallelic = n_biallelic/n_mutated) %>% View('sseee')
      select(fisher_sig, logistf_sig) %>% table
      filter(fisher_sig == 'Yes', logistf_sig == 'No') %>%
      View
      #filter(discordant_sig =='Yes') %>% View('drivers_discord')
      #filter(fisher_sig == 'Yes', logistf_sig == 'No') %>%
      #View('drivers_y')
}

## write to supplementary table 7 and 8
{
  gxb_drivers_logistf %>%  
    select(ONCOTREE_CODE, CANCER_TYPE_DETAILED, Hugo_Symbol, N=n,
           n_mut_loh = n_biallelic, n_mut_only = n_het, n_loh_only = n_loh, n_wt,
           fishers_exact_pvalue = pval, 
           tmb_coeff=CVR_TMB_SCORE_coeff, tmb_pvalue = CVR_TMB_SCORE_pvalue,
           fga_coeff:logistf_qval) %>%
    select(-logistf_selection_type) %>% 
    write.excel
}
### Compare drivers vs. vus between fishers and regression
{
  table_s8 <-fread('Table_S8_drivers_vs_vus.txt')

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
    mutate(both_at_least_10_muts = drivers_n_mutated >= 10 & vus_n_mutated >= 10) %>%
    mutate(both_sig = logistf_drivers_qval < 0.05 & logistf_vus_qval < 0.05) %>% 
    #filter(both_at_least_10_muts) %>%
    #filter(logistf_vus_qval < 0.05, logistf_drivers_qval > 0.05)
    left_join(table_s8 %>% select(Hugo_Symbol=Gene, CANCER_TYPE_DETAILED=Disease, table_s8_adjPvalue = vus_adj_p_value, table_s8_OR = vus_OR)) %>%
    mutate(sig_in_old=!is.na(table_s8_adjPvalue) & table_s8_OR > 1) %>% 
    mutate(sig_in_new = (both_at_least_10_muts & logistf_vus_qval < 0.05 & logistf_vus_coeff > 0)) %>%
    #select(sig_in_old, sig_in_new) %>% table
    #filter(!sig_in_old | !sig_in_new) %>% 
    filter(sig_in_old | sig_in_new) %>% 
    View('drivers_vs_vus_logistf')
}
    

