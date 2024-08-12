
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

