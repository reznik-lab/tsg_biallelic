#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(binom)
  library(stringr)
  library(doParallel)
})


t_vaf <- function(ref_copy_num, alt_copy_num, purity) { # adjusted for somatic variant
  ((purity * alt_copy_num)) / (( (purity * ref_copy_num) + (1-purity) ) + ( (purity * alt_copy_num) + (1-purity)))
}

# calculate zygosity for each variant in the maf
annotate_mutation_with_zygosity <- function(mut) {
  
  tcn = mut$tcn
  lcn = ifelse(tcn == 1, 0, mut$lcn)
  purity = mut$purity
  t_depth = mut$t_depth
  t_alt_count = mut$t_alt_count
  t_alt_freq = ifelse('t_alt_freq' %in% names(mut), mut$t_alt_freq, mut$t_alt_count/mut$t_depth)
  clonality = mut$clonality
  ccf_expected_copies = mut$ccf_expected_copies
  expected_alt_copies = mut$expected_alt_copies
  
  if (is.na(purity) | is.infinite(purity)) { return(c(zygosity_call = 'Indeterminate',
                                                      zygosity_call_reason = 'purity is NA')) }

  if (is.na(t_alt_freq)) { return(c(zygosity_call = 'Indeterminate',
                                    zygosity_call_reason = 't_alt_freq is NA')) }
  
  if (is.na(tcn)) { return(c(zygosity_call = 'Indeterminate',
                             zygosity_call_reason = 'tcn is NA')) }
  
  if (tcn == 0) {
    return(c(zygosity_call = 'Indeterminate',
             zygosity_call_reason = 'tcn is 0'))
  }
  
  #### Step 1: tcn = 1 (lcn is either 0 or NA). Return "biallelic" unless the mutation or VAF is outrageously subclonal (set it to INDETERMINATE, then).
  if (tcn == 1) {
    ### Future consideration: If the mutaiton is outrageously subclonal, then set the zygosity_call as 'INDETERMINATE'
    return(c(zygosity_call = 'Biallelic - Mut + LOH',
             zygosity_call_reason = 'tcn is 1'))
  }
  
  if (is.na(expected_alt_copies)) { 
    return(c(zygosity_call = 'Indeterminate', 
             zygosity_call_reason = 'expected_alt_copies is NA'))
  }
  
  #### Step 2: lcn = 0, tcn >= 1. "biallelic" if mutation is expected to be on all retaiined copies. "Heterozygous/gain-of-mutant-allele" otherwise. 
  # if (!is.na(lcn) && lcn == 0 && expected_alt_copies >= tcn) {
  #   return(c(zygosity_call = 'Biallelic - Mut + LOH',
  #            zygosity_call_reason = 'lcn = 0, mutation on all mutated copies'))
  # }
  
  if (!is.na(lcn) && lcn == 0) {
    # if expected_alt_copies >1 then, the mutation is on multiple copies. It is still possible that there are wild-type alleles (expected_alt_copies < tcn),  but those instances
    # are very rare and hard to distinguish from events were segment clonality is difficult to establish
    if (expected_alt_copies >= tcn | expected_alt_copies > 1) {
      return(c(zygosity_call = 'Biallelic - Mut + LOH',
               zygosity_call_reason = 'lcn = 0, mutation on multiple alleles'))
    }
    
    ### in instances where tcn = 2 and lcn = 0, it is possible that the mutation arose late. to make sure we do not call such events as biallelic,
    ### make sure that the observed tumor vaf is higher than the upper bound of expected VAF if there is actually a 1-copy each of ref and alt alleles in the tumor.
    expected_vaf_loh = t_vaf(1, tcn-1, purity)
    expected_vaf_loh_upper = (binom.confint(min(t_depth, expected_vaf_loh * t_depth), t_depth, methods="wilson"))$upper
    if (t_alt_freq >= expected_vaf_loh_upper) {
      return(c(zygosity_call = 'Biallelic - Mut + LOH',
               zygosity_call_reason = 't_alt_freq is consistent with multiple copies of the mutation'))
    }
  }
  
  if (is.na(lcn) && expected_alt_copies >= tcn && !is.na(ccf_expected_copies) && ccf_expected_copies > 0.75) {
    return(c(zygosity_call = 'Biallelic - Mut + LOH',
             zygosity_call_reason = 'lcn is NA, but mutation predicted to be on all mutated copies and is clonal'))
  }
  
  #### Step 3: use allele-frequency to check if the mutant allele frequency is concordant with an LCN of 0
  {
    ### determine expected VAF for num_alt_copies = tcn and num_ref_copies = 0.
    expected_vaf_loh = t_vaf(0, tcn, purity)
    expected_vaf_loh_lower = (binom.confint(min(t_depth, expected_vaf_loh * t_depth), t_depth, methods="wilson"))$lower
    
    if (t_alt_freq >= expected_vaf_loh_lower) {
      return(c(zygosity_call = 'Biallelic - Mut + LOH',
               zygosity_call_reason = 't_alt_freq consistent with complete LOH (lcn = 0) and mutation on all retained copies'))
    }
  }

  #### Step 4: lcn > 0. (includes tcn:lcn of 2:1, 3:1) distinguish between "Heterozygous/gain-of-mutant-allele" 
  if ( expected_alt_copies > tcn - expected_alt_copies) {
    return(c(zygosity_call = 'Gain-of-mutant - Mut + copy gain',
             zygosity_call_reason = 'wild-type present, n_mut_copies > n_ref_copies'))
  }

  #### Step 5: use allele-frequency to check if the mutant allele freuency is concordant with mutaion present on "mcn"  
  {
    ### determine expected VAF for num_alt_copies = tcn and num_ref_copies = 0.
    exp_mcn = ceiling(tcn/2)
    expected_vaf_loh = t_vaf(tcn - exp_mcn, exp_mcn, purity)
    expected_vaf_loh_lower = (binom.confint(min(t_depth, expected_vaf_loh * t_depth), t_depth, methods="wilson"))$lower
    if ( exp_mcn > (tcn - exp_mcn) && t_alt_freq >= expected_vaf_loh_lower) {
      return(c(zygosity_call = 'Gain-of-mutant - Mut + copy gain',
               zygosity_call_reason = 't_alt_freq consistent with mutation on mcn, with mcn > lcn'))
    }
  }
  return(c(zygosity_call = 'Heterozygous (Mutation)',
           zygosity_call_reason = ''))
}


annotate_maf_with_zygosity <- function(maf, ncores=1) {
  registerDoParallel(cores=ncores)
  maf_out = adply(maf, 1, annotate_mutation_with_zygosity, .parallel = T)
  maf_out$zygosity_call = factor(maf_out$zygosity_call, levels=c("Biallelic - Mut + LOH", "Gain-of-mutant - Mut + copy gain", "Heterozygous (Mutation)", "Heterozygous (Fusion)", "WT", "Indeterminate"))
  return(maf_out)
}
