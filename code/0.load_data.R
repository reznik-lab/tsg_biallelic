### generate biallelicEnrichmentRes table for drivers.
load("gene_zygosity_focal_10mb_lessthan_10_all_genes_June_13_2024.RData")

oncokb <-fread('oncokb.txt')
putative_germline <- read.table('putative_germline_genes_to_exclude_from_consdieration',sep='\t',header = T)

load('sysdata.rda')
non_autosomal_genes = (genes_hg19 %>% filter(grepl("X|Y|M", chrom)))$gene

table_s1 <- fread('Table_S1_clinDataFiltered.txt')

{
  gxb_cohort_combined <-
    gene_x_barcode_combined %>%
    mutate(zygosity_call = ifelse(is.na(zygosity_call), 'WT', as.character(zygosity_call))) %>% 
    filter(Tumor_Sample_Barcode %in% table_s1$SAMPLE_ID) %>%
    left_join(table_s1 %>%
                select(Tumor_Sample_Barcode = SAMPLE_ID, CANCER_TYPE, CANCER_TYPE_DETAILED, ONCOTREE_CODE)) %>%
    ### exclude specific germline  genes
    mutate(patient = gsub("\\-T.*", "", Tumor_Sample_Barcode)) %>%
    left_join(putative_germline %>% mutate(exclude_for_germline=T)) %>%
    filter(is.na(exclude_for_germline)) %>%
    ### exclude non-autosomal genes
    filter(!(Hugo_Symbol %in% non_autosomal_genes))
  
  gxb_cohort_combined <-
    rbind(gxb_cohort_combined,
          gxb_cohort_combined %>% mutate(CANCER_TYPE='all', CANCER_TYPE_DETAILED='all', ONCOTREE_CODE='all'))
  
  gxb_cohort_drivers <-
    gene_x_barcode_drivers %>%
    mutate(zygosity_call = ifelse(is.na(zygosity_call), 'WT', as.character(zygosity_call))) %>% 
    filter(Tumor_Sample_Barcode %in% table_s1$SAMPLE_ID) %>%
    left_join(table_s1 %>%
                select(Tumor_Sample_Barcode = SAMPLE_ID, CANCER_TYPE, CANCER_TYPE_DETAILED, ONCOTREE_CODE)) %>%
    ### exclude specific germline  genes
    mutate(patient = gsub("\\-T.*", "", Tumor_Sample_Barcode)) %>%
    left_join(putative_germline %>% mutate(exclude_for_germline=T)) %>%
    filter(is.na(exclude_for_germline))  %>%
    ### exclude non-autosomal genes
    filter(!(Hugo_Symbol %in% non_autosomal_genes))
  
  gxb_cohort_drivers <-
    rbind(gxb_cohort_drivers,
          gxb_cohort_drivers %>% mutate(CANCER_TYPE='all', CANCER_TYPE_DETAILED='all', ONCOTREE_CODE='all'))
    
  gxb_cohort_vus <-
    gene_x_barcode_vus %>%
    mutate(zygosity_call = ifelse(is.na(zygosity_call), 'WT', as.character(zygosity_call))) %>% 
    filter(Tumor_Sample_Barcode %in% table_s1$SAMPLE_ID) %>%
    left_join(table_s1 %>%
                select(Tumor_Sample_Barcode = SAMPLE_ID, CANCER_TYPE, CANCER_TYPE_DETAILED, ONCOTREE_CODE)) %>%
    ### exclude specific germline  genes
    mutate(patient = gsub("\\-T.*", "", Tumor_Sample_Barcode)) %>%
    left_join(putative_germline %>% mutate(exclude_for_germline=T)) %>%
    filter(is.na(exclude_for_germline))  %>%
    ### exclude non-autosomal genes
    filter(!(Hugo_Symbol %in% non_autosomal_genes))
  
  gxb_cohort_vus <-
    rbind(gxb_cohort_vus,
          gxb_cohort_vus %>% mutate(CANCER_TYPE='all', CANCER_TYPE_DETAILED='all', ONCOTREE_CODE='all'))
  
  select_diseases = c("all", (table_s1 %>% group_by(CANCER_TYPE_DETAILED) %>% summarise(total =n()) %>% filter(total >= 50))$CANCER_TYPE_DETAILED)
}




