drop_small_factor_levels_col <- function(col, min_n = 3, other_name = "Other", verbose = T) {
  
  if(class(col) != "factor") {col <- as.factor(col)}
  small_levels <- names(table(col)[table(col) < min_n])
  if(verbose==T){print(paste0("merging small levels: ", paste0(small_levels, collapse = ", ")))}
  col <- fct_lump_min(col, min = min_n, other_level = other_name)
  small_levels <- names(table(col)[table(col) < min_n])
  if(verbose==T){print(paste0("dropping small levels: ", paste0(small_levels, collapse = ", ")))}
  col <- fct_lump_min(col, min = min_n, other_level = "to_remove")
  
}

