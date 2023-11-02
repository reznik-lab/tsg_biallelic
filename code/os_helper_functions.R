pretty_HR_text <- function(model, term) {
  return(paste0("HR=",format((summary(model)$coefficients[term,"exp(coef)"]), digits=2), 
               " [95% CI: ", round((summary(model)$conf.int[term,"lower .95"]), 2), "-", round((summary(model)$conf.int[term,"upper .95"]), 2), 
               "], p-value=", format((summary(model)$coefficients[term,"Pr(>|z|)"]), digits = 2, nsmall = 2)))
}

pretty_survmedian_text <- function(surv, strata1, strata2, unit = "mo.") {
  return(paste(
    format(surv_median(surv)[surv_median(surv)$strata == strata1, "median"], digits = 2, nsmall = 1),
    unit,
    "vs.", 
    format(surv_median(surv)[surv_median(surv)$strata == strata2, "median"], digits = 2, nsmall = 1),
    unit, collapse = " "))
}