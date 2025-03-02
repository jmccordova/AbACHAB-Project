# References for imputation: https://bookdown.org/mike/data_analysis/imputation-missing-data.html

draw_missmap <- function(df, label) {
  missmap(
    df,
    main = paste("Missing Data Heatmap for", label),
    col = c("yellow", "black"),
    legend = TRUE
  )
}

draw_missmap(antibodies, "Antibodies")
draw_missmap(perc_inhibition_orig_variant, "Inhibition of Wuhan variant")
draw_missmap(perc_inhibition_delta_variant, "Inhibition of Delta variant")

# Test for MCAR missingness
# Ho: The data are missing completely at random (MCAR).
# Ha: The data are not missing completely at random (MCAR)
naniar::mcar_test(antibodies)   # p-val = 8.4e-12; reject Ho
misty::na.test(antibodies[, c('Day 0', 'Day 14', 'Day 180')])   # p-val = 0; reject Ho
# Data is either MAR or MNAR; missingness follows a pattern
# Need more advanced imputation technique

naniar::mcar_test(perc_inhibition_orig_variant)   # p-val = 0.796; accept Ho
misty::na.test(perc_inhibition_orig_variant)   # p-val = 0.804; accept Ho
# Data is MCAR; missingness is at random
# Simple imputation technique can be used

naniar::mcar_test(perc_inhibition_delta_variant)   # p-val = 0.979; accept Ho
misty::na.test(perc_inhibition_delta_variant[, c('Day 14', 'Day 30', 'Day 90', 'Day 180')])   # p-val = 0.985; accept Ho
# Data is MCAR; missingness is at random
# Simple imputation technique can be used

# Impute values
# Using missForest
antibodies.imputed.mf <- missForest(antibodies)$ximp
perc_inhibition_orig_variant.imputed.mf <- missForest(perc_inhibition_orig_variant)$ximp
perc_inhibition_delta_variant.imputed.mf <- missForest(perc_inhibition_delta_variant)$ximp

#IgG Antibodies
draw_histograms(antibodies.imputed.mf, "IgG antibodies", "IgG")
# %Inhibition original variant
draw_histograms(perc_inhibition_orig_variant.imputed.mf, "%Inhibition Original Variant", "% inhibition")
# %Inhibition delta variant
draw_histograms(perc_inhibition_delta_variant.imputed.mf, "%Inhibition Delta Variant", "% inhibition")

# Using MICE with sample method
compare_imputations <- function(df, type, method, ret = FALSE) {
    imputed_data <- mice(df, method = method, seed = 690)
    summary(imputed_data)
    # Perform linear regression on each imputed dataset
    print("----------------------")
    print("Fitted:")
    if (type != 'Delta') {
      fit <- with(imputed_data, lm(`Day 0` ~ `Day 14` + `Day 30` + `Day 90` + `Day 180`))
    } else {
      fit <- with(imputed_data, lm(`Day 14` ~ `Day 30` + `Day 90` + `Day 180`))
    }
    # View the summary of the fit
    summary(fit)
    
    # Pool the results from the multiple imputations
    print("----------------------")
    print("Pooled Results:")
    pooled_results <- pool(fit)
    # View the summary of the pooled results
    summary(pooled_results)
    
    # Visualize imputed data
 #   stripplot(imputed_data, pch = 20, cex = 1.2)
    
    if (ret) {
      return(complete(imputed_data, "stacked"))
    }
}

antibodies.imputed.mice <- compare_imputations(antibodies, "Antibodies", "norm.predict", TRUE)
perc_inhibition_orig_variant.imputed.mice <- compare_imputations(perc_inhibition_orig_variant, "Original", "norm.predict", TRUE)
perc_inhibition_delta_variant.imputed.mice <- compare_imputations(perc_inhibition_delta_variant, "Delta", "norm.predict", TRUE)

# Using Hmisc::impute with random method
antibodies.imputed.hmisc <- impute(antibodies, "random")
antibodies.imputed.hmisc <- as.data.frame(as.list(antibodies.imputed.hmisc))
colnames(antibodies.imputed.hmisc) <- colnames(antibodies)
perc_inhibition_orig_variant.imputed.hmisc <- impute(perc_inhibition_orig_variant, "random")
perc_inhibition_orig_variant.imputed.hmisc <- as.data.frame(as.list(perc_inhibition_orig_variant.imputed.hmisc))
colnames(perc_inhibition_orig_variant.imputed.hmisc) <- colnames(perc_inhibition_orig_variant)
perc_inhibition_delta_variant.imputed.hmisc <- impute(perc_inhibition_delta_variant, "random")
perc_inhibition_delta_variant.imputed.hmisc <- as.data.frame(as.list(perc_inhibition_delta_variant.imputed.hmisc))
colnames(perc_inhibition_delta_variant.imputed.hmisc) <- colnames(perc_inhibition_delta_variant)

draw_imputed <- function(col, df.orig, df.hmisc, df.mf, df.mice, title) {
  par(mfrow = c(4, 1), mar = c(2, 2, 2, 2))  # 5 rows, 4 columns
  
  # Create histograms
  hist(df.hmisc[, col], main = paste("Hmisc", title,"(", col, ")"), xlab = title, col = "yellow")
  
  if (grepl(title, 'Delta') && col == 'Day 0') {
  } else {
    hist(df.orig[, col], main = paste("Original", title, "(", col, ")"), xlab = title, col = "orange")
    #hist(df.hmisc[, col], main = paste("Hmisc", title,"(", col, ")"), xlab = title, col = "yellow")
    hist(df.mf[, col], main = paste("MissForest", title, "(", col, ")"), xlab = title, col = "green")
    hist(df.mice[, col], main = paste("MICE", title,"(", col, ")"), xlab = title, col = "blue")
  }
  
  # Reset plotting area
  par(mfrow = c(1, 1))
}

draw_imputed("Day 0", antibodies, antibodies.imputed.hmisc, antibodies.imputed.mf, antibodies.imputed.mice, "Antibodies")
draw_imputed("Day 14", antibodies, antibodies.imputed.hmisc, antibodies.imputed.mf, antibodies.imputed.mice, "Antibodies")
draw_imputed("Day 30", antibodies, antibodies.imputed.hmisc, antibodies.imputed.mf, antibodies.imputed.mice, "Antibodies")
draw_imputed("Day 90", antibodies, antibodies.imputed.hmisc, antibodies.imputed.mf, antibodies.imputed.mice, "Antibodies")
draw_imputed("Day 180", antibodies, antibodies.imputed.hmisc, antibodies.imputed.mf, antibodies.imputed.mice, "Antibodies")

draw_imputed("Day 0", perc_inhibition_orig_variant, perc_inhibition_orig_variant.imputed.hmisc, perc_inhibition_orig_variant.imputed.mf, perc_inhibition_orig_variant.imputed.mice, "%Inhibition Wuhan Variant")
draw_imputed("Day 14", perc_inhibition_orig_variant, perc_inhibition_orig_variant.imputed.hmisc, perc_inhibition_orig_variant.imputed.mf, perc_inhibition_orig_variant.imputed.mice, "%Inhibition Wuhan Variant")
draw_imputed("Day 30", perc_inhibition_orig_variant, perc_inhibition_orig_variant.imputed.hmisc, perc_inhibition_orig_variant.imputed.mf, perc_inhibition_orig_variant.imputed.mice, "%Inhibition Wuhan Variant")
draw_imputed("Day 90", perc_inhibition_orig_variant, perc_inhibition_orig_variant.imputed.hmisc, perc_inhibition_orig_variant.imputed.mf, perc_inhibition_orig_variant.imputed.mice, "%Inhibition Wuhan Variant")
draw_imputed("Day 180", perc_inhibition_orig_variant, perc_inhibition_orig_variant.imputed.hmisc, perc_inhibition_orig_variant.imputed.mf, perc_inhibition_orig_variant.imputed.mice, "%Inhibition Wuhan Variant")

draw_imputed("Day 0", perc_inhibition_delta_variant, perc_inhibition_delta_variant.imputed.hmisc, perc_inhibition_delta_variant.imputed.mf, perc_inhibition_delta_variant.imputed.mice, "%Inhibition Delta Variant")
draw_imputed("Day 14", perc_inhibition_delta_variant, perc_inhibition_delta_variant.imputed.hmisc, perc_inhibition_delta_variant.imputed.mf, perc_inhibition_delta_variant.imputed.mice, "%Inhibition Delta Variant")
draw_imputed("Day 30", perc_inhibition_delta_variant, perc_inhibition_delta_variant.imputed.hmisc, perc_inhibition_delta_variant.imputed.mf, perc_inhibition_delta_variant.imputed.mice, "%Inhibition Delta Variant")
draw_imputed("Day 90", perc_inhibition_delta_variant, perc_inhibition_delta_variant.imputed.hmisc, perc_inhibition_delta_variant.imputed.mf, perc_inhibition_delta_variant.imputed.mice, "%Inhibition Delta Variant")
draw_imputed("Day 180", perc_inhibition_delta_variant, perc_inhibition_delta_variant.imputed.hmisc, perc_inhibition_delta_variant.imputed.mf, perc_inhibition_delta_variant.imputed.mice, "%Inhibition Delta Variant")

# Choose which imputation worked for which dataset
antibodies.imputed <- antibodies.imputed.hmisc
perc_inhibition_orig_variant.imputed <- perc_inhibition_orig_variant.imputed.mf
perc_inhibition_delta_variant.imputed <- perc_inhibition_delta_variant.imputed.hmisc

# Normalize values
antibodies.imputed.norm <- log(antibodies.imputed)
perc_inhibition_orig_variant.imputed.norm <- log(perc_inhibition_orig_variant.imputed)
perc_inhibition_delta_variant.imputed.norm <- log(perc_inhibition_delta_variant.imputed)
