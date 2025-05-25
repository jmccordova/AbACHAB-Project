# Step 1: Run for Wuhan variant
for (timepoint1 in colnames(perc_inhibition_orig_variant.imputed)) {
  for (timepoint2 in colnames(perc_inhibition_orig_variant.imputed)) {
    print(paste(timepoint1, ' vs. ', timepoint2))
    # Added jitters to fix issue witn V=0
    print(wilcox.test(jitter(perc_inhibition_orig_variant.imputed[[timepoint1]]), jitter(perc_inhibition_orig_variant.imputed[[timepoint2]]), paired = TRUE))
    print("------------------")
  }
}

# Step 2: Run for Delta variant
for (timepoint1 in colnames(perc_inhibition_delta_variant.imputed)) {
  for (timepoint2 in colnames(perc_inhibition_delta_variant.imputed)) {
    print(paste(timepoint1, ' vs. ', timepoint2))
    print(wilcox.test(jitter(perc_inhibition_delta_variant.imputed[[timepoint1]]), jitter(perc_inhibition_delta_variant.imputed[[timepoint2]]), paired = TRUE))
    print("------------------")
  }
}
