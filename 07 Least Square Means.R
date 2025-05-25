# Combine all the aggregated data into one dataframe
a <- rbind(antibodies.aggregated, perc_inhibition_delta_variant.aggregated, perc_inhibition_orig_variant.aggregated)
# Factor in the source of the datatype
a$datatype <- factor(a$datatype)

# Step 1: Fit a Model
model.lmer <- lmer(values ~ datatype + time + (1 | patient), 
     data = a
)

# Step 2: Compute Least Squares Means
# provides estimated marginal means for each group at different time points
lsmeans <- emmeans(model.lmer, ~ datatype | time,
                   pbkrtest.limit = 4000,
                   lmerTest.limit = 4000)
summary(lsmeans)

# Step 3: Compute Pairwise Differences
# gives LSMeans differences between groups at each time point
pairwise_diffs <- emmeans(model.lmer, pairwise ~ datatype | time,
                          pbkrtest.limit = 4000,
                          lmerTest.limit = 4000)
summary(pairwise_diffs)

# Step 4: Extract 95% Confidence Intervals
conf_intervals <- confint(pairwise_diffs)
print(conf_intervals)
