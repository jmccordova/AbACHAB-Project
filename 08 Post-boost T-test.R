# Step 1: Log-transform the data (since geometric means are based on log-transformed values).
log_titers <- log10(antibodies.imputed)

# Step 2: Compute geometric means for each group.
gm0 <- exp(mean(log_titers$`Day 0`))
gm14 <- exp(mean(log_titers$`Day 14`))
gm30 <- exp(mean(log_titers$`Day 30`))
gm90 <- exp(mean(log_titers$`Day 90`))
gm180 <- exp(mean(log_titers$`Day 180`))

# Step 3: Perform statistical comparison (e.g., t-test or ANOVA on log-transformed values).
#a <- data.frame(
#        values = c(log_titers$'Day 0', log_titers$'Day 14', log_titers$'Day 30', log_titers$'Day 90', log_titers$'Day 180'),
#        group = rep(colnames(log_titers), each = nrow(log_titers))
#      )
#pairwise.ttest <- pairwise.t.test(a$value, a$group, p.adjust.method = "bonferroni")
#print(pairwise.ttest$p.value, digits = 10)

for (timepoint1 in colnames(log_titers)) {
  for (timepoint2 in colnames(log_titers)) {
    print(paste(timepoint1, ' vs. ', timepoint2))
    print(t.test(log_titers[[timepoint1]], log_titers[[timepoint2]]))
    print("------------------")
  }
}

# Step 4: Back-transform results to interpret in the original scale.
lm_model <- lm(log10(`Day 0`) ~ ., data = antibodies.imputed)
summary(lm_model)
print(coef(lm_model))
