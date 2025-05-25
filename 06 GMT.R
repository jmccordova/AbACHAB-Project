# Step 1: Log-transform the titers
# Since antibody titers often follow a skewed distribution, a log transformation (typically log base 10) helps normalize the data
log_titers <- log10(antibodies.imputed)

# Step 2: Perform statistical analysis
anova_model <- aov(log10(`Day 0`) ~ ., data = antibodies.imputed)
summary(anova_model)
print(coef(anova_model))

# Step 3: Take the anti-log of estimate
# After obtaining estimates (e.g., means or regression coefficients), convert them back to the original scale using the anti-log
a <- c(log_titers$'Day 0', log_titers$'Day 14', log_titers$'Day 30', log_titers$'Day 90', log_titers$'Day 180')
mean_log_titer <- mean(a)
anti_log_titer <- 10^mean_log_titer
print(paste('Mean Log Titer', mean_log_titer, ' - Anti Log Titer', anti_log_titer))
