# Part 1: Check assumptions
# Antibodies
# Normality (Shapiro-Wilk Test)
# Ho: The residuals are normally distributed.
# Ha: The residuals are not normally distributed.
shapiro.test(residuals(lm(`Day 0` ~ ., antibodies.imputed.norm)))    # p-value = 2.2 * 10^-16; reject Ho
# Violation of Normality (Kruskal)
# Ho: The distributions of expression values between patient groups are equal.
# Ha: The distributions of expression values between patient groups are not equal.
kruskal.test(`Day 0` ~ `Day 14`, antibodies.imputed.norm)   # p-value = 0.5218; accept Ho
# Homoscedasticity (Breusch-Pagan Test)
# Ho: The variances are equal.
# Ha: The variances are not equal.
bptest(lm(`Day 0` ~ ., antibodies.imputed.norm), studentize = FALSE)   # p-value = 0.3023; accept Ho

# Since normality is violated, we will switch to Repeated Measures ANOVA in the nonparametric method (Friedman Test).

# %Inhibition Wuhan variant
# Normality (Shapiro-Wilk Test)
# Ho: The residuals are normally distributed.
# Ha: The residuals are not normally distributed.
shapiro.test(residuals(lm(`Day 0` ~ ., perc_inhibition_orig_variant.imputed.norm)))    # p-value = 1.803 * 10^-8; reject Ho
# Violation of Normality (Kruskal)
# Ho: The distributions of expression values between patient groups are equal.
# Ha: The distributions of expression values between patient groups are not equal.
kruskal.test(`Day 0` ~ `Day 14`, perc_inhibition_orig_variant.imputed.norm)   # p-value = 0.3353; accept Ho
# Homoscedasticity (Breusch-Pagan Test)
# Ho: The variances are equal.
# Ha: The variances are not equal.
bptest(lm(`Day 0` ~ ., perc_inhibition_orig_variant.imputed.norm), studentize = FALSE)   # p-value = 0.000967; reject Ho

# Since normality and homoscedasticity is violated, we will switch to Repeated Measures ANOVA in the nonparametric method (Friedman Test).

# %Inhibition Wuhan variant
# Normality (Shapiro-Wilk Test)
# Ho: The residuals are normally distributed.
# Ha: The residuals are not normally distributed.
shapiro.test(residuals(lm(`Day 0` ~ ., perc_inhibition_delta_variant.imputed.norm)))    # p-value = 2.2 * 10^-16; reject Ho
# Violation of Normality (Kruskal)
# Ho: The distributions of expression values between patient groups are equal.
# Ha: The distributions of expression values between patient groups are not equal.
kruskal.test(`Day 0` ~ `Day 14`, perc_inhibition_delta_variant.imputed.norm)   # p-value = 0.3529; accept Ho
# Homoscedasticity (Breusch-Pagan Test)
# Ho: The variances are equal.
# Ha: The variances are not equal.
bptest(lm(`Day 0` ~ ., perc_inhibition_delta_variant.imputed.norm), studentize = FALSE)   # p-value = 1.065 * 10 -7; reject Ho

# Since normality and homoscedasticity is violated, we will switch to Repeated Measures ANOVA in the nonparametric method (Friedman Test).

# Part 2: Create a better dataframe to be used in analysis
build_and_friedman_test <- function (df) {
  a <- data.frame(
    patient = as.factor(rep(rownames(df), 5)),
    time = as.factor(c(rep("Day 0", nrow(df)), rep("Day 14", nrow(df)), rep("Day 30", nrow(df)), rep("Day 90", nrow(df)), rep("Day 180", nrow(df)))),
    values = c(df$`Day 0`, df$`Day 14`, df$`Day 30`, df$`Day 90`, df$`Day 180`)
  )
  
  friedman.test(y = a$values, groups = a$time, blocks = a$patients)
}

# Perform Friedman Test
build_and_friedman_test(antibodies.imputed)
build_and_friedman_test(perc_inhibition_orig_variant.imputed)
build_and_friedman_test(perc_inhibition_delta_variant.imputed)


# Part 3: Multiple regression
build_and_lm <- function (df) {
  a <- data.frame(
    patient = as.factor(rep(rownames(df), 5)),
    time = as.factor(c(rep("Day 0", nrow(df)), rep("Day 14", nrow(df)), rep("Day 30", nrow(df)), rep("Day 90", nrow(df)), rep("Day 180", nrow(df)))),
    values = c(df$`Day 0`, df$`Day 14`, df$`Day 30`, df$`Day 90`, df$`Day 180`)
  )
  
  model <- lm(values ~ time, a)
  vif(model)
}

summary(build_and_lm(antibodies.imputed))
summary(build_and_lm(perc_inhibition_orig_variant.imputed))
summary(build_and_lm(perc_inhibition_delta_variant.imputed))


# Part 2: Build dataset into one dataframe
build_dataset <- function (df) {
  return(data.frame(
    patient = as.factor(rep(df$code, 5)),
    time = as.factor(c(rep("Day 0", nrow(df)), rep("Day 14", nrow(df)), rep("Day 30", nrow(df)), rep("Day 90", nrow(df)), rep("Day 180", nrow(df)))),
    values = c(df$`Day 0`, df$`Day 14`, df$`Day 30`, df$`Day 90`, df$`Day 180`)
  ))
}

antibodies.single <- build_dataset(antibodies)
perc_inhibition_orig_variant.single <- build_dataset(perc_inhibition_orig_variant)
perc_inhibition_delta_variant.single <- build_dataset(perc_inhibition_delta_variant)

# Combine dataset of percentage and antibody count and the vaccine received
antibodies.dose <- merge(x = antibodies.single, y = data.vaccine, by = "patient", all.x = TRUE)
perc_inhibition_orig_variant.dose <- merge(x = perc_inhibition_orig_variant.single, y = data.vaccine, by = "patient", all.x = TRUE)
perc_inhibition_delta_variant.dose <- merge(x = perc_inhibition_delta_variant.single, y = data.vaccine, by = "patient", all.x = TRUE)
