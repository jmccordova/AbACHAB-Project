# Build the dataset
a <- na.omit(antibodies.dose)
'%ni%' <- Negate('%in%')
a <- upSample(x = a[,colnames(a) %ni% "time"],
                     y = a$time)
a$values <- scale(a$values)[, 1]
a$time <- factor(a$Class)
a$Class <- NULL
a$dose_1_vaccine <- factor(a$dose_1_vaccine)
a$dose_2_vaccine <- factor(a$dose_2_vaccine)

antibodies.aggregated <- a
antibodies.aggregated$datatype <- rep('antibodies', nrow(antibodies.aggregated))

# AR1 Matrix
fit.ar1 <- glmmTMB(values ~ ar1(time + 0 | patient), 
                   data = a
                   )
print(paste('Converged?', fit.ar1$sdr$pdHess))
VarCorr(fit.ar1)

# Generate Unstructured Covariance Matrix
fit.us <- glmmTMB(values ~ us(time + 0 | patient), 
                  data = a, 
                  dispformula = ~0
                  )
print(paste('Converged?', fit.us$sdr$pdHess))
VarCorr(fit.us)

# Toeplitz
fit.toep <- glmmTMB(values ~ toep(time + 0 | patient), 
                    data = a,
                    dispformula=~0
                    )
print(paste('Converged?', fit.toep$sdr$pdHess))
(vc.toep <- VarCorr(fit.toep))
vc1 <- vc.toep$cond[[1]] ## first term of var-cov for RE of conditional model
summary(diag(vc1))
summary(vc1[row(vc1)!=col(vc1)])

fit.toep.reml <- update(fit.toep, REML=TRUE)
vc1R <- VarCorr(fit.toep.reml)$cond[[1]]
summary(diag(vc1R))
summary(vc1R[row(vc1R)!=col(vc1R)])


# Compound symmetry
fit.cs <- glmmTMB(values ~ cs(time + 0 | patient), 
                  data = a, 
                  dispformula=~0
                  )
print(paste('Converged?', fit.cs$sdr$pdHess))
VarCorr(fit.cs)

# Ho: The pairwise correlations are the same.
# Ha: The pairwise correlations are not the same.
anova(fit.ar1, fit.toep, fit.us, fit.cs)

# CS is the only converged model with the lowest AIC. Thus, we can reject Ho.
# The pairwise correlations between the patient's antibody values differ every time the sample was taken from the patient.

#--------------------------

# Build the dataset
a <- na.omit(perc_inhibition_orig_variant.dose)
'%ni%' <- Negate('%in%')
a <- upSample(x = a[,colnames(a) %ni% "time"],
              y = a$time)
a$values <- scale(a$values)[, 1]
a$time <- factor(a$Class)
a$Class <- NULL
a$dose_1_vaccine <- factor(a$dose_1_vaccine)
a$dose_2_vaccine <- factor(a$dose_2_vaccine)

perc_inhibition_orig_variant.aggregated <- a
perc_inhibition_orig_variant.aggregated$datatype <- rep('inhib_orig', nrow(perc_inhibition_orig_variant.aggregated))

# AR1 Matrix
fit.ar1 <- glmmTMB(values ~ ar1(time + 0 | patient), 
                   data = a
)
print(paste('Converged?', fit.ar1$sdr$pdHess))
VarCorr(fit.ar1)

# Generate Unstructured Covariance Matrix
fit.us <- glmmTMB(values ~ us(time + 0 | patient), 
                  data = a, 
                  dispformula = ~0
)
print(paste('Converged?', fit.us$sdr$pdHess))
VarCorr(fit.us)

# Toeplitz
fit.toep <- glmmTMB(values ~ toep(time + 0 | patient), 
                    data = a,
                    dispformula=~0
)
print(paste('Converged?', fit.toep$sdr$pdHess))
(vc.toep <- VarCorr(fit.toep))
vc1 <- vc.toep$cond[[1]] ## first term of var-cov for RE of conditional model
summary(diag(vc1))
summary(vc1[row(vc1)!=col(vc1)])

fit.toep.reml <- update(fit.toep, REML=TRUE)
vc1R <- VarCorr(fit.toep.reml)$cond[[1]]
summary(diag(vc1R))
summary(vc1R[row(vc1R)!=col(vc1R)])


# Compound symmetry
fit.cs <- glmmTMB(values ~ cs(time + 0 | patient), 
                  data = a, 
                  dispformula=~0
)
print(paste('Converged?', fit.cs$sdr$pdHess))
VarCorr(fit.cs)

# Ho: The pairwise correlations are the same.
# Ha: The pairwise correlations are not the same.
anova(fit.ar1, fit.toep, fit.us, fit.cs)

# CS, TOEP, and US have all converged with US having the lowest AIC. Thus, we can reject Ho.
# The pairwise correlations between the patient's %inhibition for the Wuhan strain values differ every time the sample was taken from the patient.

#--------------------------

# Build the dataset
a <- na.omit(perc_inhibition_delta_variant.dose)
'%ni%' <- Negate('%in%')
a <- upSample(x = a[,colnames(a) %ni% "time"],
              y = a$time)
a$values <- scale(a$values)[, 1]
a$time <- factor(a$Class)
a$Class <- NULL
a$dose_1_vaccine <- factor(a$dose_1_vaccine)
a$dose_2_vaccine <- factor(a$dose_2_vaccine)

perc_inhibition_delta_variant.aggregated <- a
perc_inhibition_delta_variant.aggregated$datatype <- rep('inhib_delta', nrow(perc_inhibition_delta_variant.aggregated))

# AR1 Matrix
fit.ar1 <- glmmTMB(values ~ ar1(time + 0 | patient), 
                   data = a
)
print(paste('Converged?', fit.ar1$sdr$pdHess))
VarCorr(fit.ar1)

# Generate Unstructured Covariance Matrix
fit.us <- glmmTMB(values ~ us(time + 0 | patient), 
                  data = a, 
                  dispformula = ~0
)
print(paste('Converged?', fit.us$sdr$pdHess))
VarCorr(fit.us)

# Toeplitz
fit.toep <- glmmTMB(values ~ toep(time + 0 | patient), 
                    data = a,
                    dispformula=~0
)
print(paste('Converged?', fit.toep$sdr$pdHess))
(vc.toep <- VarCorr(fit.toep))
vc1 <- vc.toep$cond[[1]] ## first term of var-cov for RE of conditional model
summary(diag(vc1))
summary(vc1[row(vc1)!=col(vc1)])

fit.toep.reml <- update(fit.toep, REML=TRUE)
vc1R <- VarCorr(fit.toep.reml)$cond[[1]]
summary(diag(vc1R))
summary(vc1R[row(vc1R)!=col(vc1R)])


# Compound symmetry
fit.cs <- glmmTMB(values ~ cs(time + 0 | patient), 
                  data = a, 
                  dispformula=~0
)
print(paste('Converged?', fit.cs$sdr$pdHess))
VarCorr(fit.cs)

# Ho: The pairwise correlations are the same.
# Ha: The pairwise correlations are not the same.
anova(fit.ar1, fit.toep, fit.us, fit.cs)

# CS, TOEP, and US have all converged with US having the lowest AIC. Thus, we can reject Ho.
# The pairwise correlations between the patient's %inhibition for the Delta strain values differ every time the sample was taken from the patient.
