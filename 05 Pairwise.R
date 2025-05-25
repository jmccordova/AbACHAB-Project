# Combine all the aggregated data into one dataframe
a <- rbind(antibodies.aggregated, perc_inhibition_delta_variant.aggregated, perc_inhibition_orig_variant.aggregated)
# Factor in the source of the datatype
a$datatype <- factor(a$datatype)

#noise.lm <- lm(values ~ time * patient * datatype,
#               data = a
#             )

# Linear mixed-effects models
# values ~ time + datatype = This means we're modeling values as a function of time and patient (both are treated as fixed effects).
# (1 | datatype) = This part specifies that datatype is a random effect, meaning:
# Each level of datatype has its own random intercept (hence the 1).
# This accounts for variability across different datatype groups.
noise.lm <- lmer(values ~ time + datatype + (1 | patient), 
                 data = a
                 )
isSingular(noise.lm)
anova(noise.lm)
model_parameters(noise.lm)

# Pairwise comparison against different timepoints of extraction.
# These comparisons are grouped by their datatype.
emmeans.type_time <- emmeans(noise.lm, pairwise ~ time | datatype,
                      pbkrtest.limit = 4000,
                      lmerTest.limit = 4000
                      )
conf.type_time.df <- do.call(rbind, emmeans.type_time)
conf.type_time.list <- list(emmeans.type_time, emmeans.type_time)
pairs.tukey <- pairs(emmeans.type_time$emmeans, adjust="tukey")
pairs.tukey <- pairs(emmeans.type_time$contrasts, adjust="tukey")

# Pairwise comparison against different datatype.
# These comparisons are grouped by their time of extraction.
emmeans.time_type <- emmeans(noise.lm, pairwise ~ datatype | time,
        pbkrtest.limit = 4000,
        lmerTest.limit = 4000
)

conf.time_type.df <- do.call(rbind, emmeans.time_type)
conf.time_type.list <- list(emmeans.time_type, emmeans.time_type)
pairs.tukey <- pairs(emmeans.time_type$emmeans, adjust="tukey")
pairs.tukey <- pairs(emmeans.time_type$contrasts, adjust="tukey")
