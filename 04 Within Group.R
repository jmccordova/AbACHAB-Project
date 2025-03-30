# With Compound Symmetry
gls.comp <- gls(values ~ vaccineType * time,
                     correlation = corCompSymm(form = ~ patients),
                     data = a,
                     method = "REML",
                     na.action = "na.omit",
                     control = list(singular.ok = TRUE)
              )

# Unstructured
gls.un <- gls(values ~ patients:time + patients + time,
              weights = varIdent(form = ~ 1|patients),
              correlation = corCompSymm(form = ~ patients|time),
              data = a,
              method = "REML",
              na.action = "na.omit",
              control = list(singular.ok = TRUE)
            )
summary(fit.un)