# Step 1: Visualize antibodies per timepoint
hist(antibodies.imputed)
plot(antibodies.imputed)

# Step 2: Compute Spearman
for (timepoint1 in colnames(antibodies.imputed)) {
  for (timepoint2 in colnames(antibodies.imputed)) {
    print(paste(timepoint1, ' vs. ', timepoint2))
    # Since this dataset has so many ties, Spearman cannot proceed without warning.
    # Kendall can handle these ties.
    print(cor.test(antibodies.imputed[[timepoint1]], antibodies.imputed[[timepoint2]], method = "kendall"))
    print("------------------")
  }
}

ggplot(data.frame(antibodies.imputed$`Day 0`, antibodies.imputed$`Day 14`), aes(antibodies.imputed$`Day 0`, antibodies.imputed$`Day 14`)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Spearman Correlation Plot Day 0 vs. Day 14")

ggplot(data.frame(antibodies.imputed$`Day 0`, antibodies.imputed$`Day 30`), aes(antibodies.imputed$`Day 0`, antibodies.imputed$`Day 30`)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Spearman Correlation Plot Day 0 vs. Day 30")

ggplot(data.frame(antibodies.imputed$`Day 0`, antibodies.imputed$`Day 90`), aes(antibodies.imputed$`Day 0`, antibodies.imputed$`Day 90`)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Spearman Correlation Plot Day 0 vs. Day 90")

ggplot(data.frame(antibodies.imputed$`Day 0`, antibodies.imputed$`Day 180`), aes(antibodies.imputed$`Day 0`, antibodies.imputed$`Day 180`)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Spearman Correlation Plot Day 0 vs. Day 180")

ggplot(data.frame(antibodies.imputed$`Day 14`, antibodies.imputed$`Day 30`), aes(antibodies.imputed$`Day 14`, antibodies.imputed$`Day 30`)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Spearman Correlation Plot Day 14 vs. Day 30")

ggplot(data.frame(antibodies.imputed$`Day 14`, antibodies.imputed$`Day 90`), aes(antibodies.imputed$`Day 14`, antibodies.imputed$`Day 90`)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Spearman Correlation Plot Day 14 vs. Day 90")

ggplot(data.frame(antibodies.imputed$`Day 14`, antibodies.imputed$`Day 180`), aes(antibodies.imputed$`Day 14`, antibodies.imputed$`Day 180`)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Spearman Correlation Plot Day 14 vs. Day 180")

ggplot(data.frame(antibodies.imputed$`Day 30`, antibodies.imputed$`Day 90`), aes(antibodies.imputed$`Day 30`, antibodies.imputed$`Day 90`)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Spearman Correlation Plot Day 30 vs. Day 90")

ggplot(data.frame(antibodies.imputed$`Day 30`, antibodies.imputed$`Day 180`), aes(antibodies.imputed$`Day 30`, antibodies.imputed$`Day 180`)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Spearman Correlation Plot Day 30 vs. Day 180")

ggplot(data.frame(antibodies.imputed$`Day 90`, antibodies.imputed$`Day 180`), aes(antibodies.imputed$`Day 90`, antibodies.imputed$`Day 180`)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Spearman Correlation Plot Day 90 vs. Day 180")
