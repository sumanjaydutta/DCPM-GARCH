var = read.csv("dccminvar.csv")
library(ggplot2)
library(reshape2)
data_melt <- melt(var, id.vars = "Samples", variable.name = "Method", value.name = "Variance")
data_melt
ggplot(data_melt, aes(x = Samples, y = Variance, color = Method)) +
geom_line(size = 1) +
geom_point(size = 2) +
labs(title = "Variance for Different Methods vs Concentration (p/n) ratios",
x = "Concentration (p/n) Ratio",
y = "Variance") +
theme_minimal() +
theme(legend.title = element_blank(),
plot.title = element_text(hjust = 0.5, size = 16))