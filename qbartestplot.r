qbardata = read.csv("qbartest1.csv")
library(ggplot2)
library(reshape2)
data_melt <- melt(qbardata, id.vars = "Number", variable.name = "Method", value.name = "EstimationLoss")
data_melt
ggplot(data_melt, aes(x = Number, y = EstimationLoss, color = Method)) +
geom_line(size = 1) +
geom_point(size = 2) +
labs(title = "Estimation Loss for Different Methods vs Number of Simulations",
x = "Number of Simulations",
y = "Estimation Loss") +
theme_minimal() +
theme(legend.title = element_blank(),
plot.title = element_text(hjust = 0.5, size = 16))