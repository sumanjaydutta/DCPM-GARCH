library(ggplot2)
library(tidyr)
qt = read.csv("qttestggm.csv")
data_long <- gather(qt, key = "Method", value = "Value", DCCGP:DCC)
# Plot the data
ggplot(data_long, aes(x = Time, y = Value, color = Method)) +
geom_line() +
labs(title = "Frobenius Norms Over Time for DCCL1, DCCL2, DCCNL",
x = "Time Index",
y = "Frobenius Norm (Magnitude 10^-8)") +
theme_minimal() +
theme(legend.title = element_blank())