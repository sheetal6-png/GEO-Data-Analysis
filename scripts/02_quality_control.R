
#
# Description:
# Perform quality control on the downloaded expression matrix.
#############################################################

#-----------------------------
# Load Required Packages
#-----------------------------

library(ggplot2)

#-----------------------------
# Create folders if needed
#-----------------------------

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

#-----------------------------
# Load Expression Matrix
#-----------------------------

expr <- read.csv(
  "data/expression_matrix.csv",
  row.names = 1,
  check.names = FALSE
)

expr <- as.matrix(expr)

cat("Expression Matrix Loaded Successfully\n")
cat("--------------------------------------\n")
cat("Genes   :", nrow(expr), "\n")
cat("Samples :", ncol(expr), "\n")

#-----------------------------
# Missing Values
#-----------------------------

cat("--------------------------------------\n")
cat("Missing Values :", sum(is.na(expr)), "\n")

#-----------------------------
# Summary Statistics
#-----------------------------

summary(as.vector(expr))

#-----------------------------
# Boxplot
#-----------------------------

png(
  filename = "figures/01_Boxplot_Before_Normalization.png",
  width = 1800,
  height = 900
)

boxplot(
  expr,
  las = 2,
  outline = FALSE,
  col = "lightblue",
  main = "Expression Distribution Before Normalization",
  ylab = "Expression"
)

dev.off()

cat("Boxplot Saved\n")

#-----------------------------
# Density Plot
#-----------------------------

png(
  filename = "figures/02_Density_Plot.png",
  width = 1800,
  height = 900
)

plot(
  density(expr[,1]),
  main = "Density Plot",
  lwd = 2,
  col = "red"
)

for(i in 2:ncol(expr))
{
  lines(
    density(expr[,i]),
    col = rgb(0,0,1,0.15)
  )
}

dev.off()

cat("Density Plot Saved\n")

#-----------------------------
# Correlation Matrix
#-----------------------------

cor_matrix <- cor(expr)

write.csv(
  cor_matrix,
  "results/sample_correlation.csv"
)

cat("Correlation Matrix Saved\n")

#-----------------------------
# Sample Clustering
#-----------------------------

distance <- dist(t(expr))

hc <- hclust(distance)

png(
  filename = "figures/03_Sample_Clustering.png",
  width = 1800,
  height = 900
)

plot(
  hc,
  main = "Hierarchical Clustering of Samples",
  xlab = "Samples",
  sub = ""
)

dev.off()

cat("Sample Clustering Saved\n")

#-----------------------------
# Principal Variance
#-----------------------------

sample_variance <- apply(expr, 2, var)

variance_table <- data.frame(
  Sample = colnames(expr),
  Variance = sample_variance
)

write.csv(
  variance_table,
  "results/sample_variance.csv",
  row.names = FALSE
)

cat("Sample Variance Saved\n")

#-----------------------------
# Finish
#-----------------------------

cat("--------------------------------------\n")
cat("Quality Control Completed Successfully\n")
cat("--------------------------------------\n")
