#############################################################
# Project : GEO Data Analysis
# Dataset : GSE25099
# Script  : 03_preprocessing.R
# Author  : Sheetal Kumari
#
# Description:
# Log2 transformation (if required) and
# quantile normalization of GEO expression data.
#############################################################

#------------------------------------------------------------
# Load Required Packages
#------------------------------------------------------------

if (!requireNamespace("limma", quietly = TRUE))
    BiocManager::install("limma")

library(limma)

#------------------------------------------------------------
# Create folders
#------------------------------------------------------------

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

#------------------------------------------------------------
# Load Expression Matrix
#------------------------------------------------------------

expr <- read.csv(
    "data/expression_matrix.csv",
    row.names = 1,
    check.names = FALSE
)

expr <- as.matrix(expr)

cat("-------------------------------------\n")
cat("Expression Matrix Loaded\n")
cat("-------------------------------------\n")

cat("Genes   :", nrow(expr), "\n")
cat("Samples :", ncol(expr), "\n")

#------------------------------------------------------------
# Check if Log2 Transformation is Needed
#------------------------------------------------------------

qx <- quantile(expr,
               probs = c(0,0.25,0.5,0.75,0.99,1),
               na.rm = TRUE)

print(qx)

LogC <- (qx[5] > 100) ||
        (qx[6] - qx[1] > 50)

if(LogC)
{
    cat("Applying Log2 Transformation...\n")

    expr[expr <= 0] <- NA

    expr <- log2(expr)
}
else
{
    cat("Dataset is already Log2 transformed.\n")
}

#------------------------------------------------------------
# Quantile Normalization
#------------------------------------------------------------

cat("Performing Quantile Normalization...\n")

expr_norm <- normalizeBetweenArrays(
    expr,
    method = "quantile"
)

#------------------------------------------------------------
# Save Normalized Matrix
#------------------------------------------------------------

write.csv(
    expr_norm,
    "results/normalized_expression_matrix.csv"
)

cat("Normalized Matrix Saved.\n")

#------------------------------------------------------------
# Summary Statistics
#------------------------------------------------------------

summary_statistics <- data.frame(

    Minimum = min(expr_norm),

    Maximum = max(expr_norm),

    Mean = mean(expr_norm),

    Median = median(expr_norm),

    SD = sd(expr_norm)

)

write.csv(
    summary_statistics,
    "results/summary_statistics.csv",
    row.names = FALSE
)

#------------------------------------------------------------
# Boxplot After Normalization
#------------------------------------------------------------

png(
    filename = "figures/04_Boxplot_After_Normalization.png",
    width = 1800,
    height = 900
)

boxplot(
    expr_norm,
    las = 2,
    outline = FALSE,
    col = "lightgreen",
    main = "Expression After Quantile Normalization",
    ylab = "Normalized Expression"
)

dev.off()

cat("Boxplot Saved.\n")

#------------------------------------------------------------
# Density Plot After Normalization
#------------------------------------------------------------

png(
    filename = "figures/05_Density_After_Normalization.png",
    width = 1800,
    height = 900
)

plot(
    density(expr_norm[,1]),
    lwd = 2,
    col = "red",
    main = "Density Plot After Normalization"
)

for(i in 2:ncol(expr_norm))
{
    lines(
        density(expr_norm[,i]),
        col = rgb(0,0,1,0.15)
    )
}

dev.off()

cat("Density Plot Saved.\n")

#------------------------------------------------------------
# Save RData Object
#------------------------------------------------------------

save(
    expr_norm,
    file = "results/normalized_expression_matrix.RData"
)

#------------------------------------------------------------
# Finish
#------------------------------------------------------------

cat("-------------------------------------\n")
cat("Preprocessing Completed Successfully\n")
cat("-------------------------------------\n")
