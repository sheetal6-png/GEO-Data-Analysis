
#------------------------------------------------------------
# Install required packages (Run only once)
#------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c(
  "GEOquery",
  "Biobase"
)

for(pkg in packages){
  
  if(!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
  
  library(pkg, character.only = TRUE)
}

#------------------------------------------------------------
# Create Project Folders
#------------------------------------------------------------

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

#------------------------------------------------------------
# Download GEO Dataset
# Replace GSE25099 with any GEO accession if needed
#------------------------------------------------------------

geo_id <- "GSE25099"

gse <- getGEO(
  GEO = geo_id,
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)

#------------------------------------------------------------
# Display Basic Information
#------------------------------------------------------------

print(gse)

length(gse)

#------------------------------------------------------------
# Extract First ExpressionSet
#------------------------------------------------------------

eset <- gse[[1]]

#------------------------------------------------------------
# Expression Matrix
#------------------------------------------------------------

expression_matrix <- exprs(eset)

dim(expression_matrix)

head(expression_matrix)

# Save expression matrix

write.csv(
  expression_matrix,
  "data/expression_matrix.csv"
)

#------------------------------------------------------------
# Sample Metadata
#------------------------------------------------------------

sample_information <- pData(eset)

head(sample_information)

write.csv(
  sample_information,
  "data/sample_metadata.csv"
)

#------------------------------------------------------------
# Feature Annotation
#------------------------------------------------------------

feature_annotation <- fData(eset)

head(feature_annotation)

write.csv(
  feature_annotation,
  "data/feature_annotation.csv"
)

#------------------------------------------------------------
# Platform Information
#------------------------------------------------------------

platform_id <- annotation(eset)

cat("Platform:", platform_id, "\n")

platform <- getGEO(platform_id)

platform_table <- Table(platform)

write.csv(
  platform_table,
  "data/platform_annotation.csv"
)

#------------------------------------------------------------
# Save Complete ExpressionSet
#------------------------------------------------------------

save(
  eset,
  file = "data/GSE25099_ExpressionSet.RData"
)

#------------------------------------------------------------
# Summary
#------------------------------------------------------------

cat("--------------------------------------\n")
cat("Dataset Download Completed Successfully\n")
cat("--------------------------------------\n")

cat("GEO Accession :", geo_id, "\n")
cat("Platform      :", platform_id, "\n")
cat("Genes         :", nrow(expression_matrix), "\n")
cat("Samples       :", ncol(expression_matrix), "\n")

cat("--------------------------------------\n")
