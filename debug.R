metadata <- read.csv("F:/CRC_Tumour_WGS/SraRunTable.csv")
cat("Assay.Type values:", unique(metadata$Assay.Type), "\n")
meta <- metadata[metadata$Assay.Type == "WGS", ]
cat("WGS samples:", nrow(meta), "\n")

f <- "F:/CRC_Tumour_WGS/counts/SRR23821584.cds.counts.txt"
df <- read.table(f, header=TRUE, skip=1, sep="\t")
cat("Colnames:", paste(colnames(df), collapse=", "), "\n")
cat("Total counts:", sum(df[,ncol(df)]), "\n")

files <- list.files("F:/CRC_Tumour_WGS/counts", pattern=".counts.txt", full.names=TRUE)
totals <- sapply(files, function(f) {
  df <- read.table(f, header=TRUE, skip=1, sep="\t")
  sum(df[,ncol(df)])
})
names(totals) <- sub(".cds.counts.txt", "", basename(names(totals)))
print(totals)
