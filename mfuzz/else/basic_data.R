# mfuzz clustering of mouses aged 2, 6 and 12
# data DESeq2 for 2-6 and 6-12 with no log2 FC
# load the data

# Install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Mfuzz", "readxl"))

install.packages("openxlsx")


# Load the libraries
library(Mfuzz)
library(readxl)
library(openxlsx)

# --- 1. Load the data ---

file_2vs6 <- "/Users/zainabnazari/mfuzz/data/miRNA-2M_6M.xlsx"
file_6vs12 <- "/Users/zainabnazari/mfuzz/data/miRNA-6M_12M.xlsx"

# Read the data
data_2vs6 <- read.xlsx(file_2vs6)
data_6vs12 <- read.xlsx(file_6vs12)

# print

dim(data_2vs6)       # Rows and columns
dim(data_6vs12)
colnames(data_2vs6)  # Names of columns
colnames(data_2vs6)[1]
