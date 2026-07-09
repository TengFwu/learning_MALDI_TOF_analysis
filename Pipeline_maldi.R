
# Packages
library(MALDIquant)                         # MALDI-TOF data manipulation
library(MALDIquantForeign)                  # MALDI-TOF data import
library(tidyverse)                          # Data manipulation and visualization
library(patchwork)                          # Arranging ggplots
library(pheatmap)                           # Heatmap plot



setwd("D:/Pós-graduação/Mass_spectometry")


# =========================
# Variance Stabilization
# =========================

# --- 2. Exploratory analysis --- #

# Load the data

path <- "Raw MALDI Data/TLM-1/"

raw_data <- importBrukerFlex(path) # From Bruker machine

# Quality control analysis

QC(raw_data)


# ------------------------------------------- #
#        Pre-process - step by step
# ------------------------------------------- #


# This was based in default pipeline in MALDquant package

# You have to choose the best parameters in each step

# Variance stabilization
var_stab <- transformIntensity(raw_data, method = 'sqrt')

# Smoothing
smo <- smoothIntensity(var_stab, method = "SavitzkyGolay", halfWindowSize = 10)

# Baseline Correction
# Removing the baseline

baseline <- removeBaseline(smo, method = "SNIP", iterations = 100)

# Calibration and Normalization

norma <- calibrateIntensity(baseline, method = "TIC")


# Plot
all_steps_plot(raw_data, var_stab, smo, baseline, norma, choose_one_to_plot = 1)


# Considering the replicates
# We will get the mean of the replicate - same we did in PCA

samples <- factor(sapply(norma,
                         function(x)metaData(x)$sampleName)) # choose the common name of all replicates



samples
avgSpectra <- averageMassSpectra(norma, labels = samples, method = "mean")


# Plot comparing before and after
par(mfrow = c(2,2))
plot(norma[[1]], main = "")
plot(norma[[2]], main = "")
plot(norma[[3]], main = "")
plot(avgSpectra[[1]], main = "")

par(mfrow = c(1,1))
plot(avgSpectra[[1]], main = "")
lines(norma[[1]], col = "blue", main = "")
lines(norma[[2]], col = "red", main = "")
lines(norma[[3]], col = "green", main = "")
lines(avgSpectra[[1]], main = "")


# Peak detection

noise <- estimateNoise(avgSpectra[[1]])
plot(avgSpectra[[1]], xlim = c(2000, 8000), main ="")
lines(noise, col= "red")
lines(noise[, 1], noise[, 2]*2, col = 'blue')



peaks <- detectPeaks(avgSpectra, method="MAD",
                     halfWindowSize=20, SNR=2)
plot(avgSpectra[[1]], xlim=c(2000, 8000), main = "")
points(peaks[[1]], col="red", pch=4)




par(mfrow = c(1,2))

plot(avgSpectra[[1]], xlim=c(2000, 12000), main = "WBC")
points(peaks[[1]], col="red", pch=4)

plot(avgSpectra[[6]], xlim=c(2000, 15000), main = "MEL")
points(peaks[[6]], col="red", pch=4)




par(mfrow = c(1,1))
plot(avgSpectra[[1]], xlim = c(2000, 15000), type = "l", col = "blue", main ="")
lines(avgSpectra[[6]], col = "red")

legend("topright",
       legend = c("WBC", "MEL"),
       col = c("blue", "red"),
       lty = 1)

# Peak binning

peaks <- binPeaks(peaks, tolerance = 0.002)


peaks <- filterPeaks(peaks, minFrequency = 0.25)

featureMatrix <- intensityMatrix(peaks, avgSpectra)

samples <- factor(sapply(avgSpectra,
                         function(x)metaData(x)$sampleName))

rownames(featureMatrix) <- samples

head(featureMatrix[, 1:3])

featureMatrix2 <- as.data.frame(featureMatrix)

par(mfrow = c(1,1))



# ==================
# Manipulating data 
# ==================

library(tidyverse)

# Identify the rotules

featureMatrix3 <- featureMatrix2 %>%
  mutate(group = ifelse(row_number() <= 3, "WBC", "MEL"))

featureMatrix3 <- featureMatrix3 %>% 
  relocate(group)


#write.csv(featureMatrix3,
          file = "MEL_VS_WBC.csv", 
          row.names = TRUE)

# ======================
# Exploratory analysis
# =====================

pca_data <- pca_all_data(raw_data = raw_data, replicate = "sampleName") # This function do a default pre-process in the data

# Reclassifying the data
pca_data <- pca_data %>%
  mutate(Tissue = ifelse(row_number() <= 3, "DOG", "MEL"))


run_pca(pca_data, "Tissue", "DOG", "MEL")



# Heatmap

heatmap_data <- featureMatrix3 %>% 
  select(-group)


annotation_col <- data.frame(Tissue = featureMatrix3$group)
rownames(annotation_col) <- rownames(featureMatrix3)

pheatmap(t(heatmap_data),
         annotation_col = annotation_col,
         color = colorRampPalette(c("navy", "white", "red"))(30),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = FALSE,
         show_colnames = TRUE, 
         fontsize_col = 6)



# Different peaks


# Remover a coluna 1
cols <- colnames(featureMatrix3)[-1]


# Criar vetores para armazenar resultados
pvals <- numeric(length(cols))

# Loop para calcular p-value por coluna
for (i in seq_along(cols)) {
  feature <- featureMatrix3[[cols[i]]]
  
  # Teste estatístico entre grupos de Tissue
  pvals[i] <- t.test(feature ~ featureMatrix3$group)$p.value
}

# Ajuste por múltiplas comparações (FDR)
fdr_vals <- p.adjust(pvals, method = "fdr")

# Criar tabela final
results <- data.frame(
  Feature = cols,
  P_value = pvals,
  FDR = fdr_vals
)

results

results$FDR < 0.05



