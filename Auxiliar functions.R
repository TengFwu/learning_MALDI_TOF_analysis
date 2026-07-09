# Packages
library(MALDIquant)                         # MALDI-TOF data manipulation
library(MALDIquantForeign)                  # MALDI-TOF data import
library(tidyverse)                          # Data manipulation and visualization
library(patchwork)                          # Arranging ggplots
library(pheatmap)                           # Heatmap plot
library(caret)                              # Machine learning implementation
library(rpart)                              # Machine learning - Decision tree
library(rpart.plot)                         # Decision tree plots

# --- 1. Functions --- #

# Pre-processing plots

# Quality control analysis

QC <- function(x, choose_one_to_plot = 1) {
  # Inspection if any point of data are empty 
  empty <- any(sapply(x, isEmpty))
  
  # Inspect how many m/z values and how many samples
  samples <- table(sapply(x, length))
  
  # Inspect the control mass difference 
  diff <- all(sapply(x, isRegular))
  
  # Print the QC summary
  cat("# ----------------------------------------- #\n")
  cat("# ----------- Quality Control ------------- #\n")
  cat("# ----------------------------------------- #\n\n")
  cat("Empty data:", empty, "\n")
  cat("Control mass difference:", diff, "\n")
  cat("Number of samples and m/z values:\n")
  
  # Print the table 
  print(samples)
  
  # Plot one spectrum to visualize 
  plot(x[[choose_one_to_plot]], 
       main = paste("Spectrum", choose_one_to_plot, "- raw data"),
       xlab = "m/z",
       ylab = "Intensity")
}


# Exploratory PCA function --> with all data

pca_all_data <- function(raw_data, replicate = "name") {
  
  # 1. Default preprocessing 
  input_data <- transformIntensity(raw_data, method = "sqrt")
  input_data <- smoothIntensity(input_data, method = "SavitzkyGolay", halfWindowSize = 10)
  input_data <- removeBaseline(input_data, method = "SNIP", iterations = 100)
  input_data <- calibrateIntensity(input_data, method = "TIC")
  
  input_data <- alignSpectra(
    input_data,
    halfWindowSize = 20,
    SNR = 2,
    tolerance = 0.02,
    warpingMethod = "lowess"
  )
  
  # 2. Extract replicate labels
  samples <- factor(sapply(input_data, function(x) metaData(x)[[replicate]]))
  
  # 3. Average technical replicates 
  avgSpectra <- averageMassSpectra(input_data, labels = samples, method = "mean")
  
  # 4. Peak detection 
  peaks <- detectPeaks(
    avgSpectra,
    method = "MAD",
    halfWindowSize = 20,
    SNR = 2
  )
  
  # 5. Peak binning + filtering 
  peaks <- binPeaks(peaks, tolerance = 0.002)
  peaks <- filterPeaks(peaks, minFrequency = 0.25)
  
  # 6. Build feature matrix 
  featureMatrix <- intensityMatrix(peaks, avgSpectra)
  
  # Row names must match the averaged spectra order
  rownames(featureMatrix) <- names(avgSpectra)
  
  featureMatrix <- as.data.frame(featureMatrix)
  
  # 7. Return the matrix -------------------------------------------------------
  return(featureMatrix)
}


# Notice we changed 'groups' to 'group_col_name' for clarity
run_pca <- function(pca_data, group_col_name, class_1, class_2) {
  
  # Extract the group vector dynamically using the column name provided
  groups <- pca_data[[group_col_name]]
  
  # Drop the group column so only numeric data remains for PCA
  pca_number <- pca_data[, !(names(pca_data) %in% group_col_name)]
  
  # Perform PCA on the pure numeric data
  pca_res <- prcomp(pca_number, center = TRUE, scale. = TRUE)
  
  # Build Dataframe for Plotting 
  sample_ids <- rownames(pca_res$x)
  # If row names were stripped by dplyr, generate placeholder names
  if (is.null(sample_ids) || length(sample_ids) == 0) {
    sample_ids <- paste0("Sample_", seq_along(groups))
  }
  
  # Build Dataframe for Plotting
  pca_ggplot <- data.frame(
    Sample = sample_ids,
    Group = as.factor(groups), 
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2]
  )
  
  # Calculate Explained Variation
  pca.var <- pca_res$sdev^2
  pca.var.per <- round(pca.var / sum(pca.var) * 100, 1)
  
  # Identify the peaks that contribute most to PC1
  scores_pc1 <- pca_res$rotation[, 1]
  peaks_scores <- abs(scores_pc1)  
  peak_ranked <- sort(peaks_scores, decreasing = TRUE)
  
  # Text Outputs
  cat(sprintf("\n# %s #\n", paste(rep("-", 53), collapse = "")))
  cat("# ----------- Most important peaks in PC1 ------------- #\n")
  cat(sprintf("# %s #\n\n", paste(rep("-", 53), collapse = "")))
  for(i in 1:5) {
    cat(i, ":", names(peak_ranked)[i], "\n")
  }
  
  cat(sprintf("\n# %s #\n", paste(rep("-", 59), collapse = "")))
  cat("# ----------- Percentage of explained variation ------------- #\n")
  cat(sprintf("# %s #\n\n", paste(rep("-", 59), collapse = "")))
  for(i in 1:4) {
    cat(paste0("PC", i, ": ", pca.var.per[i], "%\n"))
  }
  
  # Build Plots
  cols <- setNames(c("black", "brown"), c(class_1, class_2))
  
  # PCA Plot
  p_pca <- ggplot(pca_ggplot, aes(x = PC1, y = PC2, color = Group, fill = Group)) + 
    geom_point(size = 4) +
    stat_ellipse(geom = "polygon", alpha = 0.2, show.legend = FALSE) +
    xlab(paste0("PC1 (", pca.var.per[1], "%)")) +
    ylab(paste0("PC2 (", pca.var.per[2], "%)")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    ggtitle("PCA Plot")
  
  # Scree Plot 
  scree_data <- data.frame(
    PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
    Variation = pca.var.per[1:10]
  )
  
  p_scree <- ggplot(scree_data, aes(x = PC, y = Variation)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black", alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(x = "Principal Component", y = "Percent Variation (%)", title = "Scree Plot") +
    ylim(0, max(scree_data$Variation, na.rm = TRUE) + 5)
  
  # Arrange and Print Plots
  combined_plots <- p_pca + p_scree
  print(combined_plots)
  
  # Return useful data
  return(invisible(list(
    pca_object = pca_res,
    variance_explained = pca.var.per,
    top_pc1_peaks = names(peak_ranked)[1:10],
    plot = combined_plots
  )))
}


# Step-by-step plots

all_steps_plot <- function(raw_data, var_stab, smo, baseline, norma, choose_one_to_plot = 1) {
  # GRID
  par(mfrow = c(2,3))
  
  xlim <- range(mass(raw_data[[choose_one_to_plot]]))
  # Raw data
  plot(raw_data[[choose_one_to_plot]], main = "1: raw data", sub = "", xlim = xlim)
  # Variance stabilization
  plot(var_stab[[choose_one_to_plot]], main = "2: Variance stabilization", sub = "", xlim = xlim)
  # Smoothing + baseline
  base <- estimateBaseline(smo[[choose_one_to_plot]], method = "SNIP", iterations = 100)
  plot(smo[[choose_one_to_plot]], main = "3: Smoothing", sub = "", xlim = xlim)
  lines(base, col = "red", lwd = 2)
  # Baseline removed
  plot(baseline[[choose_one_to_plot]], main = "4: Baseline correction", sub = "", xlim = xlim)
  # Peak detection
  p <- detectPeaks(baseline[[choose_one_to_plot]], method = "MAD", halfWindowSize = 20, SNR = 2)
  plot(baseline[[choose_one_to_plot]], main = "5: peak detection", sub = "", xlim = xlim)
  points(p)
  # Normalized data
  plot(norma[[choose_one_to_plot]], main = "6: Normalized", sub = "", xlim = xlim)
  par(mfrow = c(1,1))
}
