## some basic scripts for EvoGenomics course 2024


# Function to read depth information from a file
# The file is expected to have three columns: Chromosome, Position, and Depth
read_depth <- function(file) {
  # Read the file into a data frame with specified column names
  read.table(file = file,
             header = FALSE,  # No header in the file
             col.names = c("Chromosome", "Position", "Depth"))
}


# Function to plot the read depth across chromosomes
plot_depth <- function(sample) {
  # Calculate median depth for each chromosome
  medians <- sample %>%
    group_by(Chromosome) %>%
    summarize(medianDepth = median(Depth, na.rm = TRUE)) 
  
  # Plot the read depth using ggplot2 and make it interactive with ggplotly
  print(
    ggplotly(
      ggplot(sample, aes(x = Position, y = Depth)) +
        geom_line(color = "blue", linewidth = 0.1) +
        labs(title = deparse(substitute(sample)),  # Use the sample name as the title
             x = "Position in chromosome",
             y = "Per-site read depths") +  
        theme_bw() +
        facet_wrap(~ Chromosome, ncol = 1) +  # Facet by chromosome
        ylim(0, 500) +  # Set y-axis limits
        scale_x_continuous(labels = label_number()) +  # Use non-scientific labels for x-axis
        geom_hline(data = medians, aes(yintercept = medianDepth), color = "red") +  # Add median depth lines
        theme(axis.title = element_text(face = "bold"),
              axis.text.x = element_text(face = "bold")
        )
    )
  )
}


# Function to plot the histogram of read depths
hist_depth <- function(sample) {
  # Calculate median depth for each chromosome
  medians <- sample %>%
    group_by(Chromosome) %>%
    summarize(medianDepth = median(Depth, na.rm = TRUE)) 
  
  # Plot the histogram and make it interactive with ggplotly
  print(
    ggplotly(
      ggplot(sample, aes(x = Depth)) + 
        geom_histogram(fill = "lightblue", color = "black", alpha = 0.8) +
        xlim(0,200) +
        theme_bw() +
        facet_wrap(~ Chromosome, ncol = 2) +  # Facet by chromosome
        theme(axis.title = element_text(face = "bold"),
              axis.text.x = element_text(face = "bold")
        ) + 
        labs(title = deparse(substitute(sample)), x = "Depth", y = "Frequency") + 
        geom_vline(data = medians, aes(xintercept = medianDepth), color = "red") +  # Add median depth lines
        theme(axis.title = element_text(face = "bold"),
              axis.text.x = element_text(face = "bold"))
    )
  )
}


# Function to calculate the median read depth for each chromosome
calc_somy <- function(sample) {
  # Calculate median depth for each chromosome and return as a data frame
  median_depths <- sample %>%
    group_by(Chromosome) %>%
    summarize(MedianDepth = median(Depth, na.rm = TRUE)) %>%
    ungroup()
  
  return(as.data.frame(median_depths))
}


# Function to calculate local median depth within a specific chromosome
calc_local <- function(sample) {
  # Subset for a specific chromosome
  subx <- sample[sample$Chromosome == "LbrM_10_v4_47",]
  # Define windows across the chromosome
  windows <- seq(1, tail(subx$Position,n=1), 10000)
  # Initialize matrix for results
  res <- matrix(ncol = 2, nrow = length(windows))
  
  # Loop through windows to calculate median depth
  for (n in 1:length(windows)) {
    res[n,1] <- windows[n]
    res[n,2] <- median(subx[which(subx$Position > windows[n] & subx$Position <= windows[n]+9999),3], na.rm = TRUE)
  }
  res <- as.data.frame(res)
  colnames(res) <- c('Window', 'Depth')
  return(res)
}


# Function to calculate Z-score for depth values
calc_z <- function(x) (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)

plot_z <- function(x) {
  print(
    #ggplotly(
      ggplot(x, aes(x = Zscore)) + 
        geom_histogram(aes(y = ..density..), binwidth = .2, fill = "lightblue", color = "black") + 
        geom_density(color = "blue", size = 1) + 
        labs(x = "Z-score", y = "Density") +
        theme_bw() +
        theme(axis.title = element_text(face = "bold"),
              axis.text.x = element_text(face = "bold")
        ) 
    #)
  )
  
}

library(plotly)
library(ggplot2)
library(dplyr)
library(scales)

#LC2434A1
LC2434A1 <- read_depth('LC2434A1_sorted.depth')
plot_depth(LC2434A1)

hist_depth(LC2434A1)

LC2434A1.genome.median <- median(LC2434A1[,3]) 
LC2434A1.genome.median

LC2434A1.chr.median <- calc_somy(LC2434A1) 
LC2434A1.chr.median

LC2434A1.chr.median$HaploidCN <- LC2434A1.chr.median$MedianDepth/LC2434A1.genome.median 
LC2434A1.chr.median

#LC1408A1
LC1408A1 <- read_depth('LC1408A1_sorted.depth')
plot_depth(LC1408A1)

hist_depth(LC1408A1)

LC1408A1.genome.median <- median(LC1408A1[,3]) 
LC1408A1.genome.median

LC1408A1.chr.median <- calc_somy(LC1408A1) 
LC1408A1.chr.median

LC1408A1.chr.median$HaploidCN <- LC1408A1.chr.median$MedianDepth/LC1408A1.genome.median 
LC1408A1.chr.median

#PER069A1
PER069A1 <- read_depth('PER069A1_sorted.depth')
plot_depth(PER069A1)

hist_depth(PER069A1)

PER069A1.genome.median <- median(PER069A1[,3]) 
PER069A1.genome.median

PER069A1.chr.median <- calc_somy(PER069A1) 
PER069A1.chr.median

PER069A1.chr.median$HaploidCN <- PER069A1.chr.median$MedianDepth/PER069A1.genome.median 
PER069A1.chr.median

#LC2434A1 chromosome 10
LC2434A1.chr10 <- calc_local(LC2434A1)
head(LC2434A1.chr10)

LC2434A1.chr10.median <- median(LC2434A1[which(LC2434A1$Chromosome=='LbrM_10_v4_47'),3])

LC2434A1.chr10$HaploidCN <- LC2434A1.chr10$Depth/LC2434A1.chr10.median

head(LC2434A1.chr10)

#LC1408A1 chromosome 10
LC1408A1.chr10 <- calc_local(LC1408A1)
head(LC1408A1.chr10)

LC1408A1.chr10.median <- median(LC1408A1[which(LC1408A1$Chromosome=='LbrM_10_v4_47'),3])

LC1408A1.chr10$HaploidCN <- LC1408A1.chr10$Depth/LC1408A1.chr10.median

head(LC1408A1.chr10)

#PER069A1 chromosome 10
PER069A1.chr10 <- calc_local(PER069A1)
head(PER069A1.chr10)

PER069A1.chr10.median <- median(PER069A1[which(PER069A1$Chromosome=='LbrM_10_v4_47'),3])

PER069A1.chr10$HaploidCN <- PER069A1.chr10$Depth/PER069A1.chr10.median

head(PER069A1.chr10)

LC2434A1.chr10$Zscore <- calc_z(LC2434A1.chr10$HaploidCN)
plot_z(LC2434A1.chr10)

LC1408A1.chr10$Zscore <- calc_z(LC1408A1.chr10$HaploidCN)
plot_z(LC1408A1.chr10)

PER069A1.chr10$Zscore <- calc_z(PER069A1.chr10$HaploidCN)

cnv <- LC2434A1.chr10$HaploidCN - PER069A1.chr10$HaploidCN
plot(cnv, type = 'l', lwd = 2, ylab = 'scaled haploid read depth')

cnv2 <- LC2434A1.chr10$HaploidCN - LC1408A1.chr10$HaploidCN
plot(cnv2, type = 'l', lwd = 2, ylab = 'scaled haploid read depth')

