# Functions for normalized coverage analysis
### Functions



# load libraries
libraries <- c("data.table", "dplyr", "tidyr", "ggplot2", "scales", "cowplot", "stringr", "Cairo", "grid")

for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}


rm(list = ls()) # clean env.

# Function to read the files with normalized coverage data:

read.NormCov.file <- function(sample){
  file2read <- paste0(sample, "_NormalizedCoverage.txt")
  
  df.tmp <- fread(file2read)
  df.tmp$Sample <- paste0(sample)
  
  return(df.tmp)
  
}


# Function to obtained the adjusted normalized coverage:

adjust.coverage <- function(dataset.coverage, threshold.cov){
  
  threshold.cov <- as.numeric(threshold.cov)
  
  
  
  # estimate normalized coverage (coverage in each window/global median coverage)
  dataset.coverage <- dataset.coverage %>% 
    group_by(Sample) %>% 
    # Estimate a global normalized coverage for each sample
    mutate(Global.Norm.Cov = median(normalized.coverage)) %>% 
    
    
    mutate(Adj.Norm.Cov = if_else(normalized.coverage <= (Global.Norm.Cov * threshold.cov), 
                                  0.1,
                                  normalized.coverage)) %>% 
    ungroup() %>% setDT()
  
  # select a specific chromosome
  #df.cov.chr <- dataset.coverage %>% filter(chr == target.chr)
  
  
  return(dataset.coverage)  
  
}

# Function to compute log2
calculate_log2ratio <- function(data.set, ref_col, end_col) {
  
  data.set <- data.set %>%
    # Select the next columns:
    select(Sample, chr, window.start, window.end, Adj.Norm.Cov) %>% 
    # Transform the table with pivot_wider
    pivot_wider(names_from = Sample, values_from = Adj.Norm.Cov) %>% 
    # Divide each value in the indicated columns / reference column (SG200)
    mutate(across(c(!!sym(ref_col): !!sym(end_col)),
                  .fns = ~./ !!sym(ref_col),
                  .names = "{.col}_Ratio" ))
  
  ##
  ##
  ref_col_Ratio <- paste0(ref_col, "_Ratio")
  end_col_Ratio <- paste0(end_col, "_Ratio")
  
  data.set <- data.set %>% 
    select(chr, window.start, window.end, !!sym(ref_col_Ratio):!!sym(end_col_Ratio)) %>% 
    pivot_longer(cols = !!sym(ref_col_Ratio):!!sym(end_col_Ratio), 
                 names_to = "sample", values_to = "ratio") %>% 
    mutate(log2ratio = log2(ratio))
  
  data.set$sample <- gsub("_Ratio", "", data.set$sample)
  
  return(data.set)
  
}

## Function to plot the log2 ratio as heatmap
plot.heatmap.log2 <- function(data.set, target.chr){
  # Subset the data set to a specific chr
  df.AllSamples2plot.1 <- data.set %>% 
    filter(chr == target.chr)
  
  # Create a label for the x-axis 
  chr.label<-unique( paste0("\n",gsub("USMA_521_v2_", "Chromosome ", df.AllSamples2plot.1$chr), " (kb)"))
  
  # Plot with ggplot
  my.plot <- ggplot(df.AllSamples2plot.1) +
    geom_tile(aes(x = (window.end/1000), 
                  y = y.order, 
                  fill = log2ratio)) +
    #geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5), color = "white", linewidth = 2)+
    facet_grid(~ Line, space = "free", scales = "free") +
    geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5), color = "white", linewidth = 2)+
    facet_grid(~ Line, space = "free", scales = "free") +
    scale_x_continuous(limits = c(0, max(df.AllSamples2plot.1$window.end)/1000),
                       breaks = c(seq(0, max(df.AllSamples2plot.1$window.end)/1000, 150), 
                                  (max(df.AllSamples2plot.1$window.end)/1000))) +
    #scale_fill_gradient(low="white", high="blue", limits = c(-1,3)) + # GrB  
    scale_fill_gradient(low="white", high="blue", limits = c(-1, 1.75), breaks = c(0, 0.5, 1, 1.5)) + # GrB  
    #scale_fill_gradient(low="white", high="blue", limits = c(0,3.5)) + # GrB  original
    theme_classic()+
    labs(x = chr.label, y = "Colonies\n") +
    theme(
      # titles
      plot.title = element_blank(),
      axis.title.y = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 12, color = "black"),
      # axis text
      axis.text.x = element_text(size = 9, color = "black", angle = 90),
      axis.text.y = element_text(size = 9, color = "black"),
      # axis ticks
      axis.ticks.y = element_blank(),
      
      #axis.text = element_blank(),
      #axis.text.y = element_blank(),
      axis.line = element_blank(),
      # strip
      strip.background = element_blank(),
      strip.text = element_text(size = 12, color = "black"),
      legend.position = "right") +
    guides(fill = guide_colorbar(title = expression("log"["2"]*" (Ratio)"),
                                 barwidth = 0.5, ticks = T))
  
  
  return(my.plot)
  
}