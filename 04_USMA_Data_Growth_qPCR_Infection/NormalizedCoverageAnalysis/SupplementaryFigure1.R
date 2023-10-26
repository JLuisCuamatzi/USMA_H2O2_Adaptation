# Ffigure supplementary 1
# author: jcuamatzi

# load libraries
libraries <- c("data.table", "dplyr", "tidyr", "ggplot2", "scales", "cowplot", "stringr", "Cairo", "grid", "ggtext", "ggh4x")

for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}


rm(list = ls()) # clean env.

# set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read custom functions
source(file = "FunctionsForNormCov.R")

# Vector with samples to read
#samples2read <- c("2021EE01", "2021EE18", "2021EE37", "2021EE47")
samples2read <- c("2021EE01", "2021EE18")
# 2021EE01 = sequenced with DNBSeq (BGI)
# 2021EE18 = sequenced with DNBSeq (BGI)
# 2021EE37 = sequenced with Illumina NextSeq (INMEGEN)
# 2021EE47 = sequenced with Illumina NextSeq (INMEGEN)



for (sample in samples2read){
  
  # Create object with the name of the df
  df.name <- paste0("df.", sample)
  
  # Check if the file exists before reading it
  file2read <- paste0(sample, "_NormalizedCoverage.txt")
  
  if (file.exists(file2read)) {
    # Read file into a tmp df
    df.tmp <- read.NormCov.file(sample = sample)
    
    # Assign df.name and df.tmp
    assign(df.name, df.tmp)
    
    rm(df.name, df.tmp)
    
  } else {
    
    # File does not exist, print a message and continue with the next sample
    cat("File not found:", file2read, "\n")
  }
}


# List objects based on the pattern: df.2021EE
df.list <- mget(ls(pattern = "df.2021EE"))

# Merge all the data frames in the list using rbind
df.AllSamples <- do.call(rbind, df.list)

# Remove tmp files
rm(list = ls(pattern = "df.2021EE"))

# Adjusted the normalized coverage
df.AllSamples <- adjust.coverage(dataset.coverage = df.AllSamples, threshold.cov = 0.25 )


# Separe samples by sequencing technology
df.DNBSeq <- df.AllSamples %>% filter(Sample %in% c("2021EE01", "2021EE18"))
#df.NextSeq <- df.AllSamples %>% filter(Sample %in% c("2021EE37", "2021EE47"))

# Compute log2 ratio 
df.DNBSeq <- calculate_log2ratio(data.set = df.DNBSeq, ref_col = "2021EE01", end_col = "2021EE18")
#df.NextSeq <- calculate_log2ratio(data.set = df.NextSeq, ref_col = "2021EE47", end_col = "2021EE37")

# join and remove SG200 (2021EE01 & 2021EE47)
# df.2plot <- bind_rows(df.DNBSeq, df.NextSeq)
# 
# df.2plot <- df.2plot %>% 
#   filter(!sample %in% c("2021EE01", "2021EE47") ) %>% 
#   filter(chr == "USMA_521_v2_9")

# Just DNBSeq
df.2plot <- df.DNBSeq %>% 
  filter(!sample %in% c("2021EE01") ) %>% 
  filter(chr == "USMA_521_v2_9")

# extract max length of chr 9
max.chr9 <- max(df.AllSamples %>% filter(chr == "USMA_521_v2_9") %>% 
                  select(window.end) )


# plot A. (line)
plot.A <- df.2plot %>% 
  ggplot()+
  geom_line(aes(x = window.end/1000 , y = log2ratio), color = "darkred", linewidth = 0.35) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_hline(yintercept = 2, linetype = "dashed", color = "gray", linewidth = 0.4)+
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1), 
                     minor_breaks = seq(-1, 2, 0.5),
                     guide = "axis_minor") +
  scale_x_continuous(limits = c(0 , max.chr9/1000), breaks = c(seq(0, max.chr9/1000, 150), (max.chr9/1000)),
                     minor_breaks = c(seq(0, max.chr9/1000, 50), (max.chr9/1000)),
                     guide = "axis_minor") +
  
  theme_classic() +
  labs(x = "<br>Chromosome 9 (kb)", y = "Log<sub>2</sub> (Ratio)<br>") +
  theme(axis.title.y = element_markdown(size = 12, color = "black", face = "bold"),
        axis.title.x = element_markdown(size = 12, color = "black", face = "bold"),
        axis.line.x = element_blank(),
        axis.ticks.length.x = unit(0.2, "cm"),
        axis.ticks.length.y = unit(0.2, "cm"),
        axis.text = element_text(size = 11, color = "black")); plot.A
  



# plot.A <- df.2plot %>% 
#   mutate(log2Range = case_when(
#     log2ratio < -0.5 ~ "-1.0 - -0.5", # gray
#     log2ratio < 0 ~ "-0.5 - 0.0", # gray
#     log2ratio < 0.5 ~ "0.0 - 0.5", # gray
#     log2ratio < 1 ~ "0.5 - 1.0", # gray
#     log2ratio < 1.3 ~ "1.0 - 1.5",
#     log2ratio > 1.3 ~ "1.5 - 2.0")) %>% 
#   ggplot() +
#   geom_tile(aes(x = (window.end/1000), y = factor(sample, levels = c("2021EE37", "2021EE18")  ), 
#                 fill = factor(log2Range, levels = c(
#                                                     "1.5 - 2.0", 
#                                                     "1.0 - 1.5", 
#                                                     "0.5 - 1.0", 
#                                                     "0.0 - 0.5", 
#                                                     "-0.5 - 0.0",
#                                                     "-1.0 - -0.5")))) +
#   geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5), color = "white", linewidth = 2)+
#   scale_y_discrete(labels = c(
#                               "2021EE18" = "UmH<sub>2</sub>O<sub>2</sub>.R-Col",
#                               "2021EE37" = "UmH<sub>2</sub>O<sub>2</sub>.R-Col<br>+ 7 days")) +
#   
#   scale_x_continuous(limits = c(0, max.chr9/1000),
#                      breaks = c(seq(0, max.chr9/1000, 150), 
#                                 (max.chr9/1000))) +
#   theme_classic()+
#   labs(x = "\nChromosome 9 (kb)", y = "Colony\n") +
#   scale_fill_manual(values = c(
#     "#011E8F",
#     #"darkblue",
#     "#345AE7", 
#     #"steelblue",
#     "#94ABD8", 
#     #"lightblue",
#     "gray86" ,
#     "gray80", 
#     "#D49A7B"
#   ))+
#   theme(
#     # titles
#     plot.title = element_blank(),
#     #axis.title.y = element_text(size = 12, color = "black"),
#     axis.title.x = element_text(size = 12, color = "black", face = "bold"),
#     axis.title.y = element_blank(),
#     # axis text
#     axis.text.x = element_text(size = 11, color = "black", hjust = 0.5),
#     axis.text.y = element_markdown(size = 11, color = "black"),
#     # axis ticks
#     axis.ticks.y = element_blank(),
#     axis.ticks.length.x = unit(0.2, "cm"),
#     #axis.ticks.x = element_line(size = 2),
#     
#     #axis.text = element_blank(),
#     #axis.text.y = element_blank(),
#     axis.line = element_blank(),
#     # strip
#     strip.background = element_blank(),
#     strip.text = element_text(size = 12, color = "black"),
#     
#     
#     # lenged aesthetics
#     legend.position = "right",
#     # legend.direction = "vertical",
#     # legend.justification = "right",
#     legend.key.size = unit(0.8, "lines"), # Reduce the legend key size
#     legend.key.width = unit(0.6, "lines"), # Reduce key width (if needed)
#     legend.key.height = unit(0.6, "lines"), # Reduce key height (if needed)
#     legend.spacing.x = unit(0.3, "lines"), # Reduce horizontal spacing
#     legend.spacing.y = unit(0.3, "lines")  # Reduce vertical spacing
#     
#   ) +
#   guides(fill = guide_legend(title = expression(bold("log"["2"]*" (Ratio)")),
#                              barwidth = 0.5, ticks = T))
# 
# 
# plot.A 












## Plot B
## Analysis at bp
## Create data frame with coordinates of genes located between the 140 kb to the 154 kb in chr 9
df.genes <- data.frame(
  Gene = c( "UMAG_10438", "UMAG_03439", "UMAG_03440", "tr-HobS", "tr-HobS", "UMAG_10439", "UMAG_03442"),
  xstart = c(140.478, 143.381, 144.651, 147.099, 148.833, 149.932, 151.640),
  xend = c(142.766, 144.134, 145.358, 147.497, 149.295, 150.940, 153.355),
  ystart = rep(5.2, 7),
  yend = rep (5.5, 7))

df.genes <- df.genes %>% mutate(TextPositionX = (xstart + xend)/2, TextPositionY = ((ystart + yend) + 1 )/2)

df.genes$Gene <- gsub("_", "\n", df.genes$Gene)

## Read files with normalized coverage at bp resolution
samples.2read <- c("2021EE01","2021EE18")


for (sample in samples.2read){
  file.2read <- paste0(sample, ".bpNormalizedCoverage.Chr9.txt.gz")
  
  df.name <- paste0("df.", sample)
  df.tmp <- fread(file.2read)
  df.tmp <- df.tmp[Position >= 140000 & Position <=155000]
  
  assign(df.name, df.tmp)
  
  rm(df.name, df.tmp, file.2read)
}

## obtaining log2
df.log2 <- bind_rows(df.2021EE01,  df.2021EE18)

#df.log2

df.log2 <- df.log2 %>%
  select(!DepthCoverage) %>% 
  pivot_wider(names_from = Sample, values_from = normalized.coverage) %>% 
  mutate(across(c(`2021EE01`:`2021EE18`),
                .fns = ~./`2021EE01`,
                .names = "{.col}_Ratio")) %>% 
  select(!c(`2021EE01`:`2021EE18`)) %>% 
  pivot_longer(cols = `2021EE01_Ratio`:`2021EE18_Ratio`,
               names_to = "Sample", values_to = "Log2Ratio") %>% 
  filter(Sample != "2021EE01_Ratio") %>% 
  mutate(Log2Ratio = log2(Log2Ratio))

# plotting 150 kb breakpoint




# plog.log2.bp.cov <- ggplot(data = df.log2)+
plot.log2.bp.cov <- ggplot()+
  geom_line(data = df.log2, aes(x = Position/1000, y = Log2Ratio, color = Sample), linewidth = .25)+
  scale_color_manual(values = c("darkred")) +
  #geom_hline(yintercept = 0, color = "black", linewidth = 0.4)+
  geom_hline(yintercept = -1, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_hline(yintercept = 2, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_hline(yintercept = 3, linetype = "dashed", color = "gray", linewidth = 0.4) +
  #geom_hline(yintercept = 4, linetype = "dashed", color = "gray", linewidth = 0.4) +
  scale_y_continuous(limits = c(-1, 3), breaks = seq(-1, 3, 1),
                     minor_breaks = seq(-1, 3, 0.5),
                     guide = "axis_minor") +
  scale_x_continuous(breaks = c(140, 142.5, 145, 147.5, 150, 152.5, 155), 
                     minor_breaks = seq(140, 155, 1.25),
                     guide = "axis_minor")+
  geom_rect(data = df.genes, aes(xmin = xstart, xmax = xend, 
                                 ymin = ystart - 2.9, ymax = yend - 2.9),
            fill = c("gray", "gray", "gray", "darkblue", "darkblue", "gray", "gray"), alpha = 0.5, color = "black") +
  theme_classic() +
  geom_text(data = df.genes, aes(x = TextPositionX, label = Gene, y = 2.8),
            size = 3) +
  # labs(y = expression("log"["2"]*" (Ratio)"), x = "\nChromosome 9 (kb)") +
  labs(x = "<br>Chromosome 9 (kb)", y = "Log<sub>2</sub> (Ratio)<br>") +
  theme_classic() +
  geom_rect(aes(xmin = 147.099, xmax = 147.498, ymin = -1, ymax = 2.3), fill = "darkblue", alpha = 0.12)+ # tr-HobS
  geom_rect(aes(xmin = 148.883, xmax = 149.295, ymin = -1, ymax = 2.3), fill = "darkblue", alpha = 0.12 )+ # tr-HobS
  geom_rect(aes(xmin = 148.950, xmax = 149.100, ymin = -1, ymax = 2.3), fill = "darkred", alpha = 0.35 )+ #break point
  #geom_vline(xintercept = 148.950) +
  theme(
    # Hide legend box
    legend.position = "none",
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    
    # axis titles
    axis.title.x = element_markdown(face = "bold", size = 12, color = "black"),
    axis.title.y = element_markdown(face = "bold", size = 12, color = "black") ,
    
    # hide x-axis line
    axis.line.x = element_blank(),
    axis.ticks.length.x = unit(0.2, "cm"),
    axis.ticks.length.y = unit(0.2, "cm")
  ); plot.log2.bp.cov 



 
# 


figure.sup.1 <- plot_grid( plot.A, plot.log2.bp.cov, labels = c("A)", "B)"),
                           scale = 0.95, nrow = 2, align = "v")
figure.sup.1


Figure3.Paper <- paste0("Supp.Figure.1.png")

Cairo(file = Figure3.Paper, type = "png", width = 10, height = 8, units = "in", dpi = 300, bg = "white")
print(figure.sup.1)
dev.off()
getwd()

rm(list = ls())












