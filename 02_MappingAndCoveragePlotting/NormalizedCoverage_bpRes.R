## author: jcuamatzi

libraries <- c("data.table", "ggplot2", "tidyr", "cowplot", "dplyr")

for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}

rm(list = ls())

# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


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
  file.2read <- paste0("normalizedCoverageTables/bpResolution/", sample, ".bpNormalizedCoverage.Chr9.txt.gz")
  
  df.name <- paste0("df.", sample)
  df.tmp <- fread(file.2read)
  df.tmp <- df.tmp[Position >= 140000 & Position <=155000]
  
  assign(df.name, df.tmp)
  
  rm(df.name, df.tmp, file.2read)
}

## obtaining log2


df.log2 <- bind_rows(df.2021EE01,  df.2021EE18)

df.log2 <- df.log2 %>%
  select(!DepthCoverage) %>% 
  pivot_wider(names_from = Sample, values_from = normalized.coverage) %>% 
  mutate(across(c(`2021EE01`:`2021EE18`),
                .fns = ~./`2021EE01`,
                .names = "{.col}_Ratio")) %>% 
  select(!c(`2021EE01`:`2021EE18`)) %>% 
  pivot_longer(cols = `2021EE01_Ratio`:`2021EE18_Ratio`,
               names_to = "Sample", values_to = "Log2Ratio") %>% 
  filter(Sample != "2021EE01_Ratio")

# plotting 150 kb breakpoint


# plog.log2.bp.cov <- ggplot(data = df.log2)+
plot.log2.bp.cov <- ggplot()+
  geom_line(data = df.log2, aes(x = Position/1000, y = Log2Ratio, color = Sample), linewidth = .25)+
  scale_color_manual(values = c("darkred")) +
  #geom_hline(yintercept = 0, color = "black", linewidth = 0.4)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_hline(yintercept = 2, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_hline(yintercept = 3, linetype = "dashed", color = "gray", linewidth = 0.4) +
  geom_hline(yintercept = 4, linetype = "dashed", color = "gray", linewidth = 0.4) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  scale_x_continuous(breaks = c(140, 142.5, 145, 147.5, 150, 152.5, 155))+
  geom_rect(data = df.genes, aes(xmin = xstart, xmax = xend, 
                                 ymin = ystart - 0.9, ymax = yend - 0.9),
            fill = c("gray", "gray", "gray", "darkblue", "darkblue", "gray", "gray"), alpha = 0.5, color = "black") +
  theme_classic() +
  geom_text(data = df.genes, aes(x = TextPositionX, label = Gene, y = 4.9),
            size = 5) +
  # labs(y = expression("log"["2"]*" (Ratio)"), x = "\nChromosome 9 (kb)") +
  labs(y = expression(bold( "Log"["2"]*" (Ratio)") ), x = "\nChromosome 9 (kb)") +
  theme_classic() +
  geom_rect(aes(xmin = 147.099, xmax = 147.498, ymin = 0, ymax = 4.3), fill = "darkblue", alpha = 0.12)+ # tr-HobS
  geom_rect(aes(xmin = 148.883, xmax = 149.295, ymin = 0, ymax = 4.3), fill = "darkblue", alpha = 0.12 )+ # tr-HobS
  geom_rect(aes(xmin = 148.950, xmax = 149.100, ymin = 0, ymax = 4.3), fill = "darkred", alpha = 0.35 )+ #break point
  #geom_vline(xintercept = 148.950) +
  theme(
    # Hide legend box
    legend.position = "none",
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    
    # axis titles
    axis.title.x = element_text(face = "bold", size = 9, color = "black"),
    axis.title.y = element_text(face = "bold", size = 9, color = "black") ,
    
    # hide x-axis line
    #axis.line.x = element_blank()
    ); plot.log2.bp.cov 


# plot.chr9.breakpoint.log2 <- plot_grid(plot.Chr9.GC, 
#                                        plot.log2.bp.cov, 
#                                        nrow = 2, align = "hv", rel_heights = c(0.4, 1), labels = c("A)", "B)"))  
# 
# plot.chr9.breakpoint.log2
# getwd()
# ggsave(filename = "../../NewManuscript/Plots/Supp.Fig1.png", 
#        plot = plot.log2.bp.cov ,
#        width = 16, height = 11, units = "cm", dpi = 300)  

#### Plotting from 30 to 45 kb
for (sample in samples.2read){
  file.2read <- paste0("normalizedCoverageTables/bpResolution/", sample, ".bpNormalizedCoverage.Chr9.txt.gz")
  
  df.name <- paste0("df.", sample)
  df.tmp <- fread(file.2read)
  df.tmp <- df.tmp[Position >= 30000 & Position <=45000]
  
  assign(df.name, df.tmp)
  
  rm(df.name, df.tmp, file.2read)
}

# Obtaining log2
df.log2 <- bind_rows(df.2021EE01, df.2021EE18)

df.log2 <- df.log2 %>%
  select(!DepthCoverage) %>% 
  pivot_wider(names_from = Sample, values_from = normalized.coverage) %>% 
  mutate(across(c(`2021EE01`:`2021EE18`),
                .fns = ~./`2021EE01`,
                .names = "{.col}_Ratio")) %>% 
  select(!c(`2021EE01`:`2021EE18`)) %>% 
  pivot_longer(cols = `2021EE01_Ratio`:`2021EE18_Ratio`,
               names_to = "Sample", values_to = "Log2Ratio") %>% 
  filter(Sample != "2021EE01_Ratio")

#41321,43299,UMAG_03407

## Create data frame with coordinates of genes located between the 140 kb to the 154 kb in chr 9
df.genes.45kb <- data.frame(
  Gene = c( "UMAG_03403", "UMAG_03404", "UMAG_12229", "UMAG_03406", "UMAG_03407"),
  xstart = c(31.391, 36.098, 38.711, 39.995, 41.321),
  xend = c(35.074, 37.963, 39.698, 40.611, 43.299),
  ystart = rep(5.2, 5),
  yend = rep (5.5, 5))

df.genes.45kb <- df.genes.45kb %>% 
  mutate(TextPositionX = (xstart + xend)/2, TextPositionY = ((ystart + yend) + 1 )/2)

df.genes.45kb$Gene <- gsub("_", "\n", df.genes.45kb$Gene)

# GC 
# Plotting GC of chromosome 9
# GC was estimated using non-overlapping windows of 0.2 kb (200 bp)
#df.Chr9.GC <- fread("../USMA_521_v2_9_GC_100bp.txt") # windows of 50 bp
df.Chr9.GC2plot.45kb <- df.Chr9.GC %>% filter(End >= as.numeric(30000) & End <= as.numeric(45000))  

plot.Chr9.GC.45kb <- ggplot() +
  geom_line(data= df.Chr9.GC2plot.45kb, aes(x = (End/1000), y = GC ), linewidth = 0.5) +
  geom_hline(yintercept = 54.04, linetype = "dotted", color = "gray", linewidth = .9) + # mean gc for umaydis = 54.04%
  theme_classic()+
  labs(y = "\n\nGC (%)") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_text(face = "bold", color = "black", size = 8),
    axis.text.y = element_text(size = 8, color = "black")
  ) +
  scale_y_continuous(limits = c(20, 80)); plot.Chr9.GC.45kb



plot.Chr9.45kb <- ggplot()+
  geom_line(data = df.log2, aes(x = Position/1000, y = Log2Ratio, color = Sample), linewidth = .25)+
  scale_color_manual(values = c("darkgreen", "orange", "darkred")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_hline(yintercept = 2, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_hline(yintercept = 3, linetype = "dashed", color = "gray", linewidth = 0.4) +
  geom_hline(yintercept = 4, linetype = "dashed", color = "gray", linewidth = 0.4) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  scale_x_continuous(breaks = c(30, 32.5, 35, 37.5, 40, 42.5, 45)) +
  geom_rect(data = df.genes.45kb, aes(xmin = xstart, xmax = xend, 
                                  ymin = ystart - 0.9, ymax = yend - 0.9),
             fill = c("gray", "gray", "gray", "gray", "gray"), alpha = 0.5, color = "black") +
  theme_classic() +
  geom_text(data = df.genes.45kb, aes(x = TextPositionX, label = Gene, y = 4.9), size = 2.2) +
  labs(y = "\n\nLog2 - Ratio", x = "\nChromosome 9 (kb)") +
  geom_rect(aes(xmin = 38.850, xmax = 38.870, ymin = 0, ymax = 4.3), fill = "darkred" )+ #break point
  #geom_vline(xintercept = 148.950) +
  theme(
    # Hide legend box
    legend.position = "none",
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    
    # axis titles
    axis.title.x = element_text(face = "bold", size = 9, color = "black"),
    axis.title.y = element_text(face = "bold", size = 9, color = "black") ,
    
    # hide x-axis line
    axis.line.x = element_blank()
  )


plot.chr9.45kb.log2 <- plot_grid(plot.Chr9.GC.45kb, plot.Chr9.45kb, nrow = 2, align = "hv", rel_heights = c(0.4, 1), labels = c("A)", "B)"))  


getwd()
ggsave(filename = "normalizedCoverageTables/CoveragePlots/Supp.Figure.NormCov.bpResolution.log2.30to45kb.png", 
       plot = plot.chr9.45kb.log2,
       width = 16, height = 11, units = "cm", dpi = 300)  

rm(list = ls())
  
  
  
  
  
  
  
  
  
  
  
  
  
