# created by j cuamatzi
# script for figure 2

# load libraries
libraries <- c("data.table", "readxl","dplyr", "ggplot2", "tidyr", "ggpubr", 
               "cowplot", "scales", "rstatix", "ggtext", "ggh4x")

for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}

# clean env
rm(list = ls())

# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Read file with data of gene expression
df.qpcr <- read_excel(path = "Figure2/Data_2_UMAG_11067_qPCR.xlsx", sheet = "Fold_Change") # read the file
df.qpcr <- setDT(df.qpcr) # transform to data table


stat.test <- df.qpcr %>% 
  filter(Strain != "T20.LB.1") %>% 
  group_by(Condition) %>% 
  t_test(FoldChange ~ Strain) %>% 
  add_xy_position(x = "Strain")

stat.test <- stat.test[stat.test$Condition != "0 mM",]

stat.test$label <- sprintf("%.3f", stat.test$p)

# Estimate the mean of the fold change

df.mean <- df.qpcr %>% group_by(Strain, Condition) %>% 
  summarise(Mean.FC = mean(FoldChange)) %>% 
  ungroup() %>% setDT()

df.mean <- df.mean[Condition != "0 mM"]

df.qpcr <- df.qpcr[Condition != "0 mM"]
df.qpcr <- df.qpcr %>% filter(Strain != "T20.LB.1")



df.mean <- df.mean %>% mutate(Chr9.LA = case_when(
  startsWith(df.mean$Strain, "SG200") ~ "1X",
  # startsWith(df.mean$Strain, "T20.LB.1") ~ "2X",
  startsWith(df.mean$Strain, "T20.LC.1") ~ "3X"
))


df.qpcr <- df.qpcr %>% mutate(Chr9.LA = case_when(
  startsWith(df.qpcr$Strain, "SG200") ~ "1X",
  #startsWith(df.qpcr$Strain, "T20.LB.1") ~ "2X",
  startsWith(df.qpcr$Strain, "T20.LC.1") ~ "3X"
))


df.mean <- df.mean %>% filter(Strain != "T20.LB.1")


Figure.2 <- ggplot()+
  geom_col(data = df.mean, 
           aes(x = Strain, y = Mean.FC, 
               fill = Chr9.LA, color = Chr9.LA),
           alpha = 0.5, linewidth = 0.15) +
  geom_point(data = df.qpcr, aes(x = Strain, y = FoldChange), 
             size = 2.5, alpha = 0.9, color = "black")+
  theme_classic()+
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1),
                     minor_breaks = seq(0, 4, 0.25),
                     guide = "axis_minor",
                     expand = c(0,0))+
  scale_x_discrete(labels = c("SG200" = '<i>U. maydis</i> SG200<sub><span style="color: white;">2</span></sub><br> <br>Initial Strain', 
                              "T20.LC.1" = "UmH<sub>2</sub>O<sub>2</sub>-R <br> <br>Adapted Strain"))+
  labs(x = "Strain", y = "Fold Change in Catalase Gene Expression")+
  ## colors
  scale_fill_manual(values = c( 
    "white", #1X
    
    "black"  # 3X
  )) +
  scale_color_manual(values = c( 
    "black", #1X
    
    "black"  # 3X
  )) +  
  
  ##
  theme(
    # legend
    legend.position = "none",
    
    # axis title
    
    axis.title.x = element_blank(),
        axis.title.y = element_text(margin = unit(c(0,5,0,0), "mm"), face = "bold"),
        
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        
    axis.text.x = element_markdown(color = "black", size = 9),
        
        
        
        axis.text.y = element_text(size = 9, color = "black"),
    ggh4x.axis.ticks.length.minor = rel(1)) +
  stat_pvalue_manual(data = stat.test, label = "label"); Figure.2

plot.Figure2 <- paste0("Figure2/Figure_2.tiff")

ggsave(filename = plot.Figure2, plot = Figure.2,
       width = 7, height = 9, units = "cm", dpi = 300, bg = "white", device = "tiff")


rm(list = ls())


