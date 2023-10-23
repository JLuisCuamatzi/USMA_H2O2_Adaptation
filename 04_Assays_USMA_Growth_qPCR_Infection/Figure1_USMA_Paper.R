# author: jcuamatzi
# date: 2023-10-04
# Figure 1

# load libraries
libraries <- c("data.table", "dplyr", "tidyr", "ggplot2", "scales",
               "cowplot", "stringr", "Cairo", "grid", "rstatix", 
               "ggthemes", "ggpubr", "ggtext", "ggh4x")

for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
} 




rm(list = ls())

# Set as working directory, the directory in which this script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Read data of CFU
df <- fread("Umaydis_ExpEvol_GrowthData.csv")

df <- df %>% mutate(
  X.Labs = case_when(
    endsWith(df$Line, "SG200") ~ "U. maydis SG200",
    endsWith(df$Line, "C") ~ "U. maydis SG200 H2O2 resistant (UmH2O2.R)"
  )
)
# keep only CFU
df <- df[,-5]

# here we get only data for SG200 and treatment 20
df.Figure.1 <- df %>%  
  filter(TreatmentNum %in% c("SG200", "Treatment 20")) %>% filter(H2O2 %in% c("0 mM", "5 mM", "60 mM",
                                                                                   "w/o ROS 0 mM",
                                                                                   "w/o ROS 5 mM",
                                                                                   "w/o ROS 60 mM"))
# remove extrac charactets
df.Figure.1$H2O2 <- gsub("w/o ROS ", "", df.Figure.1$H2O2)

# keep only 0 mM, 5 mM and 60 mM
df.Figure.1 <- df.Figure.1 %>% 
  filter(H2O2 %in% c("0 mM", "5 mM", "60 mM"))

df.Figure.1 <- df.Figure.1 %>% mutate(CFU_Type =
                                        if_else(H2O2 != "0 mM", "CFU_Target", "CFU_Reference"))

df.Figure.1$H2O2 <- gsub(" ", "", df.Figure.1$H2O2)

# Pivot the table 
df.Figure.1 <- df.Figure.1 %>% 
  group_by(Condition, TreatmentNum, Line) %>% 
  summarise(CFU_0mM = CFU[H2O2 == "0mM"],
            CFU_5mM = CFU[H2O2 == "5mM"],
            CFU_60mM = CFU[H2O2 == "60mM"]) %>% 
  ungroup()

# estimate the percentage of cfu 
# we calculate 5 and 60 mM against 0 mM in each group 
df.Figure.1 <- df.Figure.1 %>% 
  mutate(across(c(`CFU_0mM`:`CFU_60mM`),
                .fns = ~./`CFU_0mM`,
                .names = "{.col}_Ratio"))

df.Figure.1 <- df.Figure.1 %>% 
  select(Condition, TreatmentNum, Line, CFU_0mM_Ratio: CFU_60mM_Ratio) %>% 
  pivot_longer(cols = CFU_0mM_Ratio: CFU_60mM_Ratio, names_to = "H2O2", values_to = "CFU_Ratio") %>% 
  setDT()

df.Figure.1$H2O2 <- gsub("CFU_", "", df.Figure.1$H2O2)
df.Figure.1$H2O2 <- gsub("_Ratio", "", df.Figure.1$H2O2)



# Keep just 5 mM and 60 mM
df.Figure.1 <- df.Figure.1[H2O2 %in% c("5mM", "60mM")]

df.Figure.1.Mean <- df.Figure.1 %>% 
  group_by(Condition, TreatmentNum, H2O2, Line) %>% 
  summarise(CFU_Ratio_Mean = mean(CFU_Ratio),
            CFU_Ratip_DS = sd(CFU_Ratio)) %>% 
  ungroup()

df.Figure.1 %>% 
  group_by(Condition, TreatmentNum, H2O2, Line) %>% 
  summarise(CFU_Ratio_Mean = mean(CFU_Ratio),
            CFU_Ratip_DS = sd(CFU_Ratio)) %>% 
  ungroup()


df2plot.1 <- df.Figure.1.Mean %>% 
  select(Condition, H2O2, Line, CFU_Ratio_Mean)

df2plot.1$Line <- factor(df2plot.1$Line, levels = c("SG200", "A", "B", "C", "D", "E", "F"))

df2plot.1 <- df2plot.1 %>% 
  mutate(Strip.Labs = 
           case_when(
             startsWith(as.character(df2plot.1$Line), "SG200") ~ "Initial\nStrain",
             startsWith(as.character(df2plot.1$Line), "A") ~ "Exposed Pools\n",
             startsWith(as.character(df2plot.1$Line), "B") ~ "Exposed Pools\n",
             startsWith(as.character(df2plot.1$Line), "C") ~ "Exposed Pools\n",
             startsWith(as.character(df2plot.1$Line), "D") ~ "Unexposed Pools\n",
             startsWith(as.character(df2plot.1$Line), "E") ~ "Unexposed Pools\n",
             startsWith(as.character(df2plot.1$Line), "F") ~ "Unexposed Pools\n"
  ))

df2plot.1 <- df2plot.1 %>% 
  mutate(
         X.Labs = if_else(Line == "SG200", "", as.vector(Line)))

df2plot.1$Strip.Labs <- factor(as.vector(df2plot.1$Strip.Labs), 
                             levels = c("Initial\nStrain", "Exposed Pools\n", "Unexposed Pools\n"))


df2plot.1$H2O2 <- gsub("([0-9]+)(mM)", "\\1 \\2", df2plot.1$H2O2)

df2plot.1 <- df2plot.1 %>% mutate(x.labels = case_when(
  endsWith(as.vector(df2plot.1$Line), "C") ~  "U. maydis SG200 H2O2 Resistant (UmH2O2.R)",
  endsWith(as.vector(df2plot.1$Line), "SG200") ~ "U. maydis SG200"))


plot.Figure.1 <- df2plot.1 %>% 
  filter(Line %in% c("C", "SG200") ) %>% 
  ggplot() +
  geom_col(aes(x = x.labels, y = CFU_Ratio_Mean, fill = H2O2),color = "black",
           position = position_dodge(width = 0.92), alpha = 0.7) +

  scale_y_continuous(labels = percent, limits = c(0,1), 
                     breaks = c(0, 0.20, 0.40, 0.60, 0.8, 1.0), 
                     minor_breaks = seq(0, 1.0, by = 0.1),
                     expand = c(0,0),
                     guide = "axis_minor") +
  # scale_x_discrete(labels = label_wrap(20))+
  # scale_x_discrete(labels = c("U. maydis SG200" = expression(paste(italic("U. maydis"), "SG200"))))+
  # scale_x_discrete(labels = c("U. maydis SG200" = expression( atop( italic("U. maydis"), SG200)),
  #                             "U. maydis SG200 H2O2 Resistant (UmH2O2.R)" =
  #                               expression(italic("U. maydis")*" SG200 H"["2"]*"O"["2"]*" Resistant (UmH"["2"]*"O"["2"]*".R)" )))+
  # scale_x_discrete(labels = c("U. maydis SG200" = expression(italic("U. maydis ")* "SG200"),
  #                             "U. maydis SG200 H2O2 Resistant (UmH2O2.R)" =
  #                               expression(atop( italic("U. maydis ")*" SG200", atop("H"["2"]*"O"["2"]*" Resistant", "(UmH"["2"]*"O"["2"]*".R)" ) ) )))+
  # scale_x_discrete(labels = c("U. maydis SG200" = "<i>U. maydis</i> SG200<br>Initial Strain",
  #                             "U. maydis SG200 H2O2 Resistant (UmH2O2.R)" =
  #                               "<i>U. maydis</i> SG200<br> <br>H<sub>2</sub>O<sub>2</sub> Resistant <br> <br>(UmH<sub>2</sub>O<sub>2</sub>.R)"))+
  scale_x_discrete(labels = c("U. maydis SG200" = "<i>U. maydis</i> SG200<br> <br>Initial Strain",
                              "U. maydis SG200 H2O2 Resistant (UmH2O2.R)" =
                                "UmH<sub>2</sub>O<sub>2</sub>.R<br> <br>Adapted Strain"))+

  theme_classic() +
  
  scale_fill_manual(values = c("gray60", "darkred"))+ 
  scale_color_manual(values = c("gray60", "darkred"))+ 
  labs(y = "Surviving Cells") +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    legend.position = c(0.2, 0.9),
    
    legend.key.size = unit(0.4, "lines"),   # Adjust the size of the legend key
    legend.key.height = unit(0.4, "lines"), # Adjust the height of the legend key
    legend.key.width = unit(0.6, "lines") ,  # Adjust the width of the legend key
    
    
    strip.background = element_blank(),
    
    # axis lines
    axis.line.x = element_blank(),
    #axis.line.y = element_blank(),
    
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_blank(),
    #axis.text.x = element_text(color = "black", vjust = 34.5),
    #axis.text.x = element_blank(),
    
    axis.ticks.x = element_blank(),
    #strip.text = element_text(vjust = 2, size = 11)
    strip.text = element_blank(),
    
    axis.text.y = element_text(color = "black", size = 9),
    # axis.text.x = element_text(color = "black", size = 9),
    axis.text.x = element_markdown(color = "black", size = 9),
    
    
    ggh4x.axis.ticks.length.minor = rel(1)
    ) +
  guides(fill = guide_legend(title = expression(bold( "[H"["2"]*"O"["2"]*"] Shock") ) ),
         color = "none"); plot.Figure.1







#### Panel B


## Stability of H2O2 resistance in T20.LC.1
df.halos <- fread("20230930_Halos_35GEN_H2O2_NO_H2O2.csv")

df.halos <- df.halos %>% 
  filter(Time != "200G+35G H2O2") 

df.halos <- df.halos %>% 
  mutate(x.labs = case_when(
    endsWith(df.halos$Time, "200G") ~ "UmH2O2R.Col",
    endsWith(df.halos$Time, "Initial") ~ "U. maydis SG200",
    endsWith(df.halos$Time, "35G") ~ "UmH2O2R.Col+7d"
  ) )



df.mean.InhZone <- df.halos %>% mutate(InhibitionZone = (pi* (Halo_cm/2)^2 ) ) %>% 
  group_by(x.labs) %>% 
  reframe(MeanInhibitonZone = mean(InhibitionZone)) 
  


# Statistical test


stat.test <- df.halos  %>% 
  filter(!Time %in% c("Initial" ) )%>% 
  t_test(Halo_cm ~ x.labs, paired = T) %>% add_xy_position(x = "x.labs")


stat.test$label <- sprintf("%.3f", stat.test$p)



stat.test.2 <- df.halos %>% 
  t_test(Halo_cm ~ Time) %>% add_xy_position(x = "Time")

stat.test.2 <- stat.test.2[-1,]

df.halos <- df.halos %>% filter(Time != "200G+35G H2O2")


stat.test$xmin <- stat.test$xmin +1
stat.test$xmax <- stat.test$xmax +1
stat.test.2[2,13] <- 1
stat.test.2[1,14] <- 2
stat.test.2$label <- stat.test.2$p



plot.T20LC1 <- ggplot()+
  geom_col(data = df.mean.InhZone, 
           aes(x = x.labs, y = MeanInhibitonZone,fill = x.labs, color = x.labs), 
           alpha = 0.5)+
  geom_point(data = df.halos , aes(x = x.labs, y =  pi*(Halo_cm/2)^2, color = x.labs ), size = 2.5, alpha = 0.9) +
  scale_y_continuous(limits = c(0, 14), 
                     breaks = seq(0, 14, 2),
                     expand = c(0,0),
                     minor_breaks = seq(0,14, 1),
                     guide = "axis_minor")+
  scale_x_discrete(labels = c(
    #"U. maydis SG200" = expression(italic("U. maydis")*" SG200"),
    #"U. maydis SG200" = "<i>U. maydis</i> SG200",
    "U. maydis SG200" = "Initial Strain",
    #"UmH2O2R.Col" = expression( "UmH"["2"]*"O"["2"]*".R - Col") ,
    #"UmH2O2R.Col" = "UmH<sub>2</sub>O<sub>2</sub>.R-Col",
    "UmH2O2R.Col" = "UmH<sub>2</sub>O<sub>2</sub>.R<br> <br>Treatment 20",
    # "UmH2O2R.Col+7d" = expression(atop("UmH"["2"]*"O"["2"]*".R - Col", "+ 7 days"))
    "UmH2O2R.Col+7d" = "UmH<sub>2</sub>O<sub>2</sub>.R<br> <br>After 7 days <br> <br>w/o H<sub>2</sub>O<sub>2</sub> exposure"  ) )+
  labs(
    x = "\nTimepoints", 
    y = expression( bold("Inhibition Zone by H"["2"]*"O"["2"]*" (cm"^"2"*")")) )+
  theme_classic( )+
  scale_fill_manual(values = rep("black", 6) )+
  scale_color_manual(values = rep("black", 6) )+
  theme(
    legend.position = "none",
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    #axis.text.x = element_text(color = "black", size = 9),
    axis.text.x = element_markdown(color = "black", size = 9),
    axis.title.x = element_text(color = "white"),
    
    
    #axis.text.y = element_text(color = "black", size = 9),
    axis.text.y = element_markdown(color = "black", size = 9),
    
    ggh4x.axis.ticks.length.minor = rel(1)
  ) +
  stat_pvalue_manual(data = stat.test.2, y.position = c(11.5, 12.5) )+ 
  stat_pvalue_manual(data = stat.test, y.position = c(10.5) );plot.T20LC1


plot.F1 <- plot_grid(plot.Figure.1, plot.T20LC1, labels = c("A)", "B)"),align = "h", scale = 0.95)

plot.F1

ggsave(filename = "Fig.1.png", plot = plot.F1, units = "in",
      width = 8, height = 3.8, dpi = 300, bg = "white")



rm(list = ls())
