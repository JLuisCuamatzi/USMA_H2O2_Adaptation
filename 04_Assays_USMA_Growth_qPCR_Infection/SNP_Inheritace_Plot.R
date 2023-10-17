library(ComplexHeatmap)
library(Cairo)


SNP.list <- list(
  "T20.LC.1.A" = c("USMA_521_v2_18-330272-C-A-1", "USMA_521_v2_8-464088-T-A-1", 
                 "USMA_521_v2_24-15246-G-A-1", "USMA_521_v2_1-2430690-G-C-1",
                 "USMA_521_v2_10-623352-T-A-1"),
  "T20.LC.1.B" = c("USMA_521_v2_18-330272-C-A-1", "USMA_521_v2_8-464088-T-A-1", 
                   "USMA_521_v2_24-15246-G-A-1", "USMA_521_v2_1-2430690-G-C-1"),
  "T20.LC.1.C" = c("USMA_521_v2_18-330272-C-A-1", "USMA_521_v2_8-464088-T-A-1", 
                   "USMA_521_v2_24-15246-G-A-1", "USMA_521_v2_1-2430690-G-C-1")
) 
SNP.list

SNP.list.intersection <- make_comb_mat(SNP.list)


# export
# plot and export
dir.create("../NewManuscript/Plots")

file2export <- paste0("../NewManuscript/Plots/SNP_Inheritance.png") # Create object with filename to save

CairoPNG(file = file2export, width = 3.5, height = 3, units = "in", dpi =300)
UpSet(m = SNP.list.intersection,
      top_annotation = upset_top_annotation(SNP.list.intersection,
                                            annotation_name_rot = 90,
                                            annotation_name_side = "right",
                                            axis_param = list(side = "right"),
                                            add_numbers = T,
                                            numbers_rot = 90),
      right_annotation = upset_right_annotation(SNP.list.intersection, add_numbers = T))
dev.off()

rm(list = ls())

