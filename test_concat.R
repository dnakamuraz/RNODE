source("R/concatenate.R")
setwd("~/Desktop/Doutorado/Project/B2_TEvsMOL/Other/ALL_DATA/063_Neumann2021/tnt")
concatenate(input1="063_MOL_data.tnt", input2="063_MORPH_data.tnt", output_format="tnt", output_file="063_TE_data_test", gaps=F)
