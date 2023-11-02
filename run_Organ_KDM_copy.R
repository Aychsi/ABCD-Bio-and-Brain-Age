
## This file runs the KDM_BBB_5.0_V1_Organs.qmd file to calculate the KDM for each organ

source("/Users/hansoochang/Drexel/ABCD/R Codes/BioAge_FeatureSelection_5.0_V2.Rmd")

library(BioAge)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(Metrics)
library(stats)
library(skimr)
library(pcaPP)
library(stringr)
library(stats)
library(misty)
library(emmeans)
library(jph)


### Cardiovascular ### 
jph::quarto_render_move( 
  input = "/Users/hansoochang/Drexel/ABCD/KDM_BBB_5.0_V1_Organs.qmd",
  output_file = paste0("cv", ".html"),
  execute_params = list(
    df = jsonlite::toJSON(cv),
    df_name = "cv"
  )
)

### Organ Function ### 
jph::quarto_render_move( 
  input = "/Users/hansoochang/Drexel/ABCD/KDM_BBB_5.0_V1_Organs.qmd",
  output_file = paste0("of", ".html"),
  execute_params = list(
    df = jsonlite::toJSON(of %>% select(-biospec_blood_imm_gran_abs)),
    df_name = "of"
  )
)



### Muskuloskeletal ### 
jph::quarto_render_move( 
  input = "/Users/hansoochang/Drexel/ABCD/KDM_BBB_5.0_V1_Organs.qmd",
  output_file = paste0("ms", ".html"),
  execute_params = list(
    df = jsonlite::toJSON(ms),
    df_name = "ms"
  )
)

### Immune ### 
jph::quarto_render_move( 
  input = "/Users/hansoochang/Drexel/ABCD/KDM_BBB_5.0_V1_Organs.qmd",
  output_file = paste0("im", ".html"),
  execute_params = list(
    df = jsonlite::toJSON(im),
    df_name = "im"
  )
)


### Metabolic ### 
jph::quarto_render_move( 
  input = "/Users/hansoochang/Drexel/ABCD/KDM_BBB_5.0_V1_Organs.qmd",
  output_file = paste0("mb", ".html"),
  execute_params = list(
    df = jsonlite::toJSON(mb),
    df_name = "mb"
  )
)








