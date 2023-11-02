
### File with Outcome Functions ###



## Read Data and Prepare ##
# Input: string: txt file of outcome data, string: which time point you wan
# Output: dataframe with selected parameters

read_data <- function(file, time) {
  df <- read.delim(file)
  df <- df[-1,]
  df[df == ""] <- NA
  
  df_time <- df %>% filter(eventname == time)
  
  if (nrow(df_time) == 0) {
    print("Time point not in data")
  }
  
  print(table(df$eventname))
  
  return(df_time)
}


## Histogram and Scatterplot ##
# Input: df: dataframe; outcome: string of outcome name in df
# Output: histogram and scatterplot

hist_scat <- function(df, outcome) {
  
  print(ggplot(df, aes_string(x=outcome)) + geom_histogram())
  
  print(
    ggscatter(df, x = "kdm_advance", y = outcome,
              add = "reg.line", conf.int = TRUE,
              xlab = "kdm_advance", ylab = outcome, color = "gender",
              palette = c(M = "blue", F = "pink")) +
      stat_cor(aes(color = gender)) + stat_cor(label.y.npc = 0.8)
  )
}

## Check Zero-Inflation ##
# Input: df: dataframe; outcome: string of outcome name in df
# Output: list output_list where $zero_name is name of outputs that are zero inflated
# and $reg_name is name of outputs that are not zero inflated

check_zinf <- function(df, list_col) {
  fits_reg_name <- list()
  fits_zero_name <- list()
  df$gender <- as.factor(df$gender)
  
  for (col in list_col) { 
    
    # Check if the column contains any non-NA values
    if (!any(!is.na(df[[col]]))) {
      cat("Column", col, "is entirely NA. Skipping...\n")
      next  # Skip this iteration of the loop
    }
    
    df <- df %>% mutate_at(vars(col), ~ . - min(., na.rm = T) )
    
    if (any(grepl("kdm_advance_bioage", colnames(df)))) {
      mod <- glm(data = df, paste(col, "~ kdm_advance_brainage + kdm_advance_bioage"), family = poisson, na.action = na.omit)
    } else {
    
    mod <- glm(data = df, paste(col, "~ kdm_advance*gender"), family = poisson, na.action = na.omit)
    
    }
    # print(summary(sex_lmer)$coefficients)
    if (check_zeroinflation(mod)$ratio < 0.05) {
      print(paste(col, check_zeroinflation(mod)$ratio))
      fits_zero_name <- append(fits_zero_name, col)
    } else {
      print(paste(col, check_zeroinflation(mod)$ratio))
      fits_reg_name <- append(fits_reg_name, col)
    }
  }
  output_list <- list("zero_name" = fits_zero_name, "reg_name" = fits_reg_name)
  return(output_list)
}


## Multilevel Regression (lme4) for biological age ##
#NB: Handles 3 year scores too. Adds 2 year score as a predictor.
# Input: df: dataframe; outcome: string of outcome name in df (at 2 year), exponent: if transformation needed
# Output: lmer object

lmer_bio <- function(df, outcome, exponent) {
  
  if (any(grepl("_3year", colnames(df)))) {
    df_model <- lmer(data = df, as.formula(paste(outcome, "_3year","^", exponent, 
                                                 "~ kdm_advance + sex_dummy + ", 
                                                 outcome, 
                                                 "+ (1 | site_id_l) + (1 | rel_family_id)", sep = "")), REML = F)
  } else {
    df_model <- lmer(data = df, as.formula(paste(outcome, "^", exponent, 
                                                 "~ kdm_advance + sex_dummy + (1 | site_id_l) + (1 | rel_family_id)", sep = "")), REML = F)
  }
  
  
  return(df_model)
}

## Multilevel Binomial Regression (glmer) for biological age ##
#NB: Handles 3 year scores too. Adds 2 year score as a predictor.
# Input: df: dataframe; outcome: string of outcome name in df (at 2 year), exponent: if transformation needed
# Output: lmer object

glmer_bio <- function(df, outcome) {
  
  if (any(grepl("_3year", colnames(df)))) {

    
    df_model <- glmer(data = df, as.formula(paste(outcome, "_3year", 
                                                 "~ kdm_advance + sex_dummy +", 
                                                 outcome, 
                                                 "+ (1  | site_id_l)", sep = "")), 
                      family = binomial(),
                      control=glmerControl(optimizer="bobyqa"))
  } else {
    df_model <- glmer(data = df, as.formula(paste(outcome,
                                                 "~ kdm_advance + sex_dummy + (1 | site_id_l)", sep = "")), 
                      family = binomial(),
                      control=glmerControl(optimizer="bobyqa"))
  }
  
  
  return(df_model)
}

## Multilevel Binomial Regression (glmer) for brain age ##
#NB: Handles 3 year scores too. Adds 2 year score as a predictor.
# Input: df: dataframe; outcome: string of outcome name in df (at 2 year), exponent: if transformation needed
# Output: lmer object

glmer_brainage <- function(df, outcome) {
  
  if (any(grepl("_3year", colnames(df)))) {
    df_model <- glmer(data = df, as.formula(paste(outcome, "_3year", 
                                                  "~ kdm_advance + sex_dummy + ", 
                                                  outcome, 
                                                  " + (1 | mri_info_manufacturer)", sep = "")), 
                      family = binomial(),
                      control=glmerControl(optimizer="bobyqa"))
  } else {
    df_model <- glmer(data = df, as.formula(paste(outcome,
                                                  "~ kdm_advance + sex_dummy  + (1 | mri_info_manufacturer)", sep = "")), 
                      family = binomial(),
                      control=glmerControl(optimizer="bobyqa"))
  }
  
  
  return(df_model)
}


## Multilevel Regression (lme4) for brain age ##
# Input: df: dataframe; outcome: string of outcome name in df, exponent: if transformation needed
# Output: lmer object

lmer_brainage <- function(df, outcome, exponent) {
  
  if (any(grepl("_3year", colnames(df)))) {
    df_model <- lmer(data = df, as.formula(paste(outcome, "_3year","^", exponent, 
                                                 "~ kdm_advance*sex_dummy + ", 
                                                 outcome, 
                                                 "+ (1 + sex_dummy | site_id_l) + (1 | rel_family_id:mri_info_manufacturer)", sep = "")), REML = F)
  } else {
    df_model <- lmer(data = df, as.formula(paste(outcome, "^", exponent, 
                                                 "~ kdm_advance*sex_dummy + (1 + sex_dummy | site_id_l) + (1 | rel_family_id:mri_info_manufacturer)", sep = "")), REML = F)
  }
  
  return(df_model)
}




## Cpglmm Regression for biological age ##
# Input: df: dataframe; outcome: string of outcome name in df
# Output: cpglmm object

cpglmm_bio <- function(df, form) {
  
  df$rel_family_id <- as.factor(df$rel_family_id)
  df$site_id_l <- as.factor(df$site_id_l)
  df$sex_dummy <- as.factor(ifelse(df$gender == "M", 1, 0))
  
  df_model <- cpglmm(link = "log", data = df, formula = form)
  return(df_model)
  
}



## Cpglmm Regression for brainage age ##
# Input: df: dataframe; outcome: string of outcome name in df
# Output: cpglmm object

cpglmm_brainage <- function(df, outcome) { 
  library(cplm)
  df$rel_family_id <- as.factor(df$rel_family_id)
  df$site_id_l <- as.factor(df$site_id_l)
  df$mri_info_manufacturer <- as.factor(df$mri_info_manufacturer)
  df$sex_dummy <- as.factor(ifelse(df$gender == "M", 1, 0))
  
  print(as.formula(paste(outcome, 
                         "~ kdm_advance*sex_dummy + (1 + kdm_advance + sex_dummy| site_id_l) + (1 | rel_family_id:mri_info_manufacturer)")))
  
  df_model <- cpglmm(data = df, formula = as.formula(paste(outcome, "~ kdm_advance*sex_dummy + (1 + kdm_advance + sex_dummy| site_id_l) + (1 | rel_family_id:mri_info_manufacturer)")))
  return(df_model)
  
}


## Make random effects factors and gender dummy variable ##
# Input: df: dataframe
# Output: df with correct datatypes

factor_dummy <- function(df) {
  df$rel_family_id <- as.factor(df$rel_family_id)
  df$site_id_l <- as.factor(df$site_id_l)
  df$sex_dummy <- as.factor(ifelse(df$gender == "M", 1, 0))
  
  if ("mri_info_manufacturer" %in% colnames(df)) {
    df$mri_info_manufacturer <- as.factor(df$mri_info_manufacturer)
  }
  return(df)
  
}



## Make merged 2-3 year dataframes ##
# Input: df2: dataframe at 2 years which is already merged with y2 outcomes; df3: dataframe for 3 year outcomes
# Output: merged df with all info at year 2 and outcomes at year 3 labelled (xxx_3year)

merge_23 <- function(df2, df3) {
  df3 <- df3 %>% rename_at(vars(-matches("subjectkey")), ~paste0(., "_3year"))
  df2_df3 <- merge(df2, df3, by = "subjectkey") 
  
  df2_df3[df2_df3 == ""] <- NA
  df2_df3 <- df2_df3[complete.cases(df2_df3), ]
  print(length(unique(df2_df3$subjectkey))) # 724
  
  return(df2_df3)
}





























