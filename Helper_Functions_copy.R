

# This file contains helper functions for the ABCD Analyses

# Libraries
library(BioAge)
library(ggplot2)
library(ggpubr)
library(Metrics)
library(stats)
library(skimr)
library(magrittr)
library(cplm)
library(parameters)
library(dplyr)
library(lme4)
library(performance)
library(ROSE)
library(tidyr)
library(ggpmisc)


## Add group column to distinguish Bio Age from Brain Age and get only participants with at least two timepoints
### Input 
#### df: dataframe I want to change
#### bio_or_brain: string to say whether it's biological age or brain age
### Output
#### df_filtered
add_col_at_least_two_time <- function (df, bio_or_brain) {
  
  # Check if 'group' column already exists
  if (!"group" %in% colnames(df)) {
    df$group <- bio_or_brain
  }
  
  tryCatch({
    df_filtered <- df %>%
      group_by(src_subject_id) %>%
      filter(n_distinct(eventname) >= 2) %>%
      ungroup()
    
    total_participants <- n_distinct(df$src_subject_id)
    participants_with_two <- n_distinct(df_filtered$src_subject_id)
    
    cat(paste0("Participants with at least 2 eventnames (out of ", total_participants, 
               " total participants): ", participants_with_two, "\n"))
    
  }, error = function(e) {
    cat("An error occurred:", e$message, "\n")
  })
  
  return(df_filtered)
  
}



## Add columns to indicate tertiles of kdm and kdm_advance
### Input 
#### df: dataframe I want to change
#### columns: columns I want to get tertiles of
### Output
#### df
split_into_tertiles <- function(df, columns) {
  df <- df %>%
    mutate_at(vars(columns), list(tertile = ~ ntile(., 3)))
  return(df)
}


## Clean Outcome df to be merged with aging data. Also checks to see if any columns have more than
## 50% missin data
### Input 
#### df: outcome df
#### columns: columns I want to keep
#### list of strings. If any column has any string, it is kept
### Output
#### df
clean_outcome_df <- function(df, columns = NULL, contains_list = NULL) {
  
  if (!is.null(contains_list)) {
    # Filter columns that contain any of the specified strings
    df <- df %>%
      select(contains(contains_list))
  } else if (!is.null(columns)) {
    # Select specified columns
    df <- df %>% select(columns)
  }
  
  # Convert all columns except "src_subject_id" and "eventname" to numeric
  df <- df %>%
    mutate_at(
      vars(-src_subject_id, -eventname), # Apply only to numeric columns
      as.numeric
    )
  
  # Check for columns with more than 50% NA's
  na_threshold <- nrow(df) * 0.5
  na_cols <- names(df)[sapply(df, function(col) sum(is.na(col)) > na_threshold)]
  
  if (length(na_cols) > 0) {
    cat("Columns with more than 50% NA's:", paste(na_cols, collapse = ", "), "\n")
  }
  
  return(df)
}


## Sum Columns Between Two Specified Columns
### Input
# `df`: The DataFrame containing columns to be summed.
# `start_col`: The name of the column where the summation should start (inclusive).
# `end_col`: The name of the column where the summation should end (inclusive).

### Output
# Returns a new DataFrame with an additional column containing the sum of values between `start_col` and `end_col` for each row.

sum_columns_between <- function(df, col_start, col_end, new_col_name, na_threshold = 0.5) {
  # Convert column names to column indices
  col_start_index <- match(col_start, names(df))
  col_end_index <- match(col_end, names(df))
  
  if (is.na(col_start_index) || is.na(col_end_index)) {
    stop("Column names not found in dataframe.")
  }
  
  df <- df %>%
    mutate(!!new_col_name := rowSums(across(.cols = col_start_index:col_end_index), na.rm = TRUE))
  
  # Calculate the proportion of NA values in each row
  na_proportions <- rowSums(is.na(df[, col_start_index:col_end_index])) / (col_end_index - col_start_index + 1)
  
  # Filter rows where the proportion of NA values is less than or equal to the threshold
  df <- df %>%
    filter(na_proportions <= na_threshold)
  
  return(df)
}


# Function Summary:
# Analyzes Child Behavior Checklist (CBCL) data for a specific eventname using various CBCL variables.

# Input:
#   - data: The dataframe containing CBCL data.
#   - cbcl_list: A list of CBCL variable names to analyze.
#   - eventname: The specific eventname to filter the data by.
#   - formula: The formula for the regression analysis.

# Output:
#   - Returns a list with two sets of results:
#     1. Fits for variables with low zero-inflation (fits_reg).
#     2. Fits for variables with high zero-inflation (fits_zero).
#     Also returns corresponding variable names (fits_reg_name and fits_zero_name).

# Usage:
#   result <- analyze_cbcl_data(data, cbcl_list, eventname, formula)
analyze_cbcl_data <- function(data, cbcl_list, my_formula) {
  fits_reg <- list()
  fits_reg_name <- list()
  fits_zero <- list()
  fits_zero_name <- list()
  
  
  for (cbcl in cbcl_list) { 
    
    # Check if there are any zero values in the column
    zero_check <- all(data[[cbcl]] != 0, na.rm = TRUE)
    
    # If there are no zero values, subtract the smallest value
    if (zero_check) {
      min_value <- min(data[[cbcl]], na.rm = TRUE)
      data[[cbcl]] <- data[[cbcl]] - min_value
    }
    
    sex_lmer <- glm(data = data, formula = as.formula(paste(cbcl, my_formula)), family = poisson)
    
    if (cbcl != "cbcl_scr_syn_totprob_t" && check_zeroinflation(sex_lmer)$ratio < 0.05) {
      fits_zero <- rbind(fits_zero, list(cbcl, check_zeroinflation(sex_lmer)$ratio))
      fits_zero_name <- append(fits_zero_name, cbcl)
    } else {
      fits_reg <- rbind(fits_reg, list(cbcl, check_zeroinflation(sex_lmer)$ratio))
      fits_reg_name <- append(fits_reg_name, cbcl)
    }
  }
  
  return(list(fits_reg = fits_reg, fits_reg_name = fits_reg_name, fits_zero = fits_zero, fits_zero_name = fits_zero_name))
}



# Function Summary:
# Analyzes any outcome data for a specific eventname to check which outcomes are zero-inflated.

# Input:
#   - data: The dataframe containing biomarker and outcome data.
#   - outcome_list: A list of outcome variable names to analyze.
#   - formula: The formula for the regression analysis.

# Output:
#   - Returns a list with two sets of results:
#     1. Fits for variables with low zero-inflation (fits_reg).
#     2. Fits for variables with high zero-inflation (fits_zero).
#     Also returns corresponding variable names (fits_reg_name and fits_zero_name).

# Usage:
#   result <- analyze_cbcl_data(data, cbcl_list, eventname, formula)

check_zero_inflation <- function(data, outcome_list, my_formula) {
  fits_reg <- list()
  fits_reg_name <- list()
  fits_zero <- list()
  fits_zero_name <- list()
  
  for (outcome in outcome_list) { 
    print(outcome)
    # Check if there are any zero values in the column
    zero_check <- all(data[[outcome]] != 0, na.rm = TRUE)
    
    # If there are no zero values, subtract the smallest value
    if (zero_check) {
      min_value <- min(data[[outcome]], na.rm = TRUE)
      data[[outcome]] <- data[[outcome]] - min_value
    }
    
    outcome_model <- glm(data = data, formula = as.formula(paste(outcome, "~", my_formula)), family = poisson)
    
    if (check_zeroinflation(outcome_model)$ratio < 0.05) {
      fits_zero <- rbind(fits_zero, list(outcome, check_zeroinflation(outcome_model)$ratio))
      fits_zero_name <- append(fits_zero_name, outcome)
    } else {
      fits_reg <- rbind(fits_reg, list(outcome, check_zeroinflation(outcome_model)$ratio))
      fits_reg_name <- append(fits_reg_name, outcome)
    }
  }
  
  return(list(fits_reg = fits_reg, fits_reg_name = fits_reg_name, fits_zero = fits_zero, fits_zero_name = fits_zero_name))
}


# Function Summary:
# The create_wide_dataframe function takes two dataframes, dataframe1 and dataframe2, along with 
# optional lists of event names eventnames1 and eventnames2. It combines these dataframes into a 
# single wide-format dataframe by filtering and pivoting the data based on the specified event names.
# The resulting merged dataframe contains a combination of columns from dataframe1 and dataframe2, 
# organized by subject ID.

create_wide_dataframe <- function(dataframe1, dataframe2, eventnames1 = NULL, eventnames2 = NULL, dataframe2_list) {
  
  # Determine which eventnames to use for long_df
  if (is.null(eventnames1)) {
    eventnames1 <- unique(dataframe1$eventname)
  } else {
    eventnames1 <- eventnames1
  }
  
  # Determine which eventnames to use for outcomes dataframe
  if (is.null(eventnames2)) {
    eventnames2 <- unique(dataframe2$eventname)
  } else {
    eventnames2 <- eventnames2
  }
  
  # Filter and pivot the first dataframe
  wide_df1 <- dataframe1 %>%
    filter(eventname %in% eventnames1) %>%
    pivot_wider(
      names_from = eventname,
      values_from = c(kdm_gmv, kdm_advance_gmv, kdm_bbb, kdm_advance_bbb, age_gmv, age_bbb)
    )
  
  # Filter the second dataframe
  wide_df2 <- dataframe2 %>%
    filter(eventname %in% eventnames2) %>%
    pivot_wider(
      names_from = eventname,
      values_from = dataframe2_list
    )
  
  # Merge the two dataframes
  merged_df <- merge(wide_df1, wide_df2, by = c("src_subject_id"))
  
  return(merged_df)
}



fit_zero_and_save <- function(df, outcome_list, my_formula) {
  
  print(my_formula)
  # Transform specific columns to factors
  df$rel_family_id <- as.factor(df$rel_family_id)
  df$site_id_l <- as.factor(df$site_id_l)
  
  fits_summary_bioage_zero <- list()
  fits_pvalue_bioage_zero <- list()
  
  for (name in outcome_list) {
    print(name)
    the_formula <- as.formula(paste(name, my_formula, sep = " ~ "))
    print(the_formula)
    # Fit the model
    fit <- cpglmm(as.formula(paste(name, "kdm_advance_gmv_baseline_year_1_arm_1*sex + (1 | site_id_l) + (1 + kdm_advance_gmv_baseline_year_1_arm_1 + gender| rel_family_id)", sep = " ~ ")), 
                  data = df)
    
    # Append the results to the lists
    fits_summary_bioage_zero <- rbind(fits_summary_bioage_zero, name, fit$inits$beta)
    fits_pvalue_bioage_zero <- rbind(fits_pvalue_bioage_zero, name, p_value(fit))
  }
  
  # Save the results
  return(fits_summary_bioage_zero)
  return(fits_pvalue_bioage_zero)
}



# This function fit_and_save_lmer takes the dataframe df, the formula, and the list of outcome names 
# outcome_list as inputs. It fits linear mixed-effects models using lmer for each outcome in the 
# list and stores the coefficient estimates and p-values in a list called results. 
# You can then access the results for each outcome and print or process them as needed.

fit_baseline_custom <- function(dataframe, outcome_list, formula_name, re_name, follow_up_timepoint, baseline_timepoint, significance_level = 0.05) {
  results <- list()
  
  for (name in outcome_list) { 
    tryCatch({
      follow_up_name <- paste0(name, "_",follow_up_timepoint)
      baseline_name <- paste0(name, "_",baseline_timepoint)
      
      # Check the proportion of NA values in the outcome variable for follow-up and baseline
      na_proportion_follow_up <- sum(is.na(dataframe %>% select(follow_up_name))) / nrow(dataframe %>% select(follow_up_name))
      na_proportion_baseline <- sum(is.na(dataframe %>% select(baseline_name))) / nrow(dataframe %>% select(baseline_name))
      
      if (na_proportion_follow_up <= 0.75 && na_proportion_baseline <= 0.75) {  # Skip if more than 75% NAs
        # Check if all values in the column are 0 for both follow-up and baseline
        if (all(dataframe[[follow_up_name]] == 0) || all(dataframe[[baseline_name]] == 0)) {
          print(paste(name, "contains all 0's for either follow-up or baseline. Skipping."))
          next  # Go to the next iteration
        }
        
        print(name)
        
        mod <- lmer(data = dataframe, 
                    as.formula(gsub(" ", "", paste(follow_up_name, "~", formula_name, "+", baseline_name, "+", re_name))), 
                    REML = FALSE)
        
        fit_summary <- summary(mod)$coefficients[, c("Estimate", "Pr(>|t|)")]
        
        # Check if any coefficient name contains "bbb" or "gmv" and has a significant p-value
        if (any(grepl("bbb|gmv", rownames(fit_summary), ignore.case = TRUE) &
                fit_summary[, "Pr(>|t|)"] < significance_level)) {
          results[[name]] <- fit_summary
        } else {
          print(paste(name, "has no significant coefficients containing 'bbb' or 'gmv'. Skipping."))
        }
      } else {
        print(paste("Skipping", name, "due to more than 75% NA values for either follow-up or baseline"))
      }
    }, error = function(e) {
      print(paste("Error in modeling", name))
      print(e)
    })
  }
  
  return(results)
}


# Same as above except for zero-inflated data
fit_baseline_4y_cpglmm <- function(dataframe, outcome_list, formula_name, re_name) {
  results <- list()
  
  for (name in outcome_list) { 
    print(name)
    name_base <- sub("_4_year_follow_up_y_arm_1", "_baseline_year_1_arm_1", name)
    
    mod <- cpglmm(data = dataframe, 
                  formula = as.formula(gsub(" ", "", paste(name, "~", formula_name, "+", name_base, "+", re_name,
                                                           "+ (1 +", formula_name, "| site_id_l)"))))
    
    fit_summary <- summary(mod)$coefficients[, c("Estimate", "Pr(>|t|)")]
    results[[name]] <- fit_summary
  }
  
  return(results)
}














