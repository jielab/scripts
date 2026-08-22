#Author: Linsey Jackson
#Date: 12/9/2025
#Script name: delta endpoints
#Purpose: calculate change in endpoints 

library(dplyr)
library(tidyr)
library(purrr)
library(rlang)
library(rstatix)
library(haven)
library(data.table)

#read in ADaM table
#select relevant columns 
advs <- advs %>%
  select(USUBJID, AVISIT, PARAM, AVAL, BASE, CHG, PCHG)

advs_change <- advs %>%
  # Keep only Week 72
  filter(AVISIT %in% c("Week 72 (Visit 21)")) %>%
  filter(PARAM %in% c("Weight (kg)", "Waist Circumference (cm)")) %>%
  pivot_wider(
    id_cols = c(USUBJID, AVISIT),
    names_from = PARAM,
    values_from = c(AVAL, BASE, CHG, PCHG)
  )

advs_endpoint_final <- advs_change %>%
  select(USUBJID, 'PCHG_Weight (kg)', 'CHG_Waist Circumference (cm)', 'PCHG_Waist Circumference (cm)', 'CHG_Weight (kg)')

#select height from baseline
height <- advs %>%
  filter(PARAM %in% c("Height (cm)"))
height <- height %>%
  select(USUBJID, AVAL)
height$height <- height$AVAL

#combine height and anthropometric measures into one df
df_list <- list(advs_endpoint_final, height)
combined_df <- reduce(df_list, left_join, by = "USUBJID")

#calculate change in WHtR (change in WC)/height
combined_df$WHtR_Change <- combined_df$`CHG_Waist Circumference (cm)`/combined_df$height

#filter for any missing values
final_endpoints<-combined_df %>%
  drop_na()