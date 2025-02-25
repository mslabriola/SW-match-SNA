# MATCHING OF SIGNATURE WHISTLES WITH INDIVIDUALS
# 
# Author:             Maria Silvia Labriola
# Created:            11/20/2024
# 
# R version:          R version 4.3.3
#
#
# PURPOSE
# 
# Matching bottlenose dolphin signature whistles to individuals using photo-ID 
# and PAM data from daily surveys, ensuring consistent and validated results.
# 

# Load required libraries 
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)


# User-written functions 

source("Functions/df_to_lina.R")
source("Functions/filter_repeated_ids.R")
source("Functions/match_SW_PHOTO.R")


# Load and Prepare Data ----

# Load datasets
PHOTO <- read.csv(file = "Data/PHOTO.csv", sep = ",", dec = ".") # photo-ID matrix (not in linear format)
SWlina <- read.csv(file = "Data/SWlina.csv", sep = ",", dec = ".") # SW data in linear format

# Note: "SWlina" may also be extracted directly from RavenPro selection tables following these commented steps
# (all selection tables in .txt format need to be in the Data folder)

# library(Rraven)
# rvn.dat <- imp_raven(path = "Data/RavenPro/", all.data = TRUE, name.from.file = TRUE, ext.case = "lower")
# rvn.dat <- rvn.dat %>%
#   filter(Category == "SW", Quality >= 1.5)
# SWlina <- rvn.dat  %>%
#   select(sound.files, ID_Category)
# SWlina <- SWlina %>%
#   mutate(sound.files = gsub(".*_(\\d{8})_.*", "\\1", sound.files))
# colnames(SWlina) <- c("Date", "ID")


# Reshape photo-ID data
PHOTOlina <- df_to_lina(PHOTO)

# Format dates 
PHOTOlina <- PHOTOlina %>%
  mutate(Date = substr(Date, 2, nchar(Date)))

PHOTOlina$Date <- format(as.Date(as.character(PHOTOlina$Date), format = "%Y%m%d"), "%m/%d/%Y")
SWlina$Date <- format(as.Date(as.character(SWlina$Date), format = "%Y%m%d"), "%m/%d/%Y")


# Filter data for IDs with at least two encounters
SWlina <- filter_repeated_ids(SWlina)
PHOTOlina <- filter_repeated_ids(PHOTOlina)


# Match SWs to individuals ----

# Join photo_ID and SW datasets by date
joined_by_date <- inner_join(SWlina, PHOTOlina, by = "Date", relationship = "many-to-many")

# Verify and filter combinations
result_df <- match_SW_PHOTO(joined_by_date)


# Create and Style Output Table ----
colnames(result_df) <- c("SW_ID", "PHOTO_ID")

table1 <- kable(result_df, format = "html", caption = "SWs-photo-identified individuals matchings") %>%
  kable_paper(full_width = F) %>%
  column_spec(1, width = "3cm") %>%
  column_spec(2, width = "3cm")

table1


saveRDS(result_df, "matchings.RData")
