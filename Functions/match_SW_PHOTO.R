match_SW_PHOTO <-
function(joined_data) {
  # Step 1: Verify combinations
  grouped_df <- joined_data %>%
    group_by(ID.x, ID.y) %>%
    summarise(occurrences = n(), .groups = "drop")
  
  grouped_df$verification <- sapply(1:nrow(grouped_df), function(i) {
    sound_id <- grouped_df$ID.x[i]
    individual_id <- grouped_df$ID.y[i]
    subset_xy <- joined_data %>% filter(ID.x == sound_id, ID.y == individual_id)
    subset_x <- joined_data %>% filter(ID.x == sound_id)
    missing_dates <- setdiff(unique(subset_x$Date), subset_xy$Date)
    if (length(missing_dates) == 0) {
      "All dates present"
    } else {
      paste("Missing dates:", paste(missing_dates, collapse = ", "))
    }
  })
  
  # Filter for valid combinations
  verified_data <- grouped_df %>%
    filter(verification == "All dates present") %>%
    dplyr::select(ID.x, ID.y)
  
  # Step 2: Filter unique matches
  verified_data <- verified_data %>%
    group_by(ID.x) %>%
    filter(n() == 1) %>%
    ungroup() %>%
    group_by(ID.y) %>%
    filter(n() == 1) %>%
    ungroup()
  
  return(verified_data)
}
