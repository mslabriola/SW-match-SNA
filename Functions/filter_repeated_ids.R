filter_repeated_ids <-
function(df_lina) {
  rep_ids <- df_lina %>%
    group_by(ID) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1) %>%
    pull(ID)
  
  df_lina %>%
    distinct()  %>%
    filter(ID %in% rep_ids) 
}
