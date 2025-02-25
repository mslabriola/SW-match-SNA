df_to_lina <-
function(df) {
  df %>%
    pivot_longer(cols = -ID, names_to = "Date", values_to = "Value") %>%
    filter(!is.na(Value) & Value == 1) %>%
    dplyr::select(-Value)
}
