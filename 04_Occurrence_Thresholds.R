
collapse_passes <- function(data, threshold_sec) {
  
  data %>%
    dplyr::arrange(site, species, time_sec) %>%
    dplyr::group_by(site, species, date) %>%
    dplyr::mutate(
      
      time_diff = time_sec - dplyr::lag(time_sec),
      
      new_pass = dplyr::case_when(
        is.na(time_diff) ~ TRUE,
        time_diff > threshold_sec ~ TRUE,
        TRUE ~ FALSE
      ),
      
      pass_id = cumsum(new_pass)
    ) %>%
    
    dplyr::ungroup()
}

run_thresholds <- function(data, thresholds = c(30, 60, 120, 300)) {
  
  purrr::map_dfr(thresholds, function(thresh) {
    
    collapse_passes(data, thresh) %>%
      dplyr::distinct(site, species, pass_id, date) %>%
      dplyr::count(site, species, date, name = "occurrences") %>%
      dplyr::mutate(threshold_sec = thresh)
  })
}