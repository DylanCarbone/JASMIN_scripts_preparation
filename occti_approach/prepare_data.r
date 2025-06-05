#' Prepare species Observation Data
#'
#' This function reads and preprocesses species observation data by renaming columns, formatting dates,
#' converting grid references to monads, obtaining the grid region (England, Wales, Scotland, or Northern Ireland), converting to northing and easting, and computing list lengths
#' per grid reference and date.
#'
#' @param data A dataframe containing species observation data with at least the columns `tik`, `GRIDREF`, and `lower_date`.
#'
#' @details 
#' This function performs several steps to standardise and enhance a species dataset:
#' \itemize{
#'   \item Renames `tik` to `species`, `GRIDREF` to `gridref`, and `lower_date` to `date`.
#'   \item Extracts year and day-of-year (`Year` and `yday`) using `lubridate`.
#'   \item Converts grid references to standardised monads and filters out invalid ones.
#'   \item Converts monads to northing and easting values.
#'   \item Uses reference dataframe to assign region names
#'   \item Filters out unmatched records and ensures region and spatial data are present.
#'   \item Calls `add_dates()` to add any date-related columns.
#'   \item Computes `listL`, the number of unique species per `gridref` and `date`.
#' }
#'
#' @return A cleaned and enriched data frame with spatial coordinates, region info, date parts, and list length values.
#'
#' @import dplyr
#' @importFrom lubridate yday year
#' @importFrom stats na.omit
#'
#' @export

prepare_data = function(data){

# Read and preprocess data
data <- data %>% 
  rename(species = tik, gridref = GRIDREF, date = lower_date) %>%
  mutate(date = as.Date(date),
         yday = lubridate::yday(date),
         Year = lubridate::year(date))

# Handle unique grid references to save computation time
unique_gridrefs = unique(data$gridref)

# Remove invalid formats and convert to monad
unique_gridrefs = fmt_gridref(unique_gridrefs)
unique_gridrefs = reformat_gr(unique_gridrefs, prec_out = 1000, pad_gr = FALSE)
unique_gridrefs = unique(unique_gridrefs[!is.na(unique_gridrefs)]) # unique to remove new duplicates

# Merge with the original data, whilst obtaining the country, and converting to northing and easting
easting_northing_map_df = data.frame(monad = unique_gridrefs) %>%
left_join(monad_country_df %>% select(monad, region)) %>%
mutate(region = tolower(region)) %>%
cbind(OSgrid2GB_EN(unique_gridrefs)) %>% 
rename(northing = NORTHING, easting = EASTING)

# Merge with the original data, removing unpaired rows
data = data %>%
left_join(easting_northing_map_df, join_by(gridref == monad)) %>%
filter(!is.na(region))

data = add_dates(data)

data = data %>%
group_by(gridref, date) %>%
mutate(listL = length(unique(species))) %>%
ungroup()

return(data)
}