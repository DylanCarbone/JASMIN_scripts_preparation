#' Format date information in the data
#'@param obdata Data frame containing species occurrence records with the following columns: species, Date, gridref, Year, Week, Month, and optionally covnames (see below) and listL

add_dates <- function(obdata){
                if (!requireNamespace("lubridate", quietly = TRUE)) {
                  stop("lubridate package needed for the add_dates function to work. Please install it.",
                       call. = FALSE)
                }
            obdata$date <- as.Date(obdata$date)
            obdata$Year <- lubridate::year(obdata$date)
            obdata$Month <- lubridate::month(obdata$date)
            obdata$Week <- lubridate::week(obdata$date)
            return(obdata)
}
