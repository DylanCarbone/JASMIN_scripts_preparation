#' Format date information in the data
#'@param obdata Data frame containing species occurrence records with the following columns: Species, Date, Gridref, Year, Week, Month, and optionally covnames (see below) and listL

add_dates <- function(obdata){
                if (!requireNamespace("lubridate", quietly = TRUE)) {
                  stop("lubridate package needed for the add_dates function to work. Please install it.",
                       call. = FALSE)
                }
            obdata$Date <- as.Date(obdata$Date)
            obdata$Year <- lubridate::year(obdata$Date)
            obdata$Month <- lubridate::month(obdata$Date)
            obdata$Week <- lubridate::week(obdata$Date)
            return(obdata)
}
