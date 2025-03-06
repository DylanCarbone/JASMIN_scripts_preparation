#' Add column for list lengths to the data
#' @import data.table

add_listL <- function(obdata){
                obdata <- data.table(obdata)
                obdata[, listL:=uniqueN(Species), by=.(Date, Gridref)]
                return(obdata)
                }
