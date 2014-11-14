#' @export read.gen.data
#' 
#' @title Read Genetic Data
#' @description Reads a .csv file into a \code{data.frame} formatted as genetic data. Replaces commonly
#'   used missing data values with NA and removes blank lines.
#' 
#' @param file filename of .csv file.
#' @param na.strings see \code{\link{read.table}}.
#' @param ... other arguments passed to \code{\link{read.table}}.
#' 
#' @return a \code{data.frame}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

read.gen.data <- function(file, na.strings = c(NA, "NA", "", " ", "?", "."), ...) {
  df <- read.csv(file = file, na.strings = na.strings, colClasses = "character",  stringsAsFactors = FALSE, ...)
  df <- df[rowSums(is.na(df) | df == "") != ncol(df), ]
}
