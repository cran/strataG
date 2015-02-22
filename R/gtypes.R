#' @rdname gtypes
#' @export gtypes gtypes.default gtypes.character gtypes.list is.gtypes check.gtypes stopifnot.gtypes
#' 
#' @title Create a gtypes object
#' @description Creates a coded list of stratified genotypes or haplotypes for use with 
#'   most functions in the \code{strataG} package. 
#' 
#' @param gen.data a \code{matrix}, \code{data.frame}, \code{list} of sequences, or filename where genetic data is stored.
#'   The columns designating ids and strata should come before the locus columns. The columns containing genetic data,
#'   (designated by \code{locus.col}) must be at the end.
#'   If they are diploid loci, then every two columns are considered to be two alleles of a locus. It is advised to use
#'   \code{\link{read.gen.data}} to import a .csv file for use.
#' @param id.col a number or character string specifying the column where sample IDs are located.
#' @param strata.col a number of character string specifying the column where strata designations are located
#' @param locus.col a number specifying the first column of the genetic loci. All other loci are assumed to follow this column.
#' @param dna.seq a list of DNA sequences for haploid data.
#' @param description a string naming or describing this dataset. Will be used in summaries and plots. 
#' @param g a gtypes object.
#' @param delete.missing.strata logical. Delete samples for which strata is missing (NA)?
#' @param code.start integer to start numeric recoding at. Must be >= 0.
#' @param strata.vec a character or numeric vector specifying the stratum that each sample should be assigned to.
#' @param ploidy character giving ploidy of gtypes object being tested by \code{is.gtypes}. Can be "haploid", "diploid", or "any".
#' @param x an R object to test.
#' @param show.warnings logical - show warnings for is.dna.seq describing check failures?
#' @param ... arguments to the default method 
#' 
#' @return a \code{gtypes} object which is a list containining the following elements: \cr
#' \tabular{ll}{
#'   \code{genotypes} \tab A numeric matrix where the first column ("strata") gives the stratification. 
#'     If haploid, the second column in the matrix lists haplotypes. If diploid, every two columns afterwards 
#'     are two alleles of the same locus.\cr
#'   \code{sequences} \tab If haploid, and provided on creation, a list of aligned DNA sequences. Otherwise \code{NULL}.\cr
#' }
#' 
#' @note A gtypes object has the strata and locus data coded numerically for ease of use in 
#'   analytical functions. Mapping attributes \code{strata.name} and \code{locus.name} are attached to the
#'   object. To map to the original data, one can use \code{\link{decode}}, \code{\link{as.matrix.gtypes}},
#'   or \code{\link{as.data.frame.gtypes}}, however, the resulting object from these cannot be used 
#'   where a \code{gtypes} object is required.
#' 
#' @author Eric Archer <eric.archer@@noaa.gov>
#' 
#' @seealso \code{\link{decode}}, \code{\link{as.matrix.gtypes}}, \code{\link{as.data.frame.gtypes}}
#'
#' @examples
#' data(dolph.strata)
#' data(dolph.haps)
#' data(dolph.msats)
#' 
#' # Create gtypes object for mtDNA data
#' mtdna <- gtypes(dolph.strata, id.col = 1, strata.col = 2, locus.col = 4, dna.seq = dolph.haps)
#' summary(mtdna)
#' 
#' # Create gtypes object for microsatellite data
#' # First, merge strata file and genotypes
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5, description = "msats")
#' summary(msats)

gtypes <- function(gen.data, ...) UseMethod("gtypes")
