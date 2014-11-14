#' @rdname phase.run

phase.read.pair <- function(out.file) {   
  if(!file.exists(out.file)) return(NULL)
  pair.file <- scan(file = out.file, what = "character", sep = "\n", quiet = TRUE)
  
  id.start <- grep("IND:", pair.file)
  gtype.probs <- lapply(1:length(id.start), function(i) {
    id.end <- ifelse(i == length(id.start), length(pair.file), id.start[i + 1] - 1)
    id <- sub("IND: ", "", pair.file[id.start[i]])
    t(sapply((id.start[i] + 1):id.end, function(j) {
      line.split <- unlist(strsplit(pair.file[j], " , "))
      names(line.split) <- c("hap1", "hap2", "pr")
      c(id = id, line.split)
    }))
  })             
  gtype.probs <- as.data.frame(do.call(rbind, gtype.probs), stringsAsFactors = FALSE)  
  gtype.probs$pr <- as.numeric(as.character(gtype.probs$pr))
  gtype.probs
}
