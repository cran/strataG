#' @rdname phase.run

phase.read.sample <- function(out.file, type) {
  if(!file.exists(out.file)) return(NULL)
  post.file <- scan(file = out.file, what = "character", sep = "\n", quiet = TRUE)
  iter.start <- grep(type, post.file) + 1
  lapply(iter.start, function(start) {
    num.samples <- as.integer(post.file[start - 3])
    end <- start + (num.samples * 3) - 3
    as.data.frame(t(sapply(seq(start, end, by = 3), function(i) {
      id <- strsplit(post.file[i], " ")[[1]][2]    
      hap1 <- gsub(" ", "", post.file[i + 1])
      hap2 <- gsub(" ", "", post.file[i + 2])
      c(id, hap1, hap2)
    })), stringsAsFactors = FALSE)
  })
}
