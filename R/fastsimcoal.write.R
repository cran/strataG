#' @rdname fastsimcoal.run

fastsimcoal.write <- function(num.pops, Ne, sample.size, sample.time = rep(0, num.pops),
                              growth.rate = rep(0, num.pops), mig.mat = NULL, 
                              hist.ev = NULL, num.chrom = 1, data.type = NULL,
                              locus.params = NULL, label = "fastsimcoal") {

  Ne <- rep(Ne, length.out = num.pops)
  sample.size <- rep(sample.size, length.out = num.pops)
  sample.time <- rep(sample.time, length.out = num.pops)
  hist.ev <- if(is.list(hist.ev)) do.call(rbind, hist.ev) else rbind(hist.ev)
  locus.params <- if(is.list(locus.params)) do.call(rbind, locus.params) else rbind(locus.params)
  if(nrow(locus.params) == 1 & num.chrom > 1) locus.params <- do.call(rbind, lapply(1:num.chrom, function(i) locus.params[1, ]))

  file <- paste(label, ".par", sep = "")
  write(paste("//  <<", label, ">>  (input from 'fastsimcoal.write')"), file)  
  write(paste(num.pops, "populations to sample"), file, append = T)
  
  write("//Population effective sizes", file, append = T)
  for(i in 1:length(Ne)) write(Ne[i], file, append = T)
  
  write("//Sample sizes", file, append = T)
  for(i in 1:length(sample.size)) {
    s.t <- ifelse(sample.time[i] == 0, "", paste(" ", sample.time[i], sep = ""))
    write(paste(sample.size[i], s.t, sep = ""), file, append = T)
  }
  
  write("//Growth rates", file, append = T)
  for(i in 1:length(growth.rate)) write(growth.rate[i], file, append = T)
  
  write("//Number of migration matrices", file, append = T)
  write(length(mig.mat), file, append = T)
  if(!is.null(mig.mat)) {
    for(i in 1:length(mig.mat)) {
      write("//migration matrix", file, append = T)
      for(r in 1:nrow(mig.mat[[i]])) write(mig.mat[[i]][r, ], file, append = T)
    }
  }
  
  write("//Historical events: time, source, sink, migrants, new size, growth rate, migr. matrix", file, append = T)
  write(ifelse(is.null(hist.ev), 0, nrow(hist.ev)), file, append = T)
  if(!is.null(hist.ev)) {
    for(i in 1:nrow(hist.ev)) {
      write(paste(hist.ev[i, ], collapse = " "), file, append = T)
    }
  }
  
  write("//Number of independent loci [chromosome]", file, append = T)
  write(paste(num.chrom, "0"), file, append = T)
  write("//Per chromosome: Number of linkage blocks", file, append = T)
  write("1", file, append = T)
  for(i in 1:num.chrom) {
    write("//Per block: data type, num loci, rec. rate and mut rate + optional parameters", file, append = T)
    write(paste(data.type, paste(locus.params[i, ], collapse = " ")), file, append = T)
  }
  
  file
}