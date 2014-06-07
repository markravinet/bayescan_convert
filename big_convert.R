#### Convert LARGE genepop SNP file to BayesScan input file ####
# To be used on cluster
# Mark Ravinet September 2013, University of Gothenburg
## Modified heavily from code provided by Kevin Keenan (cheers Kevin!)

big_convert <- function(infile = FALSE, outfile = FALSE){
    #require(diveRsity)

    fastScan <- function(fname) {
    # function reads genepop file in as a single vecto
    s <- file.info(fname)$size
    buf <- readChar(fname, s, useBytes = TRUE)
    return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
  }
  
  # fastscan - read in as a vector
  dat <- fastScan(fname = infile)
  
  if (length(strsplit(dat[length(dat)], split = "\\s+")[[1]]) == 
       1) {
    dat <- dat[-(length(dat))]
  }
  rm(fastScan)
  z <- gc()
  rm(z)
  popLocation <- grep("^([[:space:]]*)POP([[:space:]]*)$", 
                      toupper(dat))
  pop_pos <- c(popLocation, (length(dat) + 1))
  loci_names <- as.vector(sapply(dat[2:(pop_pos[1] - 1)], 
                                 function(x) {
                                   gsub(pattern = "\\s+", replacement = "", x)
                                 }))
  popSizes <- NULL
  for (i in 1:(length(pop_pos) - 1)) {
    popSizes[i] <- length((pop_pos[i] + 1):(pop_pos[(i + 
                                                       1)] - 1))
  }
  pops <- dat[-(c(1:(popLocation[1] - 1), popLocation))]
  popList <- lapply(seq_along(popSizes), function(i) {
    if (i == 1) {
      indx <- 1:popSizes[i]
    }
    else {
      indx <- (sum(popSizes[1:(i - 1)]) + 1):((sum(popSizes[1:(i - 
                                                                 1)])) + popSizes[i])
    }
    return(pops[indx])
  })
  npops <- length(popList)
  nloci <- length(loci_names)
  pop_sizes <- popSizes
  rm(dat, pops)
  z <- gc(reset = TRUE)
  rm(z)
  testStr <- strsplit(popList[[1]][1], split = "\\s+")[[1]]
  gpEst <- sapply(testStr, function(x) {
    if (is.character(x)) {
      nchar(x)/2
    }
    else {
      NA
    }
  })
  rm(testStr)
  gp <- as.numeric(names(sort(-table(gpEst)))[1])
  prePopList <- lapply(popList, function(x) {
    y <- array(data = NA, dim = c(length(x), (nloci + 1), 
                                  2))
    colnames(y) <- c("ind", loci_names)
    for (j in 1:length(x)) {
      data <- strsplit(x[j], split = "\\s+")[[1]]
      if (data[2] == ",") {
        data <- data[-2]
      }
      data[data == "0000"] <- NA
      data[data == "NANA"] <- NA
      data[data == "0"] <- NA
      data[data == "000000"] <- NA
      data[data == "999999"] <- NA
      data[data == "-9-9"] <- NA
      y[j, 2:(nloci + 1), 1] <- substr(data[2:(nloci + 
                                                 1)], 1, gp)
      y[j, 2:(nloci + 1), 2] <- substr(data[2:(nloci + 
                                                 1)], gp + 1, gp * 2)
      y[j, 1, 1] <- data[1]
      y[j, 1, 2] <- data[1]
    }
    return(y)
  })
  rm(popList)
  ind_names <- lapply(prePopList, function(x) {
    return(x[, 1, 1])
  })
  pop_names <- sapply(ind_names, function(x) {
    return(x[1])
  })
  
  # From here the function diverges from BigDivPart and uses
  # some code lifted from bigPreDiv
  popList <- lapply(prePopList, function(x){
    return(x[,(2:(nloci+1)),])
  })
  # count the numbers of individuals typed per population
  indtyp <- lapply(popList, function(x){
    apply(x, 2, function(y){
      length(na.omit(y[,1]))
    })
  })
  # get unique alleles per pop
  alls <- lapply(seq_along(popList), function(i){
    apply(popList[[i]], 2, function(x){
      return(unique(c(x[,1], x[,2])))
    })
  })
  # get unique alleles across pops (maybe unneccess)
  all_alleles <- lapply(1:nloci, function(i){
    alleles <- lapply(alls, function(x){
      return(x[[i]])
    })
    return(sort(unique(unlist(alleles))))
  })
  rm(alls)
  # count all observed allele numbers per population
  # (parallel is slower)
  obsAlls <- lapply(popList, function(x){
    apply(x, 2, function(y){
      als <- unique(c(na.omit(y[,1]), na.omit(y[,2])))
      counts <- sapply(als, function(z){ # takes allele, then looks for it in y
      res <- length(which(y == z))   # returns count of alleles
      names(res) <- names(z)
      return(res)
      })
    })
  })
  
  # observed alleles 
  obs_all <- lapply(1:nloci, function(i){
    loc <- matrix(nrow = length(all_alleles[[i]]),
                  ncol = npops)
    rownames(loc) <- all_alleles[[i]]
    for(j in 1:npops){
      o <- obsAlls[[j]][[i]]
      loc[names(o), j] <- o
    }
    loc[is.na(loc)] <- 0
    return(loc)
  })
  rm(obsAlls)

  # Some other necessary calculations
  lociID <- 1:nloci
  n_alls <- unlist(lapply(all_alleles, length))
  rm(all_alleles)
  
  #### Produce BayesScan file ####
  format_write <- function(info, add = TRUE){
    cat(unlist(info), file = outfile, append = add)
    write(" ", file = outfile, append = TRUE)
  }
  
  big_pop_format <- function(pop){
    total_genes_pop <- as.numeric(unlist(indtyp[[pop]]))*2
    obs_all <- sapply(1:nloci, function(x) obs_all[[x]][, pop])
    # This line is what should be written per locus to a file
    pop_info <- paste("[pop]=", pop, "\n", sep = "")
    obs_loc_pop <- sapply(1:nloci, function(j) {
      as.vector(c(lociID[j], total_genes_pop[j], n_alls[j], obs_all[[j]], "\n"))
    })
    out <- list(pop_info, obs_loc_pop)
    return(out)
  }
  
  # Write header
  header <- paste("[loci]=", nloci, "\n\n", "[populations]=", npops, "\n", sep = "")
  format_write(header, add = FALSE)
  # Write pop info
  info <- sapply(1:length(pop_sizes), big_pop_format, simplify = FALSE)
  lapply(1:length(pop_sizes), function(x) format_write(info[[x]], add = TRUE))
  
  
}












    

  