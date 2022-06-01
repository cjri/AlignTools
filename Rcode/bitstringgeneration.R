##sim code
library(MASS)
library(R.utils)
library(gtools)

##function for generating bitstrings
generatebitstring <- function(mu, omega, pos, n) {
  ##take logit of mu
  logitmu <- logit(mu)[,1]
  
  ##set up object to hold bitstrings 
  strings <- rep(NA, length = n)
  positions <- rep(NA, length = n)
  
  ##generate MVN variates on logit scale
  logittarget <- mvrnorm(n = n, logitmu, omega, tol = 1e-06, empirical = FALSE)
  logittarget <- ifelse(logittarget > 0, 1, 0)
  
  for (i in 1:n) { 
    ##apply inverse logistic function
    targetextend <- rep(0, 29900)
    targetextend[pos[,1]] <- logittarget[i,]
    
    ##save bitstring
    strings[i] <- paste(targetextend, collapse = "")
    positions[i] <- paste(pos[logittarget[i,], 1], collapse = ",")
  }
  
  ##output bitstrings
  return(list(strings,positions))
}

##function for handling uncompression and file I/O
run <- function(n) {
  ##extract locations
  locations <- list.dirs(recursive = F, full.names = F)
  
  ##extract relevant part of filename
  filename <- do.call("c", lapply(locations, function(x) {return(strsplit(x, "Week")[[1]][2])}))
  
  #for each location
  for (j in 1:length(locations)) {
    ##read in probabily vector
    mu <- read.table(paste0(getwd(), "/", locations[j], "/", "Variant_frequencies", filename[j], ".out"), header = F, sep = " ")
    
    ##unzip correlation matrix
    gunzip(filename = paste0(getwd(), "/", locations[j], "/", "Variant_correlations", filename[j], ".out.gz"), destname = paste0(getwd(), "/", locations[j], "/", "Variant_correlations", filename[j], ".out"), remove = FALSE)
    
    ##read in correlation matrix
    omega <- as.matrix(read.table(paste0(getwd(), "/", locations[j], "/", "Variant_correlations", filename[j], ".out"), header = F, sep = " "))
    if (nrow(omega) != ncol(omega)) {
      omega <- omega[,-ncol(omega)]
    }
    
    ##read in positions
    pos <- read.table(paste0(getwd(), "/", locations[j], "/", "Variant_positions", filename[j], ".out"), header = F, sep = " ")
    
    ##apply function
    bitstrings <- generatebitstring(mu, omega, pos, n)
    
    ##write out
    write.table(bitstrings[[1]], file = paste0(getwd(), "/", "bitstrings", filename[j], ".txt"), row.names = F, col.names = F, sep = "\t", eol = "\n")
    write.table(bitstrings[[2]], file = paste0(getwd(), "/", "positions", filename[j], ".txt"), row.names = F, col.names = F, sep = "\t", eol = "\n")
  }
}

##perform task
run(n = 50)

