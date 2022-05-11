##sim code
library(MASS)
library(utils)

##function for generating bitstrings
generatebitstring <- function(mu, omega, pos, n) {
  ##take logit of mu
  logitmu <- logit(mu)[,1]
  
  ##set up object to hold bitstrings 
  strings <- rep(NA, length = n)
  
  for (i in 1:n) {
    ##generate MVN variates on logit scale
    logittarget <- mvrnorm(n = 1, logitmu, omega, tol = 1e-06, empirical = FALSE)
    
    ##apply inverse logistic function
    target <- plogis(logittarget)
    targetextend <- rep(0, 29900)
    targetextend[pos[,1]] <- rbinom(rep(1, length(target)), rep(1, length(target)), target)
    
    ##save bitstring
    strings[i] <- paste(targetextend, collapse = "")
  }
  
  ##output bitstrings
  return(strings)
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
    mu <- as.vector(read.table(paste0(getwd(), "/", locations[j], "/", "Variant_frequencies", filename[j], ".out"), header = F, sep = " "))
    
    ##unzip correlation matrix
    untar(paste0(getwd(), "/", locations[j], "/", "Variant_correlations", filename[j], ".out.gz"), list = FALSE, exdir = paste0(getwd(), "/", locations[j]), compressed = T)
    
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
    write.table(bitstrings, file = paste0(getwd(), "/", "bitstrings", filename[j], ".txt"), row.names = F, col.names = F, sep = "\t", eol = "\n")
  }
}

##perform task
run(n = 50)

