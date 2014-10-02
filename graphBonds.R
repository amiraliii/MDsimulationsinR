# function to plot the hbond network
filename <- "hbonds/WTaverage_hbonds.out"
cutoff <- 0.5
plotGraph <- function(filename,
                      cutoff=0.5){
    require(Rgraphviz) # loads rgraphviz library
    networkfile <- read.table(filename,header = FALSE) # read in the file
    networkfile <- networkfile[networkfile$V5>cutoff,] # keeps only those above the cutoff value for plotting
    acceptors <- unique(networkfile$V1) # keep the acceptors
    donors <- unique(networkfile$V3) # keep the donors
    }
