# function to plot the hbond network
cutoff <- 0.5
plotGraph <- function(fn="hbonds/WTaverage_hbonds.out",
                      cutoff=0.5){
    require(Rgraphviz) # loads rgraphviz library
    networkfile <- read.table(file=fn,header = FALSE) # read in the file
    networkfile <- networkfile[networkfile$V5>cutoff,] # keeps only those above the cutoff value for plotting
    acceptors <- unique(networkfile$V1) # keep the acceptors
    donors <- unique(networkfile$V3) # keep the donors
    nodes <- paste0("Residue",unique(c(acceptors,donors))) # this is the node list
    # generate the edge list
    edgelist <- lapply(donors,function(x){w <- networkfile[networkfile$V3==x,];m <- w$V1;m <- as.vector(m);m <-paste0("Residue",m);m})
    names(edgeli)
    graph <- new("graphNEL", nodes=nodes, edgeL=edgelist, edgemode="undirected")
    graph <- layoutGraph(graph)
    }
