# function to plot the hbond network
cutoff <- 0.5
fn <- "hbonds/WTaverage_hbonds.out"
plotGraph <- function(fn="hbonds/WTaverage_hbonds.out",
                      cutoff=0.7){
    require(Rgraphviz) # loads rgraphviz library
    networkfile <- read.table(file=fn,header = FALSE) # read in the file
    networkfile <- networkfile[networkfile$V5>cutoff,] # keeps only those above the cutoff value for plotting
    acceptors <- as.vector(unique(networkfile$V1)) # keep the acceptors
    donors <- as.vector(unique(networkfile$V3)) # keep the donors
    nodes <- unique(c(acceptors,donors)) # this is the node list
    # generate the edge list
    edgelist <- lapply(donors,function(x){w <- networkfile[networkfile$V3==x,];m <- w$V1;m <- as.vector(m);m})
    names(edgelist) <- donors
    graph <- new("graphNEL", nodes=nodes, edgemode="directed")
    for(j in 1:length(edgelist)){
            graph <- addEdge(from=names(edgelist[j]),unlist(edgelist)[j],graph,1)
            }
    return(list=c(NEL=graph,Edges=list(edgelist)))
    }
        


wtGraph <- plotGraph(cutoff=0.5)
mutGraph <- plotGraph(fn="hbonds/Mutaverage_hbonds.out",cutoff=0.55)
#png("hbonds/graphviz1.png")
par(mfrow = c(2,1))
plot(wtGraph$NEL,main="wt. hydrogen bond with cutoff 0.75")
plot(mutGraph,main="Mut hydrogen bond with cutoff 0.75")
#dev.off()
