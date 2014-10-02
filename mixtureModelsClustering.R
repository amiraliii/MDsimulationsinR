library(EMCluster)

# perform EM model based clustering in R, and then identifying the centroids

# define function to identify centroid from our clustering
centroid <- function(x,y){
    # min(euclidean distance) to identify centroid
    # x = PC matrix, y= clustering result from EM mut
    stopifnot(ncol(x)==2)
    n <- nrow(y$Mu) # no of clusters
    ids <- c() # to hold ID
    pc1cor <- c() # hold the x-coordinates of the centroid
    pc2cor <- c() # hold the y-coordinates of the centroid
    colnames(x) <- c("PC1","PC2")
    targetx <- c() # hold the x-coordinates of the closest structure
    targety <- c() # hold the y-coordinates of the closest structure
    eucD <- c()
    for (i in 1:n){
        pc1 <- y$Mu[i,1]
        pc2 <- y$Mu[i,2]
        res <- which.min(((x$PC1-pc1)^2+(x$PC2-pc2)^2)^0.5) # id the closest structure
        #euc.d <- ((x$PC1-pc1)^2+(x$PC2-pc2)^2)^0.5
        ids[i] <- res+2000 # since we only considered from frame 2001 onwards, but we are reading in all frames from 1-10000
        pc1cor[i] <- pc1
        pc2cor[i] <- pc2
        target <- x[res,]
        #ucD[i] <- euc.d
        m<- target[,1]
        targetx[i] <- m
        n<- target[,2]
        targety[i] <- n
        }
    f <- data.frame(PC1=pc1cor,PC2=pc2cor,Centroid=ids,CentroidX=targetx,CentroidY=targety,Population=y$pi) # creates the data frame to be returned
    f # returns the dataframe
    }

# reads in 1 PCA data set for illustrative purpose
WT1 <- read.table("WT1projection.ppj")
WT1clustering <- exhaust.EM(WT1[,2:3],nclass=4,min.n.iter = 50) # cluster along PC1 and PC2 only (PC1 and PC2 data in columns 2 and 3 in the data)

#png("images/WT1-EMclusteringPCA.png") # for saving the image
plotem(x=WT1[,2:3],emobj = WT1clustering,xlab="PC1",ylab="PC2",main="EM-clustering of PC space (WT1)")
#dev.off() # close graphic device if the image is to be exported to a png file

WT1centers <- centroid(WT1[,2:3],WT1clustering) # to get the centroids from PC1 and PC2 data which was used to perform the clustering





