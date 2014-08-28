library(ggplot2)

# reads in an out put file from cpptraj (in this case, PCA data)
jet.colors <- colorRampPalette(c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000")) # defines the jet.color scheme

# read in the data from ptraj
WT <- read.table("WTprojection.ppj")
Mut <- read.table("Mutprojection.ppj")

# keep only PC1 and PC2
WT <- WT[,2:3]
Mut <- Mut[,2:3]
PCA <- as.data.frame(rbind(WT,Mut)) # binds them all together in 1 single data frame
colnames(PCA) <- c("PC1","PC2") # sets the column name of the data frame for easier reference
PCA$Simulation <- c(rep("WT",24000),rep("Mut",24000)) # identify which values belong to the WT and mutant respectively

# this is the actual plotting call to plot the 2D density
ggplot(PCA,aes(x=PC1,y=PC2))+stat_density2d(geom="tile",aes(fill=..density..),contour = FALSE)+scale_fill_gradientn(colours = jet.colors(7))+facet_grid(~Simulation)+theme_bw()+ggtitle("PCA clustering using Cartesian coordinates")+scale_x_continuous(breaks=seq(-50,50,10))+scale_y_continuous(breaks=seq(-20,40,10))
# in the aes(x,y) call x and y refers to the column names in the data frame where the data is being stored in


