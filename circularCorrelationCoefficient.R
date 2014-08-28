# codes to calculate CCC based on based on backbone phi and psi angles

# reads in the data first
WT.phi <- read.table("outfiles/WT_1200nsPhi.out")
WT.psi <- read.table("outfiles/WT_1200nsPsi.out")
Mut.phi <- read.table("outfiles/Mut_1200nsPhi.out")
Mut.psi <- read.table("outfiles/Mut_1200nsPsi.out")

# loads the libraries needed
library(circular)
library(reshape2)
library(ggplot2)

# identifying only data from the L3 region
WT.phi.L3 <- WT.phi[,147:161]
WT.psi.L3 <- WT.psi[,147:161]
Mut.phi.L3 <- Mut.phi[,147:161]
Mut.psi.L3 <- Mut.psi[,147:161]

# converts data to circular format first
WT.phi.L3 <- as.circular(WT.phi.L3,units="degrees")
WT.psi.L3 <- as.circular(WT.psi.L3,units="degrees")
Mut.phi.L3 <- as.circular(Mut.phi.L3,units="degrees")
WT.psi.L3 <- as.circular(WT.psi.L3,units="degrees")

# calculates the correlations
WT.L3.Cor <- cor.circular(WT.phi.L3,WT.psi.L3)
Mut.L3.Cor <- cor.circular(Mut.phi.L3,Mut.psi.L3)

# melts the data frame from wide to long format
WT.L3.Cor.melt <- melt(WT.L3.Cor)
Mut.L3.Cor.melt <- melt(Mut.L3.Cor)
colnames(Mut.L3.Cor.melt) <- colnames(WT.L3.Cor.melt) <- c("Psi","Phi","Correlation")

# identifies the mutant and the WT simulations
WT.L3.Cor.melt$Simulation <- "WT"
Mut.L3.Cor.melt$Simulation <- "Mut"

# numbering the residues
WT.L3.Cor.melt$Psi <- WT.L3.Cor.melt$Psi+235
WT.L3.Cor.melt$Phi <- WT.L3.Cor.melt$Phi+235
Mut.L3.Cor.melt$Psi <- Mut.L3.Cor.melt$Psi+235
Mut.L3.Cor.melt$Phi <- Mut.L3.Cor.melt$Phi+235

# combining both WT and mutant together to plot in a single call using facet_grid
dataplot <- rbind(WT.L3.Cor.melt,Mut.L3.Cor.melt)
dataplot$State <- unlist(lapply(dataplot$Correlation,function(x){y <- ifelse(x>0.2,x <- "red",ifelse(x<(-0.2),x <- "blue",x <- "white"));y}))

# plotting call
ggplot(dataplot,aes(y=Psi,x=Phi,fill=State))+geom_tile()+scale_fill_manual(values=c("red"="red","blue"="blue","white"="white"))+facet_grid(~Simulation)+scale_x_continuous(breaks=seq(236,250,1))+scale_y_continuous(breaks=seq(235,250,1))+xlab(expression(paste("Residue 1 " ,psi)))+ylab(expression(paste("Residue 2 ",phi)))+theme_bw()+ggtitle(expression(paste("Backbone dihedral ",phi,"/",psi," correlation at L3 loop")))+theme(axis.text.x=element_text(angle=90))
