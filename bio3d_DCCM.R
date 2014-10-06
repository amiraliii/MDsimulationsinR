# DCCM codes in Bio3d
require(bio3d)
WT.traj <- read.ncdf("WT_cat.nc")
WT.reference <- read.pdb("averagestructures/WTaverage.pdb")
WT.calpha <- atom.select(WT.reference,elety="CA")$xyz
WT.traj <- WT.traj[,WT.calpha]
WT.reference <- WT.reference$xyz[WT.calpha]
WT.traj <- fit.xyz(fixed=WT.reference,mobile=WT.traj) # performs fitting using all residue CA

Mut.traj <- read.ncdf("Mut_cat.nc")
Mut.reference <- read.pdb("averagestructures/Mutaverage.pdb")
Mut.calpha <- atom.select(Mut.reference,elety="CA")$xyz
Mut.traj <- Mut.traj[,Mut.calpha]
Mut.reference <- Mut.reference$xyz[Mut.calpha]
Mut.traj <- fit.xyz(fixed=Mut.reference,mobile=Mut.traj) # performs fitting using all residue CA

# performs the DCCM computation
WT.DCCM <- dccm(WT.traj)
Mut.DCCM <- dccm(Mut.traj)

# because I like ggplot more for graphics ...
# convert to dataframe format
WT.DCCM<- as.data.frame(WT.DCCM)
Mut.DCCM<- as.data.frame(Mut.DCCM)

# melting the DF from wide to long format
library(reshape2)
WT.DCCMmelt <- melt(WT.DCCM)
Mut.DCCMmelt <- melt(Mut.DCCM)
colnames(WT.DCCMmelt) <- colnames(Mut.DCCMmelt) <- c("Residue1","Correlation")
WT.DCCMmelt$Residue2 <- rep(91:289,199)
Mut.DCCMmelt$Residue2 <- rep(91:289,199)
WT.DCCMmelt$Residue1 <- as.numeric(WT.DCCMmelt$Residue1)+90
Mut.DCCMmelt$Residue1 <- as.numeric(Mut.DCCMmelt$Residue1)+90
WT.DCCMmelt$Simulation <- "WT"
Mut.DCCMmelt$Simulation <- "Mut"
DCCM <- rbind(WT.DCCMmelt,Mut.DCCMmelt)
# plotting the results in ggplot
library(ggplot2)
ggplot(DCCM,aes(x=Residue1,y=Residue2,fill=Correlation))+geom_tile()+scale_fill_gradientn(colours = jet.colors(7))+facet_grid(~Simulation)+theme_bw()+ggtitle(expression(paste("DCCM matrix of C-",alpha)))









