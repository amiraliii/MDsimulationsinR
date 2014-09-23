# codes for performing PCA and PCA analysis in R using bio3d

# loads libraries
library(ggplot2)
library(reshape2)
library(bio3d)


# define a function for use to export AU as b-factor values into pdb file
assignAU <- function(AU=au,pdb){
    eleno.assign <- as.vector(as.data.frame(pdb$atom)$resno)
    uniqueeleno <- unique(eleno.assign) # this is for the lapply
    # count the no of atoms per residue
    atoms.per.res <- lapply(uniqueeleno,function(x){length(which(eleno.assign==x))}) # no per atom
    bfactors <- atoms.per.res
    res <- c()
    for (i in 1:length(au)){
        k <- rep(au[i],unlist(bfactors[i]))
        res <- c(res,k)
        }
    res
    }
# for WT
WT.traj <- read.ncdf("WT_cat.nc") # read in trajectory
WT.reference <- read.pdb("averagestructures/WTaverage.pdb") # get average structure for alignment
WT.calpha <- atom.select(WT.reference,elety="CA")$xyz # get the CA coordinates index
WT.traj <- WT.traj[,WT.calpha] # keeps only CA coordinates in memory
WT.reference <- WT.reference$xyz[WT.calpha] # get reference CA coordinate
WT.traj <- fit.xyz(fixed=WT.reference,mobile=WT.traj) # performs fitting using all residue CA
WT.pca <- pca.xyz(WT.traj) # performs the PCA
#png("pca/WTpca.png")
plot(WT.pca,col=c(rep("red",8000),rep("blue",8000),rep("green",8000)),main="Principal component analysis (WT)") # plot PCA summary using base graphics
#dev.off()
WT.au <- as.data.frame(WT.pca$au) # AU is the respective contribution to each PC by each residue
colnames(WT.au) <- c(paste0("PC",1:597)) # renames column 
WT.au$Residue <- 91:289 # adds residue information
WT.au <- melt(WT.au,"Residue") # cast into long format
colnames(WT.au) <- c("Residue","Princ.Comp","AU") # renames column
WT.au <- WT.au[WT.au$Princ.Comp%in%c("PC1","PC2","PC3","PC4","PC5"),] # keeps only the top 5 PC
ggplot(WT.au,aes(x=Residue,y=AU))+geom_line()+facet_grid(~Princ.Comp)+theme_bw()+ylab("Contribution to principal component")+scale_x_continuous(breaks=seq(90,290,20))+theme(axis.text.x=element_text(angle=90))+ggtitle("Contribution of residues to PC (WT)") #plotting in ggplot with facet grid
#ggsave("pca/WT_AU.png")
# apphending au as bfactors in PDB for pymol visualization
WTpdb <- read.pdb("WT_structure.pdb")
PCs <- c("PC1","PC2","PC3","PC4","PC5")
# doing the generation of the file in a loop -- WT
for (k in 1:length(PCs)){
    cat("Doing for ",PCs[k],"\n")
    au <- WT.au[WT.au$Princ.Comp==PCs[k],]
    au <- au$AU
    au <- c(0,au,0,0)
    au <- assignAU(AU=au,pdb=WTpdb)*10
    filename <- paste0("AU/WT_AU_",PCs[k],".pdb")
    write.pdb(WTpdb,file=filename,b=au)
    cat("Saved as ",filename,"\n")
    }

# for mutant
Mut.traj <- read.ncdf("Mut_cat.nc")
Mut.reference <- read.pdb("averagestructures/Mutaverage.pdb")
Mut.calpha <- atom.select(Mut.reference,elety="CA")$xyz
Mut.traj <- Mut.traj[,Mut.calpha]
Mut.reference <- Mut.reference$xyz[Mut.calpha]
Mut.traj <- fit.xyz(fixed=Mut.reference,mobile=Mut.traj) # performs fitting using all residue CA
Mut.pca <- pca.xyz(Mut.traj) # performs the PCA
#png("pca/Mutpca.png")
plot(Mut.pca,col=c(rep("red",8000),rep("blue",8000),rep("green",8000)),main="Principal component analysis (Mut)")
#dev.off()
# to analyze residue wise contribution to PCs
Mut.au <- as.data.frame(Mut.pca$au)
colnames(Mut.au) <- c(paste0("PC",1:597))
Mut.au$Residue <- 91:289
Mut.au <- melt(Mut.au,"Residue")
colnames(Mut.au) <- c("Residue","Princ.Comp","AU")
Mut.au <- Mut.au[Mut.au$Princ.Comp%in%c("PC1","PC2","PC3","PC4","PC5"),]
ggplot(Mut.au,aes(x=Residue,y=AU))+geom_line()+facet_grid(~Princ.Comp)+theme_bw()+ylab("Contribution to principal component")+scale_x_continuous(breaks=seq(90,290,20))+theme(axis.text.x=element_text(angle=90))+ggtitle("Contribution of residues to PC (Mut)")
#ggsave("pca/Mut_AU.png")



Mutpdb <- read.pdb("Mut_structure.pdb")
for (k in 1:length(PCs)){
    cat("Doing for ",PCs[k],"\n")
    au <- Mut.au[Mut.au$Princ.Comp==PCs[k],]
    au <- au$AU
    au <- c(0,au,0,0)
    au <- assignAU(AU=au,pdb=Mutpdb)*10
    filename <- paste0("AU/Mut_AU_",PCs[k],".pdb")
    write.pdb(Mutpdb,file=filename,b=au)
    cat("Saved as ",filename,"\n")
    }


