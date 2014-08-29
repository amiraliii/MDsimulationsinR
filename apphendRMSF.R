# code request from Maryam to add export a PDB file with apphended RMSF as b-factor

library(bio3d)
# read in the pdb file
pdbobject <- read.pdb("analysis/em_clustering/centroids/WT1cluster1.pdb")

# read in the trajectory
# if you are using a netcdf file you will need additional packages ncdf
trajectory <- read.ncdf("WT_500ns.nc")

ca.index <- atom.select(pdbobject,elety = "CA")$xyz # select only the CA
trajectory.xyz <- trajectory[,ca.index] # keeps only the CA data in the trajectory
reference.ca <- pdbobject$xyz[ca.index] # reference ca coordinates
trajectory.fitted <- fit.xyz(reference.ca,trajectory.xyz) # performs alignment using CA only
ca.rmsfs <- rmsf(trajectory.fitted) # calculates the RMSF values

# to return the pdb we need the b factor column to be as long as the xyz list in the pdb
# assume that rmsf is the same for all atoms of the same residue

assignBfactors <- function(rmsf,pdb){
    eleno.assign <- as.vector(as.data.frame(pdb$atom)$resno)
    uniqueeleno <- unique(eleno.assign) # this is for the lapply
    # count the no of atoms per residue
    atoms.per.res <- lapply(uniqueeleno,function(x){length(which(eleno.assign==x))}) # no per atom
    adjustrmsf <- c(0,rmsf,0,0) # add 0 for rmsf of ACE (res1),NME (201),Zn (202) # this is specific to system, but just place 0 where there are protein residues else you will get complains
    bfactors <- atoms.per.res
    res <- c()
    for (i in 1:length(adjustrmsf)){
        k <- rep(adjustrmsf[i],unlist(bfactors[i]))
        res <- c(res,k)
        }
    res
    }

bvalues <- assignBfactors(ca.rmsfs,pdbobject) # vector of b values
write.pdb(pdbobject,file="test.pdb",b=bvalues) # exports the pdb file in which we apphend the rmsf data
# when you read the file exported pdb file into pymol you can color according to b-factors
