# function to calculate the cosine of the 2 vectors


cosineFx <- function(a1CA,a1CB,a2CA,a2CB){
    # 4 inputs are the coordinatess of the respective atoms
    # defining the vectors
    v1 <- a1CB-a1CA
    v2 <- a2CB-a2CA
    # dot product in 2D space- x1x2+y1y2+z1z2, ie the sum of the products
    dot <- apply(v1*v2,1,sum) # this is the value of the dot product
    # to get the cosine we divide the dot product by the product of the magnitudes of the vectors
    magA <- apply(v1,1,function(x){sqrt(sum(x^2))})
    magB <- apply(v2,1,function(x){sqrt(sum(x^2))})
    theta <- magA*magB
    cosv <- dot/theta # this is the cosine
    # use acos to return the angle between the 2 vectors
    angle <- acos(cosv) # returns the angle in radian
    angle <- (180*angle)/pi # converts angle to radian
    angle # returns the angle
    }

# this is a testing dataset
f <- a <- data.frame(x=c(0,0),y=c(0,0),z=c(0,0))
g <- data.frame(x=c(0,0),y=c(1,1),z=c(0,0))
b <- data.frame(x=c(1,1),y=c(0,0),z=c(0,0))
cosineFx(a,b,f,g) # should return (90,90) since both vectors are perpendicular to each other

# END

# read in the data
mutCA1 <- read.table("Mut1_Q248CA.pdb")
mutCA1 <- mutCA1[2001:10000,]
mutCA2 <- read.table("Mut2_Q248CA.pdb")
mutCA2 <- mutCA2[2001:10000,]
mutCA2 <- rbind(mutCA1,mutCA2)
mutCA2<- mutCA2[,6:8] # Q248CA in 20000 frames
mutCD1 <- read.table("Mut1_Q248CD.pdb")
mutCD1 <- mutCD1[2001:10000,]
mutCD2 <- read.table("Mut2_Q248CD.pdb")
mutCD2 <- mutCD2[2001:10000,]
mutCD2 <- rbind(mutCD1,mutCD2)
mutCD2 <- mutCD2[,6:8] #Q248CD in 20000 frames
rm(mutCA1,mutCD1)

WTCA1 <- read.table("WT2_R248CA.pdb")
WTCA1 <- WTCA1[2001:10000,]
WTCA2 <- read.table("WT3_R248CA.pdb")
WTCA2 <- WTCA2[2001:10000,]
wtCA2 <- rbind(WTCA1,WTCA2)
wtCA2 <- wtCA2[,6:8]
WTCZ1 <- read.table("WT2_R248CZ.pdb")
WTCZ1 <- WTCZ1[2001:10000,]
WTCZ2 <- read.table("WT3_R248CZ.pdb")
WTCZ2 <- WTCZ2[2001:10000,]
wtCZ2 <- rbind(WTCZ1,WTCZ2)
wtCZ2 <- wtCZ2[,6:8] #
rm(WTCA1,WTCA2,WTCZ1,WTCZ2)

# N247 data
WTCA1 <- read.table("WT2_N247CA.pdb")
WTCA1 <- WTCA1[2001:10000,]
WTCA2 <- read.table("WT3_N247CA.pdb")
WTCA2 <- WTCA2[2001:10000,]
wtCA1 <- rbind(WTCA1,WTCA2)
wtCA1 <- wtCA1[,6:8]
rm(WTCA1,WTCA2) #WT N247CA

MutCA1 <- read.table("Mut1_N247CA.pdb")
MutCA1 <- MutCA1[2001:10000,]
MutCA2 <- read.table("Mut2_N247CA.pdb")
MutCA2 <- MutCA2[2001:10000,]
mutCA1 <- rbind(MutCA1,MutCA2)
mutCA1 <- mutCA1[,6:8]
rm(MutCA1,MutCA2)

WTCG1 <- read.table("WT2_N247CG.pdb")
WTCG1 <- WTCG1[2001:10000,]
WTCG2 <- read.table("WT3_N247CG.pdb")
WTCG2 <- WTCG2[2001:10000,]
MutCG1 <- read.table("Mut1_N247CG.pdb")
MutCG1 <- MutCG1[2001:10000,]
MutCG2 <- read.table("Mut2_N247CG.pdb")
MutCG2 <- MutCG2[2001:10000,]
wtCG1 <- rbind(WTCG1,WTCG2)
wtCG1 <- wtCG1[,6:8]
mutCG1 <- rbind(MutCG1,MutCG2)
mutCG1 <- mutCG1[,6:8]
rm(WTCG1,WTCG2,MutCG1,MutCG2)
# to calculate the actual values
WTangles <- cosineFx(wtCA1,wtCG1,wtCA2,wtCZ2)
Mutangles <- cosineFx(mutCA1,mutCG1,mutCA2,mutCD2)

dataX <- data.frame(Time=rep(2001:10000*0.05,4),Angle=c(as.numeric(WTangles),as.numeric(Mutangles)),Simulation=c(rep("WT",16000),rep("Mut",16000)))
ggplot(dataX,aes(x=Angle,color=Simulation))+geom_density()+theme_bw()+scale_x_continuous(breaks=seq(0,180,30))+ggtitle("Density of angle between N247 and R/Q248 side chains")+xlab(expression(paste("Angle (",degree,")")))+theme(axis.text.x=element_text(angle=45))
ggsave("dotProduct_N247_248.png")
