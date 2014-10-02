# code to plot RMSF by trajectory, showing the error bars on ggplot

# reads in all the RMSF data
WT1RMSF <- read.table("WT1RMSF.out")
WT2RMSF <- read.table("WT2RMSF.out")
WT3RMSF <- read.table("WT3RMSF.out")
Mut1RMSF <- read.table("Mut1RMSF.out")
Mut2RMSF <- read.table("Mut2RMSF.out")
Mut3RMSF <- read.table("Mut3RMSF.out")
# bind the data into a single data frame and then rename columns for easy identification
RMSFs <- rbind(WT1RMSF,WT2RMSF,WT3RMSF,Mut1RMSF,Mut2RMSF,Mut3RMSF)
colnames(RMSFs) <- c("Residue","RMSF")

# renumbers the residues
RMSFs$Residue <- RMSFs$Residue+89

# identifies the simulations
RMSFs$Simulation <- c(rep("WT1",199),rep("WT2",199),rep("WT3",199),rep("Mut1",199),rep("Mut2",199),rep("Mut3",199))
ggplot(RMSFs,aes(x=Residue,y=RMSF,color=Simulation))+geom_line()+theme_bw()+scale_color_manual(values=c("WT1"="blue","WT2"="blue","WT3"="blue","Mut1"="red","Mut2"="red","Mut3"="red"))

# this section was taken from somewhere online...
#
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    require(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}
#

# calculate the summary statistics for the data frame that we had initially
RMSFerror <- summarySE(RMSFs,"RMSF",c("Simulation2","Residue"))

# plotting the datat, using geom_line and geom_error bars
ggplot(RMSFerror,aes(x=Residue,y=RMSF,color=Simulation2))+geom_line()+theme_bw()+scale_color_manual(values=c("Mut"="red","WT"="blue"))+scale_x_continuous(breaks=seq(90,290,10))+scale_y_continuous(breaks=seq(0,8,1))+ylab("RMSF (\uc5)")+ggtitle("RMSF of backbone CA (\uc5)")+geom_errorbar(aes(ymin=RMSF-1.96*se,ymax=RMSF+1.96*se))
