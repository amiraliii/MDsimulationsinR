# code to calculate transition frequency
# requires definition of 2 states using a cutoff condition 
# we want to calculate the frequency in which we note that there is a state change from State A to State B

calculateTransitionFrequency <- function(x,average){
    # use statistics from the WT tranjectory to get the range/deviation
    framesState2 <- which(x<average) # frames falling outside the normal range 
    framesState1 <- 1:length(x)[!framesState2]
    framesT <- framesState1+1
    frequencyTransition <- length(framesT%in%framesState2)/length(framesState1)
    return(frequencyTransition)
    } # returns a transition frequency 
