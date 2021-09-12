filepath="~/Validation/Outputs"
setwd(filepath)

mat.torus <- function(Matrix,Rad,xcord,ycord){      # Torus of Space
  
  dm    <- nrow(Matrix)                               # dimension of matrix (n x n marix, so rows = col)
  
  Crown.Pot       <-  c(seq(1,Rad,by=2)^2)-1    # arbitrarily set to 10, because that many neighbors gets silly biologically, just used for computation 
  Crowns          <-  seq(0,length(Crown.Pot))
  Crown           <-  Crowns[min(which(Crown.Pot >= (Rad-1)^2))]  #returns crown exension total (even if not whole crown)
  
  rows  <- c(  (xcord-Crown):(xcord+Crown)       )    # figure out which rows of matrix crown falls in
  cols  <- c(  (ycord-Crown):(ycord+Crown)       )    # figure out which columns of matrix crown falls in
  
  rows[which(rows<1)]    <- rows[which(rows<1)] + dm  # if crown extends to a row value less than 1, go to opposite side of torus 
  rows[which(rows>dm)]   <- rows[which(rows>dm)] - dm # if crown extends to a row value greater than dm, go to opposite side of torus
  
  cols[which(cols<1)]    <- cols[which(cols<1)] + dm  # if crown extends to a column value less than 1, go to opposite side of torus 
  cols[which(cols>dm)]   <- cols[which(cols>dm)] - dm # if crown extends to a column value greater than dm, go to opposite side of torus
  
  JC_Matrix              <- Matrix[rows,cols ]        # returns subset of matrix / trees that are in JCE zone + extras 
  
  return(JC_Matrix)
  
}
'%!in%' <- function(x,y)!('%in%'(x,y))



set.seed(181)
dm <- 275
S <- 300
S.list <- seq(1:S)
TimeSteps <- 60000 
dd <- sample(1:S, dm*dm, replace = TRUE)

Mat.S <- matrix(dd,nrow=dm,ncol=dm)
df.Props                <- data.frame( matrix(NA,ncol=S+1,nrow=(1+TimeSteps) ))
df.Props[,1]            <- seq(1:(TimeSteps+1))
df.Props[1,2:(S+1)]     <- c(table(Mat.S))/(dm*dm)

A <-  seq(4.5,4.5,length=S)
#A <- exp(-A)

set.seed(150)
Y <- rlnorm(S,mean=0,sd=.1)
names(Y) <- seq(1:S)

d <- 1
Dist.Rate <- .0025
Rad <- 7
set.seed(150)


for(mm in 1:TimeSteps){
  
  Mat.S2 <- Mat.S  
  P      <- c(table(factor(Mat.S2, levels = 1:S)))/(dm*dm)   # P goes to proportion of each species in environment
  
  Prob.Dist                             <- matrix(runif(dm*dm),ncol=dm) # Matrix that defines hte probability that each species is disurbed
  Prob.Dist[Prob.Dist >= Dist.Rate]     <- NA # Spcies with draw less than Dist.Rate are distriubed 
  
  df.Rep                      <- which(!is.na(Prob.Dist), arr.ind=TRUE) # saves the indexes of each location that is distrubed
  
  x.val  <- df.Rep[,1]   # x coordinates of disturbance
  y.val  <- df.Rep[,2]   # y coordinates of distriubance
  
  Replaceb <- length(x.val)  # total number of distrubances
  
  Replacements <- sapply(1:Replaceb, function(x){  # function that determines which speices replaces disturbed patches (apply function loops over all distrubed)
    
    Local_Species <- Mat.S2[x.val[x],y.val[x]]   # defines the locally disturbed species 
    
    
    Rad <- 7
    JC.Victims   <- mat.torus(Mat.S2,Rad,x.val[x],y.val[x])
    #JC.Victims_No <- c(unique(Mat.S2[c(Mat.S2 %!in% JC.Victims)])) # Select Species not affected by distance dependnce
    
    #JC.Victims    <- JC.Victims[JC.Victims %!in% Local_Species ] # Remove distrubed species from JC victims (bookeeping) 
    
    JC.Pred       <- c(table(factor(JC.Victims, levels = 1:S))) 
    Predation     <- exp(-A*JC.Pred)
    Seeds         <-  Y*Predation*P
    
    Total.Seeds         <- as.numeric(sum(Seeds)) # total scaled number of seeds in local patch
    Vec.Probs           <- c(Seeds)/Total.Seeds  # probability that each species wins lottery         
    
    Vec.Probs.Ordred    <- Vec.Probs[order(as.numeric(names(Vec.Probs)))]  # order previous probability values (1 to S)
    
    Vec.Sum             <- cumsum(Vec.Probs.Ordred) # creates probability intervals to determine which species wins 
    
    prob.rep <- runif(1) # draw from uniform distribution to determine the winner of the lottery
    
    Replacement <- as.numeric(names(Vec.Sum[min(which(Vec.Sum > prob.rep))])) # store winner of lottery
    return(Replacement) # return winner
  }
  )
  
  Mat.S[df.Rep] <- Replacements # put winner of lottery into correct location
  df.Props[mm+1,2:(S+1)]     <- c(table(factor(Mat.S, levels = 1:S)))/(dm*dm) # Store proportion of each species at each time step
  
} # Code that runs the simulation (notes inside)




df.PropsM <- as.matrix(df.Props)


write.csv(df.PropsM,"TS_E49_A45_Y1_AD.csv",quote=F,row.names=F)
write.csv(Mat.S,"DIST_E49_A45_Y1_AD.csv",quote=F,row.names=F)





