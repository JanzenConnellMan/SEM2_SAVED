filepath="~/Validation/Outputs"
setwd(filepath)

install.packages("raster",dependencies = TRUE, INSTALL_opts = '--no-lock')
install.packages("sp",dependencies = TRUE, INSTALL_opts = '--no-lock')

library(raster)
library(sp)


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

Distance.Func <- function(FC,Alpha,M,Yb,Xb){
  
  Ind      <- data.frame(which(M==FC, arr.ind=TRUE))
  return(Alpha*sum(exp(-R*( sqrt( (Ind$row  - Xb)^2 + (Ind$col  - Yb)^2 ) )*d_A)))
  
  
}

#Total.Neighbors <- 30

set.seed(181)
dm <- 275
S <- 300
S.list <- seq(1:S)
TimeSteps <- 25000
r2 <- raster(xmn = 0, xmx = dm, ymn = 0, ymx = dm, nrows =dm, ncols = dm)
dd <- sample(1:S, dm*dm, replace = TRUE)
r2[] <- dd

Mat.S <- matrix(dd,nrow=dm,ncol=dm)
df.Props                <- data.frame( matrix(NA,ncol=S+1,nrow=(1+TimeSteps) ))
df.Props[,1]            <- seq(1:(TimeSteps+1))
df.Props[1,2:(S+1)]     <- c(table(Mat.S))/(dm*dm)


A.Val <- .3
A1 <- seq(.5-A.Val,.5+A.Val,length=S)
A <- -log(A1)  


Y <- rep(1,S)
names(Y) <- seq(1:S)

d <- rep(.9,S)


Dist.Rate <- .0025
R <- .1
d_A <- 1/sqrt(.2)


set.seed(150)
for(mm in 1:TimeSteps){
  
  Mat.S2 <- Mat.S  
  P      <- c(table(factor(Mat.S2, levels = 1:S)))/(dm*dm)   # P goes to proportion of each species in environment
  
  Prob.Dist                             <- matrix(runif(dm*dm),ncol=dm) # Matrix that defines hte probability that each species is disturbed
  Prob.Dist[Prob.Dist >= Dist.Rate]     <- NA # Spcies with draw less than Dist.Rate are disturbed 
  
  df.Rep                      <- which(!is.na(Prob.Dist), arr.ind=TRUE) # saves the indexes of each location that is disturbed
  
  x.val  <- df.Rep[,1]   # x coordinates of disturbance
  y.val  <- df.Rep[,2]   # y coordinates of disturbance
  
  
  Replaceb <- length(x.val)  # total number of disturbances
  
  Replacements <- sapply(1:Replaceb, function(x){  # function that determines which species replaces disturbed patches (apply function loops over all disturbed)
    
    Local_Species  <- Mat.S2[x.val[x],y.val[x]]   # defines the locally disturbed species 
    NLocal_Species <- S.list[-c(Local_Species)]
    
    P_L   <- P[Local_Species]  # Proportion of local species in population (disturbed species)
    P_NL  <- P[NLocal_Species] # Proportion of non-local species in population 
    
    # Local.seeds      <-  Y[Local_Species]   * Distance.Func2(x.val[x],y.val[x],dm,Mat.S2,Local_Species,R_d,D_A,Rad2)  # Total local species pre-predation
    # NLocal.seeds     <-   Y[NLocal_Species] * sapply(1:length(P_NL),function(z) Distance.Func2(x.val[x],y.val[x],dm,Mat.S2,NLocal_Species[z],R_d,D_A,Rad2))
    
    Local.seeds      <-  Y[Local_Species]*(  (1-d[Local_Species]) + d[Local_Species]*P_L )   # Total local species pre-predation
    NLocal.seeds     <-  Y[NLocal_Species]*d[NLocal_Species]*P_NL # Total NLocal 
    
    #Distance.Func(X,Y,Dimension,Matrix,Focal_Species,Alpha,Decay,Dist_Avg)
    
    Rad <- 75
    M   <- mat.torus(Mat.S2,Rad,x.val[x],y.val[x])
    Xb   <- (Rad+1)/2
    Yb   <- (Rad+1)/2
    Sp.in <- sort(unique(c(M)))
    Predation <- rep(0,S)
    
    # Predation <- sapply(1:S,function(z) Distance.Func(x.val[x],y.val[x],dm,Mat.S2,z,A[z],R,D_A,Rad))
    
    Preds <-   sapply(Sp.in,function(z) Distance.Func(z,A[z],M,Yb,Xb))
    
    Predation[Sp.in] <- Preds
    
    # Local seeds after predation 
    Local.seeds <- Local.seeds  *exp(-1*Predation[Local_Species])
    # NLocal seeds after predation
    NLocal.seeds <- NLocal.seeds*exp(-1*Predation[NLocal_Species])
    
    # Seed probs  
    Total.Seeds         <- as.numeric(Local.seeds + sum(NLocal.seeds) ) # total scaled number of seeds in local patch
    Vec.Probs           <- c(Local.seeds,NLocal.seeds)/Total.Seeds  # probability that each species wins lottery         
    
    Vec.Probs.Ordred    <- Vec.Probs[order(as.numeric(names(Vec.Probs)))]  # order previous probability values (1 to S)
    
    Vec.Sum             <- cumsum(Vec.Probs.Ordred) # creates probability intervals to determine which species wins 
    
    prob.rep <- runif(1) # draw from uniform distribution to determine the winner of the lottery
    
    Replacement <- as.numeric(names(Vec.Sum[min(which(Vec.Sum > prob.rep))])) # store winner of lottery
    return(Replacement) # return winner
  }
  )
  
  Mat.S[df.Rep] <- Replacements  # put winner of lottery into correct location
  
  
  
  
  df.Props[mm+1,2:(S+1)]     <- c(table(factor(Mat.S, levels = 1:S)))/(dm*dm) # Store proportion of each species at each time step
  
} # Code that runs the simulation (notes inside)


df.PropsM <- as.matrix(df.Props)


write.csv(df.PropsM,"TS_v10_A3_d9.csv",quote=F,row.names=F)
write.csv(Mat.S,"DIST_v10_A3_d9.csv",quote=F,row.names=F)




