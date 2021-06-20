filepath="~/Validation/Outputs"
setwd(filepath)

install.packages("raster",dependencies = TRUE, INSTALL_opts = '--no-lock')
install.packages("sp",dependencies = TRUE, INSTALL_opts = '--no-lock')

library(raster)
library(sp)


mat.torus <- function(Matrix,Total.Neighbors,xcord,ycord){      # Torus of Space
  
  if(Total.Neighbors==0){return(Matrix[xcord,ycord])}
  
  dm    <- nrow(Matrix)                               # dimension of matrix (n x n marix, so rows = col)
  
  Crown.Pot       <-  c(seq(1,15,by=2)^2)-1    # arbitrarily set to 10, because that many neighbors gets silly biologically, just used for computation 
  Crowns          <-  seq(0,length(Crown.Pot))
  Crown           <-  Crowns[min(which(Crown.Pot >= Total.Neighbors))]  #returns crown exension total (even if not whole crown)
  Sample          <-  Total.Neighbors - (Crown.Pot[Crown])   # Returns how many to sample from the extra plots in addtion to crown
  
  rows  <- c(  (xcord-Crown):(xcord+Crown)       )    # figure out which rows of matrix crown falls in
  cols  <- c(  (ycord-Crown):(ycord+Crown)       )    # figure out which columns of matrix crown falls in
  
  rows[which(rows<1)]    <- rows[which(rows<1)] + dm  # if crown extends to a row value less than 1, go to opposite side of torus 
  rows[which(rows>dm)]   <- rows[which(rows>dm)] - dm # if crown extends to a row value greater than dm, go to opposite side of torus
  
  cols[which(cols<1)]    <- cols[which(cols<1)] + dm  # if crown extends to a column value less than 1, go to opposite side of torus 
  cols[which(cols>dm)]   <- cols[which(cols>dm)] - dm # if crown extends to a column value greater than dm, go to opposite side of torus
  
  JC_Matrix              <- Matrix[rows,cols ]        # returns subset of matrix / trees that are in JCE zone + extras 
  JC_Matrixb             <- JC_Matrix                # Copy matrix 
  
  JC_Matrix[c(1,nrow(JC_Matrix)),]    <- 0             
  JC_Matrix[,c(1,ncol(JC_Matrix))]    <- 0 
  
  vec.JC_Matrix <- c(JC_Matrix)
  vec.JC_Matrix <- vec.JC_Matrix[vec.JC_Matrix>0]
  
  
  JC_Matrixb[2:(nrow(JC_Matrix)-1),2:(ncol(JC_Matrix)-1)]    <- 0 
  
  vec.JC_Matrixb <- c(JC_Matrixb)
  vec.JC_Matrixb <- JC_Matrixb[JC_Matrixb>0]
  vec.JC_Matrixb <- sample(vec.JC_Matrixb,Sample, replace = FALSE)
  
  return(unique(c(vec.JC_Matrix,vec.JC_Matrixb )  ))
  
}
'%!in%' <- function(x,y)!('%in%'(x,y))



set.seed(181)
dm <- 275
S <- 300
S.list <- seq(1:S)
TimeSteps <- 30000 
r2 <- raster(xmn = 0, xmx = dm, ymn = 0, ymx = dm, nrows =dm, ncols = dm)
dd <- sample(1:S, dm*dm, replace = TRUE)
r2[] <- dd


Mat.S <- matrix(dd,nrow=dm,ncol=dm)
df.Props                <- data.frame( matrix(NA,ncol=S+1,nrow=(1+TimeSteps) ))
df.Props[,1]            <- seq(1:(TimeSteps+1))
df.Props[1,2:(S+1)]     <- c(table(Mat.S))/(dm*dm)

A <-  seq(1,1,length=S)
A <- exp(-A)

set.seed(150)
Y <- rlnorm(S,mean=0,sd=.45)
names(Y) <- seq(1:S)

d <- 1
Dist.Rate <- .0025
Total.Neighbors <- 49
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
    
    JC.Victims    <- mat.torus(Mat.S2,Total.Neighbors,x.val[x],y.val[x]) # Select Distance dependnt victims on Torus
    JC.Victims_No <- c(unique(Mat.S2[c(Mat.S2 %!in% JC.Victims)])) # Select Species not affected by distance dependnce
    
    JC.Victims    <- JC.Victims[JC.Victims %!in% Local_Species ] # Remove distrubed species from JC victims (bookeeping) 
    
    
    P_L   <- P[Local_Species] # Proportion of local species in population (disturbed species)
    A_L   <- A[Local_Species] # Seed rain surivial Probs 
    Y_L   <- Y[Local_Species]
    
    P.In  <- P[JC.Victims]    # Proportion of each species that are affecte by JCEs
    P.Out <- P[JC.Victims_No] # Proprotion of each species thare are not affected by JCEs
    
    A.In  <- A[JC.Victims]    # Seed rain survival probms  
    
    Y.In   <- Y[JC.Victims]
    Y.Out  <- Y[JC.Victims_No]
    
    Local.seeds         <-  Y_L*(A_L)*(  (1-d) + d*P_L )   # scaled number of seeds of local species
    NLocal.seeds_NO     <-  Y.Out*d*P.Out # scaled number of seeds of each not JCE affected non-local species
    NLocal.seeds_YES    <-  Y.In*(A.In)*d*P.In  # scaled number of seeds of each JCE affected non-local species
    
    Total.Seeds         <- as.numeric(Local.seeds + sum(NLocal.seeds_NO) + sum(NLocal.seeds_YES) ) # total scaled number of seeds in local patch
    Vec.Probs           <- c(Local.seeds,NLocal.seeds_NO,NLocal.seeds_YES)/Total.Seeds  # probability that each species wins lottery         
    
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


write.csv(df.PropsM,"TS_E49_A1_Y45.csv",quote=F,row.names=F)
write.csv(Mat.S,"DIST_E49_A1_Y45.csv",quote=F,row.names=F)





