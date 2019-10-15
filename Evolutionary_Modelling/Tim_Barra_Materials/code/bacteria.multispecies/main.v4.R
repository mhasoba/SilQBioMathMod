library(deSolve)

RateMetabolismSaturatingMM <- function(kcat, Km, S, E){
  #multiplying each row of kcat by vector S
  kcat.S <- sweep(kcat,1,S,"*")  
  ##raise to power of the cost of generalism term - added power term May 2014
  numerator <- (E^cost.generalism)*kcat.S         
  #adding vector S to each row of matrix Km
  KmplusS <- sweep(Km,1,S,"+")                    
    
  #element by element division of both arrays
  j<-numerator/KmplusS 
  
  j[lower.tri(j[,,1])] <- 0
  
  return(j)                  
}

FitnessFunction<-function(i,S,E) {
	
	##michaelis.menten bit
	michaelis.menten<-(sweep(kcat[,,i],1,S,"*")/sweep(Km[,,i],1,S,"+"))
	michaelis.menten<-michaelis.menten*natp[,,i]    ##  ADDED MAY 17th TO CORRECT MISTAKE
	michaelis.menten[!pathway[,,i]]<-0
	
	#tmp<-michaelis.menten
	#tmp[pathway[,,i]]<-michaelis.menten[pathway[,,i]]-(sum(michaelis.menten[pathway[,,i]])-michaelis.menten[pathway[,,i]])/(sum(pathway[,,i])-1)	
	
	tmp<-michaelis.menten
	E.term<-(cost.generalism*E[,,i]^(cost.generalism-1))[pathway[,,i]]
	tmp[pathway[,,i]]<-michaelis.menten[pathway[,,i]]*E.term -(sum(michaelis.menten[pathway[,,i]]*E.term)-E.term*michaelis.menten[pathway[,,i]])/(sum(pathway[,,i])-1)	

	
	query<-(tmp*kMutRate[,,i]+E[,,i])<0
	if (sum(query)>0) {
	query2<-!query&pathway[,,i]
	##add up the changes that would cause negative E, make those enzymes change just to 0, not negative
	##allocate the rest to the other enzymes in proportion to their frequency (because those with larger E better able to absorb)
	tmp[query2]<-tmp[query2]+E[,,i][query2]/sum(E[,,i][query2])*sum(tmp[query]+sum(E[,,i][query])/kMutRate[,,i][query])
	tmp[query]<--E[,,i][query]/kMutRate[,,i][query]
}
		
    return(tmp)

}

GrowthRate <- function(paramsGrowthInt){  

  attach(paramsGrowthInt, warn.conflicts = FALSE)
  w <- rep(0, n)  

  for (sp in seq(1, n)){
    w[sp]=sum(natp[,,sp]*(J[,,sp])) 
   }
   detach(paramsGrowthInt)
   return(w)
}



        
ODE <-function(times, y, parms){
	
  attach(parms, warn.conflicts = FALSE)
   
  S <- y[1:k]  
  dS.dt <- array(0,k) 
  N <-y[(k+1):(k+n)]
  dN.dt <- array(0,n)
  
    #if (times==3000) {S<-S+0.05; print("boo")}

   
   ##changed Feb2014 to make dE.dt 0 by default
  dE.dt <- array(data = 0, dim = c(k,k,n))
  E<-dE.dt
  E[pathway]<-y[(k+n+1):(k+n+sum(pathway))]
  #E <- dE.dt
  J <- RateMetabolismSaturatingMM(kcat, Km, S, E)

  
  
 ##i) growth of bacteria
  paramsGrowthInt <- list(k = k, n = n, natp = natp, J = J, E = E, S = S)
  W <- GrowthRate(paramsGrowthInt)   
  dN.dt <- W*N - dd*N
  
  
 ##ii) substrates
  J.N <-  array(data = 1, dim = c(k,k,n))
  for (species in seq(0, n)){
    J.N[,,species] <- J[,,species]*N[species]
  }

  for (subst in seq(1, k)){
    dS.dt[subst]<- D*(Q[subst] - S[subst]) - sum(J.N[subst,,]) +   
                      sum(J.N[,subst,]) - sum(J.N[,subst,][subst,])                                                
  }  
  
  ##iii) Enzymes - added Feb 2014 to skip these bits if don't want evolution
  if (evolve) {
   
  for (species in seq(1, n)) {
  
        dW.dE <- FitnessFunction(species,S,E)  
          
        #multiplying by the mutation rate
        dE.dt[,,species] <- m[,,species]*dW.dE 
  }
  }
 
  #no negative concentration of substrates    
  dS.dt[dS.dt + S < 0] <- -S[dS.dt + S < 0]
  #no negative densities of bacteria
  dN.dt[dN.dt + N < 0.000] <- -N[dN.dt + N < 0.000]
   
  detach(parms)
  
  return(list(c(dS.dt,dN.dt,dE.dt[pathway])))  
  #return(list(c(dS.dt,dN.dt)))
}

main <- function(verbose=FALSE){

  #Initial values for the differential equation solver
  y <- c(S = conc.subst, N = density.spec, E = metabolic[pathway])
  #y <- c(S = conc.subst, N = density.spec)
  #E <<- metabolic
  
  
  #parameters to be passed to lsoda
  parms <- list(k = kNumSubst, n = kNumSpec, kcat = kCat, 
            D = kSubstDilutionRate, m = kMutRate, Km = kMichaelis, Q = kInflowConc,
             natp = kAtpPerReaction, dd = kSpecDilutionRate,pathway=pathway)
 
  sol <- lsoda(y,times,ODE,parms,verbose = verbose)
 
  return(list(sol=sol,parms=parms))
}


parameters="params.R"
