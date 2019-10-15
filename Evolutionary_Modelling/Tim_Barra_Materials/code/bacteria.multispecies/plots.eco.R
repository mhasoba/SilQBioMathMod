#main plots of enzymes, species and substrates

main.plot <- function(sol){

palette(c("red","orange","yellow","green","blue","purple"))

times <- sol[,1]

  t <- length(times)
    #Separating the vectors to simplify the indexes
  #order of sol: times, S, N, E (just a reminder)
  solS <- sol[(t+1):(t*(kNumSubst + 1))]
  solN <- sol[(t * (kNumSubst+1) + 1): (((t * (kNumSubst+1) + 1)) + t*(kNumSpec)
              - 1)]
  solE <- sol[(((t * (kNumSubst+1) + 1)) + t*(kNumSpec)):length(sol)]
  
  
  wd <- getwd()
  dat <- date() 
  dat <- gsub(":",".",dat,fixed=T)
  dir.create(dat)
  setwd(dat)
  

  
  ##################Substrates##################################################
  pdf("SubstancesConcentration.pdf")
  plot(rep(times,kNumSubst), solS, type = 'n', xlab = 
       'Time', ylab = 'Substance Concentration',ylim=c(0,max(solS)))
      
  for(s in seq(1, kNumSubst))
    points(times, solS[(1+(s-1)*t):(s*t)], type = 'l', col = s)
    
  dev.off()
  

  ###################end of substrates##########################################  
  
  #####################Species##################################################
  pdf("SpeciesConcentration.pdf")
  plot(rep(times,kNumSpec), solN, type = 'n', xlab = 
       'Time', ylab = 'Species Concentration',ylim=c(0,max(solN)))
      
  for(n in seq(1, kNumSpec))
    points(times, solN[(1+(n-1)*t):(n*t)], type = 'l', col = n)

  
  dev.off()  
  
  ###################end of species#############################################  
  
  #####################Enzymes##################################################

 #par(mfrow=c(kNumSpec,1))
 
  species<-rep(1:kNumSpec,each=kNumSubst*kNumSubst)[pathway]
  enzymes<-rep(rep(1:kNumSubst,kNumSubst),kNumSpec)[pathway]
  enzyme.prod<-rep(rep(1:kNumSubst,each=kNumSubst),kNumSpec)[pathway] 

for (n in 1:kNumSpec) {

  pdf(paste("Enzymes.Species.",n,".pdf",sep=""))

  plot(rep(times,sum(pathway)), solE, type = 'n', xlab = 
       'Time', ylab = paste("Enzyme Conc Species",n,sep=" "))
  
 # for(n in seq(1, kNumSpec))
    for (e in which(species==n)){
      points(times, results$sol[,kNumSubst+kNumSpec+e+1], type = 'l', lty = 1, 
             col = enzymes[e])
    
    #text(x=max(times),y=results$sol[length(results$sol[,1]),kNumSubst+kNumSpec+e+1],labels=paste(enzymes[e],enzyme.prod[e],sep="->"))
    
   }
    dev.off() 
   }

  ###################end of enzymes#############################################
  
  #pdf("Legend.pdf")
  #plot(1,type="n",yaxt = "n",xaxt = "n",bty = "n", xlab = "",ylab = "")
  #reactions <- combn(kNumSubst,2)
  #legend("center", col = c(seq(1,kNumSpec), rep(1,length(reactions[1,]))),
  #        lty = c(rep(1,kNumSpec),seq(1,length(reactions[1,]))), legend = 
  #        c(paste("Species",seq(1,kNumSpec)), paste(reactions[1,],"→ ",reactions[2,])))
  #legend("right", lty = seq(1,n), legend = paste("Species",seq(1,n)))
  #legend("left", lty = seq(1,s), legend = paste("Substance",seq(1,s)))       
  #dev.off()
  
  
  
  sink(file="parameters",append=TRUE)
  print("kNumSubst")
  print(kNumSubst)  
  print("kNumSpec")
  print(kNumSpec)
  print("kCat")
  print(kCat)
  print("kSubstDilutionRate")
  print(kSubstDilutionRate)
  print("kMutRate")
  print(kMutRate)
  print("kMichaelis")
  print(kMichaelis)
  print("kInflowConc")
  print(kInflowConc)
  print("kAtpPerReaction")
  print(kAtpPerReaction)
   print("metabolic")
  print(metabolic)
  
  sink()
  
  setwd(wd)

} 
  
  
  

substance.plot <- function(sol, times){

  t <- length(times)
  
  #Separating the vectors to simplify the indexes
  #order of sol: times, S, N, E (just a reminder)
  solS <- sol[(t+1):(t*(kNumSubst + 1))]
  solN <- sol[(t * (kNumSubst+1) + 1): (((t * (kNumSubst+1) + 1)) + t*(kNumSpec)
              - 1)]
  solE <- sol[(((t * (kNumSubst+1) + 1)) + t*(kNumSpec)):length(sol)]
  
  #This will be a multiple plot of each substrate with the mean concentrations
  #of all enzymes that metabolise it and consume it
  pdf("SubstratePlot.pdf", quality = 100)
  par(mfrow = c(ceiling(kNumSubst/2),2))
  
  cons <- 0
  prod <- 0
      
  for(s in seq(1, kNumSubst)){
    ##plotting the substance fist, in solid black
    plot(times, solS[(1+(s-1)*t):(s*t)], main = paste("Substrate",s), type = 'l', 
         xlab = 'Time', ylab = 'Concentration', col = 1, lty = 1)
    cons <- solE[(1+t*(s-1)) : 10]
  }       
    

  dev.off()
  

  ###################end of substrates##########################################  

  
  ###################end of enzymes#############################################
  setwd(wd)

} 
  
