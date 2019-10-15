##assign values to parameters, species 1 and species 2 
N0<-c(0.1,0.1)
S0<-c(5,3)
d<-0.01
kcat<-matrix(c(0.1,0.1,0.1,0.1),nrow=2)
Km<-matrix(c(0.5,0.5,0.5,0.5),nrow=2)
c.conv<-matrix(c(1,1,1,1),nrow=2)
##set up so can have different parameter for each species X substrate
##column = species, row = substrate
E<-matrix(c(1,0,0,1),nrow=2,ncol=2)

##choose time steps and total time
time.step<-0.1
times<-seq(0,100,time.step)

##set initial values of state variables
S<-S0
N<-N0

##store results
results<-matrix(NA,nrow=length(times)+1,ncol=5)

	##store time 0 values
	results[1,1]<-0
	results[1,2:3]<-S
	results[1,4:5]<-N

##loop through the times
for (i in (1:length(times))) {
	
	##assign new values - now S1 might be used by both species, each species might grow on both substrates
	S[1] <- S[1]-(kcat[1,1]*E[1,1]*S[1]/(Km[1,1]+S[1]))*N[1]*time.step-(kcat[1,2]*E[1,2]*S[1]/(Km[1,2]+S[1]))*N[2]*time.step-d*S[1]*time.step+d*S0[1]*time.step
	S[2] <- S[2]-(kcat[2,1]*E[2,1]*S[2]/(Km[2,1]+S[2]))*N[1]*time.step-(kcat[2,2]*E[2,2]*S[2]/(Km[2,2]+S[2]))*N[2]*time.step-d*S[2]*time.step+d*S0[2]*time.step
	N[1] <- N[1]+(c.conv[1,1]*(kcat[1,1]*E[1,1]*S[1]/(Km[1,1]+S[1]))*N[1]*time.step)+(c.conv[2,1]*(kcat[2,1]*E[2,1]*S[2]/(Km[2,1]+S[2]))*N[1]*time.step)-d*N[1]*time.step
	N[2] <- N[2]+(c.conv[1,2]*(kcat[1,2]*E[1,2]*S[1]/(Km[1,2]+S[1]))*N[2]*time.step)+(c.conv[2,2]*(kcat[2,2]*E[2,2]*S[2]/(Km[2,2]+S[2]))*N[2]*time.step)-d*N[2]*time.step

	results[i+1,1]<-times[i]
	results[i+1,2:3]<-S
	results[i+1,4:5]<-N
	
}

##plot the results
##a) Densities
matplot(results[,1], results[,2:5], type="l", xlab="Time", 
        ylab="Density/Concentration",col=c("black","red","black","red"),
        lty=c(1,1,2,2))
legend(x=80,y=2,legend=c("S1","S2","N1","N2"), col=c("black","red","black","red"),lty=c(1,1,2,2))
