##A model of one bacterial species growing on 1 substrate, which assumes that:
##the per capita growth rate is proportional to the rate of metabolism
##the rate of metabolism is a Michaelis-Menten function of the concentration of substrate
##constant rate of dilution introduces new substrate, removes substrate and bacteria

##assign values to parameters
##a) starting population density of bacteria
N0<-0.1
##b) starting concentration of the substrate (resource)
S0<-3
##c) the dilution rate
d<-0.01

##c) 3 parameters for the Michaelis-Menten part, which specifies
##the rate at which substrate is metabolized, see http://en.wikipedia.org/wiki/Michaelis–Menten_kinetics

##the amount of substrate that a unit amount of enzyme can metabolise per unit time
kcat<-0.1
##the amount of enzyme per cell
E<-1
##NB Vmax = kcat*E, where Vmax is the rate of the reaction per cell when there is excess substrate

##the Michaelis constant - concentration at which the reaction rate is 1/2 its maximum
Km<-1

##the conversion rate from substrate metabolised to density increase
##i.e. per unit substrate metabolise get c.conv more cells produced
c.conv<-1

##choose time steps and total time
time.step<-0.5
times<-seq(0,500,time.step)

##set initial values of state variables
S<-S0
N<-N0

##set up matrix to store results in
results<-matrix(NA,nrow=length(times)+1,ncol=3)

	##store time 0 values
	results[1,1]<-0
	results[1,2]<-S
	results[1,3]<-N

##loop through the times
for (i in (1:length(times))) {
	
	##assign new values to S and N according to Michaelis Menten
	S <- S-(kcat*E*S/(Km+S))*N*time.step-d*S*time.step+d*S0*time.step
	N <- N+(c.conv*kcat*E*S/(Km+S))*N*time.step-d*N*time.step
	
	##store the results
	##use i+1 so as not to over-write the time 0 values
	results[i+1,1]<-times[i]
	results[i+1,2]<-S
	results[i+1,3]<-N
	
}

##plot the results on the same graph
##a) Substrate concentration
plot(results[,2]~results[,1],type="l")
##b) Density
lines(results[,3]~results[,1],col="red")