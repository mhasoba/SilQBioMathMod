##METABOLIC ASSIGNMENT FUNCTIONS

##specialists: each substrate used by just one species: if multiple products, allocation from enzyme.allocation of enzyme to each
specialists<-function() {
	metabolic <-  array(data = 0, dim = c(kNumSubst,kNumSubst,kNumSpec))
	for (i in (which(density.spec>0))) {
	metabolic[i,,i]<-enzyme.allocation[i,]}
	return(metabolic) }

##generalists, allocate enzyme to all pathways using continuous broken stick. 	
complementary.generalists<-function() {
	metabolic <-  array(data = 0, dim = c(kNumSubst,kNumSubst,kNumSpec))
	tmp<-which(whole.pathway,arr.ind=T)
	for (i in (1:kNumReac)) {
		metabolic[tmp[i,1],tmp[i,2],]<-broken.stick(kNumSpec)*enzyme.allocation[tmp[i,1],tmp[i,2]]}
	for (i in (1:kNumSpec)) {
		metabolic[,,i]<-metabolic[,,i]/sum(metabolic[,,i])}
	return(metabolic)
}

##generalists, allocate enzyme to all pathways use discrete formulation so all sums work	
complementary.generalists.discrete<-function() {
	if (alloc=="even") {num.subdiv<-prod(unique(1/enzyme.allocation[whole.pathway]))} else {num.subdiv<-10}
	metabolic <-  array(data = 0, dim = c(kNumSubst,kNumSubst,kNumSpec))
	x<-rep(which(whole.pathway,arr.ind=T)[,1],num.subdiv*enzyme.allocation[whole.pathway])
	y<-rep(which(whole.pathway,arr.ind=T)[,2],num.subdiv*enzyme.allocation[whole.pathway])
	z<-sample(rep(1:kNumSpec,num.subdiv))
	for (i in (1:(num.subdiv*kNumSpec))) {
		metabolic[x[i],y[i],z[i]]<-metabolic[x[i],y[i],z[i]]+1/num.subdiv
		}
	return(metabolic)	
}

 
##enzyme allocation: determines whether, in specialist case, splitting pathways have products 1:1 or broken stick splits
allocation<-function(x) {

enzyme.allocation<-whole.pathway

if (x=="even") {
enzyme.allocation[rowSums(whole.pathway)>1,]<-enzyme.allocation[rowSums(whole.pathway)>1,]/rowSums(whole.pathway)[rowSums(whole.pathway)>1]
return(enzyme.allocation) }

if (x=="random") {
for (i in (which(rowSums(whole.pathway)>1))) {
	enzyme.allocation[i,whole.pathway[i,]] <- broken.stick(sum(whole.pathway[i,]))}
	return(enzyme.allocation)}
	
if (x=="random.discrete") {
for (i in (which(rowSums(whole.pathway)>1))) {
	enzyme.allocation[i,whole.pathway[i,]] <- discrete.broken.stick(sum(whole.pathway[i,]))}
	return(enzyme.allocation)}
}


##Calculates predicted substrate concentrations at steady-state with specialists
specialist.predict<-function () {

S<-array(NA,kNumSubst)
N<-array(NA,kNumSpec)

metabolic <- specialists() 

whole.pathway.present<-apply(pathway&metabolic,c(1,2),any)

##positive steady-state solutions for substrates - assign zero if an end-product: use just 1 value from columns with more than 1
tmp<-(kSubstDilutionRate*kMichaelis/(kCat*kAtpPerReaction-kSubstDilutionRate))
S<-rowSums(apply(tmp*metabolic,c(1,2),sum))

##carry substrate concentration through pathway
Y<-kInflowConc

for (i in (1:kNumSubst)) {

	product<-which(whole.pathway.present[i,])
	spec<-which(metabolic[i,product[1],]>0)

	##logical to say whether positive solution is stable or not
	if ((S[i]<=0)|(S[i]>Y[i])) {
		S[i]<-Y[i]
		Y[product]<-Y[product]+0  #i.e. no substrate makes it to next step
		N[spec]<-0} else {
		#if (length(product)==1) {
		#Y[product]<-Y[product]+Y[i]-S[i]   ##i.e. of substrate arriving Y[i], S[i] remains as substrate i, rest moves on to next step
		#N[spec]<-kAtpPerReaction[i,product,spec]*(Y[i]-S[i])
		#						} else {
		Y[product]<-Y[product]+(Y[i]-S[i])*metabolic[i,product,spec]
		N[spec]<-kAtpPerReaction[i,product[1],spec]*(Y[i]-S[i])
		#}  ##i.e. E1 goes to product 1 and E2 goes to product 2 etc.
		} 	
}

##fails to find spec when that species is missing, so reassign
#N[is.na(N)]<-0

return(c(S,N))

}


##Works out reaction names in right order
p.names<-function(txt) {
pastey<-function(x) {paste(txt,paste(x,collapse="."),sep="")}
return(apply(which(pathway[,,which(density.spec>0)[1]],arr.ind=T),1,pastey))}

e.names<-function(txt) {
pastey<-function(x) {paste(txt,paste(x,collapse="."),sep="")}
return(apply(which(pathway,arr.ind=T),1,pastey))}

##generate broken stick
broken.stick<-function(n) {
		x<-runif(n-1)
		return(diff(c(0,sort(x),1)))}

discrete.broken.stick<-function(n) {
		x<-sample(1:9,n-1)
		return(diff(c(0,sort(x),10))/10)}

plot.pathway<-function() {
	
plot(c(0,24),c(1,7),type="n")	

symbols(x=c(2,6,10,14,18),y=c(3,3,3,3,3),circles=c(1,1,1,1,1),add=T,inches=F,bg=colors()[c(86,86,86,542,542,542)])#bg=c("grey","grey","grey","white","white"))
symbols(x=c(2,6,10,14,18),y=c(5,5,5,5,5),circles=c(1,1,1,1,1),add=T,inches=F,bg=colors()[c(148,148,148,43,542)])#bg=c("grey","grey","grey","white","white"))

symbols(x=22,y=4,stars=matrix(c(1,0.5,1,0.5,1,0.5,1,0.5,1,0.5,1,0.5)+0.2,nrow=1),add=T,inches=F,bg=colors()[654])

##linear
arrows(x0=c(3,7,11,15)+0.3,x1=c(5,9,13,17)-0.3,y0=c(3,3,3,3,3),length=0.1)
arrows(x0=c(3,7,11,15)+0.3,x1=c(5,9,13,17)-0.3,y0=c(5,5,5,5,5),length=0.1)

##crosslinks
arrows(x0=c(11,15)+0.3,x1=c(13,17)-0.3,y0=c(5,5)-0.3,y1=c(3,3)+0.3,length=0.1)
arrows(x0=c(11)+0.3,x1=c(13)-0.3,y0=c(3)+0.3,y1=c(5)-0.3,length=0.1)
arrows(x0=c(19,19)+0.3,x1=c(21,21)-0.3,y0=c(3.3,4.7),y1=c(3.7,4.3),length=0.1)

##curvey ones
curve(2*cos(0.3*x+1)+4,add=T,from=15,to=22)
curve(-2*cos(0.3*x+1)+4,add=T,from=15,to=22)
arrows(x0=c(21.5,21.5),x1=c(22,22),y0=c(3.213835,4.786165), y1=c(3.49748,4.50252),length=0.1)
#text(x=c(2,6,10,14,18,2,6,10,14,18,22),y=c(5,5,5,5,5,3,3,3,3,3,4),labels=c(1,3,5,7,9,2,4,6,8,10,11))
text(x=c(2,6,10,14,18,2,6,10,14,18,22),y=c(5,5,5,5,5,3,3,3,3,3,4),labels=c("sta","mo","glu","lac","pro","inu","fo","fru","ace","but","met"))

}