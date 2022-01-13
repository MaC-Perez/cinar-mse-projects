# file for Operating Model OM functions

schaefer <- function(B,C,K,r) {
  #function schaefer takes the current biomass, a catch, 
  #and the model parameters to compute next year's biomass
  res <- B + B * r * (1 - B/K) - C
  return(max(0.001,res))  # we add a constraint to prevent negative biomass
}

msschaefer <- function(B,C,K,r,alpha) {
  #function msschaefer takes a vector of current biomass, vector of catch
  #and model parameters to compute next year's biomass vector for a multispecies system
  #alpha is an nsp*nsp matrix of interaction parameters
  #nsp <- length(B)
  res <- B + B * r * (1 - B/K) - C - alpha%*%B*B
}

### Gavin's MSPROD equation
## Solves the multsipecies operating model dynamics for a single time step given parameters, a set of harvest rates, and the current biomass
# SKG would want to remove the guild stuff here and keep ony the predloss term based on alpha
# which ought to reduce to the above??
dNbydt <- function(t,N=1,parms=list(r=rep(0.4,length(N)),KGuild=rep(1,1),
                                    Ktot=10,alpha=matrix(0,nrow=1,ncol=1),
                                    Guildmembership=1,
                                    BetweenGuildComp=matrix(0,nrow=1,ncol=1),
                                    WithinGuildComp=matrix(0,nrow=1,ncol=1),hrate=0)) {
  NG <- aggregate(N,by=list(parms$Guildmembership),sum,na.rm=TRUE)$x
  NG <- t(parms$BetweenGuildComp)%*%NG
  KbyGuild <-parms$KGuild[parms$Guildmembership]
  
  predloss <-  parms$alpha%*%N*N
  cat <- parms$hrate*N
  betweenloss <- parms$r*N*NG[parms$Guildmembership]/(parms$Ktot-KbyGuild)
  withinloss <- parms$r*N*(parms$WithinGuildComp%*%N)/KbyGuild
  dN <- parms$r*N*(1-N/KbyGuild) - withinloss - betweenloss - predloss - cat
  results <- list(deriv=c(dN),catch=cat,predloss=predloss,withinloss=withinloss,betweenloss=betweenloss)
  return(results)
}

### Gavin's MSPROD get indicators function
get.Indicators <- function(Biomass=NULL,Catch=NULL,size=NULL,trophic.level=NULL,BMSY=NULL,lifespan=NULL,is.predator=NULL,is.pelagic=NULL)
{
  # revenue
  Rev <- c(2,1,2)    #c(2,2,1,1,1,1,1,3,3,3)
  Nyr <- nrow(Biomass)
  Nsp <- ncol(Biomass)
  
  # Total system biomass summed over all species
  tot.bio <- rowSums(Biomass,na.rm=TRUE)
  # Total system catch summed over all species
  tot.cat <- rowSums(Catch,na.rm=TRUE)
  # Exploitation rate
  tot.rev <- rep(0,Nyr)
  for (Iyr in 1:Nyr) tot.rev[Iyr] <- sum(Catch[Iyr,]*Rev)
  exprate <- tot.cat/tot.bio
  # mean.length <- sum(Biomass*size,na.rm=TRUE)/sum(Biomass,na.rm=TRUE)
  # Trophic level of landings
  TL.landings <- rep(NA,nrow(Biomass))
  # Trophci level of survey
  TL.survey <- rep(NA,nrow(Biomass)) 
  # Proportion of total biomass that is comprised by predatory species
  prop.predators <- rowSums(Biomass[,is.predator, drop=FALSE],na.rm=TRUE)/tot.bio
  # Pelagic demersal ratio
  pd.ratio <- rowSums(Biomass[,is.pelagic, drop=FALSE],na.rm=TRUE)/rowSums(Biomass[,-(is.pelagic), drop=FALSE],na.rm=TRUE)
  # Proportion of total biomass that is made of pelagic species
  prop.pel <- 1-(1/(pd.ratio+1))
  # Proportion of species that is overfished (less than half BMSY)
  prop.overfishedSp <- matrix(0,ncol=Nsp,nrow=Nyr)
  prop.overfished <- rep(NA,Nyr)
  for (i in 1:nrow(Biomass)) 
  {
    TL.landings[i] <- sum(trophic.level*Catch[i,],na.rm=TRUE)/sum(Catch[i,],na.rm=TRUE)
    TL.survey[i] <- sum(trophic.level*Biomass[i,],na.rm=TRUE)/sum(Biomass[i,],na.rm=TRUE)
    b.use <- Biomass[i,]
    prop.overfished[i] <- length(b.use[b.use<0.5*BMSY])/length(b.use)
    for (k in 1:Nsp)
      if (b.use[k] < 0.5*BMSY[k]) prop.overfishedSp[i,k] <- 1
  }
  # Rolling 10-year window of CVs (fill in)
  div.cv.bio <- rep(NA,nrow(Biomass))
  for (i in 10:nrow(Biomass))
    div.cv.bio[i] <- 1/(sd(tot.bio[((i-9):i)],na.rm=TRUE)/mean(tot.bio[((i-9):i)],na.rm=TRUE))
  results <- list(tot.bio=tot.bio,tot.cat=tot.cat,tot.rev=tot.rev,exprate=exprate,div.cv.bio=div.cv.bio,prop.overfished=prop.overfished,prop.pel=prop.pel,prop.predators=prop.predators,TL.landings=TL.landings,TL.survey=TL.survey,prop.overfishedSp=prop.overfishedSp)
  return(results)
}