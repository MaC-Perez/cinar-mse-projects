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