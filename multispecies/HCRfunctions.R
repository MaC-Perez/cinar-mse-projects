# for adding Harvest Control Rule HCR function adjusted for multiple species

#Constant F
control <- function(estimated.biomass, control.pars, Species) {
  control.pars$Htarg[control.pars$Species==Species]
}

#Conditional F
#control <- function(estimated.biomass, control.pars) {
#  H1 <- control.pars$H1
#  H2 <- control.pars$H2
#  Bmax <- control.pars$Bmax
#  B2 <- control.pars$B2
#  B1 <- control.pars$B1
  
#  harv <- ifelse(estimated.biomass >= B1, H1,
#                 ifelse(estimated.biomass < B2, H2,
#                        (H1-H2)/(B1-B2)*(estimated.biomass - B2) + H2))
  
#  return(harv)
  
#end function control
#}

#Plateau-slope-plateau-slop HCR  
control <- function(estimated.biomass, control.pars, Species) {
  
  H1 <- control.pars$H1[control.pars$Species==Species]
  H2 <- control.pars$H2[control.pars$Species==Species]
  H3 <- control.pars$H3[control.pars$Species==Species]
  
  B4 <- control.pars$B4[control.pars$Species==Species]
  B3 <- control.pars$B3[control.pars$Species==Species]
  B2 <- control.pars$B2[control.pars$Species==Species]
  B1 <- control.pars$B1[control.pars$Species==Species]
  
  harv <- ifelse(estimated.biomass >= B1, H1,
                 ifelse(estimated.biomass < B1 & estimated.biomass >= B2,
                        (H1-H2)/(B1-B2)*(estimated.biomass - B2) + H2,
                        ifelse(estimated.biomass < B2 & estimated.biomass >= B3, H2,
                               ifelse(estimated.biomass < B3 & estimated.biomass > B4,
                                      (H2-H3)/(B2-B3)*(estimated.biomass - B3) + H3,
                                      ifelse(estimated.biomass <= B4, H3
                                      )))))
  return(harv)
}
