# for adding Harvest Control Rule HCR function

#Constant F
control <-  function() {
  
}

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
  control <- function(estimated.biomass, control.pars) {
    H1 <- control.pars$H1
    H2 <- control.pars$H2
    H3 <- control.pars$H3
    Bmax <- control.pars$Bmax
    B3 <- control.pars$B3
    B2 <- control.pars$B2
    B1 <- control.pars$B1
    
    harv <- ifelse(estimated.biomass >= B1, H1,
                   ifelse(estimated.biomass > B1 & estimated.biomass <= B2,
                          (H1-H2)/(B1-B2)*(estimated.biomass - B2) + H2,
                          ifelse(estimated.biomass > B2 & estimated.biomass < B3, H2,
                                 ifelse(estimated.biomass > B3,
                                        (H2-H3)/(B2-B3)*(estimated.biomass - B3) + H3
                                 ))))
    return(harv)
    
    
  }
