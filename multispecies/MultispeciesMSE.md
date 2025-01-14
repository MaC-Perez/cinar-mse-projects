MultispeciesMSE
================
Sarah Gaichas
1/11/2022

Steps for multispecies MSE:

# Condition operating model on data

## Data

Simulated output of ecosystem model. Options include Norwegian-Barents
Atlantis generated dataset in
[`mskeyrun`](https://noaa-edab.github.io/ms-keyrun/) R package and
generating data from the [Georges Bank Rpath
model](https://github.com/NOAA-EDAB/GBRpath).

### mskeyrun (NOBA Atlantis) data

Which species to use in the model? Probably the ones that interact the
most. Here we use simulated survey diet composition to find the highest
average diet proportions between the 11 species in the simulated
dataset. We are assuming we only want to focus on managed fish species,
rather than interactions with lower trophic levels or protected species.

``` r
# install data package if needed
#remotes::install_github("NOAA-EDAB/ms-keyrun")
# the libraries are all in the first code chunk above now
#library(mskeyrun)

# which species pairs interact the most? 
# a too-quick analysis finding the mean proportion in diet from the diet dataset:

predprey <- mskeyrun::simSurveyDietcomp %>%
  filter(prey %in% unique(Name)) %>%
  select(Name, agecl, year, prey, value) %>%
  group_by(Name, prey) %>%
  summarise(avgpreyprop = mean(value)) %>%
  unite('Pred-Prey', Name:prey, sep = "-")
```

    ## `summarise()` has grouped output by 'Name'. You can override using the `.groups` argument.

``` r
predprey %>% arrange(desc(avgpreyprop))
```

    ## # A tibble: 34 x 2
    ##    `Pred-Prey`                  avgpreyprop
    ##    <chr>                              <dbl>
    ##  1 Redfish-Blue_whiting              0.468 
    ##  2 North_atl_cod-Capelin             0.133 
    ##  3 North_atl_cod-Haddock             0.0921
    ##  4 Redfish-Haddock                   0.0824
    ##  5 Blue_whiting-Redfish              0.0414
    ##  6 North_atl_cod-Long_rough_dab      0.0332
    ##  7 Haddock-Long_rough_dab            0.0243
    ##  8 Green_halibut-Redfish             0.0210
    ##  9 Haddock-Norwegian_ssh             0.0193
    ## 10 Haddock-North_atl_cod             0.0168
    ## # … with 24 more rows

Select species from the mskeyrun dataset and get survey biomass and
catch data for fitting. For now, pull Redfish and Blue\_whiting. Later
could add North\_atl\_cod, Capelin, Haddock as possibilities.

``` r
index.df <- mskeyrun::simSurveyIndex %>% filter(Name %in% c("Redfish", "Blue_whiting"))
harvest.df <- mskeyrun::simCatchIndex %>% filter(Name %in% c("Redfish", "Blue_whiting"))

data.years <- unique(mskeyrun::simSurveyIndex$year)
```

Plot them

``` r
# plot(data.years,index, pch=19,xlab="Year",ylab="Million tonnes (B or C)",
#      ylim=c(0,1200))
# lines(data.years,harvest,lty=2,lwd=2)

index.df %>%
  filter(variable=="biomass") %>%
  ggplot() +
  geom_point(aes(x=year, y=value, colour=survey)) +
  geom_line(data=(harvest.df %>% filter(variable=="catch")), 
            aes(x=year, y=value)) +
  theme_bw() +
  facet_wrap(~Name)
```

![](MultispeciesMSE_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Do we want to use just one of the surveys if we use this dataset? Maybe
fall because Blue whiting migrate from the system in the spring (though
maybe that is when Redfish eat them all? might be worth a look)

### Rpath data

Full dataset now in the `data` folder in this repo.

``` r
#load Rpath Georges Bank model outputs, object is called GB.data
load(here("multispecies/data/GBdata.rda"))

#units are tons/square kilometer for biomass and catch
str(GB.data)
```

    ## Classes 'data.table' and 'data.frame':   120 obs. of  5 variables:
    ##  $ Year    : int  1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 ...
    ##  $ Species : chr  "Cod" "Cod" "Cod" "Cod" ...
    ##  $ ObsBio  : num  0.442 0.574 0.552 0.578 0.586 ...
    ##  $ TotCatch: num  0.00181 0.00184 0.00186 0.00188 0.00188 ...
    ##  $ Fmort   : num  0.00409 0.00321 0.00338 0.00325 0.00321 ...
    ##  - attr(*, ".internal.selfref")=<externalptr>

``` r
#there are three spexxies
unique(GB.data$Species)
```

    ## [1] "Cod"        "Haddock"    "AtlHerring"

``` r
data.years <- unique(GB.data$Year)

GB.data %>%
  ggplot()+
  geom_point(aes(x=Year, y=ObsBio)) +
  geom_line(aes(x=Year, y=TotCatch)) +
  theme_bw() +
  facet_wrap(~Species, scales = "free")
```

![](MultispeciesMSE_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Multispecies operating model

Production model or biomass-dynamics model form with species interaction
terms. Use Gavin’s `schaefer` function as modified below and reuse as
many other functions as possible.

What interactions to model?

  - Predation mortality
  - Prey impacts on growth?

### Schaefer model (single species, directly from [my-first-mse](https://github.com/gavinfay/cinar-mse/blob/main/materials/exercises/day-01/my-first-mse.Rmd))

First, the logistic production function:

``` r
schaefer <- function(B,C,K,r) {
  #function schaefer takes the current biomass, a catch, 
  #and the model parameters to compute next year's biomass
  res <- B + B * r * (1 - B/K) - C
  return(max(0.001,res))  # we add a constraint to prevent negative biomass
}
```

Now a function to do the biomass projection:

``` r
dynamics <- function(pars,C,yrs) {
  # dynamics takes the model parameters, the time series of catch, 
  # & the yrs to do the projection over
  
  # first extract the parameters from the pars vector (we estimate K in log-space)
  K <- exp(pars[1])
  r <- exp(pars[2])
  
  # find the total number of years
  nyr <- length(C) + 1
  
  # if the vector of years was not supplied we create 
  # a default to stop the program crashing
  if (missing(yrs)) yrs <- 1:nyr
  
  #set up the biomass vector
  B <- numeric(nyr)
  
  #intialize biomass at carrying capacity
  B[1] <- K
  # project the model forward using the schaefer model
  for (y in 2:nyr) {
    B[y] <- schaefer(B[y-1],C[y-1],K,r)
  }
  
  #return the time series of biomass
  return(B[yrs])
  
  #end function dynamics
}  
```

We are going to condition the operating model by estimating the
parameters based on the historical biomass index data.

To do this we make a function that shows how well the current parameters
fit the data, we assume that the observation errors around the true
biomass are log-normally distributed.

``` r
# function to calculate the negative log-likelihood
nll <- function(pars,C,U) {  #this function takes the parameters, the catches, and the index data
  sigma <- exp(pars[3])  # additional parameter, the standard deviation of the observation error
  B <- dynamics(pars,C)  #run the biomass dynamics for this set of parameters
  Uhat <- B   #calculate the predicted biomass index - here we assume an unbiased absolute biomass estimate
  output <- -sum(dnorm(log(U),log(Uhat),sigma,log=TRUE),na.rm=TRUE)   #calculate the negative log-likelihood
  return(output)
  #end function nll
}
```

Function to perform the assessment and estimate the operating model
parameters  
(i.e. to fit the logistic model to abundance data)

``` r
assess <- function(catch,index,calc.vcov=FALSE,pars.init) {
  # assess takes catch and index data, initial values for the parameters,
  # and a flag saying whether to compute uncertainty estimates for the model parameters
  
  #fit model
  # optim runs the function nll() repeatedly with differnt values for the parameters,
  # to find the values that give the best fit to the index data
  res <- optim(pars.init,nll,C=catch,U=index,hessian=TRUE)
  
  # store the output from the model fit
  output <- list()
  output$pars <- res$par
  output$biomass <- dynamics(res$par,catch)
  output$convergence <- res$convergence
  output$nll <- res$value
  if (calc.vcov)
    output$vcov <- solve(res$hessian)
  
  return(output)
  #end function assess
}
```

## Conditioning

### Single species models?

Try fitting to each species separately first, adjusting Gavin’s code
(test with atlantis output till we get rpath):

Make separate data objects for each species

``` r
harvest.redfish <- harvest.df %>%
  filter(Name=="Redfish",
         variable=="catch") %>%
  select(value) %>%
  unlist() %>%
  unname()

index.redfish <- index.df %>%
  filter(Name=="Redfish",
         variable=="biomass", 
         survey=="BTS_fall_allbox_effic1") %>%
  select(value) %>%
  unlist() %>%
  unname()

harvest.bluewhiting <- harvest.df %>%
  filter(Name=="Blue_whiting",
         variable=="catch") %>%
  select(value) %>%
  unlist() %>%
  unname()

index.bluewhiting <- index.df %>%
  filter(Name=="Blue_whiting",
         variable=="biomass", 
         survey=="BTS_fall_allbox_effic1") %>%
  select(value) %>%
  unlist() %>%
  unname()

harvest.GBcod <- GB.data %>%
  filter(Species=="Cod") %>%
  select(TotCatch) %>%
  unlist() %>%
  unname()

index.GBcod <- GB.data %>%
  filter(Species=="Cod") %>%
  select(ObsBio) %>%
  unlist() %>%
  unname()

harvest.GBherring <- GB.data %>%
  filter(Species=="AtlHerring") %>%
  select(TotCatch) %>%
  unlist() %>%
  unname()

index.GBherring <- GB.data %>%
  filter(Species=="AtlHerring") %>%
  select(ObsBio) %>%
  unlist() %>%
  unname()

harvest.GBhaddock <- GB.data %>%
  filter(Species=="Haddock") %>%
  select(TotCatch) %>%
  unlist() %>%
  unname()

index.GBhaddock <- GB.data %>%
  filter(Species=="Haddock") %>%
  select(ObsBio) %>%
  unlist() %>%
  unname()  
```

First define initial parameter vector for each species: log(K), log(r),
log(sigma)

``` r
#ini.parms <- c(log(1200), log(0.1), log(0.3))

# use max survey B as the initial value for K
maxsvB.red <- max(index.redfish)
ini.parms.red <- c(log(maxsvB.red), log(0.1), log(0.3))

maxsvB.blue <- max(index.bluewhiting)
ini.parms.blue <- c(log(maxsvB.blue), log(0.5), log(0.3))

maxsvB.GBcod <- max(index.GBcod)
ini.parms.GBcod <- c(log(maxsvB.GBcod), log(0.2), log(0.3))

maxsvB.GBherring <- max(index.GBherring)
ini.parms.GBherring <- c(log(maxsvB.GBherring), log(0.3), log(0.3))

maxsvB.GBhaddock <- max(index.GBhaddock)
ini.parms.GBhaddock <- c(log(maxsvB.GBhaddock), log(0.3), log(0.3))
```

Fit the logistic model to data (redfish):

``` r
redfish <- assess(harvest.redfish,index.redfish,calc.vcov=TRUE,ini.parms.red)
redfish
```

    ## $pars
    ## [1] 14.7332517 -0.5727151 -2.1959659
    ## 
    ## $biomass
    ##  [1] 2503629 2503629 2503629 2503629 2503629 2503629 2503629 2503629 2503629
    ## [10] 2503629 2503629 2503629 2503629 2503629 2503629 2503629 2050018 1858715
    ## [19] 1745812 1697281 1658733 1632654 1609601 1592110 1577896 1596415 1609842
    ## [28] 1614435 1625567 1621703 1625750 1648597 1673604 1681991 1679961 1687559
    ## [37] 1682072 1680763 1681549 1687687 1705269 1896551 2038527 2130021 2187876
    ## [46] 2218672 2232304 2242516 2249058 2240322 2235494 2239173 2235118 2228290
    ## [55] 2225732 2227338 2226934 2224399 2218398 2216680 2220569 2218419 2212572
    ## [64] 2209754 2208218 2216547 2144049 2105712 2083456 2083763 2082256 2086918
    ## [73] 2081844 2086650 2081940 2077019 2076170 2075742 2078112 2078577 2074357
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $nll
    ## [1] -62.8732
    ## 
    ## $vcov
    ##               [,1]          [,2]          [,3]
    ## [1,]  3.333499e-04 -8.383304e-04 -1.935791e-06
    ## [2,] -8.383304e-04  2.700502e-03  1.126767e-05
    ## [3,] -1.935791e-06  1.126767e-05  6.162865e-03

Fit the logistic model to data (bluewhiting):

``` r
bluewhiting <- assess(harvest.bluewhiting,index.bluewhiting,calc.vcov=TRUE,ini.parms.blue)
bluewhiting
```

    ## $pars
    ## [1] 16.0585331 -0.6604354 -1.7674470
    ## 
    ## $biomass
    ##  [1] 9421766 9421766 9421766 9421766 9421766 9421766 9421766 9421766 9421766
    ## [10] 9421766 9421766 9421766 9421766 9421766 9421766 9421766 6108595 4708036
    ## [19] 3963839 3576136 3331128 3205756 3114140 3090013 3112902 3125688 3134304
    ## [28] 3173525 3184026 3212229 3243509 3273462 3318300 3354723 3394631 3441354
    ## [37] 3510169 3560318 3622509 3682686 3777356 4491775 5220294 5864792 6378247
    ## [46] 6735922 6964984 7093650 7134389 7138048 7140363 7131820 7103820 7069142
    ## [55] 7034204 7013317 7014964 6999413 6990766 6996574 6997470 7015605 7026428
    ## [64] 7019775 7041708 7042283 6608908 6351623 6210359 6102879 6066279 6071172
    ## [73] 6073240 6110343 6134802 6166780 6197304 6200742 6211462 6225605 6240488
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $nll
    ## [1] -28.22243
    ## 
    ## $vcov
    ##               [,1]          [,2]          [,3]
    ## [1,]  3.646685e-05 -6.302346e-05 -0.0000879590
    ## [2,] -6.302346e-05  1.098883e-04  0.0001511557
    ## [3,] -8.795900e-05  1.511557e-04  0.0063847211

SS model fit plots: fitted line in gray

``` r
fit.red <- data.frame(value= redfish$biomass,
                      year=c(40:120),
                      Name="Redfish")

fit.blue <- data.frame(value= bluewhiting$biomass,
                       year=c(40:120),
                       Name="Blue_whiting")

fits <- bind_rows(fit.red, fit.blue)


index.df %>%
  filter(variable=="biomass") %>%
  ggplot() +
  geom_point(aes(x=year, y=value, colour=survey)) +
  geom_line(data=(harvest.df %>% filter(variable=="catch")), 
            aes(x=year, y=value)) +
  geom_line(data=fits, aes(x=year, y=value), lwd=1.3, col=gray(0.5)) +
  theme_bw() +
  facet_wrap(~Name)
```

![](MultispeciesMSE_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Fit the logistic model to data (GBcod):

``` r
GBcod <- assess(harvest.GBcod,index.GBcod,calc.vcov=TRUE,ini.parms.GBcod)
GBcod
```

    ## $pars
    ## [1] -0.6701935 -0.2304331 -2.3782841
    ## 
    ## $biomass
    ##  [1] 0.5116096 0.5098017 0.5093893 0.5092806 0.5092458 0.5092326 0.5092276
    ##  [8] 0.5092261 0.5092265 0.5092276 0.5092290 0.5092305 0.5092318 0.5092330
    ## [15] 0.5092341 0.5092349 0.4719121 0.4636557 0.4618841 0.4618094 0.4620525
    ## [22] 0.4307199 0.4210665 0.4181843 0.4175851 0.4176187 0.3657245 0.3457196
    ## [29] 0.3370487 0.3331660 0.3312795 0.3983110 0.4386530 0.4559198 0.4612804
    ## [36] 0.4623669 0.4816092 0.4874697 0.4889433 0.4892696 0.4893518
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $nll
    ## [1] -39.3231
    ## 
    ## $vcov
    ##               [,1]          [,2]          [,3]
    ## [1,]  3.212088e-04 -1.011865e-03 -6.862185e-07
    ## [2,] -1.011865e-03  5.764018e-03  7.038378e-06
    ## [3,] -6.862185e-07  7.038378e-06  1.218913e-02

Fit the logistic model to data (GBherring):

``` r
GBherring <- assess(harvest.GBherring,index.GBherring,calc.vcov=TRUE,ini.parms.GBherring)
GBherring
```

    ## $pars
    ## [1]  2.0452567  0.2470391 -2.3721717
    ## 
    ## $biomass
    ##  [1] 7.731143 7.728940 7.729494 7.729314 7.729356 7.729341 7.729344 7.729343
    ##  [9] 7.729343 7.729343 7.729343 7.729343 7.729343 7.729343 7.729343 7.729343
    ## [17] 6.554175 6.733990 6.777934 6.790430 6.794363 5.867628 5.831486 5.868671
    ## [25] 5.899521 5.916098 4.554007 4.225116 4.103541 4.050814 4.023411 6.180407
    ## [33] 7.366130 7.368061 7.350684 7.347242 7.346002 7.346277 7.346699 7.346952
    ## [41] 7.347016
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $nll
    ## [1] -39.07345
    ## 
    ## $vcov
    ##               [,1]          [,2]          [,3]
    ## [1,]  2.698832e-04 -4.835464e-04 -7.479694e-07
    ## [2,] -4.835464e-04  1.298939e-03  7.440093e-06
    ## [3,] -7.479694e-07  7.440093e-06  1.218977e-02

Fit the logistic model to data (GBhaddock):

``` r
GBhaddock <- assess(harvest.GBhaddock,index.GBhaddock,calc.vcov=TRUE,ini.parms.GBhaddock)
GBhaddock
```

    ## $pars
    ## [1]  1.8602168  0.7813654 -2.2906425
    ## 
    ## $biomass
    ##  [1] 6.425130 6.411249 6.427585 6.408286 6.431051 6.404179 6.435881 6.398453
    ##  [9] 6.442594 6.390461 6.451915 6.379299 6.464848 6.363705 6.482765 6.341904
    ## [17] 6.417700 6.331018 6.430960 6.315595 6.448059 6.206758 6.477473 6.172147
    ## [25] 6.512577 6.127937 6.381673 6.114104 6.399022 6.092905 6.416214 6.332755
    ## [33] 6.426125 6.317310 6.441884 6.298218 6.518201 6.260053 6.559505 6.208016
    ## [41] 6.614554
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $nll
    ## [1] -35.7394
    ## 
    ## $vcov
    ##               [,1]          [,2]          [,3]
    ## [1,]  2.456998e-04 -3.163436e-06 -1.294359e-07
    ## [2,] -3.163436e-06  2.166919e-05  8.424341e-07
    ## [3,] -1.294359e-07  8.424341e-07  1.219486e-02

GB SS model fit plots: fitted line in gray

``` r
#model may not be using the last data year? misaligned. a hack here to fix

fit.cod <- data.frame(value= GBcod$biomass,
                      Year=c(1983:2023),
                      Species="Cod")

fit.herring <- data.frame(value= GBherring$biomass,
                       Year=c(1983:2023),
                       Species="AtlHerring")

fit.haddock <- data.frame(value= GBhaddock$biomass,
                       Year=c(1983:2023),
                       Species="Haddock")

fits <- bind_rows(fit.cod, fit.herring, fit.haddock)


GB.data %>%
  ggplot()+
  geom_point(aes(x=Year, y=ObsBio)) +
  geom_line(aes(x=Year, y=TotCatch)) +
  geom_line(data=fits, aes(x=Year, y=value), lwd=1.3, col=gray(0.5)) +
  theme_bw() +
  facet_wrap(~Species, scales = "free")
```

![](MultispeciesMSE_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

### Try multispecies Schaeffer model

So I like this idea, but fitting it will not be trivial.

There are options such as those described here
<https://rpubs.com/Jeet1994/Prey-predator-model>

This is a start but I don’t think we can use the nll function for
solving, may have to go with something like stan as above.

Now working with vectors and matrices

``` r
msschaefer <- function(B,C,K,r,alpha) {
  #function msschaefer takes a vector of current biomass, vector of catch
  #and model parameters to compute next year's biomass vector for a multispecies system
  #alpha is an nsp*nsp matrix of interaction parameters
  #nsp <- length(B)
  res <- B + B * r * (1 - B/K) - C - alpha%*%B*B
}
```

Input data (would be closer to GB.data if I try the Prey-predator model
approach linked above)

``` r
index.GB3 <- GB.data %>%
  select(Year, Species, ObsBio) %>%
  pivot_wider(names_from = Species, values_from = ObsBio) %>%
  select(-Year) %>%
  unname()

harvest.GB3 <- GB.data %>%
  select(Year, Species, TotCatch) %>%
  pivot_wider(names_from = Species, values_from = TotCatch)%>%
  select(-Year) %>%
  unname()
```

Modify other functions for vectors and matrices **THIS ISN’T DONE and
DOESN’T WORK**

``` r
msdynamics <- function(pars,C,yrs) {
  # dynamics takes the model parameters, the time series of catch, 
  # & the yrs to do the projection over
  
  # first extract the parameters from the pars matrix (we estimate K in log-space)
  K <- exp(pars[1,])
  r <- exp(pars[2,])
  alpha <- matrix(0,3,3)
  
  # number of species
  nsp <- dim(C)[2]
  
  # find the total number of years
  nyr <- dim(C)[1] + 1
  
  # if the vector of years was not supplied we create 
  # a default to stop the program crashing
  if (missing(yrs)) yrs <- 1:nyr
  
  #set up the biomass vector
  B <- numeric(nyr)
  
  #intialize biomass at carrying capacity
  B[1] <- K
  # project the model forward using the schaefer model
  for (y in 2:nyr) {
    B[y] <- schaefer(B[y-1],C[y-1],K,r alpha)
  }
  
  #return the time series of biomass
  return(B[yrs])
  
  #end function dynamics
}  
```

### For a 3 day project we are probably best off not trying to condition a multispecies operating model to data, or to have a multispecies assessment model

### Instead parameterize MSprod with EwE data

Use Gavin’s framework but input values from Rpath:

  - Cut input data down to 3 species model
  - Use Rpath B and C data as inputs  
  - Also use SS estimates of r and K by species as inputs

Pars from SS models above

``` r
pars.cod <- data.frame(value= GBcod$pars,
                       Par=c("logK", "logr", "logsigma"),
                       Species="Cod")

pars.herring <- data.frame(value= GBherring$pars,
                           Par=c("logK", "logr", "logsigma"), 
                           Species="AtlHerring")

pars.haddock <- data.frame(value= GBhaddock$pars,
                           Par=c("logK", "logr", "logsigma"),
                           Species="Haddock")

pars <- bind_rows(pars.cod, pars.herring, pars.haddock) %>%
  mutate(exppar = exp(value)) 
```

Reference points from single species models above for each species.

If these are going into Rpath maybe ok to leave on this scale, but I’m
wondering if we should estimate on the tons scale for input into MSprod?

``` r
SSBMSY <- pars %>%
  filter(Par=="logK") %>%
  mutate(BMSY=exppar/2) %>%
  select(Species, BMSY)

SSFmsy <- pars %>%
  filter(Par=="logr") %>%
  mutate(Fmsy=exppar/2) %>%
  select(Species, Fmsy)
```

Now build the data input (Rpath scale, t/km2)

``` r
#read in the base biological parameter values
#datfile <- "Georges.dat"

#our species order: Cod, AtlHerring, Haddock
# so guild 2 actually is pelagic as assumed below!

#Number of species
Nsp <- dim(index.GB3)[2] #scan(datfile,n=1,skip=3,quiet=T)
#guilds / functional groups
Guildmembership <- c(1,2,3) #scan(datfile,n=Nsp,skip=9,quiet=T)
NGuild = length(unique(Guildmembership))
#Initial values
Initvals <- unlist(index.GB3[1,])     #scan(datfile,n=Nsp,skip=13,quiet=T)
#carrying capacity for each guild
KGuild <-  unlist(unname(pars %>% filter(Par=="logK") %>% select(exppar)))      #scan(datfile,n=NGuild,skip=17,quiet=T)
Ktot <- sum(KGuild)
cat("Spsecies/Guilds",Nsp,NGuild,"\n")
```

    ## Spsecies/Guilds 3 3

``` r
#growth rates
r <- unlist(unname(pars %>% filter(Par=="logr") %>% select(exppar)))         #scan(datfile,n=Nsp,skip=11,quiet=T)
#interactions
BetweenGuildComp <- matrix(c(0,0,0.236, #taken from Georges.dat, Cod and Haddock compete
                             0,0,0,
                             0.236,0,0), byrow=TRUE,
                           3,3) #matrix(scan(datfile,n=NGuild^2,skip=19,quiet=T),byrow=TRUE,nrow=NGuild)
WithinGuildComp <- matrix(0,3,3)       #matrix(scan(datfile,n=Nsp^2,skip=(20+NGuild)),byrow=TRUE,nrow=Nsp)
alpha <- matrix(c(0,0,0,  # taken from Georges.dat, preds in columns, prey in rows, herring is prey of cod
                  0.00000002,0,0,
                  0,0,0), byrow=TRUE,
                3,3)       #matrix(scan(datfile,n=Nsp^2,skip=(21+NGuild+Nsp),quiet=T),byrow=TRUE,nrow=Nsp)
spatial.overlap <- matrix(1,3,3)       #matrix(scan(datfile,n=Nsp^2,skip=(22+NGuild+2*Nsp),quiet=T),byrow=TRUE,nrow=Nsp)
alpha <- alpha*spatial.overlap
# This is zero for this example
WithinGuildComp <- WithinGuildComp*spatial.overlap

# Reset harvest rates
hrate <- rep(0,Nsp)

#set values for BMSY
BMSY <- read.csv(here("multispecies/data/Bmsy.csv"),header=TRUE)
BMSY <- BMSY[c(4,5,21),]    #BMSY[c(4,5,21,22,14,23,24,6,3,7),]
BMSY[,2] <- KGuild/2
print(BMSY)
```

    ##    Species.Group      Bmsy MTL
    ## 4         GB Cod 0.2558048 4.4
    ## 5     GB Haddock 3.8655715 4.1
    ## 21       Herring 3.2125648 3.2

``` r
#initial biomass for each species
N <- Initvals
# Set projection horizon
Nyr <- 30
# Set the number of simulations
Nsims <- 5

### get historical time series of biomass and catch
NI <- as.matrix(index.GB3)  #read.table(datfile,skip=69,nrow=33,header=FALSE)
#NI <- NI[,-1] #this cuts off a year column that I don't have in the input
CI <- as.matrix(harvest.GB3)  #read.table(datfile,skip=103,nrow=33,header=FALSE)

#redefine functional groups
theguilds <- c(1,1,2)   #c(1,1,2,2,1,3,3,1,1,1)
```

Bring in the OMfunctions

``` r
source(here("multispecies/OMfunctions.R"))
```

Test the OMfunction

``` r
#inputs
# Biomass <- NI #needs to come in as.matrix or won't work
# Catch <- CI #needs to come in as.matrix or won't work
# trophic.level <- BMSY[,3]
# BMSY <- BMSY[,2]
# is.predator <- which(colSums(alpha)>0)
# is.pelagic <- which(theguilds==2)
```

Modify the MSprod run function to use the SS ref pts and HCRs below.
This is a test with Gavin’s code and our inputs, no change to HCR (aside
from dropping elasmobranchs):

``` r
############################################
# RUN MSE WITH SINGLE SPECIES ASSESSMENT
############################################

library(deSolve)

#set up a storage object to contain results for each simulation
ALL.results <- matrix(NA,nrow=Nsims,ncol=4+Nsp)

#do a bunch of simulations
for (isim in 1:Nsims)
 {
  ### calculate values for ecological indicators at start of projection
  # AEP:IS BMSY correct
  ei <- get.Indicators(Biomass=NI,Catch=CI,BMSY=BMSY[,2],trophic.level=BMSY[,3],is.predator=which(colSums(alpha)>0),is.pelagic=which(theguilds==2)) 
  ei <- as.data.frame(ei)
  ei.new <- as.numeric(ei[nrow(ei),])
  names(ei.new) = colnames(ei)
  NI.obs <- NI
  CI.obs <- CI
  
  # Output reset
  SS.results <- NULL
  
  # do projection period
  for (iyr in 2:Nyr)
   {    
    ### hrateG is for groundfish; hrateP is for pelagics; hrateE is for elasmobranches
    # WARNING These hardcode the order of species and expexct 10, adjusted here but buyer beware
    Qyr2 <- nrow(NI.obs)
    #print(NI.obs)
    # hrate for groundfish depends on cod
    if (NI.obs[Qyr2,1] > 0.5*mean(NI.obs[1:10,1]))
     hrateG <- 0.1
    else
     hrateG <-0.1*NI.obs[Qyr2,1]/mean(NI.obs[1:10,1])
    # hrate for pelagics is 0.1 unless there is not enough prey
    if (NI.obs[Qyr2,3] > 1000) #(NI.obs[Qyr2,4]+NI.obs[Qyr2,5] > 1000)
      hrateP <- 0.3
    else
      hrateP <- 0
    # # hrate for sharks and skates is 0.05 unless we are approximately overfished
    # if (NI.obs[Qyr2,7]+NI.obs[Qyr2,8] >  0.5*mean(NI.obs[1:10,7])+0.5*mean(NI.obs[1:10,8]))
    #   hrateE <- 0.05
    # else
    #   hrateE <- 0

    hrate[theguilds==1] <- hrateG
    hrate[theguilds==2] <- hrateP
    #hrate[theguilds==3] <- hrateE
    ### update operating model with new exploitation rates
    parms=list(r=r,KGuild=KGuild,Ktot=Ktot,Guildmembership=Guildmembership,BetweenGuildComp=BetweenGuildComp,WithinGuildComp=WithinGuildComp,alpha=alpha,hrate=hrate)
    x <- ode(N,seq(iyr-1,(iyr+0.5),0.5),dNbydt,parms=parms,method="rk4")

    # Extract N after one year and catches during the year; Rem is various losses
    N <- x[3,2:(Nsp+1)];N[N<=0] <- 0.01
    Cat <- 1.e-07+x[2,(Nsp+2):(2*Nsp+1)]; Cat[Cat<=0] <- 0.01
    Rem <- x[2,(2*Nsp+2):ncol(x)]

    ### store results for this time step
    SS.results <- rbind(SS.results,c(x[1,2:(Nsp+1)],x[2,(Nsp+2):ncol(x)]))
    if (iyr==Nyr) SS.results <- rbind(SS.results,c(x[3,2:(Nsp+1)],x[4,(Nsp+2):ncol(x)]))
    #generate data for this timestep and append to dataset
    Nobs <- N*exp(rnorm(10,0,0.2)-0.5*0.2*0.2)
    Cobs <- Cat
    NI.obs <- rbind(NI.obs,Nobs)
    CI.obs <- rbind(CI.obs,Cobs)
    
    #calculate ecological indicators based on new data at this time step
    # AEP:IS BMSY correct
    ei <- get.Indicators(Biomass=NI.obs,Catch=CI.obs,BMSY=BMSY[,2],trophic.level=BMSY[,3],is.predator=which(colSums(alpha)>0),is.pelagic=which(theguilds==2)) 
    ei <- as.data.frame(ei)
    ei.now <- as.numeric(ei[nrow(ei),])
    names(ei.now) = colnames(ei)
  }
 #print(ei)    
 colnames(SS.results)[1:Nsp] <- paste("Bio:",BMSY[,1],sep="")
 colnames(SS.results)[(Nsp+1):(2*Nsp)] <- paste("Catch:",BMSY[,1],sep="")
 SS.results <- SS.results[,1:(2*Nsp)]
 
 Qyr <- length(ei$tot.bio)
 ALL.results[isim,1] <- mean(ei$tot.bio[(Qyr-9):Qyr])
 ALL.results[isim,2] <- mean(ei$tot.cat[(Qyr-9):Qyr])
 ALL.results[isim,3] <- mean(ei$tot.rev[(Qyr-9):Qyr])
 ALL.results[isim,4] <- mean(ei$prop.overfished[(Qyr-9):Qyr])
 for (Isp in 1:Nsp)
   ALL.results[isim,4+Isp] <- sum(ei[(Qyr-9):Qyr,11+Isp])
 }
```

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in N * exp(rnorm(10, 0, 0.2) - 0.5 * 0.2 * 0.2): longer object length is
    ## not a multiple of shorter object length

    ## Warning in rbind(NI.obs, Nobs): number of columns of result is not a multiple of
    ## vector length (arg 2)

``` r
 print(ALL.results)
```

    ##          [,1]      [,2]      [,3] [,4] [,5] [,6] [,7]
    ## [1,] 14.54658 0.7171215 0.7547904    0    0    0    0
    ## [2,] 13.49745 0.6795621 0.7153044    0    0    0    0
    ## [3,] 14.11300 0.7136056 0.7510535    0    0    0    0
    ## [4,] 14.82354 0.7520265 0.7913274    0    0    0    0
    ## [5,] 13.38658 0.7135274 0.7509551    0    0    0    0

``` r
 cat("Expected total biomass:", mean(ALL.results[,1]),"\n")
```

    ## Expected total biomass: 14.07343

``` r
 cat("Expected total catch:", mean(ALL.results[,2]),"\n")
```

    ## Expected total catch: 0.7151686

``` r
 cat("Median revenue:", median(ALL.results[,3]),"\n")
```

    ## Median revenue: 0.7510535

``` r
 cat("Median revenue:", quantile(ALL.results[,3],prob=0.05),"\n")
```

    ## Median revenue: 0.7224346

``` r
 cat("Expected probability overfished:", mean(ALL.results[,4]),"\n")
```

    ## Expected probability overfished: 0

``` r
 for (Isp in 1:Nsp)
  cat("Probability overfished:", as.character(BMSY[Isp,1]),mean(ALL.results[,4+Isp])/10,"\n")
```

    ## Probability overfished: GB Cod 0 
    ## Probability overfished: GB Haddock 0 
    ## Probability overfished: Herring 0

# Design HCRs

Separate single species HCRs (start with Gavin’s specified as default)

# Project with HCRs

# Compare results

# Tradeoffs
