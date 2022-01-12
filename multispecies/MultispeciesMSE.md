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
    ##  $ ObsBio  : num  0.552 0.475 0.529 0.473 0.518 ...
    ##  $ TotCatch: num  0.0154 0.0154 0.0155 0.0155 0.0155 ...
    ##  $ Fmort   : num  0.0278 0.0325 0.0292 0.0327 0.0299 ...
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
  facet_wrap(~Species)
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
    ## [1] -0.6154117 -0.6558974 -2.2247195
    ## 
    ## $biomass
    ##  [1] 0.5404184 0.5250409 0.5173695 0.5133596 0.5112142 0.5100532 0.5094223
    ##  [8] 0.5090799 0.5088954 0.5087972 0.5087462 0.5087208 0.5087089 0.5087043
    ## [15] 0.5087032 0.5087038 0.5059534 0.5045932 0.5039484 0.5036558 0.5035282
    ## [22] 0.5007887 0.4994081 0.4987280 0.4983971 0.4982342 0.4929030 0.4901548
    ## [29] 0.4887389 0.4880037 0.4876125 0.4966825 0.5015735 0.5040732 0.5052920
    ## [36] 0.5058625 0.5061225 0.5062421 0.5063027 0.5063401 0.5063692
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $nll
    ## [1] -33.03797
    ## 
    ## $vcov
    ##               [,1]          [,2]          [,3]
    ## [1,]  1.423093e-03 -1.749936e-02 -4.229697e-08
    ## [2,] -1.749936e-02  2.602346e-01  8.139398e-07
    ## [3,] -4.229697e-08  8.139398e-07  1.219567e-02

Fit the logistic model to data (GBherring):

``` r
GBherring <- assess(harvest.GBherring,index.GBherring,calc.vcov=TRUE,ini.parms.GBherring)
GBherring
```

    ## $pars
    ## [1]  2.0188659  0.2364302 -2.1480052
    ## 
    ## $biomass
    ##  [1] 7.529781 7.527551 7.528080 7.527911 7.527945 7.527932 7.527935 7.527934
    ##  [9] 7.527934 7.527934 7.527934 7.527934 7.527934 7.527934 7.527934 7.527934
    ## [17] 6.350312 6.512945 6.561671 6.577174 6.582462 5.662317 5.609675 5.643432
    ## [25] 5.676540 5.696087 4.357552 4.010960 3.873795 3.808987 3.771620 5.851700
    ## [33] 7.110337 7.174143 7.146059 7.142422 7.140206 7.139688 7.139599 7.139644
    ## [41] 7.139702
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $nll
    ## [1] -29.88969
    ## 
    ## $vcov
    ##               [,1]          [,2]          [,3]
    ## [1,]  4.192576e-04 -7.178308e-04 -9.523352e-07
    ## [2,] -7.178308e-04  1.748822e-03  4.168951e-06
    ## [3,] -9.523352e-07  4.168951e-06  1.219391e-02

Fit the logistic model to data (GBhaddock):

``` r
GBhaddock <- assess(harvest.GBhaddock,index.GBhaddock,calc.vcov=TRUE,ini.parms.GBhaddock)
GBhaddock
```

    ## $pars
    ## [1]  1.8535726  0.1875406 -2.3782023
    ## 
    ## $biomass
    ##  [1] 6.382582 6.333023 6.342806 6.340995 6.341386 6.341354 6.341398 6.341422
    ##  [9] 6.341445 6.341465 6.341482 6.341497 6.341510 6.341521 6.341532 6.341542
    ## [17] 6.340236 6.340428 6.340287 6.340202 6.340116 6.338714 6.338824 6.338615
    ## [25] 6.338473 6.338338 6.335520 6.335750 6.335330 6.335046 6.334778 6.339565
    ## [33] 6.338738 6.339150 6.339351 6.339568 6.339742 6.339884 6.339996 6.340086
    ## [41] 6.340159
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $nll
    ## [1] -39.32879
    ## 
    ## $vcov
    ##               [,1]          [,2]          [,3]
    ## [1,]  2.989882e-04 -1.418830e-02 -1.006583e-07
    ## [2,] -1.418830e-02  2.186002e+00  5.287660e-06
    ## [3,] -1.006583e-07  5.287660e-06  1.219450e-02

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
  facet_wrap(~Species)
```

![](MultispeciesMSE_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

# Design HCRs

Separate single species HCRs (start with Gavin’s specified as default)

# Project with HCRs

# Compare results

# Tradeoffs
