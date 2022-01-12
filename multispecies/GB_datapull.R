#Generate data for MSE project
library(data.table); library(here); library(Rpath)

#Load GB params
load(here('multispecies', 'data-raw', 'GB_balanced_params.RData'))

#Generate balanced model
GB <- Rpath::rpath(GB.params, 'Georges Bank')

#Add variable fishing in a dynamic run
GB.scene <- rsim.scenario(GB, GB.params, years = 1983:2022)
GB.scene <- adjust.fishing(GB.scene, parameter = 'ForcedEffort', group = 'OtterTrawlLg',
                            sim.year = 1983:2017, value = c(rep(0, 15), rep(2.5, 5),
                                                            rep(5, 5), rep(10, 5),
                                                            rep(2.5, 5)))
GB.scene <- adjust.fishing(GB.scene, parameter = 'ForcedEffort', group = 'Gillnet',
                           sim.year = 1983:2012, value = c(rep(0, 10), rep(1.5, 5),
                                                           rep(2, 5), rep(4, 5),
                                                           rep(2.5, 5)))
GB.run <- rsim.run(GB.scene, years = 1983:2022)

#Plot results
GB.bio <- as.data.table(GB.run$out_Biomass)
GB.bio[, Month := 1:nrow(GB.bio)]
GB.bio <- data.table::melt(GB.bio, id.vars = 'Month', variable.name = 'Group', 
                           value.name = 'Biomass')
start <- data.table::as.data.table(GB.run$start_state$Biomass)
start[, Group := names(GB.run$start_state$Biomass)]

GB.bio <- merge(GB.bio, start, by = 'Group')
GB.bio[, Rel.biomass := Biomass / V1]

ggplot(GB.bio[Group %in% c('Cod', 'Haddock', 'AtlHerring')],
       aes(x = Month, y = Rel.biomass, col = Group)) +
  geom_line()

#Output data 
GenData <- function(Rsim.output, group, sim.year, Sigma = 0.3, bias = 1, freq = 1){
  # This routine generates data (unbiased generally) for the biological groups 
  #and returns the appropriate object
  
  node     <- extract.node(Rsim.output, group)
  TrueBio  <- node$AnnualBiomass[which(names(node$AnnualBiomass) == sim.year)]
  Observed <- c()
  for(Iyear in seq_along(sim.year)){
    if (Iyear %% freq == 0){
      Observed[Iyear] <- TrueBio[Iyear] * exp(rnorm(1, 0, Sigma) - Sigma^2/2)
    } else Observed[Iyear] <- -1
  }
  Catch <- as.numeric(node$AnnualTotalCatch[which(names(node$AnnualBiomass) == sim.year)])
  
  out <- list(ObsBio = Observed, TotCatch = Catch, Fmort = Catch / Observed)
  
  return(out)
} 

cod <- GenData(GB.run, 'Cod', 1983:2022)  
haddock <- GenData(GB.run, 'Haddock', 1983:2022)
atlherring <- GenData(GB.run, 'AtlHerring', 1983:2022)

#Merge into 1 file
GB.data <- data.table(Year = 1983:2022,
                      Species = c(rep('Cod', 40), rep('Haddock', 40), 
                                  rep('AtlHerring', 40)),
                      ObsBio = c(cod$ObsBio, haddock$ObsBio, atlherring$ObsBio),
                      TotCatch = c(cod$TotCatch, haddock$TotCatch, atlherring$TotCatch),
                      Fmort = c(cod$Fmort, haddock$Fmort, atlherring$Fmort))

save(GB.data, file = here('multispecies', 'data', 'GBdata.rda'))





