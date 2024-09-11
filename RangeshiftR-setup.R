# Parkages required to install RangeShiftR 

# install.Rtools()
# install.packages("Rcpp")
# devtools::install_github("https://github.com/RangeShifter/RangeShiftR-pkg", ref = "main", subdir="RangeShiftR")

library(Rcpp)
library(devtools)
library(Rdpack)
library(RangeShiftR)

# The basics if running RangeShiftR
# 1.1. Running a first working example

# Create a parameter master object with all the default settings and store it:
s <- RSsim()

# If you have already created a suitable directory for your simulation on your
# disc, store its path in a variable. This can either be the relative path from
# your R working directory or the absolute path.

dirpath = "Overview/"

# Create the RS folder structure, if it doesn’t yet exist:
dir.create(paste0(dirpath,"Inputs"), showWarnings = TRUE)
dir.create(paste0(dirpath,"Outputs"), showWarnings = TRUE)
dir.create(paste0(dirpath,"Output_Maps"), showWarnings = TRUE)

# With this, we are already set to run our first simulation by typing:
  
RunRS(s,dirpath)

# 2 Simulation modules
# To look at the parameter master in more detail, simply type:
s
# 2.1 Simulation
# This module is used to set general simulation parameters (e.g. simulation ID,
# number of replicates, and number of years to simulate) and to control output
# types (plus some more specific settings). For this overview, we will stick to
# the defaults:
sim <- Simulation(Simulation = 2,
                  Years = 50,
                  Replicates = 2,
                  OutIntPop = 50)
# ?Simulation

# 2.2 Landscape
# RangeShiftR can either import a map from an ASCII raster file in the ‘Inputs’ folder or generate a random map to use in the simulation.

# For each option, there is a corresponding function to create a Landscape parameter object

# ?ImportedLandscape
# ?ArtificialLandscape

land <- ArtificialLandscape(Resolution = 10,  # in meters
                            K_or_DensDep = 1500,  # ~ 15 inds/cell
                            propSuit = 0.2,
                            dimX = 129, dimY = 257, 
                            fractal = T, hurst = 0.3,
                            continuous = F)

# 2.3 Demography
demo <- Demography(Rmax = 2.2, ReproductionType = 1, PropMales = 0.45)

stg <- StageStructure(Stages = 3,
                      TransMatrix = matrix(c(0,1,0,5.7,.5,.4,3.4,0,.9),nrow = 3),
                      FecDensDep = T,
                      SurvDensDep = T)
demo <- demo + stg
# Alternatively, we define the sub-module within the Demography module:
# demo <- Demography(StageStruct = stg, ReproductionType = 1, PropMales = 0.45)

# RangeShiftR provides a number of useful functions to explore the model set-up. For example, we can plot the rates from the transition matrix:
  
plotProbs(stg)
str(stg)

#  2.4. Dispersal
# The dispersal process is modelled wih three sub-processes 
# (see the schematic figure above): Emigration(), Transfer() and Settlement().

disp <-  Dispersal(Emigration = Emigration(EmigProb = 0.2), 
                   Transfer   = DispersalKernel(Distances = 50),
                   Settlement = Settlement() )
# We can use the function plotProbs() to plot various functional relationships, 
# for example a dispersal kernel with a stage-dependent mean transfer distance:
plotProbs(DispersalKernel(Distances = matrix(c(0,1,2,70,50,30),nrow = 3), 
                          StageDep = T))
#  This is using mostly default options. For example, we can change the settlement 
# condition so that a female individual, that arrives in an unsuitable cell, 
# will wait for another time step and disperse again, while the males will die 
# if arriving in an unsuitable cell:

disp <-  disp + Settlement(SexDep = T, 
                           Settle = matrix(c(0,1,1,0), nrow = 2))
# The genetics module controls the heritability and evolution of traits and is 
# needed if inter-individual variability is enabled (IndVar = TRUE) for at least 
# one dispersal trait.
gene <- Genetics(NLoci = 3, ProbMutn = .05)

 
# 2.6 Initialisation
# In order to control the initial distribution of individuals in the landscape 
# at year 0, we set initialisation rules. We choose to initialise 3 individuals 
# per habitat cell. Additionally, since we define a stage-structured model, 
# we have to specify the initial proportion of stages:

init <- Initialise(FreeType = 0, 
                   NrCells = 2250,
                   InitDens = 2, 
                   IndsHaCell = 3, 
                   PropStages = c(0,0.7,0.3))
init

# 2.7 Parameter master
# After all settings have been made in their respective modules, we are ready 
# to combine them to a parameter master object, which is needed to run the simulation.
s <- RSsim(simul = sim, land = land, demog = demo, dispersal = disp, gene = gene, 
           init = init)

# Alternative notation:
# s <- RSsim() + land + demo + disp + sim + init + gene

# We can check the parameter master (or any single module) for potential parameter conflicts:

validateRSparams(s)

# 2.8 Run the simulation
# Once the parameter master has been defined, we can run the simulations in the specified RS directory.

RunRS(s, dirpath)

# 3 Plot results
# All results are stored in the Outputs folder. RangeShiftR provides some 
# in-built functions to access and plot these results. Here, we plot the abundance and occupancy time series:
range_df <- readRange(s, dirpath)

# ...with replicates:
par(mfrow=c(2,3))
plotAbundance(range_df)
plotOccupancy(range_df)
# ...with standard deviation:
# par(mfrow=c(1,2))
plotAbundance(range_df, sd=T, replicates = F)
plotOccupancy(range_df, sd=T, replicates = F)










