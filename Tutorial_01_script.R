# Tutorial 1: Range expansion, long-distance dispersal and environmental stochasticity
# 1 Simulating range expansions
# 1.1 Create an RS directory
# The directory in which we run RangeShiftR needs to have a certain folder structure. 
# It should contain three sub-folders named ‘Inputs’, ‘Outputs’ and ‘Output_Maps’. 
# We can create them from R:
# The following packages are required, if the function's package is not installed, 
# Use the function ´findFn("viridis")´ to find out wghich package it belongs to.



library(RangeShiftR)
library(raster)
library(RColorBrewer)
library(rasterVis) # from rasterVis package
library(latticeExtra) # from rasterVis package
library(viridis)
library(grid)
library(gridExtra)

# relative path from working directory:
dirpath = "Tutorial_01/"

# Create the folders required to run RangeShifter, If folder already exist in the directory, 
# you will be worned
dir.create(paste0(dirpath,"Inputs"), showWarnings = TRUE) 
dir.create(paste0(dirpath,"Outputs"), showWarnings = TRUE)
dir.create(paste0(dirpath,"Output_Maps"), showWarnings = TRUE)

# Copy the input files provided for exercise 1 into the ‘Inputs’ folder. The files can be downloaded here.

# 1.2 Landscape parameters

# We use a land-cover map of Great Britain at 1km resolution. Six dominant aggregated 
# habitat types were derived from LandCover Map 2007 (Morton et al. 2011).
# The map, UKmap_1km.txt, is an ASCII raster in the standard text format, 
# where each cell holds the code of its dominant habitat type. The habitat codes have 
# to be given as sequential integer numbers, starting from one:

# 1 = woodland (broadleaved and conifer)
# 2 = arable
# 3 = improved grassland
# 4 = semi-natural grassland (acid, neutral and calcareous grassland)
# 5 = heath and bog
# 6 = other (urban, water & coastal habitats)


# Furthermore, we are provided with a map that defines the (initial) species distribution, named Species_Distribution_10km.txt. It is given at a coarser resolution of 10km.

# Let’s plot the landscape map and overlay it with the initial species distribution:

UKmap <- raster(paste0(dirpath, "Inputs/UKmap_1km.txt"))
SpDist <- raster(paste0(dirpath, "Inputs/Species_Distribution_10km.txt"))
values(SpDist)[values(SpDist) < 1] <- NA # This will look for values that are less than one 
# in the raster data specified i.e. spDis and replace them with NA,  meaning no data

# plot land cover map and highlight cells with initial species distribution - option 1:
plot(UKmap, col=brewer.pal(n = 6, name = "Spectral"), axes=F)
plot(rasterToPolygons(SpDist, dissolve=F), add=T)

# For prettier mapping with a legend specifying the different land cover categories 
# and a custom colour palette, we need a little workaround:

# plot land cover map and highlight cells with initial species distribution - option 2 with categorical legend:
UKmap.f <- as.factor(UKmap)
# add the land cover classes to the raster attribute table (RAT)
rat <- levels(UKmap.f)[[1]]
rat[["landcover"]] <- c("woodland", "arable", "improved grassland", "semi-natural grassland", "heath and bog", "other")
levels(UKmap.f) <- rat

custom.pal <- c("#1A9850", "#91CF60", "#D9EF8B", "#FEE08B", "#D8B365", "#777777")
levelplot(UKmap.f, margin=F, scales=list(draw=FALSE), col.regions=custom.pal)  +
  layer(sp.polygons(rasterToPolygons(SpDist, dissolve=F), fill=NA, col='red'))


# In order to use the habitat and distribution maps in RangeShiftR, we need to set up a 
# landscape module. It takes the path to the map files along with some other parameters, check the manual

carrycap <- c(0, 0, 0, 5, 0, 0)
land <- ImportedLandscape(LandscapeFile = "UKmap_1km.txt", 
                          Resolution = 1000, 
                          Nhabitats = 6, 
                          K_or_DensDep = carrycap, 
                          SpDistFile = "Species_Distribution_10km.txt", 
                          SpDistResolution = 10000)
# 1.3 Next we specify the species parameters for the demography and dispersal modules.
# Wer choose a simple model of femals only, we only need to define maximum growth rate ´RMax´

demo <- Demography(Rmax = 1.5)

# Set the The dispersal module comprises three sub-modules representing the three phases of dispersal: 
# Emigration, Transport and Settlement.


# We assume a constant emigration probability of 0.1. 
# The transfer phase is modelled with a dispersal kernel, whose mean distance is set to 2,000m. 
# For the Settlement sub-module, we use the default options, 
# which assume that an individual will die if it arrives in an unsuitable cell and settle if it’s suitable.


disp <-  Dispersal(Emigration = Emigration(EmigProb = 0.1), 
                   Transfer = DispersalKernel(Distances = 2000), 
                   Settlement = Settlement() )

# alternative notation:
# disp <-  Dispersal() + Emigration(EmigProb = 0.1) + DispersalKernel(Distances = 2000) + Settlement()


# 1.4 Initialisation parameters

init <- Initialise(InitType = 1, # = initialisation from a loaded species distribution map
                   SpType = 0,   # = all suitable cells within all distribution presence cells
                   InitDens = 0) # = at carrying capacity

# 1.5 Simulation parameters
# Lastly, we set some basic simulation parameters, i.e. the simulation number, 
# number of replicates and number of years to be simulated. 
# Furthermore, we specify what file output will be generated. In this example, 
# we choose to output the population, occupancy and range data for every 5 years.

# Furthermore, optionally define a (Shifting) Environmental Gradient,
# Environmental Stochasticity (EnvStoch) and/or Local extinction (LocalExt).
# (These options are to be moved to separate classes in future versions.)


sim_0 <- Simulation(Simulation = 0, 
                    Replicates = 5, 
                    Years = 30,
                    OutIntPop = 5,
                    OutIntOcc = 5,
                    OutIntRange = 5)
                    # Gradient = 0, GradSteep, Optimum, f, ExtinctOptim,
                    # Shifting = FALSE, ShiftRate, ShiftStart, ShiftEnd,
                    # LocalExt = FALSE, LocalExtProb,
                    # EnvStoch = 0, EnvStochType, std, ac, minR, maxR, minK, maxK,
                    # OutIntRange = 1, OutIntOcc = 0,
                    # OutIntPop = 1, OutIntInd = 0,
                    # OutIntTraitCell = 0, OutIntTraitRow = 0,
                    # OutIntConn = 0, OutIntPaths = 0, OutIntGenetic = 0,
                    # OutGenType = 0, OutGenCrossTab = FALSE,
                    # OutStartPop = 0, OutStartInd = 0,
                    # OutStartTraitCell = 0, OutStartTraitRow = 0,
                    # OutStartConn = 0, OutStartPaths = 0, OutStartGenetic = 0,
                    # SMSHeatMap = FALSE)

# 1.6 Parameter master
# All settings we have made so far are contained in their respective module objects. 
# They have to be combined to a parameter master object, which validates all parameter settings 
# and is needed to run the simulation

s <- RSsim(land = land, demog = demo, dispersal = disp, simul = sim_0, init = init)

# alternative notation:
# s <- RSsim() + land + demo + disp + sim_0 + init

# At any stage, you can check for parameter validity of any module or the master.

validateRSparams(s)

# 1.7 Run the simulation
RunRS(s, dirpath)

# 1.8 Plot results
# For convenience and for compliance with the RangeShifter Windows GUI, the RangeShiftR package contains a few plotting functions to inspect the population time series of the replicate simulations.

# read 'range' output into a data frame
range_df <- readRange(s, dirpath)

# plot trajectories of all individual runs and overlay with mean:
par(mfrow=c(2,2))
plotAbundance(range_df)
plotOccupancy(range_df)

# or plot mean and standard deviation:
# par(mfrow=c(1,2))
plotAbundance(range_df, rep=F, sd=T)
plotOccupancy(range_df, rep=F, sd=T)


