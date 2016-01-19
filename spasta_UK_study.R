##################################################################################
## Date: 05 January 2016
## Author: Andrew Zammit-Mangion
## Title: Reproducible code for Section 5 of 'Trans-Gaussian bivariate 
## modelling with application to atmospheric trace-gas inversion' by Zammit-Mangion
## et al.
##
## Notes: Follow the instructions on https://github.com/andrewzm/atminv
## It is suggested that a server-grade machine is used since multiple MCMC chains
## are used in this program (see options for altering this below). For full
## reproducibility create a folder img/ in which to create images, tempresults/
## in which to store the results and ~/cache/Assim for temporarily storing results.
##################################################################################

library(dplyr)
library(tidyr)
library(Matrix)
library(ggplot2)        # for plotting
library(grid)           # for plotting
library(gridExtra)      # for plotting
library(atminv)         # https://github.com/andrewzm/atminv
library(hmc)            # https://github.com/andrewzm/hmc
library(linalg)         # https://github.com/andrewzm/linalg
library(foreach)        # for parallel computations
library(doRNG)          # for parallel computations
library(doMC)           # for parallel computations

#-------------------------------------------------------------------
# 1. Program setup
#------------------------------------------------------------------


rm(list=ls())

save_images   <- 0      # save images?
N             <- 12000  # number of MCMC samples per chain
adapt         <- 1000   # number of adaptation samples
burnin        <- 8000   # number of burnin samples (incl. adaptation)
dither        <- 10     # thinning factor
n_parallel    <- 10     # number of parallel chains
do_inference  <- 0      # do inference or load results?
true_inventory<- 1      # assume true inventory or wrong inventory?
real_data     <- 0      # use real data or simulated data?
model <- "Gaussian.uncorrelated" # Can be one of the following:
# Box-Cox, Lognormal, Gaussian,   Box-Cox.uncorrelated, 
# Lognormal.uncorrelated or Gaussian.uncorrelated

if(do_inference)
  print(paste0("Doing inference with Model ",model," with ", N," samples and ",n_parallel,
             " parallel chains, true_inventory = ",true_inventory, " and real data = ",real_data))

select <- dplyr::select  # ensure select is the dplyr select

md5_wrapper <- md5_cache("~/cache/Assim") # set up cache
set.seed(1) # fix seed

miss_val = -999  # missing value in datasets
ymin = 36.469    # minimum lon coord
xmin = -14.124   # minimum lat coord
dx = 0.352       # lon spacing
dy = 0.234       # lat spacing
dx_new <- dx*2   # downsampled lon spacing
dy_new <- dy*2   # downsampled lat spacing
yx <- expand.grid(y = seq(ymin, ymin + 127*dy,length=128),       # latlon grid in data
                  x = seq(xmin, xmin + 127*dx,length=128)) %>%
    as.data.frame()

yx$breaks_x <- cut(yx$x,   # super x-grid for downsampling by 2
                   seq(xmin-0.001, xmin + 127*dx + 0.001,length=64),
                   labels=FALSE)
yx$breaks_y <- cut(yx$y,   # super y-grid for downsampling by 2
                   seq(ymin-0.001, ymin + 127*dy + 0.001,length=64),
                   labels=FALSE)

### Load maps

# Load world map and subset to interesting countries
data("world",package = "atminv")
Europe <- world %>%
    subset(region %in% c("UK","France","Ireland")) %>%
    mutate(id = group, x=long,y=lat)

# Attribute each grid-cess to a country and sea/land indicator
yx <- mutate(yx, in_land= attribute_polygon(yx,Europe),
             id = in_land) %>%
    left_join(unique(dplyr::select(Europe,region,id,subregion)))

# Add UK territory (including sea areas)
UK_idx <- read.table("../data/UK_gridcells.txt")[,1]
yx$UK_territory <- 0
yx$UK_territory[UK_idx] <- 1

# Add Ireland territory (including sea areas)
Ir_idx <- read.table("../data/Ireland_gridcells.txt")[,1]
yx$Ir_territory <- 0
yx$Ir_territory[Ir_idx] <- 1

# Model Emissions field (use the one scaled by land-use statistics)
yx <- yx %>%
    mutate(z = read.table("../data/CH4_emissions_scaled_natural",skip=10)[,1], t = 1)

# Create a separate variable `Emissions` which subsamples the emissions
# Note that we sum the emissions in each grid cell and average the percentage of
# UK territory
Emissions <- group_by(yx,breaks_x,breaks_y) %>%
    summarise(y = mean(y),
              x = mean(x),
              in_land = max(in_land),
              region = region[1],
              subregion = subregion[1],
              z=sum(z),
              UK_territory =mean(UK_territory))

## Filter out the UK
Emissions_land <- subset(Emissions,
                         region %in% c("UK","Ireland"))

## Function to visualise a field of the data frame Emissions_map
Emissions_map <- function(field,ll=-2000,ul=2000) {
    Emissions_land$ul <- ul
    Emissions_land$ll <- ll
    Emissions_land$to_plot <- unlist(Emissions_land[field])
    LinePlotTheme() +
        geom_tile(data=subset(Emissions_land,x > -13 & x < 3 & y > 48 & y < 61),
                  aes(x,y,fill=pmax(pmin(to_plot,ul[1]),ll[1])))  +
        scale_fill_distiller(palette="Spectral",
                             guide = guide_legend(title="flux (g/s)")) +
        geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
        coord_map(xlim=c(-12,2),ylim=c(49,60)) +
        xlab("lon (degrees)") + ylab("lat (degrees)") +
        theme(axis.title.y = element_text(vjust = 1))
}
print(Emissions_map("z"))


yx_land <- subset(yx,region %in% c("UK","Ireland"))       # consider land cells
Dist_mat <- fields::rdist(Emissions_land[,c("x","y")])    # find distance matrix
flux_cov_fn <- function(D,scale,nu) exp(-scale*(D^nu))    # flux cov. function (Gauss scale)

## Check flux field with geoR
library(geoR)
geo_obj <- as.geodata(Emissions_land,coords.col = c("x","y"),data.col = "z")
likfit(geo_obj,ini=c(0.5,0.5),fix.lambda=FALSE,fix.nugget=T) %>%
    summary() %>% print()


#-------------------------------------------------------------------
# 2. Data preprocessing
#------------------------------------------------------------------

### Create a list `Stations` with filenames of models/coordinates and
### station coordinates
Stations <-list(MHD = list(model_filenames = c("../data/MHD_model_012014.txt",
                                               "../data/MHD_model_022014.txt",
                                               "../data/MHD_model_032014.txt"),
                           obs_filenames = c("../data/MHD_obs_012014_filtered.txt",
                                             "../data/MHD_obs_022014_filtered.txt",
                                             "../data/MHD_obs_032014_filtered.txt"),
                           x = -9.90,
                           y = 53.33),
                RGL = list(model_filenames = c("../data/RGL_model_012014.txt",
                                               "../data/RGL_model_022014.txt",
                                               "../data/RGL_model_032014.txt"),
                           obs_filenames = c("../data/RGL_obs_012014_filtered.txt",
                                             "../data/RGL_obs_022014_filtered.txt",
                                             "../data/RGL_obs_032014_filtered.txt"),
                           x = -2.54,
                           y = 52,00),
                TAC = list(model_filenames = c("../data/TAC_model_012014.txt",
                                               "../data/TAC_model_022014.txt",
                                               "../data/TAC_model_032014.txt"),
                           obs_filenames = c("../data/TAC_obs_012014_filtered.txt",
                                             "../data/TAC_obs_022014_filtered.txt",
                                             "../data/TAC_obs_032014_filtered.txt"),
                           x = 1.14,
                           y = 52.52),
                TTA = list(model_filenames = c("../data/TTA_model_012014.txt",
                                               "../data/TTA_model_022014.txt",
                                               "../data/TTA_model_032014.txt"),
                           obs_filenames = c("../data/TTA_obs_012014_filtered.txt",
                                             "../data/TTA_obs_022014_filtered.txt",
                                             "../data/TTA_obs_032014_filtered.txt"),
                           x = -2.99,
                           y = 56.56))

##We have m_obs = 4 stations:
m_obs <- length(Stations)

### Read model data from file and put into long table format and focus on UK
### I also should be multiplying by 1e9 but this produces VERY large numbers
### Below each `des` is an element of `Station`
read_process_models <- function(des,prune=T) { 
    Model_B_list <- list()                      # initialise list of SRRs
    for(i in 1:length(des$model_filenames)) {   # for each month of model output
        print(paste0("Reading ",des$model_filenames[i]))
        Model_B_list[[i]] <- read.table(des$model_filenames[i],  # read the data
                                        skip = 15) %>%
            cbind(select(yx,-t,-UK_territory,    # bind the data with some columns of yx
                         -Ir_territory,-z)) %>%
            gather(t,z,-x,-y, -in_land,          # put into long format (only put  
                   -region,-id,-subregion,       # the columns associated with time
                   -breaks_x,-breaks_y) %>%
            separate(t, into = c("V", "t"),      # now separate the V1,V2 etc. in 
                     sep = "V") %>%              # the column `t`
            mutate(t = as.numeric(t)) %>%        # now make `t` numeric and drop the `V`
            select(-V) %>%
            group_by(breaks_x,breaks_y,t) %>%    # group by the new sampling resolution
            summarise(y = mean(y), x = mean(x),  # average the SRR at this new resolution
                      in_land = max(in_land), 
                      region = region[1],
                      subregion = subregion[1],
                      z=mean(z))
        
        # If this is the second (or more) file then shift the time axis accordingly
        if(i > 1) {
            Model_B_list[[i]]$t <- Model_B_list[[i-1]]$t[nrow(Model_B_list[[i-1]])] + 
                Model_B_list[[i]]$t
        }
    }
    
    # Now append all the newly constructed data frames together
    des$Model_B <- do.call("rbind",Model_B_list)
    
    # From the long data frame des$Model_B we now want to make the big B matrix
    # To to this we cast to a data frame (from a grouped data frame), then we 
    # sort by time then lon then lat.
    # For each spatial location extract the SRR over all time as a row
    # This gives as a data frame (x,y,t1,t2,t3,...)
    # We then sort again by lon and lat so that we make sure it matches with
    # our imposed ordering convention and then drop the x,y and convert to matrix
    # We finally transport the matrix to obtain B which is now nt * ns (multiplied by
    # Yf which is of dimension ns)
    des$B <- des$Model_B %>%
        as.data.frame() %>%
        arrange(t,x,y) %>%
        plyr::ddply(c("x","y"),function(l) { l$z}) %>%
        arrange(x,y) %>%
        select(-x,-y) %>%
        as.matrix() %>%
        t()
    
    des$C <- des$B ## For now assume we have observations everywhere in space and time
    des
    
    # The 44/16 converts Nitrous Oxide  to methane mass
    # The 1e9 is to convert mol/mol to nmol/mol
    # = ppb (parts per billion)
}

### Read obs data from file and put into long table format and 
### cast as Obs object. Here `des` is an element of `Stations`
read_process_data <- function(des,t_axis) {
    
    this_x <- des$x  # Station lon
    this_y <- des$y  # Station lat
    
    # find which emissions in UK/Ireland
    id <- which(Emissions$region %in% c("UK","Ireland"))       
    
    # find the mole-fraction contributions everything outside UK/Ireland
    border_z <- as.numeric(des$B[,-id] %*% Emissions[-id,]$z)  
    
    # now we create a list of the observations
    Obs_df_list <- list()
    
    # for each observation dataset (monthly)
    for(i in 1:length(des$obs_filenames)) {
        Obs_df_list[[i]] <-
            read.table(des$obs_filenames[i],skip=1) %>%     # read data from file
            mutate(z = V1, std = V2,                        # copy columns V1 and V2
                   x = this_x, y = this_y) %>%              # add x and y station locs
            select(-V1,-V2)                                 # remove V1 and V2       
    }
    
    # create a long-format dataframe from the above
    Obs_df <- do.call("rbind",Obs_df_list) %>%  # merge above data frames into one long one
        mutate(border_z = border_z,             # add on the border contributions
               t= t_axis,                       # add on the time axis
               z = ifelse(z == miss_val,NA,z))  # use NA for missing values
    
    des$Obs_obj <- Obs(Obs_df,                    # put into Observation object
                       name=des$obs_filenames[i]) # use last filename as name
    
    x <- rep(0,length(t_axis))   # create a vector
    x[des$Obs_obj["t"]] <- 1     # which indicates which time points are observed
    des$C <- diag(x)             # create a diagonal matrix with all-zero rows for
    # unobserved time points
    
    ##  To normalise (deprecated)
    #   if(match_to_process) {
    #     des$Obs_obj["z"] <- des$Obs_obj["z"] - quantile(des$Obs_obj["z"],0.05,na.rm=T)
    #     Zpred <- log(abs(C2 %*% des$B %*% Emissions$z))
    #     z <- log(abs(des$Obs_obj["z"]))
    #
    #     des$Obs_obj["z"] <- exp(sd(Zpred)/sd(z) * (z - mean(z)) + mean(Zpred))
    #     warning("Still not adjusting observation error in transformation")
    #     #des$Obs_obj["std"] <- exp(log(des$Obs_obj["std"]) * sd(Zpred)/sd(z))
    #  }
    des
}


## Read process models and form matrices in station list
Stations <- md5_wrapper(lapply,Stations,read_process_models,prune=FALSE)

## Read data and apply to station list
t_axis <- 1:max(Stations$RGL$Model_B$t)
Stations <- lapply(Stations,read_process_data,t_axis)

## Find the 0.05 quantile at MHD = 1885.251
MHD_5 <- quantile(Stations$MHD$Obs_obj["z"],0.05,na.rm=T)

## Remove MHD_5 from all mole-fraction observations
Stations <- lapply(Stations,function(l) {l$Obs_obj["z"] <- l$Obs_obj["z"] - MHD_5; l})

## Construct a data frame of station locations with rownames as station names
Station_locs <- t(sapply(Stations,function(l) data.frame(x=l$x,y=l$y))) %>% 
    as.data.frame()
Station_locs$name <- row.names(Station_locs)

## Assign the first SRR (1 Jan 00:00) at station TTA to the data frame Emissions
Emissions$B <- Stations$TTA$B[1,]
Emissions$B_MHD1 <- Stations$MHD$B[1,]
Emissions$B_RGL1 <- Stations$RGL$B[1,]
Emissions$B_TTA1 <- Stations$TTA$B[1,]
Emissions$B_TAC1 <- Stations$TAC$B[1,]

# Plot the SRR on 01 Jan
if(save_images) {
  g_SRR <- LinePlotTheme() + geom_tile(data=subset(Emissions,x > -13 & x < 3 & y > 48 & y < 61),
                                   aes(x,y,fill=B_TTA1 + B_TAC1 + B_RGL1 + B_MHD1))  +
      #scale_fill_distiller(palette="BuPu",trans=c("sqrt"),guide = guide_legend(title="SRR (s/ng)")) +
      scale_fill_distiller(palette="Greys",trans=c("sqrt"),
                           guide = guide_legend(title=expression(bar(bold(B))[1] (s/ng))),direction=2) +
      geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
      coord_map(xlim=c(-12,2),ylim=c(49,60)) + xlab("lon (degrees)") + ylab("lat (degrees)\n") +
      geom_point(data=Station_locs,aes(x=unlist(x),y=unlist(y)),col="black",size=4) +
      geom_text(data=Station_locs,aes(x=unlist(x),y=unlist(y)+0.3,label=unlist(name)),col="black") +
      theme(axis.title.y = element_text(vjust = 1))
    ggsave("./img/SRR_01Jan.png",plot=g_SRR,width=6.7,height=6.7)
}

# Emissions_map("z") + 
#     scale_fill_distiller(palette="Greys",trans="reverse", guide = guide_legend(title="flux (g/s)")) +
#     geom_point(data=Station_locs,aes(x=unlist(x),y=unlist(y)),size=4) +
#     geom_text(data=Station_locs,aes(x=unlist(x),y=unlist(y)+0.3,label=unlist(name)),col="black")


## Now we do the pruning so we focus on the UK.
## First we remove the "background Europe" from the data
Stations <- lapply(Stations,function(l) {
    l$Obs_obj["z"] <- l$Obs_obj["z"] -
        l$Obs_obj["border_z"]; l
})

## Now we prune all matricies. First find the IDs we are going to focus on
id <- which(Emissions$region %in% c("UK","Ireland"))

## Now subset the emissions to these grid cells
Emissions_all <- Emissions
Emissions <- Emissions[id,]

## And only consider the SRRs on these grid cells
Stations <- lapply(Stations,function(l) {  l$B <- l$B[,id]; l})

#-------------------------------------------------------------------
# 3. Model + matrix construction + simulate observations
#------------------------------------------------------------------

## Now construct the big matrices we'll work with.
## First we start with B. Since we are organising our data as
## Obs1t1,Obs2t,...,Obs1t2,Obs2t2,... we have to interleave the B matrices:
## B <- gdata::interleave(as.matrix(Stations$MHD$B),as.matrix(Stations$RGL$B),...)
B <- do.call(gdata::interleave,args = lapply(Stations, function(l) as.matrix(l$B)))

## And the incidence matrices which, recall, are full diagonal matrices since
## all observations (even unobserved time points) are included, just set to NA
## We will filter out the unwanted rows after interleaving:
## Currently this is just an elaborate way of getting one big identity matrix
## but might be useful in the future when we have different observation weights
C <- do.call(gdata::interleave,args = lapply(Stations,
                                             function(l) as.matrix(diag(l$C)))) %>%
    c() %>% diag()
C <- as(C,"dgCMatrix") ## convert to sparse matrix (actually diagonal)


obs_locs <- t(sapply(Stations,                # extract station locs
                     function(l) c(l$x,l$y)))

Q_norm_zeta_fn <- function(s,t,theta_t,theta_s) {
  Q_zeta_t <- GMRF_AR(n=length(t),a=theta_t)@Q
  chol_Q_zeta_t <- chol(Q_zeta_t)

  corr_zeta_s <- atminv:::corr_s(s, theta_s)
  chol_Q_zeta_s <- chol(corr_zeta_s)
  Q_zeta_s <- chol2inv(chol_Q_zeta_s)

  Q_norm_zeta <- crossprod(kronecker(chol_Q_zeta_t, chol_Q_zeta_s))
}

if(!real_data) {

    theta_t_true <- 0.9     # temporal correlation parameter ('a' in text)
    theta_s_true <- 2.5       # spatial range parameter ('d' in text)
    sigma_zeta_true <- 10

    Q_zeta_true <- sigma_zeta_true^(-2) * Q_norm_zeta_fn(obs_locs,
                                                         t_axis,
                                                         theta_t_true,
                                                         theta_s_true)

    GMRF_true <- GMRF(Q=Q_zeta_true)

    ## Now find which observations we would keep/throw away
    s_obs <- do.call(gdata::interleave,args = lapply(Stations,
                                                     function(l) l$Obs_obj@df))
    idx <- which(s_obs$z > 0 & !is.na(s_obs$z))

    message("Setting measurement error to 1")
    Stations <- lapply(Stations,
                       function(x) {
                           x$Obs_obj["z"] <- as.numeric(x$C[which(!(colSums(x$C) == 0)),] %*%
                                                            x$B %*% Emissions$z)
                           x})

    ## Now create observations
    s_obs <- Reduce("rbind",lapply(Stations,function(l) l$Obs_obj@df)) %>%
        arrange(t) %>%
        mutate(std=1)
    s_obs$std[which(s_obs$std == -999)] <- mean(subset(s_obs,std>0)$std)

    s_obs <- mutate(s_obs,
                    dis = as.numeric(C %*% sample_GMRF(GMRF_true)),
                    obs_err =  rnorm(n = nrow(C),sd = std),
                    Ymf = z,
                    z = z + dis + obs_err)

    ## Create missing observations
    s_obs_orig <- s_obs
    C <- C[idx,]
    s_obs <- s_obs[idx,]
} else {
    ## Interleave observation data frames (recall organisation: O1t1,O2t1,...)
    s_obs <- do.call(gdata::interleave,args = lapply(Stations,
                                                     function(l) l$Obs_obj@df))

    ## Replace missing std values with the "mean error"
    message("Replacing -999 std with mean std")
    s_obs$std[which(s_obs$std == -999)] <- mean(subset(s_obs,std>0)$std)

    ## Remove all negative and NA values from both observations and incidence matrix
    message("Removing all negative values and NA")
    idx <- which(s_obs$z > 0 & !is.na(s_obs$z))
    C <- C[idx,]
    s_obs <- s_obs[idx,]
}

## Create observation precision matrix
Qobs <- Diagonal(x=1/s_obs$std^2)
Sobs <- solve(Qobs)
Z <- s_obs$z

## Create the "wrong" inventory from Europe
Emissions_all$mask <- 0
Emissions_all$mask[id + 14*63 - 15] <- 1  # Cut out England shape from France/Germany
Emissions_land$z2 <- Emissions_all$z[Emissions_all$mask==1] # Pretend this is the inventory
ggplot(Emissions_all) + geom_tile(aes(x=x,y=y,fill=mask)) +
  geom_path(data=Europe,aes(x=long,y=lat,group=group))
Emissions_map("z2")
if(true_inventory) {
  Zinv <- Emissions_land$z
} else {
  Zinv <- Emissions_land$z2
}

## Create Covariate matrix X
#X1 <- matrix(1,nrow(Dist_mat),1)
Highlands <- which(Emissions_land$y > 56.4)
X1 <- matrix(1,nrow(Dist_mat),2)
X1[Highlands,1] <- 0
X1[-Highlands,2] <- 0

p1 <- ncol(X1) # number of covariates
n1 <- nrow(X1) # number of flux field locations

if(save_images) {
  # plot simulation setup
  s_obs$obs_name <- as.character(s_obs$obs_name)
  s_obs <- s_obs %>%
    within(obs_name[grepl("MHD",obs_name)]<- "MHD") %>%
    within(obs_name[grepl("RGL",obs_name)]<- "RGL") %>%
    within(obs_name[grepl("TAC",obs_name)]<- "TAC") %>%
    within(obs_name[grepl("TTA",obs_name)]<- "TTA")
  g_stat_data <- LinePlotTheme() + geom_point(data=s_obs,
                                              aes(x=t,y=z)) +
    facet_grid(obs_name~.) +
    ylab("background-enhanced mol. fraction (ppb)\n") +  xlab("t (2 h steps)") +
    theme(axis.title.y=element_text(vjust=550))
  
  ggsave(filename = "./img/Fig3_SRR_Station_data_simmed.pdf",
         arrangeGrob(g_SRR,g_stat_data,ncol=2),
         width=17,height=8)
}
#-----------------------
# 4. MCMC algorithm
#-----------------------

## Flux field parameters
## Fix lambda based on model choice
if(grepl("Gaussian",model)) {
    message("Fixing lambda to 1")
    lambda_fix = 1
    eps <- 0.210
} else if(grepl("Lognormal",model)) {
    message("Fixing lambda to 0")
    lambda_fix = 0
    eps = 0.146
} else if(grepl("Box-Cox",model)) {
    lambda_fix = NA
    eps = 0.116
} else {
    stop("model needs to be Gaussian, Lognormal or Box-Cox")
}

## Uniformly distributed step-size for HMC
eps_gen <- function() runif(n=1,
                            min = eps,
                            max = eps+0.001)        # Step size generator for free lambda


## Set up slice sampler for the flux field parameters
if(!grepl("uncorrelated",model)) {
    nf <- 3 - !is.na(lambda_fix)
    theta_f_sampler <- slice::slice(nf)
    theta_f_samp <- matrix(0,3,N)
    theta_f_samp[,1] <- c(0.5,1,1)
    if(!is.na(lambda_fix)) theta_f_samp[3,] <- lambda_fix
} else {
    nf <- 1
    theta_f_sampler <- slice::slice(1)
    theta_f_samp <- matrix(1,1,N)
    if(!is.na(lambda_fix)) theta_f_samp[1,] <- lambda_fix
}

## HMC parameters
M <- diag(rep(1/100^2,n1))                      # scaling matrix
L <- 25L                                        # number of steps for each sample
Yf_samp <- matrix(0,nrow(M),N)                  # where to store the samples
Yf_samp[,1] <- curr_Yf_samp <- pmax(Emissions$z + 30*rnorm(n1),10)  # first 'sample'
Yf_accept <- rep(0,N)                           # initialise

## Mole-fraction field parameters
## Construct slice sampler
theta_m_sampler <- slice::slice(3)
theta_m_samp <- matrix(0,3,N)
theta_m_samp[,1] <- c(log(2000),0.7,log(2.5))
if(!real_data)
  {
     Q_zeta <- Q_zeta_true
} else  {
  Q_zeta <- (1/exp(theta_m_samp[1,1])) * Q_norm_zeta_fn(obs_locs,
                                                        t_axis,
                                                        theta_m_samp[2,1],
                                                        exp(theta_m_samp[3,1]))
}

if(do_inference) {
  print("Starting inference")

  ## Initialise parallel backend for distributing the MCMC chains (Linux only)
  if(n_parallel > 1) registerDoMC(n_parallel)

  # For each parallel MCMC chain
  All_Samps <- foreach(j = 1:n_parallel,.export='Q_norm_zeta_fn') %dorng% {

        for (i in 2:N) {  # for each sample

            
          ## Sample the flux field parameters
          if(grepl("uncorrelated",model)) {
                S_f_trans <- diag(n1)
                if(grepl("Box-Cox",model)) # if field is uncorrelated
                {
                    theta_f_samp[1,i] <- theta_f_sampler(theta_f_samp[1:nf,(i-1)],  # previous sample
                                                         log_f_theta_f,             # log-density
                                                         learn = (i <= adapt),      # learning?
                                                         Z=cbind(Zinv,Yf_samp[,i-1]), # data + inventory
                                                         D=Dist_mat,                # distance matrix
                                                         flux_cov_fn = flux_cov_fn, # flux cov. function
                                                         X=X1,                      # covariates 
                                                         lambda_fix=lambda_fix,     # lambda value
                                                         uncorrelated=TRUE)         # uncorrelated?
                }
                current_lambda <- theta_f_samp[1,i]
            } else { # if field is correlated
                theta_f_samp[1:nf,i] <- theta_f_sampler(theta_f_samp[1:nf,(i-1)],  # previous sample
                                                        log_f_theta_f,             # log density
                                                        learn = (i <= adapt),      # learning?
                                                        Z=cbind(Zinv,Yf_samps[,i-1]), # data + inventory
                                                        D=Dist_mat,                # distance matrix
                                                        flux_cov_fn = flux_cov_fn, # flux cov. function
                                                        X=X1,                      # covariates
                                                        lambda_fix=lambda_fix)     # lambda value
                current_lambda <- theta_f_samp[3,i]  # update lambda

                ## update flux correlation matrix. This is always the identity if the field is uncorrelated
                S_f_trans <- flux_cov_fn(Dist_mat,                     # distance matrix
                                         scale = theta_f_samp[1,i],    # scale parameter
                                         nu=theta_f_samp[2,i])         # smoothness parameter
            }

            
            ## Conditional distribution of the flux field
            Yf_cond <- Yf_BC_conditionals(Z = matrix(s_obs$z),      # training obs locs
                                          Zinv = matrix(Zinv),      # inventory
                                          C_m = C ,       # incidence matrix
                                          Qobs = Qobs,    # precision matrix
                                          B = B,          # SRR
                                          X = X1,         # covariates
                                          Q_zeta = Q_zeta,   # discrepancy cov matrix
                                          S_f_trans = S_f_trans/S_f_trans[1,1],  # covariance of log-flux field
                                          lambda=current_lambda,      # indices to consider (all)
                                          lambda_fix = lambda_fix)    # lambda value
                                          
            ## Sample the flux field
            Yf_sampler <- hmc_sampler(U = Yf_cond$logf,        # log density
                                      dUdq = Yf_cond$gr_logf,  # gradient of log density
                                      M = M,                # scaling matrix
                                      eps_gen = eps_gen,    # step-size generator
                                      L = L,                # number of steps
                                      lower = rep(0,ncol(B)))  # lower-bound (0)


            ## Generate flux field sample
            curr_Yf_samp <- Yf_sampler(q = curr_Yf_samp)
            Yf_samp[,i] <- curr_Yf_samp
            
            ## Every 10 samples, save the acceptance rate to file
            #print(paste0("Sample: ",i," Acceptance rate: ",(length(unique(Yf_samp[1,]))-1)/i))
            if(i%%10 == 0)
                cat(paste0("Sample: ",i," Acceptance rate: ",
                           (length(unique(Yf_samp[1,]))-1)/i,"\n"),
                file=paste0("./tempresults/",model,
                            "_Chain_",j,"TI",true_inventory,"RD",real_data),
                append=TRUE)

            if(!all(curr_Yf_samp == Yf_samp[,i-1])) Yf_accept[i] <- 1

            ## Every 10th sample with quick adaptation in beginning adjust HMC stepsize
            if(((i < 100 & i > 10) | ((i-1)%%10==0)) & i < adapt) { 
                
              # if poor acceptance reduce stepsize
                if(sum(Yf_accept[(i-10):i]) == 0) { 
                    eps = eps*0.9; print(paste0("New stepsize = ", eps))
                }
                
              # if too much acceptance  increase stepsize
                if(sum(Yf_accept[(i-10):i]) >= 9) {
                    if(grepl("uncorrelated",model))
                    {           
                      eps_max <- 1
                    } else 
                      { eps_max <- 1}
                    eps =min(eps*1.1,eps_max); print(paste0("New stepsize = ", eps))
                }
                
              # Update stepsize generator
                eps_gen <- function() runif(n=1,
                                            min = eps,         # 0.07 with original stds
                                            max = eps + 0.001)        # Step size generator for free lambda
            }


            ## Sample mole-fraction field parameters
            theta_m_samp[,i] <- theta_m_sampler(theta_m_samp[,(i-1)],   # previous samples
                                                log_f_theta_m,          # log density
                                                learn = (i <= adapt),   # learning?
                                                Yf=curr_Yf_samp,        # flux field
                                                s = obs_locs,           # observation locations
                                                t = t_axis,             # time axis
                                                C_m = C,                # incidence matrix
                                                Qobs = Qobs,            # obs. precision matrix
                                                Z = Z,                  # data
                                                B = B,                  # SRR
                                                Q_zeta_fn = Q_norm_zeta_fn)  # discr. precision matrix

            ## Update discrepancy precision matrix
            Q_zeta <- (1/exp(theta_m_samp[1,i])) * Q_norm_zeta_fn(obs_locs,           # obs. locations
                                                                  t_axis,             # time axis
                                                                  theta_m_samp[2,i],  # new samples
                                                                  exp(theta_m_samp[3,i]))

        }
    
        ## MCMC chain complete.. now dither the chain
        sub_samp <- seq(burnin,N,by=dither) 
        
        ## Save data
        if(!real_data) {
          return(list(model = model,
             i = i,
             burnin = burnin,
             n_parallel = n_parallel,
             dither = dither,
             theta_t_true = theta_t_true,
             theta_s_true = theta_s_true,
             sigma_zeta_true = sigma_zeta_true,
             theta_m_samp = theta_m_samp[,sub_samp],
             Yf_samp = Yf_samp[,sub_samp],
             theta_f_samp = theta_f_samp[,sub_samp],
             true_inventory = true_inventory))
        } else {
          if(real_data)
            return(list(model = model,
                 i = i,
                 burnin = burnin,
                 n_parallel = n_parallel,
                 dither = dither,
                 theta_m_samp = theta_m_samp[,sub_samp],
                 Yf_samp = Yf_samp[,sub_samp],
                 theta_f_samp = theta_f_samp[,sub_samp],
                 true_inventory = true_inventory))
        }

    }
    if(!real_data)
      save(All_Samps,file = paste0("./tempresults/Results_",model,"TI",true_inventory,".rda"))
    else
      save(All_Samps,file = paste0("./tempresults/Results_",model,"TI",true_inventory,"_RD.rda"))
}

#-------------------------------------------------------------------
# 5. Analyse results
#------------------------------------------------------------------

## See if there are any problematic chains (Only to be run in DEV mode)
if(0) {
  fnames <- dir("../atminv/inst/extdata/spastaMCMC/")[grepl("Chain",dir("../atminv/inst/extdata/spastaMCMC/"))]
  for (i in 1:length(fnames)) {
    #file_in <- file(fnames[i],"r")
    file_in <- file(paste0("../atminv/inst/extdata/spastaMCMC/",fnames[i]),"r")
    x <- readLines(file_in)
    temp <- gregexpr("[0-9.]+", x[length(x)])
    acc_ratio <- min(as.numeric(unique(unlist(regmatches(x[length(x)], temp)))))
    close(file_in)
    if(acc_ratio < 0.25 | acc_ratio > 0.8) print(paste0(fnames[i]," had an acceptance ratio of ",acc_ratio))
  }
}

models <- c("Box-Cox","Lognormal","Gaussian","Box-Cox.uncorrelated",
            "Lognormal.uncorrelated","Gaussian.uncorrelated")
model_names <- paste0("Model ",1:6)
e_bi <- 1 #extra burnin

## Combine the parallel chains into one long one
combine_chains <- function(l) {
  new_l <- list(Yf_samp = NULL, theta_f_samp = NULL, theta_m_samp = NULL)
  #print("Discarding Chain 9 from each model")    
  for (i in setdiff(1:length(l),NULL)) {
    new_l$Yf_samp <- cbind(new_l$Yf_samp, l[[i]]$Yf_samp[,-(1:e_bi)])
    if(is.vector( l[[i]]$theta_f_samp)) {
      new_l$theta_f_samp <- cbind(new_l$theta_f_samp, t(l[[i]]$theta_f_samp[-(1:e_bi)]))
    } else {
      new_l$theta_f_samp <- cbind(new_l$theta_f_samp, l[[i]]$theta_f_samp[,-(1:e_bi)])
    }

    new_l$theta_m_samp <- cbind(new_l$theta_m_samp, l[[i]]$theta_m_samp[,-(1:e_bi)])
  }
  list(new_l)
}

###########################
### Parameters posterior
###########################

# Load model 1
#load(paste0("Results_",models[1],"TI1.rda"))
load(system.file("extdata",paste0("spastaMCMC/Results_",models[1],"TI1.rda"), package = "atminv"))

# Posterior median of emissions for Model 1
Emissions_land$post <- apply(All_Samps[[1]]$Yf_samp,1,median)
if(save_images)
  ggsave(Emissions_map("post"),filename = "./img/median_flux_post.png")

# Number of samples
Nsamp <- (dim(All_Samps[[1]]$theta_m_samp)[2]-e_bi) * 10

# Initialise data frame
df_par <- data.frame(n = 1:Nsamp)

# For each model add six columns corresponding the the sampled parameters (on original scale)
for(i in 1:6) {
    #load(paste0("Results_",models[i],"TI1.rda"))
    load(system.file("extdata",paste0("spastaMCMC/Results_",models[i],"TI1.rda"), package = "atminv"))
    All_Samps <- combine_chains(All_Samps)
    df_par[paste0(model_names[i],"_thetam1")] <- 1/(exp(All_Samps[[1]]$theta_m_samp[1,]))
    df_par[paste0(model_names[i],"_thetam2")] <- All_Samps[[1]]$theta_m_samp[2,]
    df_par[paste0(model_names[i],"_thetam3")] <- exp(All_Samps[[1]]$theta_m_samp[3,])
    if(grepl("uncorrelated",models[i])) {
        df_par[paste0(model_names[i],"_thetaf3")] <- All_Samps[[1]]$theta_f_samp[1,]
    } else {
        df_par[paste0(model_names[i],"_thetaf1")] <- All_Samps[[1]]$theta_f_samp[1,]
        df_par[paste0(model_names[i],"_thetaf2")] <- All_Samps[[1]]$theta_f_samp[2,]
        df_par[paste0(model_names[i],"_thetaf3")] <- All_Samps[[1]]$theta_f_samp[3,]
    }
    if(grepl("Gaussian",models[i]) | grepl("Lognormal",models[i])) {
        df_par[paste0(model_names[i],"_thetaf3")] <- NULL
    }

}

## Now put the data frame into long form
df_par2 <- gather(df_par,par,sample,-n) %>%
    separate(par,c("Model","par"),sep="_") %>%
    mutate(process=ifelse(grepl("m",par),"Mole-fraction","Flux"))

## Visalise some of the chains to make sure they make sense
ggplot(subset(df_par2, Model=="Model 5")) + geom_point(aes(x=n,y=sample)) + facet_grid(par~.,scales="free")

## Labeller function for labelling the facets
my.label_bquote <- function(variable, value) {
    value <- as.character(value)
    if(variable == 'par')
        lapply(value,function(x)
            if(x == "thetam1") {
                bquote(tau[2])
            } else if (x == "thetam2") {
                bquote(a)
            } else if (x == "thetam3") {
                bquote(d)
            } else if (x == "thetaf1") {
                bquote(theta[11])
            } else if (x == "thetaf2") {
                bquote(theta[12])
            } else if (x == "thetaf3") {
                bquote(lambda)
            })
}

## Plot the posterior distributions of the flux parameters
g1 <- LinePlotTheme() + 
    geom_density(data=subset(df_par2,process=="Flux"),aes(x=sample,linetype=Model),fill="black",alpha=0.1,adjust=2) +
    geom_density(data=subset(df_par2,process=="Flux"),aes(x=sample,linetype=Model),adjust=2) +
  facet_grid(.~par,scales="free",labeller =my.label_bquote) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  scale_linetype_manual(values = c(1:4)) + xlab("")

## Plot the posterior distributions of the mole-fraction parameters
g2 <- LinePlotTheme() + geom_boxplot(data=subset(df_par2,process=="Mole-fraction"),aes(x=Model,y=sample)) +
  facet_grid(par~.,scales="free",labeller =my.label_bquote) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("")

## Save the image
g_all <- arrangeGrob(g1,g2,ncol=2)
if(save_images)
  ggsave(filename = "./img/Fig4_theta_post.pdf",
         arrangeGrob(g1,g2,ncol=2),
         width = 17,height=8)


## Find the posterior mean of lambda in Model 1
#load(paste0("Results_",models[1],"TI1.rda"))
load(system.file("extdata",paste0("spastaMCMC/Results_",models[1],"TI1.rda"), package = "atminv"))
All_Samps <- combine_chains(All_Samps)
mean(All_Samps[[1]]$theta_f_samp[3,])

###########################
### Flux field posterior
###########################

## Initialise data frame
df_flux <- data.frame(n = 1:Nsamp)

## Initialise data frame for MAE, CRPS quantiles
MAE_flux <- crps_flux <- data.frame(Model = 1:6,LQ = 0, Median = 0, UQ =0)

## Initialise data frame for errors 
error_flux <- data.frame(Model = 1:6, MAPE = 0, RMSPE = 0,MCRPS=0)

## for both the true inventory and wrong inventory models
for(TI in c(1,0))
  for(i in 1:6) {
      #if(file.exists(paste0("Results_",models[i],"TI",TI,".rda"))) {
      if(file.exists(system.file("extdata",paste0("spastaMCMC/Results_",models[i],"TI",TI,".rda"), 
                                 package = "atminv"))) {
          #load(paste0("Results_",models[i],"TI",TI,".rda"))
          
          # load file 
          load(system.file("extdata",paste0("spastaMCMC/Results_",models[i],"TI",TI,".rda"), package = "atminv"))
          
          # combine chains
          All_Samps <- combine_chains(All_Samps)
          
          # extract flux samples
          Yf_samp <- All_Samps[[1]]$Yf_samp
          
          # model number (accounting for false inventory model, Model 7)
          ii <- i + (!TI)*6
          
          # initialise crps, cdf and error fields
          Emissions_land[paste0("crps",ii)] <- 0
          Emissions_land[paste0("cdf", ii)] <- 0
          Emissions_land[paste0("error",ii)] <- apply(Yf_samp,1,mean) - Emissions_land$z
          
          # for each flux field location
          for (j in 1:nrow(Emissions_land)){
            
              # find the ECDF at this location
              Emissions_land[paste0("cdf",ii)][j,] <- ecdf(Yf_samp[j,1000:(i-1)])(Emissions_land$z[j])
              # compute the CRPS
              Emissions_land[paste0("crps",ii)][j,] <- sum((ecdf(Yf_samp[j,])(0:5000) -
                                                               as.numeric((0:5000)>Emissions_land$z[j]))^2) * 1 #unit grid width
          }
          ## add CRPS to CRPS data frame
          crps_flux[ii,] <- c(ii, quantile(unlist(Emissions_land[paste0("crps",i)]),
                                         c(0.25,0.5,0.75)))
          
          ## add MAE to MAE data frame
          MAE_flux[ii*6,] <- c(i + (!TI)*6, quantile(mean(abs(unlist(Emissions_land[paste0("error",i)]))),
                                                     c(0.25,0.5,0.75)))
          
          ## add errors to error data frame
          error_flux[ii,] <- c(ii, mean(abs(unlist(Emissions_land[paste0("error",ii)]))),
                              sqrt(mean(abs(unlist(Emissions_land[paste0("error",ii)]))^2)),
                              mean(unlist(Emissions_land[paste0("crps",ii)])))
          }
  }

## print LaTeX table of flux diagnostics
print(xtable::xtable(error_flux,digits=c(NA,1,1,1,1)),
      include.rownames=FALSE)


## Compare the Gaussian vs Box-Cox vs Lognormal
Emissions_land$wtf <- (abs(Emissions_land$error1)  < abs(Emissions_land$error3))
g1 <- Emissions_map("wtf") + scale_fill_distiller(breaks=c(0,1), guide = guide_legend(title="")) +
    ggtitle("Model 1 vs. Model 3")
Emissions_land$wtf <- (abs(Emissions_land$error2)  < abs(Emissions_land$error1))
g2 <- Emissions_map("wtf")+ scale_fill_distiller(breaks=c(0,1), guide = guide_legend(title=""))  +
    ggtitle("Model 2 vs. Model 1")
if(save_images) {
  ggsave(g1,filename = "./img/error_map1.png",width = 8,height=10)
  ggsave(g2,filename = "./img/error_map2.png",width = 8,height=10)

  png(filename="./img/qqplots.png",width = 4000,height=1600,res = 300)
  par(mfrow=c(1,3))
  qqnorm(Emissions_land$z,ylim=c(-600,2500),main="Identity"); qqline(Emissions$z)
  qqnorm(Emissions_land$z^0.5,main="Square root"); qqline(Emissions$z^0.5)
  qqnorm(Emissions_land$z^0.2,main="Fifth root"); qqline(Emissions$z^0.2)
  dev.off()
}

## some plots of errors of Model A vs errors of Model B (not included in paper)
par(mfrow=c(1,2))
plot(Emissions_land$error1 - Emissions_land$error2,
     ylim = c(-300,200),
     ylab = "Difference in prediction errors (g/s)",
     xlab = "Location index");
lines(c(0,124),c(0,0),lty=2)
plot(Emissions_land$error1 - Emissions_land$error3,
     ylim = c(-300,200),
     ylab = "Difference in prediction errors (g/s)",
     xlab = "Location index");  lines(c(0,124),c(0,0),lty=2)
dev.off()

## VIOLIN PLOT
## Load Box-Cox results
#load(paste0("Results_",models[1],"TI1.rda"))
load(system.file("extdata",paste0("spastaMCMC/Results_",models[1],"TI1.rda"), package = "atminv"))
All_Samps <- combine_chains(All_Samps)

## Extract flux field parameters
Yf_samp <- All_Samps[[1]]$Yf_samp

## Arrange samples into a long data frame
Q <- as.data.frame(t(Yf_samp)) %>%
    tidyr::gather(t,z) %>%
    separate(t, into=c("V","s"),sep="V") %>%
    select(-V) %>%
    mutate(x = Emissions$x[as.numeric(s)],
           y = Emissions$y[as.numeric(s)]) %>%
    select(-s)

## Fuse Q above with the Emissions data frame (after replacing "z" with "NAEI")
Q2 <- left_join(Q,select(mutate(Emissions,NAEI=z),
                         x,y,NAEI,region,subregion)) %>%
        group_by(x,y) %>%
        mutate(med_z = median(z)) %>%
        as.data.frame()

## facet labeller
dec_labeller <- function(variable, value) {
    lapply(as.numeric(value),function(x) round(x,1))
}

## Plot the violin plot
if(save_images) {
    violin_plot <- ggplot(mutate(Q2,y=-y)) +
        geom_rect(data = mutate(Emissions_land,y=-y),
                  aes(xmin = 0.5,xmax = 1.5,
                  ymin = 0,ymax = 3000,fill = z)) +
        geom_violin(data=mutate(Q2,y=-y),aes(x=1,y=z,fill=med_z)) +
        geom_point(aes(x=1,y=NAEI),col="black",shape=4) +
        facet_grid(y~x,labeller = dec_labeller) +
        scale_fill_distiller(palette="Spectral",
                             guide = guide_legend(title="flux (g/s)"))  +
        ylim(0,3000) + scale_x_continuous(breaks=NULL) + xlab("") + ylab("flux g/s") +
        theme(panel.background = element_rect(fill='white', colour='grey'),
              panel.grid=element_blank(),axis.ticks=element_blank(),
              panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_blank())


    ggsave(violin_plot,filename = "./img/Fig5_violin_plot.png",height=12,width=12)
}

################
## ADMIN FLUX
################

## Admin boundaries obtained from Tthe GADM database (www.gadm.org), version 2.8, November 2015. They can be used for non-commercial purposes only.  It is not allowed to redistribute these data, or use them for commercial purposes, without prior consent. See the website for more information.

## id: 0: England
##     1: Northern Ireland
##     2: Scotland
##     3: Wales
data(UK_regions)

Emissions_land$admin <- "Ireland"
for(i in c("0.1","1.1","2.1","3.1")) {
    this_admin <- filter(UK_regions,group==i)
    idx <- which(   sp::point.in.polygon(Emissions_land$x,Emissions_land$y,
                                         this_admin$long,this_admin$lat)  |
                        sp::point.in.polygon(Emissions_land$x + dx_new/3,Emissions_land$y + dy_new/3,
                                             this_admin$long,this_admin$lat)  |
                        sp::point.in.polygon(Emissions_land$x - dx_new/3,Emissions_land$y + dy_new/3,
                                             this_admin$long,this_admin$lat) |
                        sp::point.in.polygon(Emissions_land$x + dx_new/3,Emissions_land$y - dy_new/3,
                                             this_admin$long,this_admin$lat) |
                        sp::point.in.polygon(Emissions_land$x - dx_new/3,Emissions_land$y - dy_new/3,
                                             this_admin$long,this_admin$lat)  == 1)
    Emissions_land$admin[idx] <- switch(i, 
                                        "0.1" = "England",
                                        "1.1" = "Northern Ireland",
                                        "2.1" = "Scotland",
                                        "3.1" = "Wales")
}
ggplot(Emissions_land) + geom_tile(aes(x=x,y=y,fill=admin))


## Initialise data frame
admin_flux <- data.frame(Model = NULL,samp = NULL)

load(system.file("extdata",paste0("spastaMCMC/Results_Box-CoxTI1.rda"), package = "atminv"))
All_Samps <- combine_chains(All_Samps)

## For each model compute the total flux in each sample
for(i in c("England","Northern Ireland","Scotland","Wales","Ireland")) {
    idx <- which(Emissions_land$admin == i)
    admin_flux <- rbind(admin_flux,
                        data.frame(admin = i, 
                                   samp = apply(All_Samps[[1]]$Yf_samp[idx,],2,sum)*3600*24*365/10^12))
}

## Table of aggregates
Emissions_admin <- group_by(Emissions_land,admin) %>% 
    summarise(flux = sum(z)*3600*24*365/10^12)

g_admin <- LinePlotTheme() + geom_boxplot(data=admin_flux,aes(x=admin,y=samp)) +
    geom_point(data=Emissions_admin,aes(x=admin,y=flux),shape=2,size=3) +
    xlab("\n Administrative region (land only)") + 
    ylab("Total flux (Tg/yr) \n")

###############
## Total flux
###############

## Initialise data frame
tot_flux <- data.frame(Model = NULL,samp = NULL)

## For each model compute the total flux in each sample
for(j in 1:6) {
    #load(paste0("Results_",models[j],"TI1.rda"))
    load(system.file("extdata",paste0("spastaMCMC/Results_",models[j],"TI1.rda"), package = "atminv"))
    All_Samps <- combine_chains(All_Samps)
    tot_flux <- rbind(tot_flux,
                      data.frame(Model = paste0("Model ",j), 
                                 samp = apply(All_Samps[[1]]$Yf_samp,2,sum)*3600*24*365/10^12))
}

## Add on Model 1* (Model 7 in this program)
#load(paste0("Results_",models[1],"TI0.rda"))
load(system.file("extdata",paste0("spastaMCMC/Results_",models[1],"TI0.rda"), package = "atminv"))
All_Samps <- combine_chains(All_Samps)
tot_flux <- rbind(tot_flux,
                  data.frame(Model = "Model 1*", 
                             samp = apply(All_Samps[[1]]$Yf_samp,2,sum)*3600*24*365/10^12))

## Table of aggregates
group_by(tot_flux,Model) %>%
    summarise(mean = mean(samp),
              lq = quantile(samp,0.25),
              uq = quantile(samp,0.75),
              s2 = mean(samp) - 2*sd(samp))

## Total flux in UK + Ireland in the correct inventory
tot_UK_Ir <- sum(Emissions_land$z*3600*24*365/10^12)

## Total flux in the wrong (Europe) inventory
tot_Europe <- sum(Emissions_land$z2*3600*24*365/10^12)

## Plot the bar chart of aggregates
g_total <- LinePlotTheme() + geom_boxplot(data=tot_flux,aes(x=Model,y=samp)) +
    geom_point(data=data.frame(y=rep(tot_UK_Ir,7),x=1:7),aes(x,y),shape=2,size=3) +
    ylab("Total flux (Tg/yr)\n") + xlab("\n Model")
if(save_images)
  ggsave(arrangeGrob(g_admin,g_total,ncol=2),file="./img/Fig6_tot_flux.pdf",height=9,width=18)

### Mole-fraction
LinePlotTheme() + geom_point(data=s_obs,aes(t,z)) + facet_grid(obs_name~.)

#############################################
### Posterior distribution of mole fraction
#############################################

## Initialise list
All_Samps_list <- list()

# Put all the samples (for each model) into one big list
for(j in 1:6){
    #load(paste0("Results_",models[j],"TI1.rda"))
    load(system.file("extdata",paste0("spastaMCMC/Results_",models[j],"TI1.rda"), package = "atminv"))
    All_Samps <- combine_chains(All_Samps)
    All_Samps_list[[j]] <- All_Samps
}

## Add on Model 1*
#load(paste0("Results_Box-CoxTI0.rda"))
load(system.file("extdata",paste0("spastaMCMC/Results_Box-CoxTI0.rda"), package = "atminv"))
All_Samps_list[[7]] <- combine_chains(All_Samps)

## To compute this set to if(1). 
if(0) {
  
  # parallelise using 4 cores
  registerDoMC(4)
  print("Sampling mole fractions, this will take a while...")
  
  ## for each model
  All_MFs <- foreach(All_Samps = All_Samps_list,.export='Q_norm_zeta_fn') %dorng% {
      
      # extract the relevant quantities
      Yf_samp <- All_Samps[[1]]$Yf_samp
      theta_m_samp <- All_Samps[[1]]$theta_m_samp
      theta_f_samp <- All_Samps[[1]]$theta_f_samp

      # use only 1000 samples from what we have (should be sufficient)
      use_samps <- round(seq(1,Nsamp,length=1000))
      Ns <- length(use_samps)
      
      # initialise mole-fraction sample matrix
      mf_samp <- matrix(0,nrow(B),Ns)         
      cnt <- 1
      
      # for each sample
      for (i in use_samps){
          # construct discrepancy precision matrix
          Q_zeta <- (1/exp(theta_m_samp[1,i])) * Q_norm_zeta_fn(obs_locs,
                                                                t_axis,
                                                                theta_m_samp[2,i],
                                                                exp(theta_m_samp[3,i]))
          
          # posterior precision matrix
          Q_post <- t(C) %*% Qobs %*% C + Q_zeta
          X <- cholPermute(Q_post)
          
          # standard Gaussian sampling
          SQB <- linalg:::cholsolve(Q_post,Q_zeta %*% B,cholQp = X$Qpermchol, P = X$P)
          StCQoz <- linalg:::cholsolve(Q_post,(t(C) %*% Qobs %*% Z),
                                       cholQp = X$Qpermchol, P = X$P)
          mu_post <- as.vector(SQB %*% Yf_samp[,i] + StCQoz)
          this_GMRF <- GMRF(mu=matrix(mu_post),Q = Q_post)
          mf_samp[,cnt] <- sample_GMRF(this_GMRF,L=X$Qpermchol,P=X$P)
          cnt <- cnt + 1
      }
      
      # save into list
      list(mf_samp = mf_samp,
           mf_mean = apply(mf_samp,1,mean),
           mf_lq = apply(mf_samp,1,quantile,0.25),
           mf_uq = apply(mf_samp,1,quantile,0.75))
  }
  
  # plot mole fractions and observations for Model 1 (sanity check)
  s_obs_orig$mf_mean <- All_MFs[[1]]$mf_mean
  LinePlotTheme() + geom_point(data=s_obs_orig,aes(t,z)) +
    geom_line(data=s_obs_orig,aes(t,mf_mean),col="red") +
    facet_grid(obs_name~.)
  
  # create the diagnostic mole-fraction data frame and initialise
  MF_df <- data.frame(MAPE=0,MSPE=0,CRPS=0)
  crps_axis <- seq(-100,200,by=0.1)
  
  # for each model
  for (j in 1:7){
    
    # extract the mean
    All_MFs[[j]]$mf_mean
    
    # initialise count
    this_crps <- cnt <- 0
    
    # for every unobserved location
    for (i in setdiff(1:nrow(s_obs_orig),idx)){
      
      # find the CRPS at this location
      cnt <- cnt + 1
      this_crps[cnt] <- sum((ecdf(All_MFs[[j]]$mf_samp[i,])(crps_axis) -
                               as.numeric((crps_axis)>s_obs_orig$Ymf[i]))^2) * 0.1 #unit grid width
    }
    
    # compute the diagnostics for all locations
    MF_df[j,] <- c(mean(abs(All_MFs[[j]]$mf_mean[-idx] - s_obs_orig$Ymf[-idx])),
                   sqrt(mean((All_MFs[[j]]$mf_mean[-idx] - s_obs_orig$Ymf[-idx])^2)),
                   mean(this_crps))
  }

  # save results for quick load
  save(MF_df,file = "./tempresults/MF_df.rda")
} else {
  #load("MF_df.rda")
  load(system.file("extdata",paste0("spastaMCMC/MF_df.rda"), package = "atminv"))
}

## print both the flux and mole-fraction diagnostics side-by-side
print(xtable::xtable(cbind(error_flux,MF_df),
                     align = c("c","|c", "|c", "c", "c", "|c","c","c|"),
                     digits=c(NA,1,1,1,1,1,1,1)),
      include.rownames=FALSE)

