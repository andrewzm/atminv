#---------------------
# 1. Model Setup
#---------------------

## Important: set wkdir to base path of package atminv

library(dplyr)
library(tidyr)
library(Matrix)
library(devtools)
library(ggplot2)
library(grid)
library(gridExtra)
library(sp)      # Needed for variogram modelling
library(gstat)   # Needed for variogram modelling
library(atminv)  # https://github.com/andrewzm/atminv
library(hmc)     # https://github.com/andrewzm/hmc
library(fitdistrplus)

rm(list=ls())
cache_results <- 1
load_results  <- 0
save_images   <- 1

md5_wrapper <- md5_cache("~/cache/Assim") # set up cache
set.seed(1)

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

### Load maps

# Load world map and subset to interesting countries
Europe <- map_data("world") %>%  
    subset(region %in% c("UK","France","Ireland")) %>%
    mutate(id = group, x=long,y=lat)

# Attribute each grid-cess to a country and sea/land indicator
yx <- mutate(yx, in_land= attribute_polygon(yx,Europe),
             id = in_land) %>%
    left_join(unique(select(Europe,region,id,subregion)))

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

# Plot Emissions map
#ggplot() + geom_tile(data=Emissions,aes(x,y,fill=pmin(z,2000))) + coord_fixed()

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

# Plot the emissions in the UK/Ireland regions
g <- LinePlotTheme() + 
    geom_tile(data=subset(Emissions,x > -13 & x < 3 & y > 48 & y < 61),
                                 aes(x,y,fill=pmax(pmin(z,2000),-2000)))  +
    scale_fill_distiller(palette="Spectral", 
                         trans="reverse", 
                         guide = guide_legend(title="flux (g/s)")) +
    geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
    coord_map(xlim=c(-12,2),ylim=c(49,60)) + 
    xlab("lon (degrees)") + ylab("lat (degrees)") +
    theme(axis.title.y = element_text(vjust = 1))
if(save_images)
    ggsave(filename = "./img/NAEI.png",plot = g,width=6.4,height=6.7)

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

### Deprecated -- Find the important locations to keep (areas of influence)
### -------------------------------------------------------
if(0) {
    # Stations <- md5_wrapper(lapply,Stations,read_process_models,prune=F)
    # max_detect <- lapply(Stations, function(des) {
    #   des$tot_influence <- group_by(des$Model_B,x,y) %>%
    #     summarise(max_detect = max(z))
    #   des$tot_influence$tot <- des$tot_influence$max_detect * Emissions$z
    # })
    # Emissions$max_detect <- apply(Reduce("cbind",max_detect),1,max)
    # yx$max_detect <- Emissions$max_detect
    # Stations <- lapply(Stations,function(des) { des$Model_B$max_detect = Emissions$max_detect; des })
    #ggplot() + geom_tile(data=subset(Emissions,!(in_land==0)),aes(x,y,fill=max_detect > max(max_detect)/100)) + coord_fixed()
}

### Deprecated -- Emulation
### -------------------------------------------------------
if(0) {
    ### Emulation example
    subregion <- function(des) {
        des$Model_C$s0x <- des$x
        des$Model_C$s0y <- des$y
        subset(des$Model_C,t==70 & x < 10 & x > -10 & y < 58 & y > 45)
    }
    C_subs <- lapply(Stations,subregion)
    
    ### Covariance approach
    C_subs_ccat <- Reduce("rbind",C_subs) %>% subset(z > 0)
    Cov_fun <- function(x1,x2) {
        sc = 1
        x1 <- mutate(x1,s0x = s0x/sc, s0y =s0y/sc)
        x2 <- mutate(x2,s0x = s0x/sc, s0y =s0y/sc)
        
        R1 <- rdist(select(x1,s0x,s0y),select(x2,s0x,s0y))
        R2 <- rdist(select(x1,h1,h2),select(x2,h1,h2))
        R <- rdist(select(x1,s0x,s0y,h1,h2),select(x2,s0y,h1,h2))
        R <- 0.05*R1 + R2
        C <- MVST::Matern(R,kappa=1,nu=1/2,var=0.001)
        C
    }
    C_subs_ccat <- C_subs_ccat %>% mutate(h1 = x - s0x, h2 = y - s0y)
    Cov <- Cov_fun(x1 = select(C_subs_ccat,s0x,s0y,h1,h2),
                   x2 = select(C_subs_ccat,s0x,s0y,h1,h2)) #+ 0.01 *diag(rep(1,nrow(C_subs_ccat)))
    C_predict <- C_subs[[1]] %>% mutate(s0x = -2, s0y = 54.00,
                                        h1 = x - s0x, h2 = y - s0y)
    Cov2 <- Cov_fun(x1 = C_predict,
                    x2 = C_subs_ccat)
    C_predict$zhat <- Cov2 %*% chol2inv(chol(Cov)) %*% C_subs_ccat$z
    ggplot(C_predict) + geom_tile(aes(x,y,fill=zhat)) + geom_point(aes(x=s0x,y=s0y),colour="white")
    ggplot(C_subs[[2]]) + geom_tile(aes(x,y,fill=z)) + geom_point(aes(x=s0x,y=s0y),colour="white")
    C_predict$z <- as.numeric(C_predict$zhat)
    XXX <- rbind(select(C_predict,x,y,z),select(C_subs_ccat,x,y,z)) %>%
        group_by(x,y) %>%
        summarise(ztot = sum(z))
    ggplot(XXX) + geom_tile(aes(x,y,fill=pmax(ztot,0))) + scale_fill_gradient(low="white",high="red")
    
    ### IDW approach
    subregion <- function(des) {
        s0x <- des$Model_C$x[which.min(abs(des$Model_C$x - des$x))]
        s0y <- des$Model_C$y[which.min(abs(des$Model_C$y - des$y))]
        des$Model_C %>%
            subset(t==300) %>%
            mutate(s0x = s0x,
                   s0y = s0y,
                   h1 = x - s0x,
                   h1rnd = round(h1,1),
                   h2 = y - s0y,
                   h2rnd = round(h2,1),
                   r = sqrt(h1^2 + h2^2),
                   theta = atan2(h2,h1))
        
    }
    C_subs <- lapply(Stations,subregion)
    for(i in 1:length(C_subs)) {C_subs[[i]]$number = i; C_subs[[i]]$idx <- 1:nrow(C_subs[[i]])}
    C_subs_ccat <- Reduce("rbind",C_subs)
    
    # Prediction data frame
    spred <- C_subs[[1]] %>% mutate(s0xp = -2, s0yp = 54,
                                    xp = x, yp = y) %>%
        select(s0xp,s0yp,xp,yp) %>%
        mutate(h1p = xp -s0xp, h2p = yp - s0yp,
               rp = sqrt(h1p^2 + h2p^2), thetap = atan2(h2p,h1p))
    
    # Group long data frame by stations
    Station_groups <- C_subs_ccat %>%
        group_by(s0x,s0y)
    
    # for each prediction location, find the "closest" index for each station
    for(i in 1:nrow(spred)) {
        cl_idx <- summarise(Station_groups,
                            idx = idx[which.min(Mod(r*exp(1i*theta) - spred$rp[i]*exp(1i*spred$thetap[i])))],
                            number = number[1])
        for(j in 1 : nrow(cl_idx)) spred[i,paste0("Station.",cl_idx$number[j])] <- cl_idx$idx[j]
        print(i)
    }
    
    idw <- function(s0xp,s0yp,s0x,s0y,z) {
        D <- rdist(cbind(s0x,s0y),cbind(s0xp,s0yp))
        lambda <- 1/(D)^2
        if(all(lambda==0)) browser()
        sum(lambda * z)/sum(lambda)
    }
    
    spred2 <- gather(spred,Station,idx,Station.1,Station.2,Station.3,Station.4) %>%
        separate(Station, into = c("Station", "number"), sep = "\\.") %>%
        select(-Station)
    spred2$number <- as.integer(spred2$number)
    
    ## Only this needs to be run for each time point
    spred3 <- left_join(spred2,select(C_subs_ccat,number,idx,z,s0x,s0y)) %>%
        group_by(xp,yp) %>%
        summarise(z = idw(s0xp[1],s0yp[1],s0x,s0y,z))
    spred3 <- mutate(spred, x = xp, y = yp)
    XXX <- rbind(select(spred3,x,y,z),select(C_subs_ccat,x,y,z)) %>%
        group_by(x,y) %>%
        summarise(ztot = sum(z))
    ggplot(subset(XXX,ztot > 0.001)) + geom_tile(aes(x,y,fill=ztot)) +
        scale_fill_gradient(low="white",high="red")
    
    
    
    idw_all <- function(X,Xp) {
        XX <- X %>%
            group_by(s0x,s0y) %>%
            summarise(z = z[which.min(Mod(r*exp(1i*theta) - Xp$rp*exp(1i*Xp$thetap)))])
        XXp <- select(Xp, s0xp,s0yp)
        D <- rdist(select(XX,s0x,s0y),XXp)
        lambda <- 1/(D)^2
        sum(lambda * XX$z)/sum(lambda)
    }
    library(doMC)
    registerDoMC(6)
    spred$z <- foreach(i = 1:nrow(spred),.combine = c) %dopar% {idw_all(C_subs_ccat,spred[i,])}
    spred <- mutate(spred, x = xp, y = yp)
    XXX <- rbind(select(spred,x,y,z),select(C_subs_ccat,x,y,z)) %>%
        group_by(x,y) %>%
        summarise(ztot = sum(z))
    ggplot(subset(XXX,ztot > 0.001)) + geom_tile(aes(x,y,fill=ztot)) +
        scale_fill_gradient(low="white",high="red")
    #+
    #  geom_point(data=subset(spred,z > 0.001),aes(x,y,colour=z),size=3) +
    #  scale_colour_gradient(low="white",high="blue")
    
}

### Preprocess
### -------------------------------------------------------

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

## Plot the SRR at this time point
g <- LinePlotTheme() + geom_tile(data=subset(Emissions,x > -13 & x < 3 & y > 48 & y < 61),
                                 aes(x,y,fill=B))  +
    scale_fill_distiller(palette="BuPu",trans="sqrt",guide = guide_legend(title="SRR (s/ng)")) +
    geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
    coord_map(xlim=c(-12,2),ylim=c(49,60)) + xlab("lon (degrees)") + ylab("lat (degrees)") +
    geom_point(data=Station_locs,aes(x=unlist(x),y=unlist(y)),col="red",size=4) +
    geom_text(data=Station_locs,aes(x=unlist(x),y=unlist(y)+0.3,label=unlist(name)),col="red") +
    theme(axis.title.y = element_text(vjust = 1))
if(save_images)
   ggsave(filename = "./img/TTA_01_01.png",plot = g,width=6.4,height=6.7)

## Plot the SRR at this time point
g <- LinePlotTheme() + geom_point(data=subset(Stations$TAC$Obs_obj@df,t < 372),aes(x=t,y=z)) +
    geom_line(data=data.frame(t=1:372,z=(Stations$TAC$B[1:372,] %*% Emissions$z)[,1]),aes(x=t,y=z),col="red") +
    ylab("mol. fraction (ppb)") +
    xlab("t (2 h steps)")
if(save_images)
    ggsave(filename = "./img/TAC_01.png",plot = g,width=10,height=3.7)

## Plot the inventory-predicted and the observations at TAC. Note that this still
## include all the background SRR and background (Europeat) emissions. Nothing is
## truncated yet
g <- LinePlotTheme() + geom_point(data=subset(Stations$TAC$Obs_obj@df,t < 372),
                                  aes(x=t,y=z)) +
    geom_line(data=data.frame(t=1:372,z=(Stations$TAC$B[1:372,] %*% Emissions$z)[,1]),
              aes(x=t,y=z),col="red") + 
    ylab("mol. fraction (ppb)") +  xlab("t (2 h steps)")
if(save_images)
    ggsave(filename = "./img/TTA_01.png",plot = g,width=10,height=3.7)

## Now we do the pruning so we focus on the UK.
## First we remove the "background Europe" from the data
Stations <- lapply(Stations,function(l) {  
    l$Obs_obj["z"] <- l$Obs_obj["z"] - 
                     l$Obs_obj["border_z"]; l
    })

## Now we prune all matricies. First find the IDs we are going to focus on
id <- which(Emissions$region %in% c("UK","Ireland"))

## Now subset the emissions to these grid cells
Emissions <- Emissions[id,]

## And only consider the SRRs on these grid cells
Stations <- lapply(Stations,function(l) {  l$B <- l$B[,id]; l})

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

#-----------------------
# 2. Variogram modelling
#-----------------------

## Variogram modelling of emissions in the UK. Here we calibrate the Emissions prior model 
## to a combination of NAEI, EDGAR etc. (Anita's flux map)
Emissions_sp <- as.data.frame(Emissions)          # Create copy data frame
Emissions_sp$z <- Emissions_sp$z /(dx_new*dy_new) # Convert to density -- shift in log space
coordinates(Emissions_sp) <- ~x+y                 # Convert to sp object
Emissions.vgm = variogram(log(z)~1,               # Create variogram model in log space
                          subset(Emissions_sp,
                                 region %in% c("UK","Ireland")),
                          width=0.4)
Emissions.fit = fit.variogram(Emissions.vgm,      # Fit spherical variogram
                              model = vgm(1,model= "Sph",1,1))
Emissions.fit2 = fit.variogram(Emissions.vgm,     # Fit exponential variogram
                               model = vgm(1,model= "Exp",1,1))
Emissions.fit3 = fit.variogram(Emissions.vgm,     # Fit Gaussian variogram
                               model = vgm(1,model= "Gau",1,1))

## The attr(Emissions.fit,"SSErr") of each one shows that the Spherical model 
## has the lowest sum of squared errors
## attr(Emissions.fit,"SSErr")
## [1] 2.212792
## attr(Emissions.fit2,"SSErr")
## [1] 2.438081
## attr(Emissions.fit3,"SSErr")
## [1] 2.268855

if(save_images) {
    ## Plot histogram of fluxes in the UK/Ireland
    png("./img/histUK.png", width=4, height=4, units="in", res=300)
    hist(Emissions_sp$z,main="",xlab=expression(flux (g/s/"degree"^2)),ylab="frequency",8); dev.off()

    ## Plot best variogram fit
    png("./img/variogram_est.png", width=4, height=4, units="in", res=300)
    plot(Emissions.vgm,Emissions.fit,col="black",xlab=("h (degrees)")); dev.off()
}
    
## Extract the covariance matrix of the field (in log space) from the variogram
Sigma_log <- variogramLine(Emissions.fit,
                           dist_vector = spDists(Emissions_sp,
                                                 Emissions_sp),
                           covariance = TRUE)

## Extract the mean (in log space)
log_mu <- mean(log(Emissions$z))

## and put into vector form
mu_log <- matrix(log_mu,nrow(Emissions),1)


#-------------------------------------------------------------------
# 3. Deal with missing observations (and simulate new ones if desired)
#------------------------------------------------------------------

obs_locs <- t(sapply(Stations,                # extract station locs
                     function(l) c(l$x,l$y)))


if(0) {
    corr_s_mat <- function(theta_s) corr_s(s = obs_locs,theta_s = theta_s)
    corr_t_mat <- function(theta_t) corr_t(t = t_axis,theta_t = theta_t)
    d_corr_s_mat <- function(theta_s) d_corr_s(s = obs_locs,theta_s = theta_s)
    d_corr_t_mat <- function(theta_t) d_corr_t(t = t_axis,theta_t = theta_t)
    
    theta_t_true <- 0.8
    theta_s_true <- 1.5
    sigma_zeta_true <- 1
    S_zeta_true <- sigma_zeta_true^2 * kronecker(corr_t_mat(theta_t_true),corr_s_mat(theta_s_true))
    
    Emissions$sim <- exp(mu_log + t(chol(Sigma_log)) %*% rnorm(n=ncol(B)))
    warning("CHANGING OBSERVATIONS TO EDGAR PREDICTIVE ONES")
    warning("CHANGING OBSERVATION ERRORS TO HAVE EXACT ERROR AND DISCREPANCY CHARACTERISTICS")
    warning("CHANGING OBSERVATION ERRORS TO 1")
    Stations <- lapply(Stations,function(x) { x$Obs_obj["z"] <- as.numeric(x$C[which(!(colSums(x$C) == 0)),] %*% x$B %*% Emissions$z) ; x})
    s_obs <- Reduce("rbind",lapply(Stations,function(l) l$Obs_obj@df)) %>%
        arrange(t) %>%
        mutate(std=1) %>%
        mutate(dis = as.numeric(C %*% t(chol(S_zeta_true)) %*% rnorm(n = nrow(B))),
               obs_err =  rnorm(n = nrow(C),sd = std),
               z = z + dis + obs_err)
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
Qobs <- as(diag(1/s_obs$std^2),"dgCMatrix")


#-----------------------
# 4. EM algorithm
#-----------------------

mode <- c(exp(mu_log),          # Yf
          B %*% Emissions$z)    # Ym
n_EM <- 25L                    # number of EM iterations

s_obs_all <- s_obs              # copy for validation/testing

## For reproducibility we save the indices that were used for testing
## These were generated using test_idx <- sort(sample(nrow(s_obs),round(0.8*nrow(s_obs))))
load(system.file("extdata","train_idx.rda", package = "atminv"))

## Create the training matrices
C_train <- C[train_idx,]
Qobs_train <- Qobs[train_idx,train_idx]
s_obs_train <- s_obs[train_idx,]

## Now set up the EM algorithm. Comments in-line
EM_alg <- EM(s_obs = s_obs_train, # use training data observations
             C_m = C_train,       # incidence matrix
             Qobs = Qobs_train,   # precision matrix
             B = B,               # SRR
             S_f_log = Sigma_log, # log of flux-field covariance matrix
             mu_f_log = mu_log,   # log of flux-field mean
             Yf_thresh = 1e-8,    # threshold when to assume convergence
             Y_init =  mode,      # initial value
             theta_init = c(10^2,0.5,0.5), # initial value
             n_EM = n_EM,         # number of EM iterations
             t_mol = t_axis,      # time axis on which to do inference
             s_mol = obs_locs)    # spatial locs on which to do inference

## Now that EM_alg is setup we iterate through it
if(!load_results) {
    for(i in 1:(n_EM-2)) {
        X <-  EM_alg(max_E_it = 1e6,  # max E steps (very high)
                     max_M_it = 10,   # max M steps (comp. intensive, keep low)
                     fine_tune_E = (i %in% c(1,2,8)))  # fine-tune these steps 
                                                       # (needed for positive definiteness of Hessian)
    }
} else {
    load(system.file("extdata","UK_results_X.rda", package = "atminv"))
}

if(cache_results) {
    ## Remove some unwanted large objects and save
    X$lap_approx$Ym <- X$lap_approx$S_mm <- X$lap_approx$S_mf <- X$lap_approx$S_fm <- NULL
    save(X,file="inst/extdata/UK_results_X.rda")
}


#-----------------------
# 5. HMC
#-----------------------

## Find the discrepancy term covariance matrix
corr_zeta <- corr_zeta_fn(obs_locs,t_axis,
                          X$theta[2,n_EM],
                          X$theta[3,n_EM])
S_zeta <- X$theta[1,n_EM] * corr_zeta
Q_zeta <- chol2inv(chol(S_zeta))

## Now find the derivative functions (same as Laplace equations)
lap2 <- Yf_marg_approx_fns(s = s_obs_train,      # training ob locs
                           C_m = C_train,        # incidence matrix
                           Qobs = Qobs_train,    # precision matrix
                           B = B,                # SRR
                           S_zeta = S_zeta,      # discrepancy cov matrix
                           mu_f_log = mu_log,    # mean of log-flux field
                           S_f_log = Sigma_log,  # covariance of log-flux field
                           ind=X$ind)            # indices to consider (all)

## Scaling matrix (base on Laplace approximation E-step)
M <- diag(1/pmax(X$lap_approx$Yf,100)^2)

## Step size generator
eps_gen <- function() runif(n=1,min = 0.07, max = 0.071)

## Number of steps for each sample
L <- 20L

## Now set up the sampler
sampler <- hmc_sampler(U = lap2$logf,        # log density
                       dUdq = lap2$gr_logf,  # gradient of log density
                       M = M,                # scaling matrix
                       eps_gen = eps_gen,    # step-size generator
                       L = L,                # number of steps
                       lower = rep(0,length(X$ind)))  # lower-bound (0)

## More HMC parameters and variables
N <- 10000                          # number of HMC steps
q <- matrix(0,nrow(M),N)            # where to store the samples
qsamp <- pmax(X$lap_approx$Yf,100)  # first 'sample'
dither <- 1                         # no dithering
i <- count <- 1                     # count variables

## Now run the sampler
if(!load_results) {
    while(i < N) {
        qsamp <- sampler(q = qsamp)
        if(count  == dither) {
            q[,i] <- qsamp
            count = 1
            i <- i + 1
            print(paste0("Sample: ",i," Acceptance rate: ",(nrow(unique(t(q)))-1)/i))
        } else {
            count <- count + 1
        }
    }
} else {
    load(system.file("extdata","UK_results_q.rda", package = "atminv"))
}

if(cache_results) {
    save(q,i,file="inst/extdata/UK_results_q.rda")
}

#----------------------------
# 6A. Results and plots: FLUX
#----------------------------

## Arrange samples into a long data frame
Q <- as.data.frame(t(q[1:length(X$ind),1:(i-1)])) %>%
    tidyr::gather(t,z) %>%
    separate(t, into=c("V","s"),sep="V") %>%
    select(-V) %>%
    mutate(x = Emissions$x[X$ind[as.numeric(s)]],
           y = Emissions$y[X$ind[as.numeric(s)]]) %>%
    select(-s)

## Fuse Q aboce with the Emissions data frame (after replacing "z" with "NAEI")
Q2 <- left_join(Q,select(mutate(Emissions,NAEI=z),
                         x,y,NAEI,region,subregion))

## Violin Plot of Ireland
ggplot(mutate(subset(Q2,region=="Ireland"),y=-y)) + geom_boxplot(aes(x=1,y=z)) +
    geom_point(aes(x=1,y=NAEI),col="red") +
    facet_grid(y~x) + ylim(0,5000) +
    theme(panel.background = element_rect(fill='white', colour='white'),panel.grid=element_blank(),axis.ticks=element_blank(),
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_blank())

## Median emissions map: UK and Ireland (land only)
Emissions$x_mean <- apply(q[,1000:(i-1)],1,median)
Em50 <- LinePlotTheme() + geom_tile(data=Emissions,aes(x,y,fill=pmax(pmin(x_mean,2000),-2000)))  +
    scale_fill_distiller(palette="Spectral",trans="reverse",guide=guide_legend(title = "flux (g/s)")) +
    geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
    coord_map(xlim=c(-12,2),ylim=c(49,60)) + xlab("lon (degrees)") + ylab("lat (degrees)") +
    theme(axis.title.y = element_text(vjust = 1))
if(save_images)
    ggsave(filename = "./img/Em50.png",plot = Em50,width=6.4,height=6.7)

## 95-percentile emissions map: UK and Ireland (land only)
Emissions$x_95 <- apply(q[,1000:(i-1)],1,quantile,0.95)
Em95 <- LinePlotTheme() + geom_tile(data=Emissions,aes(x,y,fill=pmax(pmin(x_95,2000),-2000)))  +
    scale_fill_distiller(palette="Spectral",trans="reverse",guide=guide_legend(title = "flux (g/s)")) +
    geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
    coord_map(xlim=c(-12,2),ylim=c(49,60)) + xlab("lon (degrees)") + ylab("lat (degrees)") +
    theme(axis.title.y = element_text(vjust = 1))
if(save_images)
    ggsave(filename = "./img/Em95.png",plot = Em95,width=6.4,height=6.7)

## Difference between median and flux inventory
Emissions$x_mean <- apply(q[,1000:(i-1)],1,median)
Em50_diff <- LinePlotTheme() + geom_tile(data=Emissions,aes(x,y,fill=x_mean - z))  +
    scale_fill_distiller(palette="Spectral",trans="reverse",guide=guide_legend(title = "est - inventory (g/s)")) +   geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
    coord_map(xlim=c(-12,2),ylim=c(49,60)) + xlab("lon") + ylab("lat")

## Function to plot histogram and overlayed Laplace approximation of flux field
comp_density <- function(j,bw,xu,xl=0) {
    x <- seq(xl,xu,by=1)
    plot(density(q[j,1:(i-1)],bw = bw),
         xlim=c(xl,xu),
         xlab=c("flux (g/s)"),
         ylab="[Yf | Zm]",main="")
    lines(x,dnorm(x,mean=X$lap_approx$Yf[j],sd=sqrt(X$lap_approx$S_ff[j,j])),lty=2)
}

## Comparison between Laplace approximation and histogram at some locations
if(save_images) {
    png("./img/density40.png", width=4, height=4, units="in", res=300)
    comp_density(40,78,2000); dev.off()
    png("./img/density76.png", width=4, height=4, units="in", res=300)
    comp_density(76,110,3000,0); dev.off()
}

## Find the emissions contribution by Ireland that are in the sea
Sea_Ir <-  subset(yx,Ir_territory==1 & (is.na(region) | !(region == "UK" | region == "Ireland")))
Ir <- which(Emissions$region == "Ireland")
Sea_emissions_Ir <- sum(Sea_Ir$z*3600*24*365)  # convert to total flux in a year

## Find the emissions contribution by Ireland that are in the UK
Sea_UK <-  subset(yx,UK_territory==1 & (is.na(region) | !(region == "UK" | region == "Ireland")))
Sea_emissions_UK <- sum(Sea_UK$z*3600*24*365)  # convert to total flux in a year

options("scipen"=-100, "digits"=4)
print(paste0("Mean UK: ", mean(apply(q[-Ir,1000:(i-1)],2,sum))*3600*24*365 + Sea_emissions_UK))
print(paste0("Sd UK: ", sd(apply(q[-Ir,1000:(i-1)],2,sum))*3600*24*365))
print(paste0("Inventory UK: ", sum(subset(yx,UK_territory == 1)$z) *3600*24*365))

print(paste0("Mean Ireland: ", mean(apply(q[Ir,1000:(i-1)],2,sum))*3600*24*365 + Sea_emissions_Ir))
print(paste0("Sd Ireland: ", sd(apply(q[Ir,1000:(i-1)],2,sum))*3600*24*365))
print(paste0("Inventory Ireland: ", sum(subset(yx,Ir_territory == 1)$z) *3600*24*365))

## London
print(paste0("Median London: ", median(q[109,1000:(i-1)])*3600*24*365))
print(paste0("Inventory London: ", Emissions$z[109]*3600*24*365))


#--------------------------------------
# 6B. Results and plots: MOLE FRACTION
#--------------------------------------

## Sample new mole fractions
mf_samp <- matrix(0,nrow(B),(i-1))                          # initialise mole-fraction sample matrix
mu_post <- matrix(0,nrow(B),(i-1))                          # initialise mole-fraction sample matrix
Q_post <- t(C_train) %*% Qobs_train %*% C_train + Q_zeta    # Compute posterior precision
S_post <- as.matrix(chol2inv(chol(Q_post)))                 # Compute posterior covariance
L_S_post <- t(chol(S_post))                                 # Cholesky
SQB <- S_post %*% Q_zeta %*% B[,X$ind]                      # Samples MF  
StCQoz <- S_post %*% (t(C_train) %*% Qobs_train %*% s_obs_train$z)
for (j in 1:(i-1)){
    mu_post[,j] <- as.vector(SQB %*% q[,j] + StCQoz)
    mf_samp[,j] <- as.vector(mu_post[,j] + L_S_post %*% rnorm(nrow(B)))
}

## Plot mole fraction at stations in January
plot_mf_Jan <- function(i) {
    # Get out January samples for station i
    stat1 <- seq(i,1440,by=4)  
    
    # Find certain statistics of this station
    mu <- apply(mf_samp[stat1,-c(1:1000,i)],1,median)
    uq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.75)
    lq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.25)
    uuq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.95)
    llq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.05)
    
    # Put into data frame
    df <- data.frame(mean = mu, uq=uq,lq=lq,uuq=uuq,llq=llq,t = 1:length(mu))
    
    # Plot the median and credibility intervals
    g <- LinePlotTheme() + geom_ribbon(data=df,aes(x=t,ymax=uq,ymin=lq),alpha=0.6) +
        geom_ribbon(data=df,aes(x=t,ymax=uuq,ymin=llq),alpha=0.3) +
        geom_point(data=subset(s_obs[-train_idx,],obs_name==s_obs[i,]$obs_name & t < 361),
                   aes(x=t,y = z),colour='red',size=3,shape=4) +
        geom_point(data=subset(s_obs_train,obs_name==s_obs[i,]$obs_name & t < 361),
                   aes(x=t,y = z),colour='black',size=3,shape=4)+ ylab("Ym (ppb)") +xlab("")
    
}

## Do the plots at each of the stations and save
g_MHD <- plot_mf_Jan(1) + ggtitle("MHD")
g_RGL <- plot_mf_Jan(2) + ggtitle("RGL")
g_TAC <- plot_mf_Jan(3) + ggtitle("TAC")
g_TTA <- plot_mf_Jan(4) + ggtitle("TTA") + xlab("t (2 h steps)")
g_all <- arrangeGrob(g_MHD,g_RGL,g_TAC,g_TTA,ncol=1)
if(save_images)
    ggsave(filename = "./img/MF_Stations.png",plot = g_all,width=13,height=14)

## Is lognormal spatial process better fit?
pred_log <- NULL
for(i in 1:nrow(Emissions)) {
    Emissions.vgm_i = variogram(log(z)~1, subset(Emissions_sp[-i,],region %in% c("UK","Ireland")),width=0.4)
    Emissions.fit_i = fit.variogram(Emissions.vgm, model = vgm(0.8,model= "Sph",3.33,0.005))
    pred_log_i <- gstat(NULL,"log.z",log(z)~1,Emissions_sp[-i,],model=Emissions.fit_i) %>%
        predict(newdata=Emissions_sp[i,])
    pred_log <- rbind(pred_log,cbind(pred_log_i@data,
                                     log.z = log(Emissions_sp[i,]$z),
                                     z = Emissions_sp[i,]$z,
                                     z.pred = exp(pred_log_i$log.z.pred + 0.5*pred_log_i$log.z.var),
                                     z.var = (exp(pred_log_i$log.z.pred + 0.5*pred_log_i$log.z.var))^2*
                                         (exp(pred_log_i$log.z.var) - 1)))
}

res_norm <- (pred_log["log.z"] - pred_log["log.z.pred"])/sqrt(pred_log["log.z.var"])
if(save_images) {
    png("./img/norm_resid_log_flux.png", width=10, height=8, units="in", res=300)
    plotdist(res_norm[,1],"norm",para=list(mean=0,sd=1)); dev.off()
}

#--------------------------------------
# 7. Other exploratory plotting fns
#--------------------------------------

### The below were used to explore the data. However they are currently inoperational with the
### the new format. They are left here as they can be fixed with a little bit of effort.

## PLOT ONE MODEL OUTPUT
plot_field <- function(g,field,i) {
    g + geom_tile(data=filter(field,t==i & z > 1e6),aes(x=x,y=y,fill=z))
}

## PLOT ALL MODEL OUTPUTS
plot_time <- function(i) {
    g <- LinePlotTheme() %>%
        plot_field(Model_C_MHD,i) %>%
        plot_field(Model_C_RGL,i) %>%
        plot_field(Model_C_TAC,i) %>%
        plot_field(Model_C_TTA,i) +
        scale_fill_gradient(low="light yellow",high="purple") +
        geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
        coord_fixed(xlim=c(-14,6),ylim=c(45,65))
}

## PLOT MODEL PREDICTED VS. OBSERVATIONS
plot_Obs <- function(des) {
    Zmodel <- des$C2 %*% Emissions$z
    Zpred <- des$C2 %*% Emissions$x_mode
    plot(des$Obs_obj["z"],ylim = c(0,70))
    lines(Zpred,col='red')
    lines(Zmodel)
    message(sum((Zmodel - des$Obs_obj["z"])^2))
    message(sum((Zpred - des$Obs_obj["z"])^2))
}

if(0) {
    plot_Obs(Stations$RGL)

    Zpred_RGL <- scale(C[["RGL"]] %*% Emissions$z)
    plot(Obs_RGL$t,scale(Obs[["RGL"]]$z)); lines(t_axis,Zpred_RGL)
    
    Zpred_TAC <- scale(C[["TAC"]] %*% Emissions$z)
    plot(Obs_TAC$t,scale(Obs[["TAC"]]$z)); lines(t_axis,Zpred_TAC)
    
    Zpred_TTA <- scale(C[["TTA"]] %*% Emissions$z)
    plot(Obs_TTA$t,scale(Obs[["TTA"]])); lines(t_axis,Zpred_TTA)
}

## DO ANIMATION
if(0) {
    library(animation)
    oopt <- animation::ani.options(interval = 0.1)
    
    FUN2 <- function() {
        t = 1:372
        lapply(t, function(i) {
            print(plot_time(i))
            ani.pause()
        })
    }
    
    FUN2()
    saveHTML(FUN2(), autoplay = FALSE, loop = FALSE, verbose = FALSE, outdir = ".",
             single.opts = "'controls': ['first', 'previous', 'play', 'next', 'last', 'loop', 'speed'], 'delayMin': 0")
}

