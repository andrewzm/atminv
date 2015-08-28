library(dplyr)
library(tidyr)
library(Matrix)
library(devtools)
library(ggplot2)
library(grid)
library(gridExtra)
load_all("../../../../PostDoc Bristol/R Code/pkg/MVST/")
load_all("../../../pkgs/hmc/")
load_all("../../../pkgs/atminv")

rm(list=ls())
md5_wrapper <- MVST::md5_cache("~/cache/Assim")

prune = FALSE  ## deprecated
match_to_process = TRUE

miss_val = -999
ymin = 36.469
xmin = -14.124
xco <- 100
yco <- 0
dx = 0.352
dy = 0.234
yx <- as.data.frame(cbind(expand.grid(y = seq(ymin, ymin + 127*dy,length=128),
                                      x = seq(xmin, xmin + 127*dx,length=128))))

yx$breaks_x <- cut(yx$x,seq(xmin-0.001, xmin + 127*dx + 0.001,length=64),labels=FALSE)
yx$breaks_y <- cut(yx$y,seq(ymin-0.001, ymin + 127*dy + 0.001,length=64),labels=FALSE)


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

## LOAD MAPS
Europe <- map_data("world") %>%
  subset(region %in% c("UK","France","Ireland")) %>%
  mutate(id = group, x=long,y=lat)

yx <- mutate(yx, in_land= attribute_polygon(yx,Europe),
             id = in_land) %>%
  left_join(unique(select(Europe,region,id,subregion)))

## Add UK territory
UK_idx <- read.table("../data/UK_gridcells.txt")[,1]
yx$UK_territory <- 0
yx$UK_territory[UK_idx] <- 1

## Add Ireland territory
Ir_idx <- read.table("../data/Ireland_gridcells.txt")[,1]
yx$Ir_territory <- 0
yx$Ir_territory[Ir_idx] <- 1


## Model Emissions field
yx <- yx %>%
  #mutate(z = read.table("../data/CH4_emissions_g_per_s",skip=10)[,1], t = 1)
  mutate(z = read.table("../data/CH4_emissions_scaled_natural",skip=10)[,1], t = 1)
#ggplot() + geom_tile(data=Emissions,aes(x,y,fill=pmin(z,2000))) + coord_fixed()



Emissions <- group_by(yx,breaks_x,breaks_y) %>%
              summarise(y = mean(y), x = mean(x), in_land = max(in_land), region = region[1],
                        subregion = subregion[1],z=sum(z),UK_territory =mean(UK_territory))

dx_new <- dx*2
dy_new <- dy*2

g <- LinePlotTheme() + geom_tile(data=subset(Emissions,x > -13 & x < 3 & y > 48 & y < 61),
                                 aes(x,y,fill=pmax(pmin(z,2000),-2000)))  +
  scale_fill_distiller(palette="Spectral", trans="reverse", guide = guide_legend(title="flux (g/s)")) +
  geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
  coord_map(xlim=c(-12,2),ylim=c(49,60)) + xlab("lon (degrees)") + ylab("lat (degrees)") +
    theme(axis.title.y = element_text(vjust = 1))
ggsave(filename = "../NAEI.png",plot = g,width=6.4,height=6.7)

conv_B_df_to_mat <- function(B_df) {
  B_df %>%
    select(-y,-x,-in_land,-id,-region,-subregion) %>%
    as.matrix() %>%
    t() %>%
    as("dgCMatrix")

}

### Read model data from file and put into long table format and focus on UK
### I also should be multiplying by 1e9 but this produces VERY large numbers
read_process_models <- function(des,prune=T) {
  Model_B_list <- list()
  for(i in 1:length(des$model_filenames)) {
      browser()
    Model_B_list[[i]] <- read.table(des$model_filenames[i],skip = 15) %>%
      cbind(yx) %>%
      gather(t,z,-x,-y, -in_land,-region,-id,-subregion,-breaks_x,-breaks_y) %>%
      separate(t, into = c("V", "t"), sep = "V") %>%
      mutate(t = as.numeric(t)) %>%
      select(-V) %>%
      group_by(breaks_x,breaks_y,t) %>%
      summarise(y = mean(y), x = mean(x), in_land = max(in_land), region = region[1],
                subregion = subregion[1],z=mean(z))
    if(i > 1) {
      Model_B_list[[i]]$t <- Model_B_list[[i-1]]$t[nrow(Model_B_list[[i-1]])] + Model_B_list[[i]]$t
    }
  }

  des$Model_B <- do.call("rbind",Model_B_list)

  #if(prune) Model_B <- subset(Model_B,(max_detect > max(max_detect)/20))# & !(in_land==0))

  #des$B <- conv_B_df_to_mat(Model_B)
  #des$C <- des$B ## For now assume we have observations everywhere in space and time

#   des$Model_B <- Model_B %>%
#     gather(t,z,-x,-y, -in_land,-region,-id,-subregion,-breaks_x,-breaks_y) %>%
#     separate(t, into = c("V", "t"), sep = "V") %>%
#     mutate(t = as.numeric(t)) %>%
#     select(-V) %>%
#     group_by(breaks_x,breaks_y,t) %>%
#     summarise(y = mean(y), x = mean(x), in_land = max(in_land), region = region[1],
#               subregion = subregion[1],z=mean(z))

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

### Read obs data from file and put into long table format and MVST::Obs object
read_process_data <- function(des) {
  this_x <- des$x
  this_y <- des$y

  id <- which(Emissions$region %in% c("UK","Ireland"))
  border_z <- as.numeric(des$B[,-id] %*% Emissions[-id,]$z)

  Obs_df_list <- list()
  for(i in 1:length(des$obs_filenames)) {
    Obs_df_list[[i]] <-
      read.table(des$obs_filenames[i],skip=1) %>%
      mutate(z = V1, std = V2, x = this_x, y = this_y) %>%
      select(-V1,-V2)
  }
  Obs_df <- do.call("rbind",Obs_df_list) %>%
            mutate(border_z = border_z,
                   t= t_axis,
                   z = ifelse(z == miss_val,NA,z))

  des$Obs_obj <- Obs(Obs_df,name=des$obs_filenames[i])


#   des$Obs_obj <-
#     read.table(des$obs_filename,skip=1) %>%
#     mutate(z = V1, std = V2, x = this_x, y = this_y, t = t_axis) %>%
#     select(-V1,-V2) %>%
#     mutate(border_z = border_z) %>%
#     subset(!(z==miss_val)) %>%
#     Obs(name=des$obs_filename)

  x <- rep(0,length(t_axis))
  x[des$Obs_obj["t"]] <- 1
  des$C <- diag(x)

  C2 <- des$C[des$Obs_obj["t"],]
  ## I'm going to normalise... see notes

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

### Find the important locations to keep (areas of influence)
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

# Read process models and form matrices
Stations <- md5_wrapper(lapply,Stations,read_process_models,prune=FALSE)

# Read data
t_axis <- 1:max(Stations$RGL$Model_B$t)
Stations <- lapply(Stations,read_process_data)

# Remove the 0.05 quantile
MHD_5 <- quantile(Stations$MHD$Obs_obj["z"],0.05,na.rm=T)
#warning("Fixing MHD_5 to 1893.091")
#MHD_5 <- 1893.091
Stations <- lapply(Stations,function(l) {l$Obs_obj["z"] <- l$Obs_obj["z"] - MHD_5; l})

Station_locs <- t(sapply(Stations,function(l) data.frame(x=l$x,y=l$y))) %>% as.data.frame()
Station_locs$name <- row.names(Station_locs)
Emissions$B <- Stations$TTA$B[1,]
g <- LinePlotTheme() + geom_tile(data=subset(Emissions,x > -13 & x < 3 & y > 48 & y < 61),
                                 aes(x,y,fill=B))  +
    scale_fill_distiller(palette="BuPu",trans="sqrt",guide = guide_legend(title="SRR (s/ng)")) +
    geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
    coord_map(xlim=c(-12,2),ylim=c(49,60)) + xlab("lon (degrees)") + ylab("lat (degrees)") +
    geom_point(data=Station_locs,aes(x=unlist(x),y=unlist(y)),col="red",size=4) +
    geom_text(data=Station_locs,aes(x=unlist(x),y=unlist(y)+0.3,label=unlist(name)),col="red") +
    theme(axis.title.y = element_text(vjust = 1))
ggsave(filename = "../TTA_01_01.png",plot = g,width=6.4,height=6.7)

g <- LinePlotTheme() + geom_point(data=subset(Stations$TAC$Obs_obj@df,t < 372),aes(x=t,y=z)) +
    geom_line(data=data.frame(t=1:372,z=(Stations$TAC$B[1:372,] %*% Emissions$z)[,1]),aes(x=t,y=z),col="red") +
    ylab("mol. fraction (ppb)") +
    xlab("t (2 h steps)")
ggsave(filename = "../TAC_01.png",plot = g,width=10,height=3.7)

g <- LinePlotTheme() + geom_point(data=subset(Stations$TAC$Obs_obj@df,t < 372),aes(x=t,y=z)) +
                  geom_line(data=data.frame(t=1:372,z=(Stations$TAC$B[1:372,] %*% Emissions$z)[,1]),aes(x=t,y=z),col="red") + ylab("mol. fraction (ppb)") +  xlab("t (2 h steps)")
ggsave(filename = "../TTA_01.png",plot = g,width=10,height=3.7)

# # Remove offset from data
# offset <- quantile(Stations$RGL$Obs_obj["z"],0.5)
# Stations <- lapply(Stations,function(l) {  l$Obs_obj["z"] <- l$Obs_obj["z"] - offset;  l})

# Now prune problem to the UK.
# First we remove the "background Europe" from the data
Stations <- lapply(Stations,function(l) {  l$Obs_obj["z"] <- l$Obs_obj["z"] - l$Obs_obj["border_z"]; l})

# Now we prune all matricies
id <- which(Emissions$region %in% c("UK","Ireland"))
Emissions <- Emissions[id,]
Stations <- lapply(Stations,function(l) {  l$B <- l$B[,id]; l})

# Now construct the required objects
m_obs <- length(Stations)
# B <- gdata::interleave(as.matrix(Stations$MHD$B),
#                        as.matrix(Stations$RGL$B),
#                        as.matrix(Stations$TAC$B),
#                        as.matrix(Stations$TTA$B))
B <- do.call(gdata::interleave,args = lapply(Stations, function(l) as.matrix(l$B)))

# C <- gdata::interleave(as.matrix(diag(Stations$MHD$C)),
#                        as.matrix(diag(Stations$RGL$C)),
#                        as.matrix(diag(Stations$TAC$C)),
#                        as.matrix(diag(Stations$TTA$C))) %>%
#   c() %>% diag()
C <- do.call(gdata::interleave,args = lapply(Stations, function(l) as.matrix(diag(l$C)))) %>%
  c() %>% diag()
C <- as(C,"dgCMatrix")

breakme
## Calibrate Emissions prior model to EDGAR
library(sp)
library(gstat)
Emissions_sp <- as.data.frame(Emissions)
Emissions_sp$z <- Emissions_sp$z /(dx_new*dy_new) # Convert to density
log_mu_density <- mean(log(Emissions_sp$z))
coordinates(Emissions_sp) <- ~x+y
Emissions.vgm = variogram(log(z)~1, subset(Emissions_sp,region %in% c("UK","Ireland")),width=0.4)
Emissions.fit = fit.variogram(Emissions.vgm, model = vgm(1,model= "Sph",1,1))

png("../histUK.png", width=4, height=4, units="in", res=300)
hist(Emissions_sp$z,main="",xlab=expression(flux (g/s/"degree"^2)),ylab="frequency",8); dev.off()
png("../variogram_est.png", width=4, height=4, units="in", res=300)
plot(Emissions.vgm,Emissions.fit,col="black",xlab=("h (degrees)")); dev.off()


Sigma_log <- variogramLine(Emissions.fit,
                        dist_vector = spDists(Emissions_sp,
                                       Emissions_sp),
                        covariance = TRUE)

log_mu <- mean(log(Emissions$z))
mu_log <- matrix(log_mu,nrow(Emissions),1)

# D <- as.matrix(dist(cbind(Emissions$x,Emissions$y)))
# Sigma_log <- Matern(r = D,nu = 1.5,var = log_var,kappa = kappa_from_l(5,nu = 1.5))
#mu_Em <- matrix(Emissions$z,nrow(Emissions),1)

### New EM algorithm
n_EM = 10
mode <- c(exp(mu_log),B %*% Emissions$z)

obs_locs <- t(sapply(Stations,function(l) c(l$x,l$y)))


# warning("CHANGING OBSERVATIONS TO SIMULATION ONES")
# Stations <- lapply(Stations,function(x) { x$Obs_obj["z"] <- as.numeric(x$C[which(!(colSums(x$C) == 0)),] %*% x$B %*% Emissions$sim) ; x})

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
  #s_obs <- Reduce("rbind",lapply(Stations,function(l) l$Obs_obj@df))
  s_obs <- do.call(gdata::interleave,args = lapply(Stations, function(l) l$Obs_obj@df))
  warning("Replacing -999 std with mean std")
  s_obs$std[which(s_obs$std == -999)] <- mean(subset(s_obs,std>0)$std)
  warning("Removing all negative values and NA")
  idx <- which(s_obs$z > 0 & !is.na(s_obs$z))
  C <- C[idx,]
  s_obs <- s_obs[idx,]
}


Qobs <- as(diag(1/s_obs$std^2),"dgCMatrix")

#mu_log <- matrix(log(abs(Emissions$z)))
mode <- c(exp(mu_log),B %*% Emissions$z)
n_EM <- 100L

s_obs_all <- s_obs
test_idx <- sort(sample(nrow(s_obs),round(0.8*nrow(s_obs))))
C_test <- C[test_idx,]
Qobs_test <- Qobs[test_idx,test_idx]
s_obs_test <- s_obs[test_idx,]

breakme
EM_alg <- EM(s_obs = s_obs_test,
        C_m = C_test,
        Qobs = Qobs_test,
        B = B,
        S_f_log = Sigma_log,
        mu_f_log = mu_log,
        Yf_thresh = 1e-8,
        Y_init =  mode,
        theta_init = c(10^2,0.5,0.5),
        n_EM = n_EM,
        t_mol = t_axis,
        s_mol = obs_locs)
for(i in 1:99) {
  X <-  EM_alg(max_E_it = 1e6,
               max_M_it = 10,
               fine_tune_E = (i %in% c(1,2,8)))
}

save(X,file="~/Desktop/UK_results4.dat")
i <- 55
S_zeta <- X$theta[1,i] * kronecker(corr_t(t_axis,X$theta[2,i]),corr_s(s = obs_locs,theta_s = X$theta[3,i]))

lap2 <- Yf_marg_approx_fns(s = s_obs_test,
                           C_m = C_test,
                           Qobs = Qobs_test,
                           B = B,
                           S_zeta = S_zeta, #sigma_zeta_true^2 * diag(rep(1,nrow(s_obs)))
                           mu_f_log = mu_log,
                           S_f_log = Sigma_log,
                           ind=X$ind)

### HMC SAMPLING
M <- diag(1/pmax(X$lap_approx$Yf,100)^2)
eps_gen <- function() runif(n=1,min = 0.07, max = 0.071)
L <- 20L

### GOOD SETTINGS L = 100, eps = 0.03--0.032
sampler <- hmc_sampler(U = lap2$logf,dUdq = lap2$gr_logf,
                       M = M,
                       eps_gen = eps_gen,
                       L = L,
                       lower = rep(0,length(X$ind)))
N <- 10000
q <- matrix(0,nrow(M),N)
qsamp <- pmax(X$lap_approx$Yf,100)
dither <- 1
i <- count <- 1

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
save(list=ls(),file="~/Desktop/UK_results_all4.dat")

# mode <- optim(X$lap_approx$Yf,
#                fn = lap2$logf,
#                gr = lap2$gr_logf,
#                method = "BFGS",
#                control=list(trace=4,reltol=1e-8,maxit=1e6,))$par

if(0) {
  ### SLICE SAMPLING
  slice_sampler <- slice::oaatslice(d = length(X$ind), scale=X$lap_approx$Yf)
  log_f <- function(x) -lap2$logf(x)

  N2 <- 1000
  q2 <- matrix(0,length(X$ind),N)
  qsamp2 <- pmax(X$lap_approx$Yf,10)
  dither2 <- 1
  i2 <- count2 <- 1

  while(i2 < N) {
    qsamp2 <- slice_sampler(qsamp2,log_f,learn = i2 <= 200)


    if(count2  == dither2) {
      q2[,i] <- qsamp2
      count2 = 1
      i2 <- i2 + 1
      print(i2)
    } else {
      count2 <- count2 + 1
    }
  }
}

Q <- as.data.frame(t(q[1:length(X$ind),1:(i-1)])) %>%
  tidyr::gather(t,z) %>%
  separate(t, into=c("V","s"),sep="V") %>%
  select(-V) %>%
  mutate(x = Emissions$x[X$ind[as.numeric(s)]],
         y = Emissions$y[X$ind[as.numeric(s)]]) %>%
  select(-s)

# Q2 <- rbind(Q,data.frame(x = Emissions$x[setdiff(1:nrow(Emissions),X$ind)],y = Emissions$y[setdiff(1:nrow(Emissions),X$ind)], z = 0))
Q2 <- left_join(Q,select(mutate(Emissions,NAEI=z),x,y,NAEI,region,subregion))

ggplot(mutate(subset(Q2,region=="Ireland"),y=-y)) + geom_boxplot(aes(x=1,y=z)) +
  geom_point(aes(x=1,y=NAEI),col="red") +
  facet_grid(y~x) + ylim(0,5000) +
  theme(panel.background = element_rect(fill='white', colour='white'),panel.grid=element_blank(),axis.ticks=element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_blank())
#   geom_point(data=subset(st_grid),aes(x=as.factor(s),y = Yf),colour='red',size=3) +
#   geom_point(data=subset(s_obs,t==1),aes(x=as.factor(s),y = 0),colour='blue',size=3)
# geom_point(data=data.frame(s=s_axis[X$ind],y=mode[1:length(X$ind)],sd=sqrt(diag(gr2)[1:length(X$ind)])),aes(x=factor(s),y = y),colour='blue',size=3) +
#   geom_point(data=data.frame(s=s_axis[X$ind],y=mode[1:length(X$ind)],sd=sqrt(diag(gr2)[1:length(X$ind)])),aes(x=factor(s),y = y + sd),colour='blue',size=3) +
#   geom_point(data=data.frame(s=s_axis[X$ind],y=mode[1:length(X$ind)],sd=sqrt(diag(gr2)[1:length(X$ind)])),aes(x=factor(s),y = y-sd),colour='blue',size=3)
if(0) {
  ### BAYES LINEAR
  mu_Y <- exp(mu_log)
  Sigma_exp <- exp(Sigma_log[1,1]) * (exp(Sigma_log) - 1)
  Cov_YZ <- Sigma_exp %*% t(C %*% B)
  Cov_Z <- (C %*% B) %*% Sigma_exp %*% t(C %*% B) + solve(Qobs)
  Emissions$x_mean <- as.numeric(mu_Y + Cov_YZ %*% solve(Cov_Z) %*% (s_obs$z - C %*% B %*% mu_Y) )
  Emissions$std <- sqrt(diag(Sigma_exp - Cov_YZ %*% solve(Cov_Z) %*% t(Cov_YZ)))
}

### PLOTTING
Emissions$x_mean <- apply(q[,1000:(i-1)],1,median)
Em50 <- LinePlotTheme() + geom_tile(data=Emissions,aes(x,y,fill=pmax(pmin(x_mean,2000),-2000)))  +
  scale_fill_distiller(palette="Spectral",trans="reverse",guide=guide_legend(title = "flux (g/s)")) +
  geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
  coord_map(xlim=c(-12,2),ylim=c(49,60)) + xlab("lon (degrees)") + ylab("lat (degrees)") +
    theme(axis.title.y = element_text(vjust = 1))
ggsave(filename = "../Em50.png",plot = Em50,width=6.4,height=6.7)

Emissions$x_95 <- apply(q[,1000:(i-1)],1,quantile,0.95)
Em95 <- LinePlotTheme() + geom_tile(data=Emissions,aes(x,y,fill=pmax(pmin(x_95,2000),-2000)))  +
  scale_fill_distiller(palette="Spectral",trans="reverse",guide=guide_legend(title = "flux (g/s)",)) +
  geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
  coord_map(xlim=c(-12,2),ylim=c(49,60)) + xlab("lon (degrees)") + ylab("lat (degrees)") +
    theme(axis.title.y = element_text(vjust = 1))
ggsave(filename = "../Em95.png",plot = Em95,width=6.4,height=6.7)

Emissions$x_mean <- apply(q[,1000:(i-1)],1,median)
Em50_diff <- LinePlotTheme() + geom_tile(data=Emissions,aes(x,y,fill=x_mean - z))  +
    scale_fill_distiller(palette="Spectral",trans="reverse",guide=guide_legend(title = "est - inventory (g/s)")) +   geom_path(data=Europe,aes(x=long,y=lat,group=group)) +
    coord_map(xlim=c(-12,2),ylim=c(49,60)) + xlab("lon") + ylab("lat")

comp_density <- function(j,bw,xu,xl=0) {
  x <- seq(xl,xu,by=1)
  plot(density(q[j,1:(i-1)],bw = bw),xlim=c(xl,xu),xlab=c("flux (g/s)"),ylab="[Yf | Zm]",main="")
  lines(x,dnorm(x,mean=X$lap_approx$Yf[j],sd=sqrt(X$lap_approx$S_ff[j,j])),lty=2)
}

png("../density40.png", width=4, height=4, units="in", res=300)
comp_density(40,78,2000); dev.off()
png("../density76.png", width=4, height=4, units="in", res=300)
comp_density(76,110,3000,0); dev.off()

Sea_Ir <-  subset(yx,Ir_territory==1 & (is.na(region) | !(region == "UK" | region == "Ireland")))
Ir <- which(Emissions$region == "Ireland")
Sea_emissions_Ir <- sum(Sea_Ir$z*3600*24*365)

Sea_UK <-  subset(yx,UK_territory==1 & (is.na(region) | !(region == "UK" | region == "Ireland")))
Sea_emissions_UK <- sum(Sea_UK$z*3600*24*365)

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

### MOLE FRACTION AT STATIONS

corr_zeta <- corr_zeta_fn(s_obs$s[1:(m_obs)],t_axis,X$theta[2,n_EM],X$theta[3,n_EM])
S_zeta <- X$theta[1,n_EM] * corr_zeta
Q_zeta <- chol2inv(chol(S_zeta))

mf_samp <- matrix(0,nrow(B),(i-1))
mu_post <- matrix(0,nrow(B),(i-1))
Q_post <- t(C_test) %*% Qobs_test %*% C_test + Q_zeta
S_post <- as.matrix(chol2inv(chol(Q_post)))
L_S_post <- t(chol(S_post))
SQB <- S_post %*% Q_zeta %*% B[,X$ind]
StCQoz <- S_post %*% (t(C_test) %*% Qobs_test %*% s_obs_test$z)
for (j in 1:(i-1)){
    mu_post[,j] <- as.vector(SQB %*% q[,j] + StCQoz)
    mf_samp[,j] <- as.vector(mu_post[,j] + L_S_post %*% rnorm(nrow(B)))
}

### Plot mole fraction at stations in January
plot_mf_Jan <- function(i) {
    stat1 <- seq(i,1440,by=4)
    mu <- apply(mf_samp[stat1,-c(1:1000,i)],1,median)
    uq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.75)
    lq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.25)
    uuq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.95)
    llq <- apply(mf_samp[stat1,-c(1:1000,i)],1,quantile,0.05)
    df <- data.frame(mean = mu, uq=uq,lq=lq,uuq=uuq,llq=llq,t = 1:length(mu))

    g <- LinePlotTheme() + geom_ribbon(data=df,aes(x=t,ymax=uq,ymin=lq),alpha=0.6) +
                           geom_ribbon(data=df,aes(x=t,ymax=uuq,ymin=llq),alpha=0.3) +
        geom_point(data=subset(s_obs[-test_idx,],obs_name==s_obs[i,]$obs_name & t < 361),aes(x=t,y = z),colour='red',size=3,shape=4) +
        geom_point(data=subset(s_obs_test,obs_name==s_obs[i,]$obs_name & t < 361),aes(x=t,y = z),colour='black',size=3,shape=4)+ ylab("Ym (ppb)") +xlab("")


#     g <- LinePlotTheme() + geom_boxplot(data=MF2,aes(x=t,y=z))  +
#         scale_x_discrete(breaks = t_axis[seq(1,length(t_axis),length=28)]) +
#         #geom_point(data=subset(st_grid, s < 4.4),aes(x=as.factor(s),y = Yf),colour='red',size=3) +
#         geom_point(data=subset(s_obs,obs_name==s_obs[i,]$obs_name & t < 361),aes(x=as.factor(t),y = z),colour='red',size=3,shape=4)+ ylab("Ym")
#
#     MF <- as.data.frame(t(mf_samp[stat1,-c(1:1000,i)])) %>%
#         tidyr::gather(s,z) %>%
#         separate(s, into=c("V","t"),sep="V") %>%
#         select(-V) %>%
#         mutate(t = t_axis[as.numeric(t)])
#
#     MF2 <- MF %>%
#         mutate(t = as.factor(t))
#     g <- LinePlotTheme() + geom_boxplot(data=MF2,aes(x=t,y=z))  +
#         scale_x_discrete(breaks = t_axis[seq(1,length(t_axis),length=28)]) +
#         #geom_point(data=subset(st_grid, s < 4.4),aes(x=as.factor(s),y = Yf),colour='red',size=3) +
#         geom_point(data=subset(s_obs,obs_name==s_obs[i,]$obs_name & t < 361),aes(x=as.factor(t),y = z),colour='red',size=3,shape=4)+ ylab("Ym")
#     return(g)
}
options("scipen"=0, "digits"=4)
g_MHD <- plot_mf_Jan(1) + ggtitle("MHD")
g_RGL <- plot_mf_Jan(2) + ggtitle("RGL")
g_TAC <- plot_mf_Jan(3) + ggtitle("TAC")
g_TTA <- plot_mf_Jan(4) + ggtitle("TTA") + xlab("t (2 h steps)")
g_all <- arrangeGrob(g_MHD,g_RGL,g_TAC,g_TTA,ncol=1)
ggsave(filename = "../MF_Stations.png",plot = g_all,width=13,height=14)

# ### MY WRONG BAYES LINEAR
# mu_Y <- exp(diag(Sigma)/2)
# Sigma_exp <- exp(Sigma[1,1]) * (exp(Sigma) - 1)
# Qtot <- CtQoC + solve(Sigma_exp)
# R <- chol(Qtot)
# ybar = CtQoz + solve(Sigma_exp) %*% mu_Y
# x_mean <- matrix(backsolve(R,backsolve(R,ybar,transpose=T)),ncol=1)
# S <- chol2inv(R)
# Emissions$x_mean <- x_mean
# Emissions$std <- sqrt(diag(S))


### LOGORMAL PRIOR NEW
mu = 0
f  <- function(x) {
  -(-0.5 * crossprod(chol_Qobs %*% (y_tot$z - C_full %*% exp(x))) - 0.5*crossprod(chol_Q %*% (x - mu)))
}

gradf <- function(x) {
  diagexp_x <- diag(exp(x))
  -(-Q_full %*% (x - mu) + diagexp_x %*% CtQoz - diagexp_x %*% CtQoC %*% exp(x))
}

grad2f <- function(x) {
  diagexp_x <- diag(exp(x))
  -(-Q_full + diag(as.numeric(CtQoz)) %*% diagexp_x - diagexp_x %*%(CtQoC %*% diagexp_x + diag(CtQoC %*% exp(x))))
}

mode <- optim(Emissions$z*0,fn = f,gr = gradf,method = "BFGS",
              control=list(maxit=500,trace=4,parscale=log(abs(Emissions$z))))$par
Emissions$x_var <- diag(-(solve(grad2f(mode))))
Emissions$x_mean <- exp(mode + Emissions$x_var/2) ## Need to add variance
Emissions$x_mode <- exp(mode)

#### LOGNORMAL PRIOR OLD
mu <- log(abs(Emissions$z))
sigma2 <- mu^2
f1 <- function(x) log(x)
f2 <- function(x) 1/x
f <- function(x) {
  if (any(x < 0))
  {
    Inf
  } else {
    -(-0.5 * crossprod(chol_Qobs %*% (y_tot$z - C_full %*% x)) - sum(f1(x)) - sum(0.5*(f1(x) - mu)^2/sigma2))
  }
}

gradf <- function(x) {
  -(-CtQoC %*% x  + CtQoz - f2(x) + (mu * f2(x) - f1(x) * f2(x)) / sigma2)
}
mode <- optim(abs(Emissions$z),fn = f,gr = gradf,method = "BFGS",control=list(trace=4,parscale=abs(Emissions$z)))$par
#mode <- nloptr::stogo(abs(Emissions$z),fn = f,gr = gradf,lower=abs(Emissions$z)/10,upper=abs(Emissions$z)*10)$par
Emissions$x_mean <- mode

#### Metropolis-Hastings
log_pold = -1e99
xold <- abs(Emissions$z)
N <- 5000000
thin <- 100
count <- 1
xkeep <- matrix(0,length(mu),N/thin)
for(i in 1:N) {
  xnew <- xold + rnorm(n = length(mu),mean  =0, sd = sqrt(sigma2)/30)
  log_pnew <- -f(xnew)
  accept = min(1, exp(log_pnew - log_pold)) > runif(1)
  if(accept) {
    log_pold <- log_pnew
    xold <- xnew
  }
  if (i %% thin == 0) {xkeep[,count] <- xold; print(i); count = count + 1}
  print(accept)
}







## PLOT ONE MODEL OUTPUT
plot_field <- function(g,field,i) {
  g + geom_tile(data=filter(field,t==i & z > 1e6),aes(x=x,y=y,fill=z))
}

## PLOT ALL MODEL OUTPUTS
plot_time <- function(i) {
  g <- MVST::LinePlotTheme() %>%
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
  Zmodel <- des$C %*% Emissions$z
  Zpred <- des$C %*% Emissions$x_mode
  plot(des$Obs_obj["z"],ylim = c(0,70))
  lines(Zpred,col='red')
  lines(Zmodel)
  message(sum((Zmodel - des$Obs_obj["z"])^2))
  message(sum((Zpred - des$Obs_obj["z"])^2))
}
plot_Obs(Stations$RGL)

Zpred_RGL <- scale(C[["RGL"]] %*% Emissions$z)
plot(Obs_RGL$t,scale(Obs[["RGL"]]$z)); lines(t_axis,Zpred_RGL)

Zpred_TAC <- scale(C[["TAC"]] %*% Emissions$z)
plot(Obs_TAC$t,scale(Obs[["TAC"]]$z)); lines(t_axis,Zpred_TAC)

Zpred_TTA <- scale(C[["TTA"]] %*% Emissions$z)
plot(Obs_TTA$t,scale(Obs[["TTA"]])); lines(t_axis,Zpred_TTA)

## DO ANIMATION
if(F) {
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

