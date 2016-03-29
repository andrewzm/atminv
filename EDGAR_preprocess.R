library(RNetCDF)
library(dplyr)
library(ggplot2)
library(geoR)

#fid <- open.nc("../data/2014.emissions.CH4.surface.MERRA.monthly.global.nc")
fid <- open.nc("~/Desktop/v42_FT2010_CH4_2010_TOT.0.1x0.1.nc")
#print.nc(fid)
dat <- read.nc(fid)
lon <- dat[["lon"]]
lon <- data.frame(lon_idx = 1:length(lon), lon = as.numeric(lon))  %>%
    mutate(lon = lon - 360*(lon > 180))
lat <- dat[["lat"]]
lat <- data.frame(lat_idx = 1:length(lat), lat = as.numeric(lat))

#z <- dat[["emissions"]] %>% reshape2::melt()
#names(z) <- c("lon_idx","lat_idx","month","CH4flux")
z <- dat[["emi_ch4"]] %>% reshape2::melt()
names(z) <- c("lon_idx","lat_idx","CH4flux")

z <- left_join(z, lon,by="lon_idx") %>%
    left_join(lat,by="lat_idx") %>%
    dplyr::select(-lon_idx, -lat_idx)

(ggplot(z %>% filter(lon > -12 & lat > 49 & lon < 2 & lat < 61)) + geom_tile(aes(lon,lat,fill=CH4flux)) + scale_fill_distiller(palette="Spectral")) %>% FRK:::draw_world() + xlim(c(-12,2)) + ylim(c(49,61))

## Subsample
dlon_new = 0.352*2       # lon spacing
dlat_new = 0.234*2       # lat spacing

### GRID INTO OUR ONE
new_lon <- seq(-180, 180,by = dlon_new)
new_lat <- seq(-90, 90,by = dlat_new)
z$breaks_lon <- cut(z$lon,   # super x-grid for downsampling by 2
                  new_lon,
                  labels=FALSE)
z$breaks_lat <- cut(z$lat,   # super y-grid for downsampling by 2
                  new_lat,
                  labels=FALSE)
Emissions <- group_by(z,breaks_lon,breaks_lat) %>%
    summarise(lat = new_lat[breaks_lat[1]],
              lon = new_lon[breaks_lon[1]],
              CH4flux=mean(CH4flux)) ## Units are in kg m-2 s-1

library(geosphere)
Emissions$area_m2 <- apply(Emissions,1,function(df) {
    polygon <- matrix(c(df["lon"] - dlon_new/2,df["lat"] - dlat_new/2,
                        df["lon"] + dlon_new/2,df["lat"] - dlat_new/2,
                        df["lon"] + dlon_new/2,df["lat"] + dlat_new/2,
                        df["lon"] - dlon_new/2,df["lat"] + dlat_new/2),4,2,byrow = TRUE)
    return(areaPolygon(polygon))
})

Emissions <- mutate(Emissions,
                    CH4flux_tot = CH4flux * 1000 * area_m2) # convert to g then find total flux

## World emissions
gworld <- (LinePlotTheme() + geom_tile(data= Emissions, aes(lon,lat,fill=pmin(CH4flux_tot,2000))) +     scale_fill_distiller(palette="Spectral", name= "flux (g/s)\n")) %>% FRK:::draw_world() + xlab("lon (deg)") + ylab("lat (deg)")
ggsave(gworld,filename = "./img/EDGAR_world.pdf",width = 12,height=7)

## UK emissions
(ggplot(Emissions %>% filter(lon > -12 & lat > 49 & lon < 2 & lat < 61)) + geom_tile(aes(lon,lat,fill=CH4flux_tot)) + scale_fill_distiller(palette="Spectral")) %>% FRK:::draw_world() + xlim(c(-12,2)) + ylim(c(49,61)) + coord_map()

## US midwest emissions
(ggplot(Emissions %>% filter(lon > -105 & lat > 38.5 & lon < -93 & lat < 47)) + geom_tile(aes(lon,lat,fill=pmin(CH4flux_tot,2000))) + scale_fill_distiller(palette="Spectral")) %>% FRK:::draw_world() + xlim(c(-105,-61)) + ylim(c(17.5,47)) + coord_map()

#Emissions_Midwest <- filter(Emissions,lon > -105 & lat > 38.5 & lon < -93 & lat < 47)
Emissions_Midwest <- filter(Emissions,lon > -105 & lat > 17.5 & lon < -61 & lat < 47)

## Canada emissions
Emissions_Canada <- filter(Emissions,lon > -115 & lat > 40 & lon < -70.5 & lat < 69.5)
(ggplot(Emissions_Canada) + geom_tile(aes(lon,lat,fill=CH4flux_tot)) + scale_fill_distiller(palette="Spectral",limits=c(0,100))) %>% FRK:::draw_world() + xlim(c(-115,-70.5)) + ylim(c(40,69.5)) + coord_map()



## Australia emissions
(ggplot(Emissions %>% filter(lon > 127 & lat > -55.25 & lon < 171 & lat < -14.75)) + geom_tile(aes(lon,lat,fill=pmin(CH4flux_tot,2000))) + scale_fill_distiller(palette="Spectral")) %>% FRK:::draw_world() + xlim(c(-105,-61)) + ylim(c(17.5,47)) + coord_map()


Emissions_Australia <- filter(Emissions,lon > 112 & lat > -40.25 & lon < 156.5 & lat < -10.75)
gAustralia <- (ggplot() + geom_tile(data=Emissions_Australia,aes(lon,lat,fill=pmin(CH4flux_tot,200))) + scale_fill_distiller(palette="Spectral", name = "flux (g/s)\n")) %>% FRK:::draw_world()+ xlim(c(112,156.5)) + ylim(c(-40.25,-10.75)) + xlab("lon (deg)") + ylab("lat (deg)")
ggsave(gAustralia,filename = "./img/EDGAR_Australia.pdf",width = 10,height=7)

dfdsfds
save(Emissions_Midwest,file="inst/extdata/Emissions_Midwest.rda")
save(Emissions_Australia,file="inst/extdata/Emissions_Australia.rda")
save(Emissions_Canada,file="inst/extdata/Emissions_Canada.rda")


## Now divide into more boxes
### GRID INTO OUR ONE
trans_lon <- seq(-180.001, 180.001,by = 15)
trans_lat <- seq(-90.001, 90.001,by = 12)
Emissions$trans_lon <- cut(Emissions$lon,   
                    trans_lon,
                    labels=FALSE)
Emissions$trans_lat <- cut(Emissions$lat,   
                    trans_lat,
                    labels=FALSE)


# lambda_map <- plyr::ddply(Emissions,c("trans_lon","trans_lat"), function(df) {
#     df <- filter(df, CH4flux_tot > 0)
#     if(nrow(df) > 100) {
#         #df$CH4flux <- df$CH4flux/(1e-16)
#         geo_obj <- as.geodata(df,coords.col = c("lon","lat"),data.col = "CH4flux_tot")
#         model.fit <- tryCatch(likfit(geo_obj,ini=c(0.5,0.5),fix.lambda=FALSE,fix.nugget=T),
#                               error=function(e) list(lambda=-9999))
#         lambda <- model.fit$lambda
#     } else {
#         lambda <- NA
#     }
#     lambda
# })
load("./inst/extdata/lambda_map.rda")

lambda_map$lambda <- lambda_map$V1
lambda_map$lon <- trans_lon[lambda_map$trans_lon]
lambda_map$lat <- trans_lat[lambda_map$trans_lat]
(ggplot(lambda_map) + geom_tile(aes(lon,lat,fill=lambda))) %>% FRK::draw_world()

