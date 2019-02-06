## ----- import data
library(readr)  
library(readxl)
library(raster)
library(rgdal)
library(plyr)
library(reshape)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(textclean)
library(ade4)
library(sjPlot)

##Paris human population (year 2015)
pop2015paris <- read.csv("data/env_data/pop2015paris.csv")
head(pop2015paris)  
##Paris human flow (based on entrances in metro station)
humanflow <- read.csv("data/env_data/humanflow_metro_paris2013.csv", sep=";", h = TRUE)
head(humanflow)
##Paris commerces activities
commercesparis <- read_excel("data/env_data/commercesparis.xlsx")
commercesparis <- as.data.frame(commercesparis)
head(commercesparis)
##Paris food markets locations, and opening hours' count
parismarkets.shp <- readOGR(dsn="data/env_data",layer="liste_des_marches_de_quartier_a_paris")
hrscount_markets<- read.csv("data/env_data/localmarketsparis_hourscount.csv", sep=",", h = TRUE)
head(hrscount_markets)
##Air pollution
pollutNO2paris<- read.csv("data/env_data/airpollutionNO2paris.csv",  sep=",", h = TRUE)
head(pollutNO2paris)
pollutPM10paris<- read.csv("data/env_data/airpollutionPM10paris.csv", sep=",", h = TRUE)
head(pollutPM10paris)
indices_pollut <- read.csv("data/env_data/indices_QA_commune_IDF_2014.csv", sep=",", h = TRUE)
indices_pollut_paris <- subset(indices_pollut,ninsee %in% c(75101:75120))
# pollutO3paris<- read.csv("data/env_data/airpollutionO3paris.csv", sep=";", h = TRUE)
# head(pollutO3paris)
triris_stationair<- read.csv("data/env_data/triris_stationairpollution.csv", sep=",", h = TRUE)
head(triris_stationair)
triris_stationair$TRIRIS <- as.factor(formatC(as.numeric(triris_stationair$TRIRIS), width=3, flag='0', format="d")) #format with 'total field width'' = 3 for all codes, and using "0" pads leading zeros
##Paris waste 
recyclwaste <- read.csv("data/env_data/waste_recycl_paris2011.csv",  sep=",", h = TRUE)
head(recyclwaste)
nonrecyclwaste <- read.csv("data/env_data/waste_nonrecycl_paris2011.csv",  sep=",", h = TRUE)
head(nonrecyclwaste)

# marketsparis <- read_excel("data/env_data/localmarketsparis.xlsx")
# View(marketsparis)
##IRIS contours
# dpt_france.shp <- readOGR(dsn="data/departements_france_2014",layer="departements-20140306-100m")
# iris.shp <- readOGR(dsn="data/gis_iris_lamb93",layer="iris-2013-01-01")
## select IRIS within Paris (<=> departement 75)
# paris.shp <- subset(dpt_france.shp, code_insee == 75)
# iris_paris.shp <- intersect(iris.shp,paris.shp)
# writeOGR(iris_paris.shp,dsn="data/gis_iris_lamb93",layer="iris_paris",driver="ESRI Shapefile")
# plot(iris_paris.shp)
iris_paris.shp <- readOGR(dsn="data/gis_iris_lamb93",layer="iris_paris")
##metro stations locations
loc_ratp.shp <- readOGR(dsn="data/gis_ratp_paris",layer="positions-geographiques-des-stations-du-reseau-ratp")
plot(loc_ratp.shp)
stationsmetro_triris <- read.csv("data/env_data/stations_metro_triris.csv",h=T)


## ----- overlay metro stations and IRIS contours
ov <- over(loc_ratp.shp,iris_paris.shp)
loc_ratp.shp@data$IRIS <- ov$IRIS
loc_ratp.df <- loc_ratp.shp@data

## ----- subset only food markets and overlay with IRIS
foodmarkets.shp <- subset(parismarkets.shp, type=="Alimentaire découvert" | type=="Alimentaire couvert" | type=="Biologique")
ov_markets_iris <- over(foodmarkets.shp, iris_paris.shp)
foodmarkets.shp@data$IRIS <- ov_markets_iris$IRIS
foodmarkets.shp@data <- foodmarkets.shp@data[,!names(foodmarkets.shp@data) %in% c("lundi","mardi","mercredi","jeudi","vendredi","samedi","dimanche","jour_ferie")]
foodmarkets.shp@data <- merge(foodmarkets.shp@data,hrscount_markets,by.x="marche",by.y="Marche")


## ----- summarize data by TRIRIS
## need first to summarize by IRIS for some variables before summarizing by TRIRIS

## create a table with IRIS-TRIRIS-DISTRICT correspondance
iris_triris <- pop2015paris %>%
  group_by(IRIS,TRIRIS) %>%
  select(IRIS, TRIRIS)
iris_district <- commercesparis %>%
  group_by(IRIS,ARRONDISSEMENT) %>%
  select(IRIS, ARRONDISSEMENT)
iris_triris_district <- merge(iris_triris,iris_district,by="IRIS")
iris_triris_district <- iris_triris_district[!duplicated(iris_triris_district$IRIS),]
iris_triris_district$iris_simplif <- substr(iris_triris_district$IRIS, start = 6, stop = 9)
iris_triris_district$sect_triris <- substr(iris_triris_district$TRIRIS, start = 3, stop = 5)
iris_triris_district <- plyr::rename(iris_triris_district, c("ARRONDISSEMENT" = "district", "IRIS" = "sect_iris"))
#table with only district-triris corresp.
district_triris <- iris_triris_district %>%
  select(sect_triris,district,TRIRIS) %>%
  distinct(sect_triris,.keep_all=T)

detach("package:plyr", unload=TRUE)

## human pop / TRIRIS
humpop_triris <- pop2015paris %>%
  select(TRIRIS, P15_POP) %>%
  group_by(TRIRIS) %>%
  summarize(humpop_2015 = sum(P15_POP))
humpop_triris <- merge(humpop_triris,district_triris,by="TRIRIS")
humpop_triris <- humpop_triris %>%
  select(-TRIRIS)
humpop_triris <- as.data.frame(humpop_triris)


## calculate human flow / TRIRIS
humanflow <- humanflow %>%
  select(Station, Trafic,Arrondissement.pour.Paris)
#merge with triris
humanflow_triris <- merge(humanflow,stationsmetro_triris,by.x="Station",by.y="STATION")
#summarize by triris 
humanflow_triris <- humanflow_triris %>%
  group_by(TRIRIS) %>% 
  summarise(nb_passers = sum(Trafic))
humanflow_triris <- as.data.frame(humanflow_triris)
humanflow_triris$TRIRIS <- as.factor(formatC(as.numeric(humanflow_triris$TRIRIS), width=3, flag='0', format="d")) #format with 'total field width'' = 3 for all codes, and using "0" pads leading zeros
humanflow_triris <- merge(district_triris,humanflow_triris,by.x="sect_triris",by.y="TRIRIS",all.x=T)
humanflow_triris <- humanflow_triris %>%
  select(-c(district,TRIRIS))

## calculate number of commerces (bakeries, haircut, fabric manufacture)
commercesparis <- commercesparis %>%
  filter(`LIBELLE ACTIVITE` %in% c("Boulangerie - Boulangerie PÃ¢tisserie","PÃ¢tisserie","Coiffure","Tissus - Textile - Mercerie","Commerce de gros fabrication textile"))
commercesparis <- merge(commercesparis,iris_triris_district,by.x="IRIS",by.y="sect_iris", all.x=T)
commerces_triris <- commercesparis %>%
  mutate(activity = recode(`LIBELLE ACTIVITE`,"Boulangerie - Boulangerie PÃ¢tisserie"="bakeries","PÃ¢tisserie"="bakeries","Coiffure"="haircut","Tissus - Textile - Mercerie"="fabric_manuf","Commerce de gros fabrication textile"="fabric_manuf")) %>%
  select(-`LIBELLE ACTIVITE`)  %>%
  group_by(sect_triris, activity) %>%
  summarise(nb_commerces = n()) 
nbcommerces_triris <- cast(commerces_triris, sect_triris ~ activity)
nbcommerces_triris[is.na(nbcommerces_triris)] <- 0

## calculate number of food markets / TRIRIS
markets_iris <- merge(foodmarkets.shp@data,iris_triris_district,by.x="IRIS",by.y="iris_simplif", all.x=T)
markets_triris <- markets_iris %>%
  group_by(sect_triris,type) %>%
  summarise(nb_markets = n())
markets_triris_count <- cast(markets_triris,sect_triris~type, sum)
colnames(markets_triris_count) <- c("sect_triris","nb_closedmarkets","nb_openmarkets","nb_biomarkets")
markets_triris_count$nb_foodmarkets <- rowSums(markets_triris_count[,c(2:4)])
markets_triris_count$openmarkets <- ifelse(markets_triris_count$nb_openmarkets>0,1,0)
markets_triris_count$foodmarkets <- ifelse(markets_triris_count$nb_foodmarkets>0,1,0)
hrs_markets_triris <- markets_iris %>%
  group_by(sect_triris,type) %>%
  summarise(nb_hours = sum(TOTAL))
hrs_markets_triris <- cast(hrs_markets_triris,sect_triris~type, sum)
colnames(hrs_markets_triris) <- c("sect_triris","hrs_closedmarkets","hrs_openmarkets","hrs_biomarkets")
hrs_markets_triris$tot_openhours <- rowSums(hrs_markets_triris[,c(2:4)])
data_markets_triris <- merge(markets_triris_count,hrs_markets_triris,by="sect_triris")

## air pollution / TRIRIS

#merge station-triris with triris-district to get station-district
triris_stationair <- plyr::rename(triris_stationair,c("TRIRIS"="sect_triris"))
stationair_triris_district <- merge(triris_stationair,district_triris,by="sect_triris",all.x=T)
stationair_district <- stationair_triris_district %>%
  select(STATION,district) %>%
  distinct()
#indices air quality
#calculate the mean for each poluutant / district(ninsee)
mean_pollut <- indices_pollut_paris %>% 
  select(-date) %>%
  group_by(ninsee) %>%
  summarise_all(mean)
district <- mean_pollut$ninsee %>%
  str_replace_all("751","750")
#run pca to reduce to a single principal components
library(FactoMineR)
library(factoextra)
meanpollut.pca <- PCA(mean_pollut[,-1], graph=F)
summary(meanpollut.pca)
# plot PCA
pdf(file="figs/PCA.pdf",width=8.5,height=11)
plot.PCA(meanpollut.pca,choix="var")
dev.off()
# sqrt of eigenvalues
pca.eigen <- meanpollut.pca$eig
write.csv(pca.eigen,"outputs/airqual_pca_eigenvalues.csv")
# correlations between variables and PCs
pca.coord <- meanpollut.pca$var$coord
write.csv(pca.coord,"outputs/airqual_pca_coord.csv")
# PCs (aka scores)
pca.scores <- meanpollut.pca$ind$coord
write.csv(pca.scores,"outputs/airqual_pca_scores.csv")
airqual_district <- data.frame(district=district,airqual_pca1 = pca.scores[,1])
#merge with triris
airqual_triris <- merge(district_triris,airqual_district,by="district")
airqual_triris <- airqual_triris %>%
  select(-c(district,TRIRIS))

# #NO2 (raw measures) !!!! there is something wrong with those values - probably an issue with the allocation of the closest station to each triris!!!
# max_pollutNO2paris <- pollutNO2paris %>%
#   group_by(date) %>%
#   select(-heure)  %>%
#   summarise_all(max)
# max_pollutNO2paris <- as.data.frame(max_pollutNO2paris)
# #reshape
# max_pollutNO2paris <- melt(max_pollutNO2paris,id="date")
# colnames(max_pollutNO2paris) <-c("date","station","maxNO2")
# #ad district
# pollutNO2paris_district <- merge(max_pollutNO2paris, stationair_district, by.x = "station", by.y="STATION")
# pollutNO2paris_district$maxNO2 <- as.numeric(as.character(pollutNO2paris_district$maxNO2))
# meanpollutNO2paris_district <- pollutNO2paris_district %>%
#   na.omit() %>%
#   select(district,maxNO2) %>%
#   group_by(district) %>%
#   summarise(NO2 = mean(maxNO2))


## waste
library(tidyverse)
recyclwaste <- recyclwaste %>%
  select(Granularité,TOT_2011)
recyclwaste <-  plyr::rename(recyclwaste, c("Granularité"= "district","TOT_2011"="recyclwaste"))
nonrecyclwaste <- nonrecyclwaste %>%
  select(Granularité,TOT_2011)
nonrecyclwaste <-  plyr::rename(nonrecyclwaste, c("Granularité"= "district","TOT_2011"="nonrecyclwaste"))
waste_triris <- Reduce(function(x,y) merge(x = x, y = y, by = "district",all.x=T), 
       list(district_triris, recyclwaste, nonrecyclwaste))
waste_triris$waste <- rowSums(waste_triris[,4:5])
waste_triris <- waste_triris %>%
  select(-c(district,TRIRIS))

## group environmental data / TRIRIS
envdata_triris <- Reduce(function(x,y) merge(x = x, y = y, by = "sect_triris",all.x=T), 
                       list(humpop_triris, humanflow_triris,nbcommerces_triris,data_markets_triris,airqual_triris, waste_triris))
# envdata_triris <- plyr::rename(envdata_triris,c("district.x"="district"))
write.csv(envdata_triris,"data/env_data/allenvdata_triris.csv")

## combine with pigeon data
data2018 <- read.csv("data/data2018.csv", h=T,sep=",")
head(data2018)
colnames(data2018) <- c("district","sect_triris","morph","colour","colour_simplif","extleft_fing","medleft_fing","intleft_fing","backleft_fing","extright_fing","medright_fing","intright_fing","backright_fing","left_tot","right_tot","total","mutilation","age", "weather","nbindiv_checked",	"nbindiv_mutil","abund_sect","area","humpop_2008","humpop_dens","nb_passers","flow_passers", "openmarkets","bakeries", "bakeries_dens","greenspace_area", "greenspace_dens", "dovecote","noise_pollut","PM10","GHG","waste","waste_dens")
data2018 <- data2018 %>%
  select(c("district","sect_triris","morph","colour","colour_simplif","extleft_fing","medleft_fing","intleft_fing","backleft_fing","extright_fing","medright_fing","intright_fing","backright_fing","left_tot","right_tot","total","mutilation","age", "weather","nbindiv_checked",	"nbindiv_mutil","abund_sect","area","greenspace_area", "greenspace_dens","dovecote","noise_pollut"))
data2018$sect_triris <- as.factor(formatC(as.numeric(data2018$sect_triris), width=3, flag='0', format="d")) 
data2018_update <- merge(data2018,envdata_triris,by="sect_triris",all.x=T)
summary(data2018_update)
head(data2018_update)
data2018_update <- replace_na(data2018_update,list("nb_closedmarkets"=0,"nb_openmarkets"=0,"nb_biomarkets"=0,"nb_foodmarkets"=0,"openmarkets" =0,"foodmarkets" =0, "hrs_closedmarkets"=0, "hrs_openmarkets"=0,   "hrs_biomarkets"=0,  "tot_openhours"=0))
data2018_update <- data2018_update %>%
  select(-district.y)
data2018_update$humpop_dens <- data2018_update$humpop_2015/data2018_update$area
data2018_update$human_flow <- data2018_update$nb_passers/data2018_update$area
data2018_update$openmarkets_dens <- data2018_update$nb_openmarkets/data2018_update$area
data2018_update$foodmarkets_dens <- data2018_update$nb_foodmarkets/data2018_update$area
data2018_update$hrsarea_openmarkets <- data2018_update$hrs_openmarkets/data2018_update$area
data2018_update$hrsarea_totopenhrs <- data2018_update$tot_openhours/data2018_update$area
data2018_update$bakeries_dens <- data2018_update$bakeries/data2018_update$area
data2018_update$fabricmanuf_dens <- data2018_update$fabric_manuf/data2018_update$area
data2018_update$hairdres_dens <- data2018_update$haircut/data2018_update$area
data2018_update$waste_dens <- data2018_update$waste/data2018_update$area
head(data2018_update)
write.csv(data2018_update,"data/data2018_update.csv")


