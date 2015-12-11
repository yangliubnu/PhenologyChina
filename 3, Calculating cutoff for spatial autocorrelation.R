#### Code in this file uses initial occurrence records of bees as an example
#### The same procedures were conducted for every phenophase and every taxon/taxon-pair

library(nlme)
library(spdep)
library(ggplot2)
library(grid)
library(scales)

#### Setting directory

#### Setting up function for small-sample AIC (AICc)
AICc <- function(LL, n, k) 
{ -2*LL + 2*k*n/(n-k-1) } 

#### Using breeding season initial sighting data as an example (initial sighting is represented as "arrival" in data set)

### Reading in data file
arrival.breeding.date  = read.csv('BreedingArrival.csv', header=T)
arrival.breeding.date$SiteID = as.factor(arrival.breeding.date $SiteID)
arrival.breeding.date$Year = as.numeric(arrival.breeding.date $Year)
head(arrival.breeding.date )

arrival.breeding.date.5_8 <- arrival.breeding.date[arrival.breeding.date$m == 5 & arrival.breeding.date$n == 8,]

### Calculating how many sites there are for each AnalysisGroup
arrival.breeding.date.5_8_Sites <- data.frame(arrival.breeding.date.5_8$AnalysisGroup, arrival.breeding.date.5_8$SiteID, rep(1,nrow(arrival.breeding.date.5_8)))
names(arrival.breeding.date.5_8_Sites) <- c("AnalysisGroup", "SiteID", "Site")
arrival.breeding.date.5_8_Sites_NoDup <- arrival.breeding.date.5_8_Sites[!duplicated(arrival.breeding.date.5_8_Sites[,1:2]),]
aggregate(Site ~ AnalysisGroup, FUN="sum", data=arrival.breeding.date.5_8_Sites_NoDup)
## All AnalysisGroups had >=10 sites except for hoopoe and swift 
## Do not run autocorrelation tests for n<10. 

species = list('Bee','Cicada','Cricket','Cuckoo',
               'Frog','Goose','Hoopoe',
               'Swallow','Swift',
               'Toad')

### Analyses of the trend of initial occurrence of individual taxonomic units
arrival = arrival.breeding.date.5_8[arrival.breeding.date.5_8$AnalysisGroup==species[[1]],]

d1 <- data.frame(lon=arrival$Site_long, lat=arrival$Site_Lat)
d <- d1[!duplicated(d1),]    ## longlat of sites

## Running random slope + intercept models (allows for correlation between random slope and intercept)
m1 <- lme(scale(Arrival)~scale(Year), random = ~scale(Year) | SiteID, data = na.omit(arrival))
#m1$coefficients
#m1$coefficients$random$SiteID[,1]  ##level-2 resid for intercept random effect
#m1$coefficients$random$SiteID[,2]  ##level-2 resid for slope random effect
AICc(LL=logLik(m1), n=nrow(na.omit(arrival)), k=6)

## Setting up data frame to store Moran's I values for each 200km lag
Moran_Bee<-matrix(NA, nrow=10, ncol=4) ## data table
	
## 0-200km lag
dnn <- dnearneigh(coordinates(d), 0, 200, longlat=TRUE)
d.listw <- nb2listw(dnn, style="W", zero.policy=TRUE)

set.seed(190)# need to set seed so the program samples the same slope and intercept residuals
moran_int<-moran.mc(m1$coefficients$random$SiteID[,1], listw=d.listw, zero.policy=TRUE, nsim=999)	
set.seed(190)
moran_slope<-moran.mc(m1$coefficients$random$SiteID[,2], listw=d.listw, zero.policy=TRUE, nsim=999)
	
Moran_Bee[1,]<-c(moran_int$stat,
		 moran_int$p.value,
		 moran_slope$stat,
		 moran_slope$p.value)

## 200-400km lag
dnn <- dnearneigh(coordinates(d), 200, 400, longlat=TRUE)
d.listw <- nb2listw(dnn, style="W", zero.policy=TRUE)

set.seed(910)# need to set seed so the program samples the same slope and intercept residuals
moran_int<-moran.mc(m1$coefficients$random$SiteID[,1], listw=d.listw, zero.policy=TRUE, nsim=999)	
set.seed(910)
moran_slope<-moran.mc(m1$coefficients$random$SiteID[,2], listw=d.listw, zero.policy=TRUE, nsim=999)
	
Moran_Bee[2,]<-c(moran_int$stat,
		 moran_int$p.value,
		 moran_slope$stat,
		 moran_slope$p.value)

## 400-600km lag
dnn <- dnearneigh(coordinates(d), 400, 600, longlat=TRUE)
d.listw <- nb2listw(dnn, style="W", zero.policy=TRUE)

set.seed(565)# need to set seed so the program samples the same slope and intercept residuals
moran_int<-moran.mc(m1$coefficients$random$SiteID[,1], listw=d.listw, zero.policy=TRUE, nsim=999)	
set.seed(565)
moran_slope<-moran.mc(m1$coefficients$random$SiteID[,2], listw=d.listw, zero.policy=TRUE, nsim=999)
	
Moran_Bee[3,]<-c(moran_int$stat,
	 	 moran_int$p.value,
		 moran_slope$stat,
		 moran_slope$p.value)

## 600-800km lag
dnn <- dnearneigh(coordinates(d), 600, 800, longlat=TRUE)
d.listw <- nb2listw(dnn, style="W", zero.policy=TRUE)

set.seed(566)# need to set seed so the program samples the same slope and intercept residuals (account for random effects covariance structure)
moran_int<-moran.mc(m1$coefficients$random$SiteID[,1], listw=d.listw, zero.policy=TRUE, nsim=999)	
set.seed(566)
moran_slope<-moran.mc(m1$coefficients$random$SiteID[,2], listw=d.listw, zero.policy=TRUE, nsim=999)
	
Moran_Bee[4,]<-c(moran_int$stat,
		 moran_int$p.value,
		 moran_slope$stat,
		 moran_slope$p.value)

## 800-1000km lag
dnn <- dnearneigh(coordinates(d), 800, 1000, longlat=TRUE)
d.listw <- nb2listw(dnn, style="W", zero.policy=TRUE)

set.seed(677)# need to set seed so the program samples the same slope and intercept residuals
moran_int<-moran.mc(m1$coefficients$random$SiteID[,1], listw=d.listw, zero.policy=TRUE, nsim=999)	
set.seed(677)
moran_slope<-moran.mc(m1$coefficients$random$SiteID[,2], listw=d.listw, zero.policy=TRUE, nsim=999)
	
Moran_Bee[5,]<-c(moran_int$stat,
		 moran_int$p.value,
		 moran_slope$stat,
		 moran_slope$p.value)
				
## 1000-1200km lag
dnn <- dnearneigh(coordinates(d), 1000, 1200, longlat=TRUE)
d.listw <- nb2listw(dnn, style="W", zero.policy=TRUE)

set.seed(776)# need to set seed so the program samples the same slope and intercept residuals
moran_int<-moran.mc(m1$coefficients$random$SiteID[,1], listw=d.listw, zero.policy=TRUE, nsim=999)	
set.seed(776)
moran_slope<-moran.mc(m1$coefficients$random$SiteID[,2], listw=d.listw, zero.policy=TRUE, nsim=999)
	
Moran_Bee[6,]<-c(moran_int$stat,
		 moran_int$p.value,
		 moran_slope$stat,
		 moran_slope$p.value)

## 1200-1400km lag
dnn <- dnearneigh(coordinates(d), 1200, 1400, longlat=TRUE)
d.listw <- nb2listw(dnn, style="W", zero.policy=TRUE)

set.seed(145)# need to set seed so the program samples the same slope and intercept residuals
moran_int<-moran.mc(m1$coefficients$random$SiteID[,1], listw=d.listw, zero.policy=TRUE, nsim=999)	
set.seed(145)
moran_slope<-moran.mc(m1$coefficients$random$SiteID[,2], listw=d.listw, zero.policy=TRUE, nsim=999)
	
Moran_Bee[7,]<-c(moran_int$stat,
		 moran_int$p.value,
		 moran_slope$stat,
		 moran_slope$p.value)
				
## 1400-1600km lag
dnn <- dnearneigh(coordinates(d), 1400, 1600, longlat=TRUE)
d.listw <- nb2listw(dnn, style="W", zero.policy=TRUE)

set.seed(148)# need to set seed so the program samples the same slope and intercept residuals
moran_int<-moran.mc(m1$coefficients$random$SiteID[,1], listw=d.listw, zero.policy=TRUE, nsim=999)	
set.seed(148)
moran_slope<-moran.mc(m1$coefficients$random$SiteID[,2], listw=d.listw, zero.policy=TRUE, nsim=999)
	
Moran_Bee[8,]<-c(moran_int$stat,
		 moran_int$p.value,
		 moran_slope$stat,
		 moran_slope$p.value)		
			
## 1600-1800km lag
dnn <- dnearneigh(coordinates(d), 1600, 1800, longlat=TRUE)
d.listw <- nb2listw(dnn, style="W", zero.policy=TRUE)

set.seed(569)# need to set seed so the program samples the same slope and intercept residuals
moran_int<-moran.mc(m1$coefficients$random$SiteID[,1], listw=d.listw, zero.policy=TRUE, nsim=999)	
set.seed(569)
moran_slope<-moran.mc(m1$coefficients$random$SiteID[,2], listw=d.listw, zero.policy=TRUE, nsim=999)
	
Moran_Bee[9,]<-c(moran_int$stat,
		 moran_int$p.value,
		 moran_slope$stat,
		 moran_slope$p.value)
		
## 1800-2000km lag
dnn <- dnearneigh(coordinates(d), 1800, 2000, longlat=TRUE)
d.listw <- nb2listw(dnn, style="W", zero.policy=TRUE)

set.seed(831)# need to set seed so the program samples the same slope and intercept residuals
moran_int<-moran.mc(m1$coefficients$random$SiteID[,1], listw=d.listw, zero.policy=TRUE, nsim=999)	
set.seed(831)
moran_slope<-moran.mc(m1$coefficients$random$SiteID[,2], listw=d.listw, zero.policy=TRUE, nsim=999)
	
Moran_Bee[10,]<-c(moran_int$stat,
		  moran_int$p.value,
		  moran_slope$stat,
		  moran_slope$p.value)
Moran_Bee <- data.frame(Moran_Bee); names(Moran_Bee) <- c("Int_I", "Int_I_pval", "Slope_I", "Slope_I_pval")
Moran_Bee$Lags <- c(1,2,3,4,5,6,7,8,9,10) ## 1=1-200 km, 2=201-400 km etc...
#write.csv(Moran_Bee, "Moran_Bee_Arrival.csv")
