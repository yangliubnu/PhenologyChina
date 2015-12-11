#### Linear mixed models for analysis pooling time series across the country

library(nlme)

#### Analysing temporal trend in phenophase

### Initial occurrence during breeding season (initial occurrence represented as "arrival" in data set)
## Reading in data file
site.data = read.csv('arrival.csv', header=T)
site.data$SiteID = as.factor(site.data$SiteID)
site.data$Year = as.numeric(site.data$Year)
head(site.data)

criteria.all=c(5,6,10,11)

## Setting up list to store results
arrival.trend = list(rep(NA,4))
intercept = data.frame(matrix(NA,ncol=9,nrow=1))
names(intercept)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")
slope = data.frame(matrix(NA,ncol=9,nrow=1))
names(slope)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")

## Looping through analysis
for (l in 1:4){
  print(l)
  m=criteria.all[l]
  site.data.sub = site.data[site.data$m==m,]
  species.all = as.character(unique(site.data.sub$AnalysisGroup))
  species.number = length(species.all)
  arrival.trend[[l]]=list(rep(NA,species.number))
  
  for (i in 1:species.number){
    print(i)
    arrival.trend[[l]][[i]]=NA
    arrival = site.data.sub[site.data.sub$AnalysisGroup==species.all[[i]],]
    
    s.y=sd(arrival$Arrival)
    s.x=sd(arrival$Year)
    mean.y = mean(arrival$Arrival)
    mean.x = mean(arrival$Year)
    
    intercept.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(intercept.temporary)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")
    intercept.temporary[1,1]=criteria.all[l]
    intercept.temporary[1,2]=species.all[[i]]
    
    slope.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(slope.temporary)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")
    slope.temporary[1,1]=l
    slope.temporary[1,2]=species.all[[i]]
    
    tryCatch({
      #trend = lme(scale(Arrival)~scale(Year), random = ~scale(Year) | SiteID, data = na.omit(arrival))
      trend = lme(scale(Arrival)~scale(Year), random = ~1 | SiteID, data = na.omit(arrival))
      arrival.trend[[l]][[i]] = summary(trend)
      
      intercept.temporary[1,3:7]=summary(trend)[[20]][1,]
      intercept.temporary[1,8]=summary(trend)[[20]][1,1]*s.y + mean.y - summary(trend)[[20]][2,1]*s.y/s.x*mean.x
      intercept.temporary[1,9]=summary(trend)[[20]][1,2]*s.y + mean.y - summary(trend)[[20]][2,2]*s.y/s.x*mean.x
      
      slope.temporary[1,3:7]=summary(trend)[[20]][2,]
      slope.temporary[1,8]=summary(trend)[[20]][2,1]*s.y/s.x
      slope.temporary[1,9]=summary(trend)[[20]][2,2]*s.y/s.x
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})      
    
    intercept = rbind(intercept, intercept.temporary)
    slope = rbind(slope, slope.temporary)
  }
}
arrival.trend
write.csv(intercept,"intercept.csv")
write.csv(slope,"slope.csv")

### Disappearance during breeding season (disappearance represented as "departure" in data set)
## Reading in data file
site.data = read.csv('departure.csv', header=T)
site.data$SiteID = as.factor(site.data$SiteID)
site.data$Year = as.numeric(site.data$Year)
head(site.data)

criteria.all=c(5,6,10,11)

## Setting up list to store results
departure.trend = list(rep(NA,4))
intercept = data.frame(matrix(NA,ncol=9,nrow=1))
names(intercept)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")
slope = data.frame(matrix(NA,ncol=9,nrow=1))
names(slope)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")

## Looping through analysis
for (l in 1:4){
  print(l)
  m=criteria.all[l]
  site.data.sub = site.data[site.data$m==m,]
  species.all = as.character(unique(site.data.sub$AnalysisGroup))
  species.number = length(species.all)
  departure.trend[[l]]=list(rep(NA,species.number))
  
  for (i in 1:species.number){
    print(i)
    departure.trend[[l]][[i]]=NA
    departure = site.data.sub[site.data.sub$AnalysisGroup==species.all[[i]],]
    
    s.y=sd(departure$Departure)
    s.x=sd(departure$Year)
    mean.y = mean(departure$Departure)
    mean.x = mean(departure$Year)
    
    intercept.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(intercept.temporary)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")
    intercept.temporary[1,1]=criteria.all[l]
    intercept.temporary[1,2]=species.all[[i]]
    
    slope.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(slope.temporary)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")
    slope.temporary[1,1]=criteria.all[l]
    slope.temporary[1,2]=species.all[[i]]
    
    tryCatch({
      #trend = lme(scale(Departure)~scale(Year), random = ~scale(Year) | SiteID, data = na.omit(departure))
      trend = lme(scale(Departure)~scale(Year), random = ~1 | SiteID, data = na.omit(departure))
      departure.trend[[l]][[i]] = summary(trend)
      
      intercept.temporary[1,3:7]=summary(trend)[[20]][1,]
      intercept.temporary[1,8]=summary(trend)[[20]][1,1]*s.y + mean.y - summary(trend)[[20]][2,1]*s.y/s.x*mean.x
      intercept.temporary[1,9]=summary(trend)[[20]][1,2]*s.y + mean.y - summary(trend)[[20]][2,2]*s.y/s.x*mean.x
      
      slope.temporary[1,3:7]=summary(trend)[[20]][2,]
      slope.temporary[1,8]=summary(trend)[[20]][2,1]*s.y/s.x
      slope.temporary[1,9]=summary(trend)[[20]][2,2]*s.y/s.x
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})      
    
    intercept = rbind(intercept, intercept.temporary)
    slope = rbind(slope, slope.temporary)
  }
}
departure.trend
write.csv(intercept,"intercept.csv")
write.csv(slope,"slope.csv")

### Length of temporal occurrence span during breeding season (length of temporal occurrence span represented as "span" in data set)
## Reading in data file
site.data = read.csv('stay.csv', header=T)
site.data$SiteID = as.factor(site.data$SiteID)
site.data$Year = as.numeric(site.data$Year)
head(site.data)

criteria.all=c(5,6,10,11)

## Setting up list to store results
stay.trend = list(rep(NA,4))
intercept = data.frame(matrix(NA,ncol=9,nrow=1))
names(intercept)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")
slope = data.frame(matrix(NA,ncol=9,nrow=1))
names(slope)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")

## Looping through analysis
for (l in 1:4){
  print(l)
  m=criteria.all[l]
  site.data.sub = site.data[site.data$m==m,]
  species.all = as.character(unique(site.data.sub$AnalysisGroup))
  species.number = length(species.all)
  stay.trend[[l]]=list(rep(NA,species.number))
  
  for (i in 1:species.number){
    print(i)
    stay.trend[[l]][[i]]=NA
    stay = site.data.sub[site.data.sub$AnalysisGroup==species.all[[i]],]
    
    s.y=sd(stay$Stay)
    s.x=sd(stay$Year)
    mean.y = mean(stay$Stay)
    mean.x = mean(stay$Year)
    
    intercept.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(intercept.temporary)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")
    intercept.temporary[1,1]=criteria.all[l]
    intercept.temporary[1,2]=species.all[[i]]
    
    slope.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(slope.temporary)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")
    slope.temporary[1,1]=criteria.all[l]
    slope.temporary[1,2]=species.all[[i]]
    
    tryCatch({
      #trend = lme(scale(Stay)~scale(Year), random = ~scale(Year) | SiteID, data = na.omit(stay))
      trend = lme(scale(Stay)~scale(Year), random = ~1 | SiteID, data = na.omit(stay))
      stay.trend[[l]][[i]] = summary(trend)
      
      intercept.temporary[1,3:7]=summary(trend)[[20]][1,]
      intercept.temporary[1,8]=summary(trend)[[20]][1,1]*s.y + mean.y - summary(trend)[[20]][2,1]*s.y/s.x*mean.x
      intercept.temporary[1,9]=summary(trend)[[20]][1,2]*s.y + mean.y - summary(trend)[[20]][2,2]*s.y/s.x*mean.x
      
      slope.temporary[1,3:7]=summary(trend)[[20]][2,]
      slope.temporary[1,8]=summary(trend)[[20]][2,1]*s.y/s.x
      slope.temporary[1,9]=summary(trend)[[20]][2,2]*s.y/s.x
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})      
    
    intercept = rbind(intercept, intercept.temporary)
    slope = rbind(slope, slope.temporary)
  }
}
stay.trend
write.csv(intercept,"intercept.csv")
write.csv(slope,"slope.csv")

### Temporal overlap between taxa
## Reading in data file
site.data = read.csv('CommunityOverlap.csv', header=T)
site.data$SiteID = as.factor(site.data$SiteID)
site.data$Year = as.numeric(site.data$Year)
head(site.data)

criteria.all=c(5,6,10,11)

## Setting up list to store results
overlap.trend = list(rep(NA,4))
intercept = data.frame(matrix(NA,ncol=9,nrow=1))
names(intercept)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")
slope = data.frame(matrix(NA,ncol=9,nrow=1))
names(slope)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")

## Looping through analysis
for (l in 1:4){
  print(l)
  m=criteria.all[l]
  site.data.sub = site.data[site.data$m==m,]
  species.all = as.character(unique(site.data.sub$TaxonPair))
  species.number = length(species.all)
  overlap.trend[[l]]=list(rep(NA,species.number))
  
  for (i in 1:species.number){
    print(i)
    overlap.trend[[l]][[i]]=NA
    overlap = site.data.sub[site.data.sub$TaxonPair==species.all[[i]],]
    
    s.y=sd(overlap$Overlap)
    s.x=sd(overlap$Year)
    mean.y = mean(overlap$Overlap)
    mean.x = mean(overlap$Year)
    
    intercept.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(intercept.temporary)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")
    intercept.temporary[1,1]=criteria.all[l]
    intercept.temporary[1,2]=species.all[[i]]
    
    slope.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(slope.temporary)=c("Criterion","Taxon","beta.std","SE.std","DF","t-value","p-value","beta","SE")
    slope.temporary[1,1]=criteria.all[l]
    slope.temporary[1,2]=species.all[[i]]
    
    tryCatch({
      #trend = lme(scale(Overlap)~scale(Year), random = ~scale(Year) | SiteID, data = na.omit(overlap))
      trend = lme(scale(Overlap)~scale(Year), random = ~1 | SiteID, data = na.omit(overlap))
      overlap.trend[[l]][[i]] = summary(trend)
      
      intercept.temporary[1,3:7]=summary(trend)[[20]][1,]
      intercept.temporary[1,8]=summary(trend)[[20]][1,1]*s.y + mean.y - summary(trend)[[20]][2,1]*s.y/s.x*mean.x
      intercept.temporary[1,9]=summary(trend)[[20]][1,2]*s.y + mean.y - summary(trend)[[20]][2,2]*s.y/s.x*mean.x
      
      slope.temporary[1,3:7]=summary(trend)[[20]][2,]
      slope.temporary[1,8]=summary(trend)[[20]][2,1]*s.y/s.x
      slope.temporary[1,9]=summary(trend)[[20]][2,2]*s.y/s.x
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})      
    
    intercept = rbind(intercept, intercept.temporary)
    slope = rbind(slope, slope.temporary)
  }
}
overlap.trend
write.csv(intercept,"intercept.csv")
write.csv(slope,"slope.csv")
