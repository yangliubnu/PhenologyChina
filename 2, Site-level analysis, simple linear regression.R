#### All analyses in this code file are on individual time series

library(nlme)

#### set directory

#### Taxon-level analysis

### Backtransformation (because data has previously been centered and scaled) is only conducted for simple linear regression
### i.e., values fed into quadratic regression have not been backtransformed
### All conclusions about phenological responses are based on simple linear regression

### Temporal trend of phenophase: breeding taxa
##  Read in data files
arrival.breeding.date = read.table('breeding arrival.txt', header=T)
head(arrival.breeding.date)
departure.breeding.date = read.table('breeding departure.txt', header=T)
head(departure.breeding.date)
stay.breeding.date = read.table('breeding stay.txt', header=T)
head(stay.breeding.date)

species = list('Bee','Cicada','Cricket','Cuckoo',
               'Frog','Goose','Hoopoe',
               'Swallow','Swift',
               'Toad')

## define different data sets (of different trimming criteria)
N_Year = c(5,6,11,10)#criterion on the minimum number of time steps in a time series
Span_Year = c(8,10,11,14)#criterion on the minimum span of time series

## Analysis of initial occurrence data (initial occurrence represented as "arrival" in data set)
arrival.trend = list(rep(NA,4))
write.table(NA, file = "trend.csv", sep = " ")
for (a in 1:4){
  M=N_Year[a]
  N=Span_Year[a]
  data = arrival.breeding.date[arrival.breeding.date$m==M,]
  arrival.trend[[a]]=list(rep(NA,10))
  for (i in 1:10){
    print(i)
    arrival = data[data$AnalysisGroup==species[[i]],]  
    sites = as.character(unique(arrival$SiteID))
    site.number = length(sites)
    arrival.trend[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=23))
    for (j in 1:site.number){
      print(j)
      arrival.subset = na.omit(arrival[arrival$SiteID==sites[j],])
      if (nrow(arrival.subset)!=0){
        trend = lm(scale(Arrival)~scale(Year),
                   data = arrival.subset)
        
        s.y=sd(arrival.subset$Arrival)
        s.x=sd(arrival.subset$Year)
        
        residual = resid(trend)
        predicted = predict(trend)
        test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
        test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
        test= min(test.2)# the smallest p-value from test.2
        arrival.trend[[a]][[i]][j,1]=species[[i]]
        arrival.trend[[a]][[i]][j,2]=sites[j]
        arrival.trend[[a]][[i]][j,19]=M
        arrival.trend[[a]][[i]][j,20]=N
        arrival.trend[[a]][[i]][j,22]=summary(trend)[[4]][2]*s.y/s.x
        arrival.trend[[a]][[i]][j,23]=summary(trend)[[4]][4]*s.y/s.x
        
        if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
          arrival.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,4,6,8)] # summary of linear regression outputs
          arrival.trend[[a]][[i]][j,11]=summary(trend)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          arrival.trend[[a]][[i]][j,13]=summary(trend)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          arrival.trend[[a]][[i]][j,15]="linear"
          arrival.trend[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
        } else{#violation of the linear assmption: use quadratic form of the equation instead
          trend = lm(scale(Arrival)~scale(Year)+I(scale(Year)^2),
                     data = arrival.subset)
          residual = resid(trend)
          predicted = predict(trend)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          arrival.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
          arrival.trend[[a]][[i]][j,c(7:10)]=summary(trend)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
          arrival.trend[[a]][[i]][j,11]=summary(trend)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          arrival.trend[[a]][[i]][j,13]=summary(trend)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          arrival.trend[[a]][[i]][j,12]=summary(trend)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          arrival.trend[[a]][[i]][j,14]=summary(trend)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          arrival.trend[[a]][[i]][j,15]="quadratic"
          arrival.trend[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
          if (min(predicted)==head(predicted,n=1)){
            if(max(predicted)==tail(predicted,n=1)){
              arrival.trend[[a]][[i]][j,21]="mono_pos" # to check whether within the range of x variable, y variable shows a monotonic positive trend
            }
          }
          if (min(predicted)==tail(predicted,n=1)){
            if(max(predicted)==head(predicted,n=1)){
              arrival.trend[[a]][[i]][j,21]="mono_neg" # to check whether within the range of x variable, y variable shows a monotonic positive trend
            }
          }      
        }
      }
    }
    write.table(arrival.trend[[a]][[i]], file = "trend.csv", sep = " ",append=TRUE)
  }
}

## Analysis of disappearance data (disappearance represented as "departure" in data set)
departure.trend = list(rep(NA,4))
write.table(NA, file = "trend.csv", sep = " ")
for (a in 1:4){
  M=N_Year[a]
  N=Span_Year[a]
  data = departure.breeding.date[departure.breeding.date$m==M,]
  departure.trend[[a]]=list(rep(NA,10))
  for (i in 1:10){
    print(i)
    departure = data[data$AnalysisGroup==species[[i]],]  
    sites = as.character(unique(departure$SiteID))
    site.number = length(sites)
    departure.trend[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=22))
    for (j in 1:site.number){
      print(j)
      departure.subset = na.omit(departure[departure$SiteID==sites[j],])
      if (nrow(departure.subset)!=0){
        trend = lm(scale(Departure)~scale(Year),
                   data = departure.subset)
        
        s.y=sd(departure.subset$Departure)
        s.x=sd(departure.subset$Year)
        
        residual = resid(trend)
        predicted = predict(trend)
        test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
        test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
        test= min(test.2)# the smallest p-value from test.2
        departure.trend[[a]][[i]][j,1]=species[[i]]
        departure.trend[[a]][[i]][j,2]=sites[j]
        departure.trend[[a]][[i]][j,19]=M
        departure.trend[[a]][[i]][j,20]=N
        departure.trend[[a]][[i]][j,21]=summary(trend)[[4]][2]*s.y/s.x
        departure.trend[[a]][[i]][j,22]=summary(trend)[[4]][4]*s.y/s.x
        
        if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
          departure.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,4,6,8)] # summary of linear regression outputs
          departure.trend[[a]][[i]][j,11]=summary(trend)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          departure.trend[[a]][[i]][j,13]=summary(trend)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          departure.trend[[a]][[i]][j,15]="linear"
          departure.trend[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
        } else{#violation of the linear assmption: use quadratic form of the equation instead
          trend = lm(scale(Departure)~scale(Year)+I(scale(Year)^2),
                     data = departure.subset)
          residual = resid(trend)
          predicted = predict(trend)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          departure.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
          departure.trend[[a]][[i]][j,c(7:10)]=summary(trend)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
          departure.trend[[a]][[i]][j,11]=summary(trend)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          departure.trend[[a]][[i]][j,13]=summary(trend)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          departure.trend[[a]][[i]][j,12]=summary(trend)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          departure.trend[[a]][[i]][j,14]=summary(trend)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          departure.trend[[a]][[i]][j,15]="quadratic"
          departure.trend[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
        }
      }
    }
    write.table(departure.trend[[a]][[i]], file = "trend.csv", sep = " ",append=TRUE)
  }
}

## Analysis of data on the length of temporal occurrence windows (length of temporal occurrence windows represented as "stay" in data set)
stay.trend = list(rep(NA,4))
write.table(NA, file = "trend.csv", sep = " ")
for (a in 1:4){
  M=N_Year[a]
  N=Span_Year[a]
  data = stay.breeding.date[stay.breeding.date$m==M,]
  stay.trend[[a]]=list(rep(NA,10))
  for (i in 1:10){
    print(i)
    stay = data[data$AnalysisGroup==species[[i]],]  
    sites = as.character(unique(stay$SiteID))
    site.number = length(sites)
    stay.trend[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=22))
    for (j in 1:site.number){
      print(j)
      stay.subset = na.omit(stay[stay$SiteID==sites[j],])
      if (nrow(stay.subset)!=0){
        trend = lm(scale(Stay)~scale(Year),
                   data = stay.subset)
        
        s.y=sd(stay.subset$Stay)
        s.x=sd(stay.subset$Year)
        
        residual = resid(trend)
        predicted = predict(trend)
        test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
        test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
        test= min(test.2)# the smallest p-value from test.2
        stay.trend[[a]][[i]][j,1]=species[[i]]
        stay.trend[[a]][[i]][j,2]=sites[j]
        stay.trend[[a]][[i]][j,19]=M
        stay.trend[[a]][[i]][j,20]=N
        stay.trend[[a]][[i]][j,21]=summary(trend)[[4]][2]*s.y/s.x
        stay.trend[[a]][[i]][j,22]=summary(trend)[[4]][4]*s.y/s.x
        
        if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
          stay.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,4,6,8)] # summary of linear regression outputs
          stay.trend[[a]][[i]][j,11]=summary(trend)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          stay.trend[[a]][[i]][j,13]=summary(trend)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          stay.trend[[a]][[i]][j,15]="linear"
          stay.trend[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
        } else{#violation of the linear assmption: use quadratic form of the equation instead
          trend = lm(scale(Stay)~scale(Year)+I(scale(Year)^2),
                     data = stay.subset)
          residual = resid(trend)
          predicted = predict(trend)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          stay.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
          stay.trend[[a]][[i]][j,c(7:10)]=summary(trend)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
          stay.trend[[a]][[i]][j,11]=summary(trend)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          stay.trend[[a]][[i]][j,13]=summary(trend)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          stay.trend[[a]][[i]][j,12]=summary(trend)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          stay.trend[[a]][[i]][j,14]=summary(trend)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          stay.trend[[a]][[i]][j,15]="quadratic"
          stay.trend[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
        }
      }
    }
    write.table(stay.trend[[a]][[i]], file = "trend.csv", sep = " ",append=TRUE)
  }
}

### Temporal trend of phenophase: nonbreeding taxa
##  Read in data files
arrival.wintering.date = read.table('wintering arrival.txt', header=T)
head(arrival.wintering.date)
departure.wintering.date = read.table('wintering departure.txt', header=T)
head(departure.wintering.date)
stay.wintering.date = read.table('wintering stay.txt', header=T)
head(stay.wintering.date)

species = list('Goose','OtherBirds')

## define different data sets (of different trimming criteria)
N_Year = c(5,6,11,10)#criterion on the minimum number of time steps in a time series
Span_Year = c(8,10,11,14)#criterion on the minimum span of time series

## Analysis of initial occurrence data (initial occurrence represented as "arrival" in data set)
arrival.trend = list(rep(NA,4))
write.table(NA, file = "trend.csv", sep = " ")
for (a in 1:4){
  M=N_Year[a]
  N=Span_Year[a]
  data = arrival.wintering.date[arrival.wintering.date$m==M,]
  arrival.trend[[a]]=list(rep(NA,2))
  for (i in 1:2){
    print(i)
    arrival = data[data$AnalysisGroup==species[[i]],]  
    sites = as.character(unique(arrival$SiteID))
    site.number = length(sites)
    arrival.trend[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=22))
    for (j in 1:site.number){
      print(j)
      arrival.subset = na.omit(arrival[arrival$SiteID==sites[j],])
      if (nrow(arrival.subset)!=0){
        trend = lm(scale(Arrival)~scale(Year),
                   data = arrival.subset)
        
        s.y=sd(arrival.subset$Arrival)
        s.x=sd(arrival.subset$Year)
        
        residual = resid(trend)
        predicted = predict(trend)
        test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
        test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
        test= min(test.2)# the smallest p-value from test.2
        arrival.trend[[a]][[i]][j,1]=species[[i]]
        arrival.trend[[a]][[i]][j,2]=sites[j]
        arrival.trend[[a]][[i]][j,19]=M
        arrival.trend[[a]][[i]][j,20]=N
        arrival.trend[[a]][[i]][j,21]=summary(trend)[[4]][2]*s.y/s.x
        arrival.trend[[a]][[i]][j,22]=summary(trend)[[4]][4]*s.y/s.x
        
        if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
          arrival.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,4,6,8)] # summary of linear regression outputs
          arrival.trend[[a]][[i]][j,11]=summary(trend)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          arrival.trend[[a]][[i]][j,13]=summary(trend)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          arrival.trend[[a]][[i]][j,15]="linear"
          arrival.trend[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
        } else{#violation of the linear assmption: use quadratic form of the equation instead
          trend = lm(scale(Arrival)~scale(Year)+I(scale(Year)^2),
                     data = arrival.subset)
          residual = resid(trend)
          predicted = predict(trend)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          arrival.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
          arrival.trend[[a]][[i]][j,c(7:10)]=summary(trend)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
          arrival.trend[[a]][[i]][j,11]=summary(trend)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          arrival.trend[[a]][[i]][j,13]=summary(trend)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          arrival.trend[[a]][[i]][j,12]=summary(trend)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          arrival.trend[[a]][[i]][j,14]=summary(trend)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          arrival.trend[[a]][[i]][j,15]="quadratic"
          arrival.trend[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive      
        }
      }
    }
    write.table(arrival.trend[[a]][[i]], file = "trend.csv", sep = " ",append=TRUE)
  }
}

## Analysis of disappearance data (disappearance represented as "departure" in data set)
departure.trend = list(rep(NA,4))
write.table(NA, file = "trend.csv", sep = " ")
for (a in 1:4){
  M=N_Year[a]
  N=Span_Year[a]
  data = departure.wintering.date[departure.wintering.date$m==M,]
  departure.trend[[a]]=list(rep(NA,2))
  for (i in 1:2){
    print(i)
    departure = data[data$AnalysisGroup==species[[i]],]  
    sites = as.character(unique(departure$SiteID))
    site.number = length(sites)
    departure.trend[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=22))
    for (j in 1:site.number){
      print(j)
      departure.subset = na.omit(departure[departure$SiteID==sites[j],])
      if (nrow(departure.subset)!=0){
        trend = lm(scale(Departure)~scale(Year),
                   data = departure.subset)
        
        s.y=sd(departure.subset$Departure)
        s.x=sd(departure.subset$Year)
        
        residual = resid(trend)
        predicted = predict(trend)
        test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
        test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
        test= min(test.2)# the smallest p-value from test.2
        departure.trend[[a]][[i]][j,1]=species[[i]]
        departure.trend[[a]][[i]][j,2]=sites[j]
        departure.trend[[a]][[i]][j,19]=M
        departure.trend[[a]][[i]][j,20]=N
        departure.trend[[a]][[i]][j,21]=summary(trend)[[4]][2]*s.y/s.x
        departure.trend[[a]][[i]][j,22]=summary(trend)[[4]][4]*s.y/s.x
        
        if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
          departure.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,4,6,8)] # summary of linear regression outputs
          departure.trend[[a]][[i]][j,11]=summary(trend)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          departure.trend[[a]][[i]][j,13]=summary(trend)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          departure.trend[[a]][[i]][j,15]="linear"
          departure.trend[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
        } else{#violation of the linear assmption: use quadratic form of the equation instead
          trend = lm(scale(Departure)~scale(Year)+I(scale(Year)^2),
                     data = departure.subset)
          residual = resid(trend)
          predicted = predict(trend)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          departure.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
          departure.trend[[a]][[i]][j,c(7:10)]=summary(trend)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
          departure.trend[[a]][[i]][j,11]=summary(trend)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          departure.trend[[a]][[i]][j,13]=summary(trend)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          departure.trend[[a]][[i]][j,12]=summary(trend)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          departure.trend[[a]][[i]][j,14]=summary(trend)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          departure.trend[[a]][[i]][j,15]="quadratic"
          departure.trend[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
        }
      }
    }
    write.table(departure.trend[[a]][[i]], file = "trend.csv", sep = " ",append=TRUE)
  }
}

## Analysis of data on the length of temporal occurrence windows (length of temporal occurrence windows represented as "stay" in data set)
stay.trend = list(rep(NA,4))
write.table(NA, file = "trend.csv", sep = " ")
for (a in 1:4){
  M=N_Year[a]
  N=Span_Year[a]
  data = stay.wintering.date[stay.wintering.date$m==M,]
  stay.trend[[a]]=list(rep(NA,2))
  for (i in 1:2){
    print(i)
    stay = data[data$AnalysisGroup==species[[i]],]  
    sites = as.character(unique(stay$SiteID))
    site.number = length(sites)
    stay.trend[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=22))
    for (j in 1:site.number){
      print(j)
      stay.subset = na.omit(stay[stay$SiteID==sites[j],])
      if (nrow(stay.subset)!=0){
        trend = lm(scale(Stay)~scale(Year),
                   data = stay.subset)
        
        s.y=sd(stay.subset$Stay)
        s.x=sd(stay.subset$Year)
        
        residual = resid(trend)
        predicted = predict(trend)
        test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
        test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
        test= min(test.2)# the smallest p-value from test.2
        stay.trend[[a]][[i]][j,1]=species[[i]]
        stay.trend[[a]][[i]][j,2]=sites[j]
        stay.trend[[a]][[i]][j,19]=M
        stay.trend[[a]][[i]][j,20]=N
        stay.trend[[a]][[i]][j,21]=summary(trend)[[4]][2]*s.y/s.x
        stay.trend[[a]][[i]][j,22]=summary(trend)[[4]][4]*s.y/s.x
        
        if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
          stay.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,4,6,8)] # summary of linear regression outputs
          stay.trend[[a]][[i]][j,11]=summary(trend)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          stay.trend[[a]][[i]][j,13]=summary(trend)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          stay.trend[[a]][[i]][j,15]="linear"
          stay.trend[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
        } else{#violation of the linear assmption: use quadratic form of the equation instead
          trend = lm(scale(Stay)~scale(Year)+I(scale(Year)^2),
                     data = stay.subset)
          residual = resid(trend)
          predicted = predict(trend)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          stay.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
          stay.trend[[a]][[i]][j,c(7:10)]=summary(trend)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
          stay.trend[[a]][[i]][j,11]=summary(trend)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          stay.trend[[a]][[i]][j,13]=summary(trend)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          stay.trend[[a]][[i]][j,12]=summary(trend)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          stay.trend[[a]][[i]][j,14]=summary(trend)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          stay.trend[[a]][[i]][j,15]="quadratic"
          stay.trend[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
        }
      }
    }
    write.table(stay.trend[[a]][[i]], file = "trend.csv", sep = " ",append=TRUE)
  }
}

### Relationship of phenophase to temperature: breeding taxa
## read in data file
arrival.breeding.date = read.table('breeding arrival w temp.txt', header=T)
head(arrival.breeding.date)
departure.breeding.date = read.table('breeding departure w temp.txt', header=T)
head(departure.breeding.date)
stay.breeding.date = read.table('breeding stay w temp.txt', header=T)
head(stay.breeding.date)

species = list('Bee','Cicada','Cricket','Cuckoo',
               'Frog','Goose','Hoopoe',
               'Swallow','Swift',
               'Toad')

## define different data sets (of different trimming criteria)
N_Year = c(5,6,11,10)#criterion on the minimum number of time steps in a time series
Span_Year = c(8,10,11,14)#criterion on the minimum span of time series

## Analysis of initial occurrence data (initial occurrence represented as "arrival" in data set)
arrival.rel = list(rep(NA,4))
write.table(NA, file = "relation.csv", sep = " ")
for (a in 1:4){
  M=N_Year[a]
  N=Span_Year[a]
  data = arrival.breeding.date[arrival.breeding.date$m==M,]
  arrival.rel[[a]]=list(rep(NA,10))
  for (i in 1:10){
    print(i)
    arrival = data[data$AnalysisGroup==species[[i]],]  
    sites = as.character(unique(arrival$SiteID))
    site.number = length(sites)
    arrival.rel[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=22))
    for (j in 1:site.number){
      print(j)
      arrival.subset = na.omit(arrival[arrival$SiteID==sites[j],])
      if (nrow(arrival.subset)!=0){
        relation = lm(scale(Arrival)~scale(SpringTemp),
                   data = arrival.subset)
        
        s.y=sd(arrival.subset$Arrival)
        s.x=sd(arrival.subset$SpringTemp)
        
        residual = resid(relation)
        predicted = predict(relation)
        test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
        test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
        test= min(test.2)# the smallest p-value from test.2
        arrival.rel[[a]][[i]][j,1]=species[[i]]
        arrival.rel[[a]][[i]][j,2]=sites[j]
        arrival.rel[[a]][[i]][j,19]=M
        arrival.rel[[a]][[i]][j,20]=N
        arrival.rel[[a]][[i]][j,21]=summary(relation)[[4]][2]*s.y/s.x
        arrival.rel[[a]][[i]][j,22]=summary(relation)[[4]][4]*s.y/s.x
        
        if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
          arrival.rel[[a]][[i]][j,c(3:6)]=summary(relation)[[4]][c(2,4,6,8)] # summary of linear regression outputs
          arrival.rel[[a]][[i]][j,11]=summary(relation)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          arrival.rel[[a]][[i]][j,13]=summary(relation)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          arrival.rel[[a]][[i]][j,15]="linear"
          arrival.rel[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
        } else{#violation of the linear assmption: use quadratic form of the equation instead
          relation = lm(scale(Arrival)~scale(SpringTemp)+I(scale(SpringTemp)^2),
                     data = arrival.subset)
          residual = resid(relation)
          predicted = predict(relation)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          arrival.rel[[a]][[i]][j,c(3:6)]=summary(relation)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
          arrival.rel[[a]][[i]][j,c(7:10)]=summary(relation)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
          arrival.rel[[a]][[i]][j,11]=summary(relation)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          arrival.rel[[a]][[i]][j,13]=summary(relation)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          arrival.rel[[a]][[i]][j,12]=summary(relation)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          arrival.rel[[a]][[i]][j,14]=summary(relation)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          arrival.rel[[a]][[i]][j,15]="quadratic"
          arrival.rel[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
        }
      }
    }
    write.table(arrival.rel[[a]][[i]], file = "relation.csv", sep = " ",append=TRUE)
  }
}

## Analysis of disappearance data (disappearance represented as "departure" in data set)
departure.rel = list(rep(NA,4))
write.table(NA, file = "relation.csv", sep = " ")
for (a in 1:4){
  M=N_Year[a]
  N=Span_Year[a]
  data = departure.breeding.date[departure.breeding.date$m==M,]
  departure.rel[[a]]=list(rep(NA,10))
  for (i in 1:10){
    print(i)
    departure = data[data$AnalysisGroup==species[[i]],]  
    sites = as.character(unique(departure$SiteID))
    site.number = length(sites)
    departure.rel[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=22))
    for (j in 1:site.number){
      print(j)
      departure.subset = na.omit(departure[departure$SiteID==sites[j],])
      if (nrow(departure.subset)!=0){
        relation = lm(scale(Departure)~scale(FallTemp),
                   data = departure.subset)
        
        s.y=sd(departure.subset$Departure)
        s.x=sd(departure.subset$FallTemp)
        
        residual = resid(relation)
        predicted = predict(relation)
        test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
        test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
        test= min(test.2)# the smallest p-value from test.2
        departure.rel[[a]][[i]][j,1]=species[[i]]
        departure.rel[[a]][[i]][j,2]=sites[j]
        departure.rel[[a]][[i]][j,19]=M
        departure.rel[[a]][[i]][j,20]=N
        departure.rel[[a]][[i]][j,21]=summary(relation)[[4]][2]*s.y/s.x
        departure.rel[[a]][[i]][j,22]=summary(relation)[[4]][4]*s.y/s.x
        
        if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
          departure.rel[[a]][[i]][j,c(3:6)]=summary(relation)[[4]][c(2,4,6,8)] # summary of linear regression outputs
          departure.rel[[a]][[i]][j,11]=summary(relation)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          departure.rel[[a]][[i]][j,13]=summary(relation)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          departure.rel[[a]][[i]][j,15]="linear"
          departure.rel[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
        } else{#violation of the linear assmption: use quadratic form of the equation instead
          relation = lm(scale(Departure)~scale(FallTemp)+I(scale(FallTemp)^2),
                     data = departure.subset)
          residual = resid(relation)
          predicted = predict(relation)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          departure.rel[[a]][[i]][j,c(3:6)]=summary(relation)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
          departure.rel[[a]][[i]][j,c(7:10)]=summary(relation)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
          departure.rel[[a]][[i]][j,11]=summary(relation)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          departure.rel[[a]][[i]][j,13]=summary(relation)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          departure.rel[[a]][[i]][j,12]=summary(relation)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          departure.rel[[a]][[i]][j,14]=summary(relation)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          departure.rel[[a]][[i]][j,15]="quadratic"
          departure.rel[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
        }
      }
    }
    write.table(departure.rel[[a]][[i]], file = "relation.csv", sep = " ",append=TRUE)
  }
}

## Analysis of data on the length of temporal occurrence windows (length of temporal occurrence windows represented as "stay" in data set)
stay.rel = list(rep(NA,4))
write.table(NA, file = "relation.csv", sep = " ")
for (a in 1:4){
  M=N_Year[a]
  N=Span_Year[a]
  data = stay.breeding.date[stay.breeding.date$m==M,]
  stay.rel[[a]]=list(rep(NA,10))
  for (i in 1:10){
    print(i)
    stay = data[data$AnalysisGroup==species[[i]],]  
    sites = as.character(unique(stay$SiteID))
    site.number = length(sites)
    stay.rel[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=22))
    for (j in 1:site.number){
      print(j)
      stay.subset = na.omit(stay[stay$SiteID==sites[j],])
      if (nrow(stay.subset)!=0){
        relation = lm(scale(Stay)~scale(SpringTemp),
                   data = stay.subset)
        
        s.y=sd(stay.subset$Stay)
        s.x=sd(stay.subset$SpringTemp)
        
        residual = resid(relation)
        predicted = predict(relation)
        test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
        test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
        test= min(test.2)# the smallest p-value from test.2
        stay.rel[[a]][[i]][j,1]=species[[i]]
        stay.rel[[a]][[i]][j,2]=sites[j]
        stay.rel[[a]][[i]][j,19]=M
        stay.rel[[a]][[i]][j,20]=N
        stay.rel[[a]][[i]][j,21]=summary(relation)[[4]][2]*s.y/s.x
        stay.rel[[a]][[i]][j,22]=summary(relation)[[4]][4]*s.y/s.x
        
        if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
          stay.rel[[a]][[i]][j,c(3:6)]=summary(relation)[[4]][c(2,4,6,8)] # summary of linear regression outputs
          stay.rel[[a]][[i]][j,11]=summary(relation)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          stay.rel[[a]][[i]][j,13]=summary(relation)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          stay.rel[[a]][[i]][j,15]="linear"
          stay.rel[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
        } else{#violation of the linear assmption: use quadratic form of the equation instead
          relation = lm(scale(Stay)~scale(SpringTemp)+I(scale(SpringTemp)^2),
                     data = stay.subset)
          residual = resid(relation)
          predicted = predict(relation)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          stay.rel[[a]][[i]][j,c(3:6)]=summary(relation)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
          stay.rel[[a]][[i]][j,c(7:10)]=summary(relation)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
          stay.rel[[a]][[i]][j,11]=summary(relation)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          stay.rel[[a]][[i]][j,13]=summary(relation)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          stay.rel[[a]][[i]][j,12]=summary(relation)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          stay.rel[[a]][[i]][j,14]=summary(relation)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          stay.rel[[a]][[i]][j,15]="quadratic"
          stay.rel[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
        }
      }
    }
    write.table(stay.rel[[a]][[i]], file = "relation.csv", sep = " ",append=TRUE)
  }
}

### Relationship of phenophase to temperature: nonbreeding taxa
## read in data file
arrival.wintering.date = read.table('wintering arrival w temp.txt', header=T)
head(arrival.wintering.date)
departure.wintering.date = read.table('wintering departure w temp.txt', header=T)
head(departure.wintering.date)
stay.wintering.date = read.table('wintering stay w temp.txt', header=T)
head(stay.wintering.date)

species = list('Goose','OtherBirds')

## define different data sets (of different trimming criteria)
N_Year = c(5,6,11,10)#criterion on the minimum number of time steps in a time series
Span_Year = c(8,10,11,14)#criterion on the minimum span of time series

## Analysis of initial occurrence data (initial occurrence represented as "arrival" in data set)
arrival.rel = list(rep(NA,4))
write.table(NA, file = "relation.csv", sep = " ")
for (a in 1:4){
  M=N_Year[a]
  N=Span_Year[a]
  data = arrival.wintering.date[arrival.wintering.date$m==M,]
  arrival.rel[[a]]=list(rep(NA,2))
  for (i in 1:2){
    print(i)
    arrival = data[data$AnalysisGroup==species[[i]],]  
    sites = as.character(unique(arrival$SiteID))
    site.number = length(sites)
    arrival.rel[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=22))
    for (j in 1:site.number){
      print(j)
      arrival.subset = na.omit(arrival[arrival$SiteID==sites[j],])
      if (nrow(arrival.subset)!=0){
        relation = lm(scale(Arrival)~scale(FallTemp),
                      data = arrival.subset)
        
        s.y=sd(arrival.subset$Arrival)
        s.x=sd(arrival.subset$FallTemp)
        
        residual = resid(relation)
        predicted = predict(relation)
        test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
        test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
        test= min(test.2)# the smallest p-value from test.2
        arrival.rel[[a]][[i]][j,1]=species[[i]]
        arrival.rel[[a]][[i]][j,2]=sites[j]
        arrival.rel[[a]][[i]][j,19]=M
        arrival.rel[[a]][[i]][j,20]=N
        arrival.rel[[a]][[i]][j,21]=summary(relation)[[4]][2]*s.y/s.x
        arrival.rel[[a]][[i]][j,22]=summary(relation)[[4]][4]*s.y/s.x
        
        if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
          arrival.rel[[a]][[i]][j,c(3:6)]=summary(relation)[[4]][c(2,4,6,8)] # summary of linear regression outputs
          arrival.rel[[a]][[i]][j,11]=summary(relation)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          arrival.rel[[a]][[i]][j,13]=summary(relation)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          arrival.rel[[a]][[i]][j,15]="linear"
          arrival.rel[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
        } else{#violation of the linear assmption: use quadratic form of the equation instead
          relation = lm(scale(Arrival)~scale(FallTemp)+I(scale(FallTemp)^2),
                        data = arrival.subset)
          residual = resid(relation)
          predicted = predict(relation)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          arrival.rel[[a]][[i]][j,c(3:6)]=summary(relation)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
          arrival.rel[[a]][[i]][j,c(7:10)]=summary(relation)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
          arrival.rel[[a]][[i]][j,11]=summary(relation)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          arrival.rel[[a]][[i]][j,13]=summary(relation)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          arrival.rel[[a]][[i]][j,12]=summary(relation)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          arrival.rel[[a]][[i]][j,14]=summary(relation)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          arrival.rel[[a]][[i]][j,15]="quadratic"
          arrival.rel[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
        }
      }
    }
    write.table(arrival.rel[[a]][[i]], file = "relation.csv", sep = " ",append=TRUE)
  }
}


## Analysis of disappearance data (disappearance represented as "departure" in data set)
departure.rel = list(rep(NA,4))
write.table(NA, file = "relation.csv", sep = " ")
for (a in 1:4){
  M=N_Year[a]
  N=Span_Year[a]
  data = departure.wintering.date[departure.wintering.date$m==M,]
  departure.rel[[a]]=list(rep(NA,2))
  for (i in 1:2){
    print(i)
    departure = data[data$AnalysisGroup==species[[i]],]  
    sites = as.character(unique(departure$SiteID))
    site.number = length(sites)
    departure.rel[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=22))
    for (j in 1:site.number){
      print(j)
      departure.subset = na.omit(departure[departure$SiteID==sites[j],])
      if (nrow(departure.subset)!=0){
        relation = lm(scale(Departure)~scale(SpringTemp),
                      data = departure.subset)
        
        s.y=sd(departure.subset$Departure)
        s.x=sd(departure.subset$SpringTemp)
        
        residual = resid(relation)
        predicted = predict(relation)
        test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
        test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
        test= min(test.2)# the smallest p-value from test.2
        departure.rel[[a]][[i]][j,1]=species[[i]]
        departure.rel[[a]][[i]][j,2]=sites[j]
        departure.rel[[a]][[i]][j,19]=M
        departure.rel[[a]][[i]][j,20]=N
        departure.rel[[a]][[i]][j,21]=summary(relation)[[4]][2]*s.y/s.x
        departure.rel[[a]][[i]][j,22]=summary(relation)[[4]][4]*s.y/s.x
        
        if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
          departure.rel[[a]][[i]][j,c(3:6)]=summary(relation)[[4]][c(2,4,6,8)] # summary of linear regression outputs
          departure.rel[[a]][[i]][j,11]=summary(relation)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          departure.rel[[a]][[i]][j,13]=summary(relation)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          departure.rel[[a]][[i]][j,15]="linear"
          departure.rel[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
        } else{#violation of the linear assmption: use quadratic form of the equation instead
          relation = lm(scale(Departure)~scale(SpringTemp)+I(scale(SpringTemp)^2),
                        data = departure.subset)
          residual = resid(relation)
          predicted = predict(relation)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          departure.rel[[a]][[i]][j,c(3:6)]=summary(relation)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
          departure.rel[[a]][[i]][j,c(7:10)]=summary(relation)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
          departure.rel[[a]][[i]][j,11]=summary(relation)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          departure.rel[[a]][[i]][j,13]=summary(relation)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          departure.rel[[a]][[i]][j,12]=summary(relation)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          departure.rel[[a]][[i]][j,14]=summary(relation)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          departure.rel[[a]][[i]][j,15]="quadratic"
          departure.rel[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
        }
      }
    }
    write.table(departure.rel[[a]][[i]], file = "relation.csv", sep = " ",append=TRUE)
  }
}

## Analysis of data on the length of temporal occurrence windows (length of temporal occurrence windows represented as "stay" in data set)
stay.rel = list(rep(NA,4))
write.table(NA, file = "relation.csv", sep = " ")
for (a in 1:4){
  M=N_Year[a]
  N=Span_Year[a]
  data = stay.wintering.date[stay.wintering.date$m==M,]
  stay.rel[[a]]=list(rep(NA,2))
  for (i in 1:2){
    print(i)
    stay = data[data$AnalysisGroup==species[[i]],]  
    sites = as.character(unique(stay$SiteID))
    site.number = length(sites)
    stay.rel[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=22))
    for (j in 1:site.number){
      print(j)
      stay.subset = na.omit(stay[stay$SiteID==sites[j],])
      if (nrow(stay.subset)!=0){
        relation = lm(scale(Stay)~scale(SpringTemp),
                      data = stay.subset)
        
        s.y=sd(stay.subset$Stay)
        s.x=sd(stay.subset$SpringTemp)
        
        residual = resid(relation)
        predicted = predict(relation)
        test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
        test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
        test= min(test.2)# the smallest p-value from test.2
        stay.rel[[a]][[i]][j,1]=species[[i]]
        stay.rel[[a]][[i]][j,2]=sites[j]
        stay.rel[[a]][[i]][j,19]=M
        stay.rel[[a]][[i]][j,20]=N
        stay.rel[[a]][[i]][j,21]=summary(relation)[[4]][2]*s.y/s.x
        stay.rel[[a]][[i]][j,22]=summary(relation)[[4]][4]*s.y/s.x
        
        if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
          stay.rel[[a]][[i]][j,c(3:6)]=summary(relation)[[4]][c(2,4,6,8)] # summary of linear regression outputs
          stay.rel[[a]][[i]][j,11]=summary(relation)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          stay.rel[[a]][[i]][j,13]=summary(relation)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          stay.rel[[a]][[i]][j,15]="linear"
          stay.rel[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
        } else{#violation of the linear assmption: use quadratic form of the equation instead
          relation = lm(scale(Stay)~scale(SpringTemp)+I(scale(SpringTemp)^2),
                        data = stay.subset)
          residual = resid(relation)
          predicted = predict(relation)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          stay.rel[[a]][[i]][j,c(3:6)]=summary(relation)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
          stay.rel[[a]][[i]][j,c(7:10)]=summary(relation)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
          stay.rel[[a]][[i]][j,11]=summary(relation)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          stay.rel[[a]][[i]][j,13]=summary(relation)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          stay.rel[[a]][[i]][j,12]=summary(relation)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
          stay.rel[[a]][[i]][j,14]=summary(relation)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
          stay.rel[[a]][[i]][j,15]="quadratic"
          stay.rel[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
        }
      }
    }
    write.table(stay.rel[[a]][[i]], file = "relation.csv", sep = " ",append=TRUE)
  }
}

### Temporal trend of temperature
## read in data file
temp = read.table('Site temperature.txt', header=TRUE)
temp$Site=as.factor(temp$Site)
head(temp,n=15)
spring = temp[,c(1:4,12)]
head(spring)
fall = temp[,c(1:4,13)]
head(fall)

## correlation between spring and fall temp
cor.test(temp$SpringTemp,temp$FallTemp)

## Analysis of spring temperature
site = as.character(unique(spring$Site))
site.number=length(site)
spring.trend=as.data.frame(matrix(NA, nrow=site.number,ncol=17))

for (j in 1:site.number){
  print(j)
  data = spring[spring$Site==site[[j]],]  
  if (nrow(data)!=0){
    trend = lm(scale(SpringTemp)~scale(Year),
               data = data)
    residual = resid(trend)
    predicted = predict(trend)
    test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
    test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
    test= min(test.2)# the smallest p-value from test.2
    spring.trend[j,1]=site[j]
    if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
      spring.trend[j,c(2:5)]=summary(trend)[[4]][c(2,4,6,8)] # summary of linear regression outputs
      spring.trend[j,10]=summary(trend)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
      spring.trend[j,12]=summary(trend)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
      spring.trend[j,14]="linear"
      spring.trend[j,c(15:17)]=c(test.1,test.2)-0.05
    } else{#violation of the linear assmption: use quadratic form of the equation instead
      trend = lm(scale(SpringTemp)~scale(Year)+I(scale(Year)^2),
                 data = data)
      residual = resid(trend)
      predicted = predict(trend)
      test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
      test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
      spring.trend[j,c(2:5)]=summary(trend)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
      spring.trend[j,c(6:9)]=summary(trend)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
      spring.trend[j,10]=summary(trend)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
      spring.trend[j,12]=summary(trend)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
      spring.trend[j,11]=summary(trend)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
      spring.trend[j,13]=summary(trend)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
      spring.trend[j,14]="quadratic"
      spring.trend[j,c(15:17)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
    }
  }  
}
write.table(spring.trend, file = "trend.csv", sep = " ")


#### Community-level analysis

community = read.table('overlap.txt', header=T)
community$SiteID = as.factor(community$SiteID)
community$TaxonPair = as.character(community$TaxonPair)
head(community)

pair = as.character(unique(community$TaxonPair))
length(pair)#27 pairs

## define different data sets (of different trimming criteria)
N_Year = c(5,6,11,10)#criterion on the minimum number of time steps in a time series
Span_Year = c(8,10,11,14)#criterion on the minimum span of time series

## overlap analysis
overlap.trend = list(rep(NA,4))
# in running, realized the following time series have all 0's
# i=2(j=2),i=7(j=6)

write.table(NA, file = "trend.csv", sep = " ")
for (a in 1:4){
  print(a)
  M=N_Year[a]
  N=Span_Year[a]
  data = community[community$m==M,]
  overlap.trend[[a]]=list(rep(NA,27))
  for (i in c(1,3:6,8:27)){#skip pairs i=2 and i=7 for now
    print(i)
    overlap = data[data$TaxonPair==pair[i],]  
    sites = as.character(unique(overlap$SiteID))
    site.number = length(sites)
    overlap.trend[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=20))
    if (nrow(overlap)>0){
      for (j in 1:site.number){
        print(j)
        overlap.subset = na.omit(overlap[overlap$SiteID==sites[j],])
        if (nrow(overlap.subset)!=0){
          trend = lm(scale(Overlap)~scale(Year),
                        data = overlap.subset)
          residual = resid(trend)
          predicted = predict(trend)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          test= min(test.2)# the smallest p-value from test.2
          overlap.trend[[a]][[i]][j,1]=pair[i]
          overlap.trend[[a]][[i]][j,2]=sites[j]
          overlap.trend[[a]][[i]][j,19]=M
          overlap.trend[[a]][[i]][j,20]=N
          if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
            overlap.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,4,6,8)] # summary of linear regression outputs
            overlap.trend[[a]][[i]][j,11]=summary(trend)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,13]=summary(trend)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,15]="linear"
            overlap.trend[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
          } else{#violation of the linear assmption: use quadratic form of the equation instead
            trend = lm(scale(Overlap)~scale(Year)+I(scale(Year)^2),
                          data = overlap.subset)
            residual = resid(trend)
            predicted = predict(trend)
            test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
            test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
            overlap.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
            overlap.trend[[a]][[i]][j,c(7:10)]=summary(trend)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
            overlap.trend[[a]][[i]][j,11]=summary(trend)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,13]=summary(trend)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,12]=summary(trend)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,14]=summary(trend)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,15]="quadratic"
            overlap.trend[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
          }
        }
      }
      write.table(overlap.trend[[a]][[i]], file = "trend.csv", sep = " ",append=TRUE)
    }
  }
  
  for (i in 2:2){
    print(i)
    overlap = data[data$TaxonPair==pair[i],]  
    overlap.trend[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=20))
    if (nrow(overlap)>0){
      sites = as.character(unique(overlap$SiteID))
      site.number = length(sites)
      for (j in c(1,3:site.number)){
        print(j)
        overlap.subset = na.omit(overlap[overlap$SiteID==sites[j],])
        if (nrow(overlap.subset)!=0){
          trend = lm(scale(Overlap)~scale(Year),
                     data = overlap.subset)
          residual = resid(trend)
          predicted = predict(trend)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          test= min(test.2)# the smallest p-value from test.2
          overlap.trend[[a]][[i]][j,1]=pair[i]
          overlap.trend[[a]][[i]][j,2]=sites[j]
          overlap.trend[[a]][[i]][j,19]=M
          overlap.trend[[a]][[i]][j,20]=N
          if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
            overlap.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,4,6,8)] # summary of linear regression outputs
            overlap.trend[[a]][[i]][j,11]=summary(trend)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,13]=summary(trend)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,15]="linear"
            overlap.trend[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
          } else{#violation of the linear assmption: use quadratic form of the equation instead
            trend = lm(scale(Overlap)~scale(Year)+I(scale(Year)^2),
                       data = overlap.subset)
            residual = resid(trend)
            predicted = predict(trend)
            test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
            test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
            overlap.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
            overlap.trend[[a]][[i]][j,c(7:10)]=summary(trend)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
            overlap.trend[[a]][[i]][j,11]=summary(trend)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,13]=summary(trend)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,12]=summary(trend)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,14]=summary(trend)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,15]="quadratic"
            overlap.trend[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
          }
        }
      }
      write.table(overlap.trend[[a]][[i]], file = "trend.csv", sep = " ",append=TRUE)
    }    
  }
  
  for (i in 7:7){
    print(i)
    overlap = data[data$TaxonPair==pair[i],]  
    overlap.trend[[a]][[i]]=as.data.frame(matrix(NA, nrow=site.number,ncol=20))
    if (nrow(overlap)>0){
      sites = as.character(unique(overlap$SiteID))
      site.number = length(sites)
      for (j in c(1:5,7:site.number)){
        print(j)
        overlap.subset = na.omit(overlap[overlap$SiteID==sites[j],])
        if (nrow(overlap.subset)!=0){
          trend = lm(scale(Overlap)~scale(Year),
                     data = overlap.subset)
          residual = resid(trend)
          predicted = predict(trend)
          test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
          test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
          test= min(test.2)# the smallest p-value from test.2
          overlap.trend[[a]][[i]][j,1]=pair[i]
          overlap.trend[[a]][[i]][j,2]=sites[j]
          overlap.trend[[a]][[i]][j,19]=M
          overlap.trend[[a]][[i]][j,20]=N
          if (test>0.05){#no violation of linearity assumption; results from linear regression fed directly into results table
            overlap.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,4,6,8)] # summary of linear regression outputs
            overlap.trend[[a]][[i]][j,11]=summary(trend)[[4]][8]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,13]=summary(trend)[[4]][8]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,15]="linear"
            overlap.trend[[a]][[i]][j,16:18]=c(test.1,test.2)-0.05
          } else{#violation of the linear assmption: use quadratic form of the equation instead
            trend = lm(scale(Overlap)~scale(Year)+I(scale(Year)^2),
                       data = overlap.subset)
            residual = resid(trend)
            predicted = predict(trend)
            test.1 = summary(lm(residual~predicted))[[4]][8]#p value for the linear fit of test.1
            test.2 = summary(lm(residual ~ predicted + I(predicted^2)))[[4]][c(11,12)] #p value for the quadralinear fit of test.2
            overlap.trend[[a]][[i]][j,c(3:6)]=summary(trend)[[4]][c(2,5,8,11)] # summary of regression outputs: linear term
            overlap.trend[[a]][[i]][j,c(7:10)]=summary(trend)[[4]][c(3,6,9,12)] # summary of regression outputs: quadratic term
            overlap.trend[[a]][[i]][j,11]=summary(trend)[[4]][11]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,13]=summary(trend)[[4]][11]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,12]=summary(trend)[[4]][12]-0.05 # comparing p value of slope to 0.05 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,14]=summary(trend)[[4]][12]-0.1 # comparing p value of slope to 0.1 -> significant ones would be equal or smaller than 0
            overlap.trend[[a]][[i]][j,15]="quadratic"
            overlap.trend[[a]][[i]][j,c(16:18)]=c(test.1,test.2)-0.05 #assumption check result: different between p values and 0.05 -> nonsignificant ones would be positive
          }
        }
      }
      write.table(overlap.trend[[a]][[i]], file = "trend.csv", sep = " ",append=TRUE)
    }    
  }
}
