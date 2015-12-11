library(nlme)

#### arrival

### change directory to the arrival w year storage folder

taxa.all = list(c("Cicada","Cricket","Cuckoo","Frog","Goose","Swallow","Toad"),
                c("Cicada","Cricket","Frog","Goose","Swallow","Toad"),
                c("Cicada","Frog","Goose","Swallow","Toad"),
                c("Cicada","Frog","Goose","Swallow","Toad"))

m.all=c(5,6,10,11)
resample=1000

intercept = data.frame(matrix(NA,ncol=11,nrow=1))
names(intercept)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")
slope = data.frame(matrix(NA,ncol=11,nrow=1))
names(slope)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")

### running GLMMs
for (l in 1:4){
  print(l)
  
  taxa = taxa.all[[l]]
  taxa.number=length(taxa)
  
  for (t in 1:taxa.number){
    print(l)
    print(t)
    
    for (i in 1:resample){
      print(i)
      
      filename = paste(l,taxa[[t]],i, "subsample.csv")
      arrival = read.csv(filename, header=TRUE)
      arrival = arrival[arrival$whether_independent==1,]
      
      s.y=sd(arrival$Arrival)
      s.x=sd(arrival$Year)
      mean.y = mean(arrival$Arrival)
      mean.x = mean(arrival$Year)
      
      intercept.temporary = data.frame(matrix(NA,ncol=11,nrow=1))
      names(intercept.temporary)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")
      slope.temporary = data.frame(matrix(NA,ncol=11,nrow=1))
      names(slope.temporary)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")
      
      intercept.temporary[1,1]=m.all[l]
      intercept.temporary[1,2]=taxa[[t]]
      intercept.temporary[1,3]=i
      intercept.temporary[1,4]=length(unique(as.character(arrival$SiteID)))
      
      slope.temporary[1,1]=m.all[l]
      slope.temporary[1,2]=taxa[[t]]
      slope.temporary[1,3]=i
      slope.temporary[1,4]=length(unique(as.character(arrival$SiteID)))
      
      tryCatch({
        ## run random slope + intercept models (allows for correlation between random slope and intercept)
        m1 = lme(scale(Arrival)~scale(Year), random = ~scale(Year) | SiteID, data = na.omit(arrival))
        #m1 = lme(scale(Arrival)~scale(Year), random = ~scale(Year) | draw, data = na.omit(arrival))
        
        intercept.temporary[1,5:9]=summary(m1)[[20]][1,]
        intercept.temporary[1,10]=summary(m1)[[20]][1,1]*s.y + mean.y - summary(m1)[[20]][2,1]*s.y/s.x*mean.x
        intercept.temporary[1,11]=summary(m1)[[20]][1,2]*s.y + mean.y - summary(m1)[[20]][2,2]*s.y/s.x*mean.x
        
        slope.temporary[1,5:9]=summary(m1)[[20]][2,]
        slope.temporary[1,10]=summary(m1)[[20]][2,1]*s.y/s.x
        slope.temporary[1,11]=summary(m1)[[20]][2,2]*s.y/s.x
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      intercept = rbind(intercept, intercept.temporary)
      slope = rbind(slope, slope.temporary)
    }
    
    write.csv(intercept,"intercept arrival w year.csv")
    write.csv(slope,"slope arrival w year.csv")
  }
}

### calculate 95% intervals
## clear environment

taxa.all = list(c("Cicada","Cricket","Cuckoo","Frog","Goose","Swallow","Toad"),
                c("Cicada","Cricket","Frog","Goose","Swallow","Toad"),
                c("Cicada","Frog","Goose","Swallow","Toad"),
                c("Cicada","Frog","Goose","Swallow","Toad"))
m.all=c(5,6,10,11)
resample=1000

all.intercept = data.frame(matrix(NA,ncol=9,nrow=1))
names(all.intercept)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")
all.slope = data.frame(matrix(NA,ncol=9,nrow=1))
names(all.slope)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")

intercept=read.csv("intercept arrival w year.csv", header=TRUE)
slope=read.csv("slope arrival w year.csv", header=TRUE)

for (l in 1:4){
  print(l)
  
  taxa = taxa.all[[l]]
  taxa.number=length(taxa)

  for (t in 1:taxa.number){
    print(t)
        
    intercept.subset = na.omit(intercept[intercept$Criterion==m.all[[l]] & intercept$Taxon==taxa[t],])
    slope.subset = na.omit(slope[slope$Criterion==m.all[[l]] & slope$Taxon==taxa[t],])
    subsample.number = nrow(intercept.subset)
    
    all.intercept.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(all.intercept.temporary)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")
    all.slope.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(all.slope.temporary)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")
    
    all.intercept.temporary[1,1]=m.all[l]
    all.intercept.temporary[1,2]=taxa[[t]]
    all.intercept.temporary[1,3]=subsample.number
    all.intercept.temporary[1,4]=min(intercept.subset$sites)
    all.intercept.temporary[1,5]=max(intercept.subset$sites)
    all.slope.temporary[1,1]=m.all[l]
    all.slope.temporary[1,2]=taxa[[t]]
    all.slope.temporary[1,3]=subsample.number
    all.slope.temporary[1,4]=min(slope.subset$sites)
    all.slope.temporary[1,5]=max(slope.subset$sites)
    
    upper95=NA
    lower95=NA
    crosszero=NA
    test.1=rep(NA,subsample.number)
    test.2=rep(NA,subsample.number)
    test.1.update=NA
    test.2.update=NA
    
    tryCatch({
      for (a in 1:subsample.number){
        test.1[a] = sum(intercept.subset$beta[a]>intercept.subset$beta)/subsample.number#how much percentage of subsamples has values of intercept below the current subsample
        test.2[a] = sum(intercept.subset$beta[a]<intercept.subset$beta)/subsample.number#how much percentage of subsamples has values of intercept above the current subsample
      }# test.1: a vector of length subsample.number, of the above percentage for each subsample
    
      #test.1.update = test.1[which(test.1>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.1.update = test.1[which(test.1>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #intercept.beta.1 = intercept.subset$beta[which(test.1>0.9499)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      intercept.beta.1 = intercept.subset$beta[which(test.1>0.8999)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      if (length(intercept.beta.1)>0){
        upper95 = intercept.beta.1[which.min(test.1.update)]
      }else{
        upper95=max(intercept.subset$beta)
      }
            
      #test.2.update = test.2[which(test.2>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.2.update = test.2[which(test.2>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #intercept.beta.2 = intercept.subset$beta[which(test.2>0.9499)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      intercept.beta.2 = intercept.subset$beta[which(test.2>0.8999)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      if (length(intercept.beta.2)>0){
        lower95 = intercept.beta.2[which.min(test.2.update)]
      }else{
        lower95=min(intercept.subset$beta)
      }
      
      if(upper95>0){
        if(lower95<0){
          crosszero=1
        }else{
          crosszero=0
        }
      }else{
        crosszero=0
      }
      ## above: if 95% interval just touches zero, considered not crossing zero.
      
      all.intercept.temporary[1,6]=mean(intercept.subset$beta)
      all.intercept.temporary[1,7]=upper95
      all.intercept.temporary[1,8]=lower95
      all.intercept.temporary[1,9]=crosszero      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})      
    
    all.intercept=rbind(all.intercept,all.intercept.temporary)      
      
    upper95=NA
    lower95=NA
    crosszero=NA
    test.1=rep(NA,subsample.number)
    test.2=rep(NA,subsample.number)
    test.1.update=NA
    test.2.update=NA
    
    tryCatch({
      for (a in 1:subsample.number){
        test.1[a] = sum(slope.subset$beta[a]>slope.subset$beta)/subsample.number#how much percentage of subsamples has values of slope below the current subsample
        test.2[a] = sum(slope.subset$beta[a]<slope.subset$beta)/subsample.number#how much percentage of subsamples has values of slope above the current subsample
      }# test.1: a vector of length subsample.number, of the above percentage for each subsample
      #test.1.update = test.1[which(test.1>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.1.update = test.1[which(test.1>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #slope.beta.1 = slope.subset$beta[which(test.1>0.9499)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      slope.beta.1 = slope.subset$beta[which(test.1>0.8999)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      if (length(slope.beta.1)>0){
        upper95 = slope.beta.1[which.min(test.1.update)]
      }else{
        upper95=max(slope.subset$beta)
      }
      
      #test.2.update = test.2[which(test.2>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.2.update = test.2[which(test.2>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #slope.beta.2 = slope.subset$beta[which(test.2>0.9499)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      slope.beta.2 = slope.subset$beta[which(test.2>0.8999)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      if (length(slope.beta.2)>0){
        lower95 = slope.beta.2[which.min(test.2.update)]
      }else{
        lower95=min(slope.subset$beta)
      }
      
      if(upper95>0){
        if(lower95<0){
          crosszero=1
        }else{
          crosszero=0
        }
      }else{
        crosszero=0
      }
      
      all.slope.temporary[1,6]=mean(slope.subset$beta)
      all.slope.temporary[1,7]=upper95
      all.slope.temporary[1,8]=lower95
      all.slope.temporary[1,9]=crosszero 
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})    
      
    all.slope=rbind(all.slope,all.slope.temporary)  
  }
}

write.csv(all.intercept,"overall arrival intercept.csv")
write.csv(all.slope,"overall arrival slope.csv")



### departure

## change directory to the departure w year folder

taxa.all = list(c("Cicada","Frog","Goose","Swallow"),
                c("Cicada","Frog","Swallow"),
                c("Goose","Swallow"),
                c("Cicada","Frog","Swallow"))

m.all=c(5,6,10,11)
resample=1000

intercept = data.frame(matrix(NA,ncol=11,nrow=1))
names(intercept)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")
slope = data.frame(matrix(NA,ncol=11,nrow=1))
names(slope)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")

### running GLMMs
for (l in 1:4){
  print(l)
  
  taxa = taxa.all[[l]]
  taxa.number=length(taxa)
  
  for (t in 1:taxa.number){
    print(l)
    print(t)
    
    for (i in 1:resample){
      print(i)
      
      filename = paste(l,taxa[[t]],i, "subsample.csv")
      departure = read.csv(filename, header=TRUE)
      departure = departure[departure$whether_independent==1,]
      
      s.y=sd(departure$Departure)
      s.x=sd(departure$Year)
      mean.y = mean(departure$Departure)
      mean.x = mean(departure$Year)
      
      intercept.temporary = data.frame(matrix(NA,ncol=11,nrow=1))
      names(intercept.temporary)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")
      slope.temporary = data.frame(matrix(NA,ncol=11,nrow=1))
      names(slope.temporary)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")
      
      intercept.temporary[1,1]=m.all[l]
      intercept.temporary[1,2]=taxa[[t]]
      intercept.temporary[1,3]=i
      intercept.temporary[1,4]=length(unique(as.character(departure$SiteID)))
      
      slope.temporary[1,1]=m.all[l]
      slope.temporary[1,2]=taxa[[t]]
      slope.temporary[1,3]=i
      slope.temporary[1,4]=length(unique(as.character(departure$SiteID)))
      
      tryCatch({
        ## run random slope + intercept models (allows for correlation between random slope and intercept)
        m1 = lme(scale(Departure)~scale(Year), random = ~scale(Year) | SiteID, data = na.omit(departure))
        #m1 = lme(scale(departure)~scale(Year), random = ~scale(Year) | draw, data = na.omit(departure))
        
        intercept.temporary[1,5:9]=summary(m1)[[20]][1,]
        intercept.temporary[1,10]=summary(m1)[[20]][1,1]*s.y + mean.y - summary(m1)[[20]][2,1]*s.y/s.x*mean.x
        intercept.temporary[1,11]=summary(m1)[[20]][1,2]*s.y + mean.y - summary(m1)[[20]][2,2]*s.y/s.x*mean.x
        
        slope.temporary[1,5:9]=summary(m1)[[20]][2,]
        slope.temporary[1,10]=summary(m1)[[20]][2,1]*s.y/s.x
        slope.temporary[1,11]=summary(m1)[[20]][2,2]*s.y/s.x
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      intercept = rbind(intercept, intercept.temporary)
      slope = rbind(slope, slope.temporary)
    }
    
    write.csv(intercept,"intercept departure w year.csv")
    write.csv(slope,"slope departure w year.csv")
  }
}

### calculate 95% intervals
## clear environment

taxa.all = list(c("Cicada","Frog","Goose","Swallow"),
                c("Cicada","Frog","Swallow"),
                c("Goose","Swallow"),
                c("Cicada","Frog","Swallow"))

m.all=c(5,6,10,11)
resample=1000

all.intercept = data.frame(matrix(NA,ncol=9,nrow=1))
names(all.intercept)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")
all.slope = data.frame(matrix(NA,ncol=9,nrow=1))
names(all.slope)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")

intercept=read.csv("intercept departure w year.csv", header=TRUE)
slope=read.csv("slope departure w year.csv", header=TRUE)

for (l in 1:4){
  print(l)
  
  taxa = taxa.all[[l]]
  taxa.number=length(taxa)
  
  for (t in 1:taxa.number){
    print(t)
    
    intercept.subset = na.omit(intercept[intercept$Criterion==m.all[[l]] & intercept$Taxon==taxa[t],])
    slope.subset = na.omit(slope[slope$Criterion==m.all[[l]] & slope$Taxon==taxa[t],])
    subsample.number = nrow(intercept.subset)
    
    all.intercept.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(all.intercept.temporary)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")
    all.slope.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(all.slope.temporary)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")
    
    all.intercept.temporary[1,1]=m.all[l]
    all.intercept.temporary[1,2]=taxa[[t]]
    all.intercept.temporary[1,3]=subsample.number
    all.intercept.temporary[1,4]=min(intercept.subset$sites)
    all.intercept.temporary[1,5]=max(intercept.subset$sites)
    all.slope.temporary[1,1]=m.all[l]
    all.slope.temporary[1,2]=taxa[[t]]
    all.slope.temporary[1,3]=subsample.number
    all.slope.temporary[1,4]=min(slope.subset$sites)
    all.slope.temporary[1,5]=max(slope.subset$sites)
    
    upper95=NA
    lower95=NA
    crosszero=NA
    test.1=rep(NA,subsample.number)
    test.2=rep(NA,subsample.number)
    test.1.update=NA
    test.2.update=NA
    
    tryCatch({
      for (a in 1:subsample.number){
        test.1[a] = sum(intercept.subset$beta[a]>intercept.subset$beta)/subsample.number#how much percentage of subsamples has values of intercept below the current subsample
        test.2[a] = sum(intercept.subset$beta[a]<intercept.subset$beta)/subsample.number#how much percentage of subsamples has values of intercept above the current subsample
      }# test.1: a vector of length subsample.number, of the above percentage for each subsample
      
      #test.1.update = test.1[which(test.1>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.1.update = test.1[which(test.1>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #intercept.beta.1 = intercept.subset$beta[which(test.1>0.9499)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      intercept.beta.1 = intercept.subset$beta[which(test.1>0.8999)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      if (length(intercept.beta.1)>0){
        upper95 = intercept.beta.1[which.min(test.1.update)]
      }else{
        upper95=max(intercept.subset$beta)
      }
      
      #test.2.update = test.2[which(test.2>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.2.update = test.2[which(test.2>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #intercept.beta.2 = intercept.subset$beta[which(test.2>0.9499)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      intercept.beta.2 = intercept.subset$beta[which(test.2>0.8999)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      if (length(intercept.beta.2)>0){
        lower95 = intercept.beta.2[which.min(test.2.update)]
      }else{
        lower95=min(intercept.subset$beta)
      }
      
      if(upper95>0){
        if(lower95<0){
          crosszero=1
        }else{
          crosszero=0
        }
      }else{
        crosszero=0
      }
      ## above: if 95% interval just touches zero, considered not crossing zero.
      
      all.intercept.temporary[1,6]=mean(intercept.subset$beta)
      all.intercept.temporary[1,7]=upper95
      all.intercept.temporary[1,8]=lower95
      all.intercept.temporary[1,9]=crosszero      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})      
    
    all.intercept=rbind(all.intercept,all.intercept.temporary)      
    
    upper95=NA
    lower95=NA
    crosszero=NA
    test.1=rep(NA,subsample.number)
    test.2=rep(NA,subsample.number)
    test.1.update=NA
    test.2.update=NA
    
    tryCatch({
      for (a in 1:subsample.number){
        test.1[a] = sum(slope.subset$beta[a]>slope.subset$beta)/subsample.number#how much percentage of subsamples has values of slope below the current subsample
        test.2[a] = sum(slope.subset$beta[a]<slope.subset$beta)/subsample.number#how much percentage of subsamples has values of slope above the current subsample
      }# test.1: a vector of length subsample.number, of the above percentage for each subsample
      #test.1.update = test.1[which(test.1>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.1.update = test.1[which(test.1>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #slope.beta.1 = slope.subset$beta[which(test.1>0.9499)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      slope.beta.1 = slope.subset$beta[which(test.1>0.8999)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      if (length(slope.beta.1)>0){
        upper95 = slope.beta.1[which.min(test.1.update)]
      }else{
        upper95=max(slope.subset$beta)
      }
      
      #test.2.update = test.2[which(test.2>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.2.update = test.2[which(test.2>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #slope.beta.2 = slope.subset$beta[which(test.2>0.9499)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      slope.beta.2 = slope.subset$beta[which(test.2>0.8999)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      if (length(slope.beta.2)>0){
        lower95 = slope.beta.2[which.min(test.2.update)]
      }else{
        lower95=min(slope.subset$beta)
      }
      
      if(upper95>0){
        if(lower95<0){
          crosszero=1
        }else{
          crosszero=0
        }
      }else{
        crosszero=0
      }
      
      all.slope.temporary[1,6]=mean(slope.subset$beta)
      all.slope.temporary[1,7]=upper95
      all.slope.temporary[1,8]=lower95
      all.slope.temporary[1,9]=crosszero 
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})    
    
    all.slope=rbind(all.slope,all.slope.temporary)  
  }
}

write.csv(all.intercept,"overall departure intercept.csv")
write.csv(all.slope,"overall departure slope.csv")


### stay

## change directory to the stay w year folder

taxa.all = list(c("Cicada","Frog","Goose","Swallow"),
                c("Cicada","Cricket","Frog","Goose","Swallow"),
                c("Goose","Swallow","Toad"),
                c("Frog","Goose","Swallow"))

m.all=c(5,6,10,11)
resample=1000

intercept = data.frame(matrix(NA,ncol=11,nrow=1))
names(intercept)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")
slope = data.frame(matrix(NA,ncol=11,nrow=1))
names(slope)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")

### running GLMMs
for (l in 1:4){
  print(l)
  
  taxa = taxa.all[[l]]
  taxa.number=length(taxa)
  
  for (t in 1:taxa.number){
    print(l)
    print(t)
    
    for (i in 1:resample){
      print(i)
      
      filename = paste(l,taxa[[t]],i, "subsample.csv")
      stay = read.csv(filename, header=TRUE)
      stay = stay[stay$whether_independent==1,]
      
      s.y=sd(stay$Stay)
      s.x=sd(stay$Year)
      mean.y = mean(stay$Stay)
      mean.x = mean(stay$Year)
      
      intercept.temporary = data.frame(matrix(NA,ncol=11,nrow=1))
      names(intercept.temporary)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")
      slope.temporary = data.frame(matrix(NA,ncol=11,nrow=1))
      names(slope.temporary)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")
      
      intercept.temporary[1,1]=m.all[l]
      intercept.temporary[1,2]=taxa[[t]]
      intercept.temporary[1,3]=i
      intercept.temporary[1,4]=length(unique(as.character(stay$SiteID)))
      
      slope.temporary[1,1]=m.all[l]
      slope.temporary[1,2]=taxa[[t]]
      slope.temporary[1,3]=i
      slope.temporary[1,4]=length(unique(as.character(stay$SiteID)))
      
      tryCatch({
        ## run random slope + intercept models (allows for correlation between random slope and intercept)
        m1 = lme(scale(Stay)~scale(Year), random = ~scale(Year) | SiteID, data = na.omit(stay))
        #m1 = lme(scale(stay)~scale(Year), random = ~scale(Year) | draw, data = na.omit(stay))
        
        intercept.temporary[1,5:9]=summary(m1)[[20]][1,]
        intercept.temporary[1,10]=summary(m1)[[20]][1,1]*s.y + mean.y - summary(m1)[[20]][2,1]*s.y/s.x*mean.x
        intercept.temporary[1,11]=summary(m1)[[20]][1,2]*s.y + mean.y - summary(m1)[[20]][2,2]*s.y/s.x*mean.x
        
        slope.temporary[1,5:9]=summary(m1)[[20]][2,]
        slope.temporary[1,10]=summary(m1)[[20]][2,1]*s.y/s.x
        slope.temporary[1,11]=summary(m1)[[20]][2,2]*s.y/s.x
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      intercept = rbind(intercept, intercept.temporary)
      slope = rbind(slope, slope.temporary)
    }
    
    write.csv(intercept,"intercept stay w year.csv")
    write.csv(slope,"slope stay w year.csv")
  }
}

### calculate 95% intervals
## clear environment

taxa.all = list(c("Cicada","Frog","Goose","Swallow"),
                c("Cicada","Cricket","Frog","Goose","Swallow"),
                c("Goose","Swallow","Toad"),
                c("Frog","Goose","Swallow"))

m.all=c(5,6,10,11)
resample=1000

all.intercept = data.frame(matrix(NA,ncol=9,nrow=1))
names(all.intercept)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")
all.slope = data.frame(matrix(NA,ncol=9,nrow=1))
names(all.slope)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")

intercept=read.csv("intercept stay w year.csv", header=TRUE)
slope=read.csv("slope stay w year.csv", header=TRUE)

for (l in 1:4){
  print(l)
  
  taxa = taxa.all[[l]]
  taxa.number=length(taxa)
  
  for (t in 1:taxa.number){
    print(t)
    
    intercept.subset = na.omit(intercept[intercept$Criterion==m.all[[l]] & intercept$Taxon==taxa[t],])
    slope.subset = na.omit(slope[slope$Criterion==m.all[[l]] & slope$Taxon==taxa[t],])
    subsample.number = nrow(intercept.subset)
    
    all.intercept.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(all.intercept.temporary)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")
    all.slope.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(all.slope.temporary)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")
    
    all.intercept.temporary[1,1]=m.all[l]
    all.intercept.temporary[1,2]=taxa[[t]]
    all.intercept.temporary[1,3]=subsample.number
    all.intercept.temporary[1,4]=min(intercept.subset$sites)
    all.intercept.temporary[1,5]=max(intercept.subset$sites)
    
    all.slope.temporary[1,1]=m.all[l]
    all.slope.temporary[1,2]=taxa[[t]]
    all.slope.temporary[1,3]=subsample.number
    all.slope.temporary[1,4]=min(slope.subset$sites)
    all.slope.temporary[1,5]=max(slope.subset$sites)
    
    upper95=NA
    lower95=NA
    crosszero=NA
    test.1=rep(NA,subsample.number)
    test.2=rep(NA,subsample.number)
    test.1.update=NA
    test.2.update=NA
    
    tryCatch({
      for (a in 1:subsample.number){
        test.1[a] = sum(intercept.subset$beta[a]>intercept.subset$beta)/subsample.number#how much percentage of subsamples has values of intercept below the current subsample
        test.2[a] = sum(intercept.subset$beta[a]<intercept.subset$beta)/subsample.number#how much percentage of subsamples has values of intercept above the current subsample
      }# test.1: a vector of length subsample.number, of the above percentage for each subsample
      
      #test.1.update = test.1[which(test.1>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.1.update = test.1[which(test.1>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #intercept.beta.1 = intercept.subset$beta[which(test.1>0.9499)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      intercept.beta.1 = intercept.subset$beta[which(test.1>0.8999)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      if (length(intercept.beta.1)>0){
        upper95 = intercept.beta.1[which.min(test.1.update)]
      }else{
        upper95=max(intercept.subset$beta)
      }
      
      #test.2.update = test.2[which(test.2>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.2.update = test.2[which(test.2>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #intercept.beta.2 = intercept.subset$beta[which(test.2>0.9499)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      intercept.beta.2 = intercept.subset$beta[which(test.2>0.8999)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      if (length(intercept.beta.2)>0){
        lower95 = intercept.beta.2[which.min(test.2.update)]
      }else{
        lower95=min(intercept.subset$beta)
      }
      
      if(upper95>0){
        if(lower95<0){
          crosszero=1
        }else{
          crosszero=0
        }
      }else{
        crosszero=0
      }
      ## above: if 95% interval just touches zero, considered not crossing zero.
      
      all.intercept.temporary[1,6]=mean(intercept.subset$beta)
      all.intercept.temporary[1,7]=upper95
      all.intercept.temporary[1,8]=lower95
      all.intercept.temporary[1,9]=crosszero      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})      
    
    all.intercept=rbind(all.intercept,all.intercept.temporary)      
    
    upper95=NA
    lower95=NA
    crosszero=NA
    test.1=rep(NA,subsample.number)
    test.2=rep(NA,subsample.number)
    test.1.update=NA
    test.2.update=NA
    
    tryCatch({
      for (a in 1:subsample.number){
        test.1[a] = sum(slope.subset$beta[a]>slope.subset$beta)/subsample.number#how much percentage of subsamples has values of slope below the current subsample
        test.2[a] = sum(slope.subset$beta[a]<slope.subset$beta)/subsample.number#how much percentage of subsamples has values of slope above the current subsample
      }# test.1: a vector of length subsample.number, of the above percentage for each subsample
      #test.1.update = test.1[which(test.1>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.1.update = test.1[which(test.1>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #slope.beta.1 = slope.subset$beta[which(test.1>0.9499)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      slope.beta.1 = slope.subset$beta[which(test.1>0.8999)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      if (length(slope.beta.1)>0){
        upper95 = slope.beta.1[which.min(test.1.update)]
      }else{
        upper95=max(slope.subset$beta)
      }
      
      #test.2.update = test.2[which(test.2>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.2.update = test.2[which(test.2>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #slope.beta.2 = slope.subset$beta[which(test.2>0.9499)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      slope.beta.2 = slope.subset$beta[which(test.2>0.8999)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      if (length(slope.beta.2)>0){
        lower95 = slope.beta.2[which.min(test.2.update)]
      }else{
        lower95=min(slope.subset$beta)
      }
      
      if(upper95>0){
        if(lower95<0){
          crosszero=1
        }else{
          crosszero=0
        }
      }else{
        crosszero=0
      }
      
      all.slope.temporary[1,6]=mean(slope.subset$beta)
      all.slope.temporary[1,7]=upper95
      all.slope.temporary[1,8]=lower95
      all.slope.temporary[1,9]=crosszero 
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})    
    
    all.slope=rbind(all.slope,all.slope.temporary)  
  }
}

write.csv(all.intercept,"overall stay intercept.csv")
write.csv(all.slope,"overall stay slope.csv")


### overlap

## change directory to the overlap w year folder

taxa.all = list(c("Cicada_Cricket","Cicada_Cuckoo","Cicada_Swallow","Frog_Swallow","Frog_Toad","Goose_Swallow"),
                c("Cicada_Swallow","Cricket_Frog","Frog_Swallow","Frog_Toad","Goose_Swallow"),
                c("Cicada_Swallow","Frog_Swallow"))

m.all=c(5,6,11)
resample=1000

intercept = data.frame(matrix(NA,ncol=11,nrow=1))
names(intercept)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")
slope = data.frame(matrix(NA,ncol=11,nrow=1))
names(slope)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")

### running GLMMs
for (l in 1:3){
  print(l)
  
  taxa = taxa.all[[l]]
  taxa.number=length(taxa)
  
  for (t in 1:taxa.number){
    print(l)
    print(t)
    
    for (i in 1:resample){
      print(i)
      
      filename = paste(l,taxa[[t]],i, "subsample.csv")
      overlap = read.csv(filename, header=TRUE)
      overlap = overlap[overlap$whether_independent==1,]
      
      s.y=sd(overlap$Overlap)
      s.x=sd(overlap$Year)
      mean.y = mean(overlap$Overlap)
      mean.x = mean(overlap$Year)
      
      intercept.temporary = data.frame(matrix(NA,ncol=11,nrow=1))
      names(intercept.temporary)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")
      slope.temporary = data.frame(matrix(NA,ncol=11,nrow=1))
      names(slope.temporary)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")
      
      intercept.temporary[1,1]=m.all[l]
      intercept.temporary[1,2]=taxa[[t]]
      intercept.temporary[1,3]=i
      intercept.temporary[1,4]=length(unique(as.character(overlap$SiteID)))
      
      slope.temporary[1,1]=m.all[l]
      slope.temporary[1,2]=taxa[[t]]
      slope.temporary[1,3]=i
      slope.temporary[1,4]=length(unique(as.character(overlap$SiteID)))
      
      tryCatch({
        ## run random slope + intercept models (allows for correlation between random slope and intercept)
        m1 = lme(scale(Overlap)~scale(Year), random = ~scale(Year) | SiteID, data = na.omit(overlap))
        #m1 = lme(scale(overlap)~scale(Year), random = ~scale(Year) | draw, data = na.omit(overlap))
        
        intercept.temporary[1,5:9]=summary(m1)[[20]][1,]
        intercept.temporary[1,10]=summary(m1)[[20]][1,1]*s.y + mean.y - summary(m1)[[20]][2,1]*s.y/s.x*mean.x
        intercept.temporary[1,11]=summary(m1)[[20]][1,2]*s.y + mean.y - summary(m1)[[20]][2,2]*s.y/s.x*mean.x
        
        slope.temporary[1,5:9]=summary(m1)[[20]][2,]
        slope.temporary[1,10]=summary(m1)[[20]][2,1]*s.y/s.x
        slope.temporary[1,11]=summary(m1)[[20]][2,2]*s.y/s.x
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      intercept = rbind(intercept, intercept.temporary)
      slope = rbind(slope, slope.temporary)
    }
    
    write.csv(intercept,"intercept overlap w year.csv")
    write.csv(slope,"slope overlap w year.csv")
  }
}

### calculate 95% intervals
## clear environment

taxa.all = list(c("Cicada_Cricket","Cicada_Cuckoo","Cicada_Swallow","Frog_Swallow","Frog_Toad","Goose_Swallow"),
                c("Cicada_Swallow","Cricket_Frog","Frog_Swallow","Frog_Toad","Goose_Swallow"),
                c("Cicada_Swallow","Frog_Swallow"))

m.all=c(5,6,11)
resample=1000

all.intercept = data.frame(matrix(NA,ncol=9,nrow=1))
names(all.intercept)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")
all.slope = data.frame(matrix(NA,ncol=9,nrow=1))
names(all.slope)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")

intercept=read.csv("intercept overlap w year.csv", header=TRUE)
slope=read.csv("slope overlap w year.csv", header=TRUE)

for (l in 1:3){
  print(l)
  
  taxa = taxa.all[[l]]
  taxa.number=length(taxa)
  
  for (t in 1:taxa.number){
    print(t)
    
    intercept.subset = na.omit(intercept[intercept$Criterion==m.all[[l]] & intercept$Taxon==taxa[t],])
    slope.subset = na.omit(slope[slope$Criterion==m.all[[l]] & slope$Taxon==taxa[t],])
    subsample.number = nrow(intercept.subset)
    
    all.intercept.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(all.intercept.temporary)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")
    all.slope.temporary = data.frame(matrix(NA,ncol=9,nrow=1))
    names(all.slope.temporary)=c("Criterion","Taxon","subsamples","min_sites","max_sites","beta","upper95","lower95","crosszero")
    
    all.intercept.temporary[1,1]=m.all[l]
    all.intercept.temporary[1,2]=taxa[[t]]
    all.intercept.temporary[1,3]=subsample.number
    all.intercept.temporary[1,4]=min(intercept.subset$sites)
    all.intercept.temporary[1,5]=max(intercept.subset$sites)
    
    all.slope.temporary[1,1]=m.all[l]
    all.slope.temporary[1,2]=taxa[[t]]
    all.slope.temporary[1,3]=subsample.number
    all.slope.temporary[1,4]=min(slope.subset$sites)
    all.slope.temporary[1,5]=max(slope.subset$sites)
    
    upper95=NA
    lower95=NA
    crosszero=NA
    test.1=rep(NA,subsample.number)
    test.2=rep(NA,subsample.number)
    test.1.update=NA
    test.2.update=NA
    
    tryCatch({
      for (a in 1:subsample.number){
        test.1[a] = sum(intercept.subset$beta[a]>intercept.subset$beta)/subsample.number#how much percentage of subsamples has values of intercept below the current subsample
        test.2[a] = sum(intercept.subset$beta[a]<intercept.subset$beta)/subsample.number#how much percentage of subsamples has values of intercept above the current subsample
      }# test.1: a vector of length subsample.number, of the above percentage for each subsample
      
      #test.1.update = test.1[which(test.1>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.1.update = test.1[which(test.1>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #intercept.beta.1 = intercept.subset$beta[which(test.1>0.9499)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      intercept.beta.1 = intercept.subset$beta[which(test.1>0.8999)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      if (length(intercept.beta.1)>0){
        upper95 = intercept.beta.1[which.min(test.1.update)]
      }else{
        upper95=max(intercept.subset$beta)
      }
      
      #test.2.update = test.2[which(test.2>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.2.update = test.2[which(test.2>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #intercept.beta.2 = intercept.subset$beta[which(test.2>0.9499)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      intercept.beta.2 = intercept.subset$beta[which(test.2>0.8999)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      if (length(intercept.beta.2)>0){
        lower95 = intercept.beta.2[which.min(test.2.update)]
      }else{
        lower95=min(intercept.subset$beta)
      }
      
      if(upper95>0){
        if(lower95<0){
          crosszero=1
        }else{
          crosszero=0
        }
      }else{
        crosszero=0
      }
      ## above: if 95% interval just touches zero, considered not crossing zero.
      
      all.intercept.temporary[1,6]=mean(intercept.subset$beta)
      all.intercept.temporary[1,7]=upper95
      all.intercept.temporary[1,8]=lower95
      all.intercept.temporary[1,9]=crosszero      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})      
    
    all.intercept=rbind(all.intercept,all.intercept.temporary)      
    
    upper95=NA
    lower95=NA
    crosszero=NA
    test.1=rep(NA,subsample.number)
    test.2=rep(NA,subsample.number)
    test.1.update=NA
    test.2.update=NA
    
    tryCatch({
      for (a in 1:subsample.number){
        test.1[a] = sum(slope.subset$beta[a]>slope.subset$beta)/subsample.number#how much percentage of subsamples has values of slope below the current subsample
        test.2[a] = sum(slope.subset$beta[a]<slope.subset$beta)/subsample.number#how much percentage of subsamples has values of slope above the current subsample
      }# test.1: a vector of length subsample.number, of the above percentage for each subsample
      #test.1.update = test.1[which(test.1>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.1.update = test.1[which(test.1>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #slope.beta.1 = slope.subset$beta[which(test.1>0.9499)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      slope.beta.1 = slope.subset$beta[which(test.1>0.8999)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      if (length(slope.beta.1)>0){
        upper95 = slope.beta.1[which.min(test.1.update)]
      }else{
        upper95=max(slope.subset$beta)
      }
      
      #test.2.update = test.2[which(test.2>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      test.2.update = test.2[which(test.2>0.8999)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      #slope.beta.2 = slope.subset$beta[which(test.2>0.9499)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      slope.beta.2 = slope.subset$beta[which(test.2>0.8999)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      if (length(slope.beta.2)>0){
        lower95 = slope.beta.2[which.min(test.2.update)]
      }else{
        lower95=min(slope.subset$beta)
      }
      
      if(upper95>0){
        if(lower95<0){
          crosszero=1
        }else{
          crosszero=0
        }
      }else{
        crosszero=0
      }
      
      all.slope.temporary[1,6]=mean(slope.subset$beta)
      all.slope.temporary[1,7]=upper95
      all.slope.temporary[1,8]=lower95
      all.slope.temporary[1,9]=crosszero 
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})    
    
    all.slope=rbind(all.slope,all.slope.temporary)  
  }
}

write.csv(all.intercept,"overall overlap intercept.csv")
write.csv(all.slope,"overall overlap slope.csv")





#### old version, not run: combining synthesis in the same loop as individual models: 
for (l in 1:4){
  print(l)
  
  taxa = taxa.all[[l]]
  taxa.number=length(taxa)
  
  for (t in 1:taxa.number){
    print(t)
    
    beginning = r+1 # indexing: the row in the intercept and slope file in which a taxon begins
    
    for (i in 1:resample){
      print(i)
      
      filename = paste(l,taxa[[t]],i, "subsample.csv")
      arrival = read.csv(filename, header=TRUE)
      arrival = arrival[arrival$whether_independent==1,]
      
      s.y=sd(arrival$Arrival)
      s.x=sd(arrival$Year)
      mean.y = mean(arrival$Arrival)
      mean.x = mean(arrival$Year)
      
      intercept.temporary = data.frame(matrix(NA,ncol=11,nrow=1))
      names(intercept.temporary)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")
      slope.temporary = data.frame(matrix(NA,ncol=11,nrow=1))
      names(slope.temporary)=c("Criterion","Taxon","subsample","sites","beta.std","SE.std","DF","t-value","p-value","beta","SE")
      
      intercept.temporary[1,1]=m.all[l]
      intercept.temporary[1,2]=taxa[[t]]
      intercept.temporary[1,3]=i
      intercept.temporary[1,4]=length(unique(as.character(arrival$SiteID)))
      
      slope.temporary[1,1]=m.all[l]
      slope.temporary[1,2]=taxa[[t]]
      slope.temporary[1,3]=i
      slope.temporary[1,4]=length(unique(as.character(arrival$SiteID)))
      
      tryCatch({
        ## run random slope + intercept models (allows for correlation between random slope and intercept)
        m1 = lme(scale(Arrival)~scale(Year), random = ~scale(Year) | SiteID, data = na.omit(arrival))
        #m1 = lme(scale(Arrival)~scale(Year), random = ~scale(Year) | draw, data = na.omit(arrival))
        
        intercept.temporary[1,5:9]=summary(m1)[[20]][1,]
        intercept.temporary[1,10]=summary(m1)[[20]][1,1]*s.y + mean.y - summary(m1)[[20]][2,1]*s.y/s.x*mean.x
        intercept.temporary[1,11]=summary(m1)[[20]][1,2]*s.y + mean.y - summary(m1)[[20]][2,2]*s.y/s.x*mean.x
        
        slope.temporary[1,5:9]=summary(m1)[[20]][2,]
        slope.temporary[1,10]=summary(m1)[[20]][2,1]*s.y/s.x
        slope.temporary[1,11]=summary(m1)[[20]][2,2]*s.y/s.x
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      intercept = rbind(intercept, intercept.temporary)
      slope = rbind(slope, slope.temporary)
      
      r=r+1 # indexing: the number of runs the loop has been executing
    }
    
    write.csv(intercept,"intercept arrival w year.csv")
    write.csv(slope,"slope arrival w year.csv")
    
    end = r # indexing: the row in the intercept and slope file in which a taxon ends
    
    intercept.subset = na.omit(intercept[beginning:end,])
    slope.subset = na.omit(slope[beginning:end,])
    subsample.number = nrow(intercept.subset)
    
    all.intercept.temporary = data.frame(matrix(NA,ncol=7,nrow=1))
    names(all.intercept.temporary)=c("Criterion","Taxon","subsamples","beta","upper95","lower95","crosszero")
    all.slope.temporary = data.frame(matrix(NA,ncol=7,nrow=1))
    names(all.slope.temporary)=c("Criterion","Taxon","subsamples","beta","upper95","lower95","crosszero")
    
    all.intercept.temporary[1,1]=m.all[l]
    all.intercept.temporary[1,2]=taxa[[t]]
    all.intercept.temporary[1,3]=subsample.number
    all.slope.temporary[1,1]=m.all[l]
    all.slope.temporary[1,2]=taxa[[t]]
    all.slope.temporary[1,3]=subsample.number
    
    upper95=NA
    lower95=NA
    crosszero=NA
    test.1=rep(NA,subsample.number)
    test.2=rep(NA,subsample.number)
    
    tryCatch({
      for (a in 1:subsample.number){
        test.1[a] = sum(intercept.subset$beta[a]>intercept.subset$beta)/subsample.number#how much percentage of subsamples has values of intercept below the current subsample
        test.2[a] = sum(intercept.subset$beta[a]<intercept.subset$beta)/subsample.number#how much percentage of subsamples has values of intercept above the current subsample
      }# test.1: a vector of length subsample.number, of the above percentage for each subsample
      
      test.1.update = test.1[which(test.1>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      intercept.beta.1 = intercept.subset$beta[which(test.1>0.9499)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      upper95 = intercept.beta.1[which.min(test.1.update)]
      
      test.2.update = test.2[which(test.2>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      intercept.beta.2 = intercept.subset$beta[which(test.2>0.9499)]#a reduced vector of intercept values that satisfy the above 0.9499 criterion
      lower95 = intercept.beta.2[which.min(test.2.update)]
      
      if(upper95>0){
        if(lower95<0){
          crosszero=1
        }else{
          crosszero=0
        }
      }else{
        crosszero=0
      }
      
      all.intercept.temporary[1,4]=mean(intercept.subset$beta)
      all.intercept.temporary[1,5]=upper95
      all.intercept.temporary[1,6]=lower95
      all.intercept.temporary[1,7]=crosszero      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})      
    
    all.intercept=rbind(all.intercept,all.intercept.temporary)      
    
    upper95=NA
    lower95=NA
    crosszero=NA
    test.1=rep(NA,subsample.number)
    test.2=rep(NA,subsample.number)
    
    tryCatch({
      for (a in 1:subsample.number){
        test.1[a] = sum(slope.subset$beta[a]>slope.subset$beta)/subsample.number#how much percentage of subsamples has values of slope below the current subsample
        test.2[a] = sum(slope.subset$beta[a]<slope.subset$beta)/subsample.number#how much percentage of subsamples has values of slope above the current subsample
      }# test.1: a vector of length subsample.number, of the above percentage for each subsample
      test.1.update = test.1[which(test.1>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      slope.beta.1 = slope.subset$beta[which(test.1>0.9499)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      upper95 = slope.beta.1[which.min(test.1.update)]
      
      test.2.update = test.2[which(test.2>0.9499)]#a reduced vector from test.1 that includes only those components that are above 0.9499
      slope.beta.2 = slope.subset$beta[which(test.2>0.9499)]#a reduced vector of slope values that satisfy the above 0.9499 criterion
      lower95 = slope.beta.2[which.min(test.2.update)]
      
      if(upper95>0){
        if(lower95<0){
          crosszero=1
        }else{
          crosszero=0
        }
      }else{
        crosszero=0
      }
      all.slope.temporary[1,4]=mean(slope.subset$beta)
      all.slope.temporary[1,5]=upper95
      all.slope.temporary[1,6]=lower95
      all.slope.temporary[1,7]=crosszero 
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})    
    
    all.slope=rbind(all.slope,all.slope.temporary)  
    
    write.csv(all.intercept,"intercept w year, across subsamples.csv")
    write.csv(all.slope,"slope w year, across subsamples.csv")
  }
}
