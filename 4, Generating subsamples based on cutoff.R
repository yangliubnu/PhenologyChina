#### For data with spatial autocorrelation issues, generating subsamples based on cut-off distance

library(sp)
library(boot)
library(nlme)

### For initial occurrence data during breeding season (initial occurrence represented as "arrival" in data set)
arrival.breeding.date  = read.csv('Breedingarrival.csv', header=T)
arrival.breeding.date$SiteID = as.factor(arrival.breeding.date $SiteID)
arrival.breeding.date$Year = as.numeric(arrival.breeding.date $Year)
head(arrival.breeding.date)

resample=1000
taxa.all = list(c("Cicada","Cricket","Cuckoo","Frog","Goose","Swallow","Toad"),
                c("Cicada","Cricket","Frog","Goose","Swallow","Toad"),
                c("Cicada","Frog","Goose","Swallow","Toad"),
                c("Cicada","Frog","Goose","Swallow","Toad"))
cutoff.all=list(c(800,200,200,1200,600,1400,200),
                c(800,200,1200,600,1400,200),
                c(200,1000,600,1400,200),
                c(800,1000,600,1400,200))
m.all=c(5,6,10,11)
subsample=list(NA)

book.keeping=as.data.frame(matrix(NA,ncol=8,nrow=1))
names(book.keeping)=c("Criterion","Taxon","OriginalSampleSize","Subsample","IndependentSampleSize","AllSampleSize_post","SatisfyMinimumDistance","SatisfyDataCompletion")

for (l in 1:4){# the large loop of criteria
  print(l)
  taxa = taxa.all[[l]]
  cutoff.suball = cutoff.all[[l]]
  m=m.all[l]
  taxa.number=length(taxa)
  subsample[[l]]=list(NA)
  
  for (t in 1:taxa.number){#looping through each taxa
    print(t)
    subsample[[l]][[t]]=list(NA)
    cutoff=cutoff.suball[[t]]
    
    # information about the original dataset that is not going to be updated in subsequent loops
    data=arrival.breeding.date[arrival.breeding.date$AnalysisGroup==taxa[[t]] & arrival.breeding.date$m==m,]
    site.list=unique(data$SiteID)#the list of viable site list - to be updated in each j loop
    site.number=length(site.list) # the number of all sites - remains unchanged
    site.geo=as.matrix(data[!duplicated(data[,5:6]),5:6]) #the geo location of all sites, keeping only lat long
    geo.dist=spDists(x=site.geo,longlat=TRUE)#a matrix of distance between all point pairs, to be referenced to in checking distanc
    
    for (i in 1:resample){
      print(i)
      
      pick=rep(NA,site.number)#an empty vector to store all the picks
      
      # the first draw
      pick[1]=sample(x=1:site.number,size=1)#generate 1st "resample" draw - to be updated in each j loop
      subsample[[l]][[t]][[i]]=data[data$SiteID==site.list[pick[1]],]#put data from the first draw into sub - to be udpated in each j loop
      # add an identifier to note this is data from draws
      identifier = data.frame(matrix(NA,ncol=3,nrow=nrow(subsample[[l]][[t]][[i]])));names(identifier)=c("draw","whether_independent","cluster")
      identifier[,1]=1
      identifier[,2]=1
      identifier[,3]=site.list[pick[1]]
      subsample[[l]][[t]][[i]]=cbind(subsample[[l]][[t]][[i]],identifier)
      
      distance=geo.dist[pick[1],]#the distance of all sites from this first draw - to be updated in each j loop
      
      # sites dropped from this first draw
      pool.drop=site.list[which(distance<=cutoff & distance>0)]# the list of sites viable for further draw - to be updated in each j loop
      pool.size.drop = length(pool.drop)
      
      if(pool.size.drop>0){
        data.drop=data[data$SiteID==pool.drop[1],]# a stub of data.drop to be potentially added on
        if (pool.size.drop>1){
          for(p in 1:(pool.size.drop-1)){
            data.drop=rbind(data.drop,data[data$SiteID==pool.drop[1+p],])
          }
        }else{
          data.drop=data.drop
        }
          
        identifier.drop = data.frame(matrix(NA,ncol=3,nrow=nrow(data.drop)));names(identifier.drop)=c("draw","whether_independent","cluster")
        identifier.drop[,1]=1
        identifier.drop[,2]=0
        identifier.drop[,3]=site.list[pick[1]]
        
        data.drop=cbind(data.drop,identifier.drop)
        
        # site picked and the ones dropped
        subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],data.drop)
      }else{
        subsample[[l]][[t]][[i]]=subsample[[l]][[t]][[i]]
      }
      
      # viable sites after the first draw
      pool=site.list[which(distance>cutoff)]# the list of sites viable for further draw - to be updated in each j loop
      pool.size=length(pool)# the number of sites viable for further draw - to be updated in each j loop
      pool.dist=geo.dist[which(distance>cutoff),which(distance>cutoff)]# the distance matrix only for sites kept after applying the cutoff - to be updated in each j loop
            
      # subsequent draws based on these viable sites - viable sites are to be updated after each draw
      for(j in 1:(site.number-1)){#the max number of samples that could be drawn into the subsample
        print(j)
        if (pool.size > 1) {
          # a new draw
          pick[1+j]=sample(x=1:pool.size,size=1)#generate "resample"-number of start sites (in the form of its numbering postion)
          subsample.temp=data[data$SiteID==pool[pick[1+j]],]
          identifier = data.frame(matrix(NA,ncol=3,nrow=nrow(subsample.temp)));names(identifier)=c("draw","whether_independent","cluster")
          identifier[,1]=1+j
          identifier[,2]=1
          identifier[,3]=pool[pick[1+j]]
          subsample.temp=cbind(subsample.temp,identifier)
          
          subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],subsample.temp)
                    
          distance=pool.dist[pick[1+j],]
          
          # sites dropped from this previous draw
          pool.drop=pool[which(distance<=cutoff & distance>0)]# the list of sites viable for further draw - to be updated in each j loop
          pool.size.drop = length(pool.drop)
          
          if(pool.size.drop>0){
            data.drop=data[data$SiteID==pool.drop[1],]# a stub of data.drop to be potentially added on
            if (pool.size.drop>1){
              for(p in 1:(pool.size.drop-1)){
                data.drop=rbind(data.drop,data[data$SiteID==pool.drop[1+p],])
              }
            }else{
              data.drop=data.drop
            }
            identifier.drop = data.frame(matrix(NA,ncol=3,nrow=nrow(data.drop)));names(identifier.drop)=c("draw","whether_independent","cluster")
            identifier.drop[,1]=1+j
            identifier.drop[,2]=0
            identifier.drop[,3]=pool[pick[1+j]]
            
            data.drop=cbind(data.drop,identifier.drop)
            
            # site picked and the ones dropped
            subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],data.drop)
          }else{
            subsample[[l]][[t]][[i]]=subsample[[l]][[t]][[i]]
          }
          
          # viable sites after the previous draw
          pool=pool[which(distance>cutoff)]
          pool.size=length(pool)# the number of sites viable for further draw - to be updated in each j loop
          pool.dist=pool.dist[which(distance>cutoff),which(distance>cutoff)]# the distance matrix only for sites kept after applying the cutoff - to be updated in each j loop
        }else{
          if(pool.size==1){
            # the only possible new draw
            pick[1+j]=1
            
            subsample.temp=data[data$SiteID==pool[pick[1+j]],]
            identifier = data.frame(matrix(NA,ncol=3,nrow=nrow(subsample.temp)));names(identifier)=c("draw","whether_independent","cluster")
            identifier[,1]=1+j
            identifier[,2]=1
            identifier[,3]=pool[pick[1+j]]
            subsample.temp=cbind(subsample.temp,identifier)
            
            subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],subsample.temp)
            
            pool.size=0
          }else{
            subsample[[l]][[t]][[i]]=subsample[[l]][[t]][[i]]
          }
        }  
      }
      
      filename = paste(l,taxa[[t]],i, "subsample.csv")
      write.csv(subsample[[l]][[t]][[i]], filename)
      
      # double check if sites too close to each other have been elimiated
      independent = subsample[[l]][[t]][[i]][subsample[[l]][[t]][[i]]$whether_independent==1,]
      independent.number=length(unique(independent$SiteID))
      new.site.geo=as.matrix(independent[!duplicated(independent[,5:6]),5:6]) #the geo location of all sites, keeping only lat long
      new.geo.dist=spDists(x=new.site.geo,longlat=TRUE)#a matrix of distance between all point pairs, to be referenced to in checking distanc
      
      if(sum(new.geo.dist>cutoff)==(dim(new.geo.dist)[1]^2-dim(new.geo.dist)[1])){#a logical test whether all distances except for the diagonals are greater than cutoff)
        check.independent=1
      }else{
        check.independent=0
      }
      
      # double check if the produced data set is of the same size as the original data set
      if(length(unique(data$SiteID))==length(unique(subsample[[l]][[t]][[i]]$SiteID))){
        check.complete=1
      }else{
        check.complete=0
      }
      
      book.keeping.temporary=as.data.frame(matrix(NA,ncol=8,nrow=1))
      book.keeping.temporary[1]=m
      book.keeping.temporary[2]=taxa[t]
      book.keeping.temporary[3]=site.number
      book.keeping.temporary[4]=i
      book.keeping.temporary[5]=independent.number
      book.keeping.temporary[6]=length(unique(subsample[[l]][[t]][[i]]$SiteID))
      book.keeping.temporary[7]=check.independent
      book.keeping.temporary[8]=check.complete
      
      names(book.keeping.temporary)=c("Criterion","Taxon","OriginalSampleSize","Subsample","IndependentSampleSize","AllSampleSize_post","SatisfyMinimumDistance","SatisfyDataCompletion")
      
      book.keeping=rbind(book.keeping,book.keeping.temporary)
    }
  }
}

## double check quality of subsamples
book.keeping[book.keeping$FilteredSampleSize>book.keeping$OriginalSampleSize,]
# returns only the NA row, which is the stub/first row -> no problem
unique(book.keeping$SatisfyMinimumDistance)
# returns only 1 and NA -> all good
# write book.keeping data frame to .csv
write.csv(book.keeping, "bookkeeping.csv")

### For disapperance data during breeding season (disappearance represented as "departure" in data set)
departure.breeding.date  = read.csv('Breedingdeparture.csv', header=T)
departure.breeding.date$SiteID = as.factor(departure.breeding.date $SiteID)
departure.breeding.date$Year = as.numeric(departure.breeding.date $Year)
head(departure.breeding.date)

resample=1000
taxa.all = list(c("Cicada","Frog","Goose","Swallow"),
                c("Cicada","Frog","Swallow"),
                c("Goose","Swallow"),
                c("Cicada","Frog","Swallow"))
cutoff.all=list(c(400,1000,600,800),
                c(200,1000,800),
                c(600,800),
                c(800,1000,800))
m.all=c(5,6,10,11)
subsample=list(NA)

book.keeping=as.data.frame(matrix(NA,ncol=8,nrow=1))
names(book.keeping)=c("Criterion","Taxon","OriginalSampleSize","Subsample","IndependentSampleSize","AllSampleSize_post","SatisfyMinimumDistance","SatisfyDataCompletion")

for (l in 1:4){# the large loop of criteria
  print(l)
  taxa = taxa.all[[l]]
  cutoff.suball = cutoff.all[[l]]
  m=m.all[l]
  taxa.number=length(taxa)
  subsample[[l]]=list(NA)
  
  for (t in 1:taxa.number){#looping through each taxa
    print(t)
    subsample[[l]][[t]]=list(NA)
    cutoff=cutoff.suball[[t]]
    
    # information about the original dataset that is not going to be updated in subsequent loops
    data=departure.breeding.date[departure.breeding.date$AnalysisGroup==taxa[[t]] & departure.breeding.date$m==m,]
    site.list=unique(data$SiteID)#the list of viable site list - to be updated in each j loop
    site.number=length(site.list) # the number of all sites - remains unchanged
    site.geo=as.matrix(data[!duplicated(data[,5:6]),5:6]) #the geo location of all sites, keeping only lat long
    geo.dist=spDists(x=site.geo,longlat=TRUE)#a matrix of distance between all point pairs, to be referenced to in checking distanc
    
    for (i in 1:resample){
      print(i)
      
      pick=rep(NA,site.number)#an empty vector to store all the picks
      
      # the first draw
      pick[1]=sample(x=1:site.number,size=1)#generate 1st "resample" draw - to be updated in each j loop
      subsample[[l]][[t]][[i]]=data[data$SiteID==site.list[pick[1]],]#put data from the first draw into sub - to be udpated in each j loop
      # add an identifier to note this is data from draws
      identifier = data.frame(matrix(NA,ncol=3,nrow=nrow(subsample[[l]][[t]][[i]])));names(identifier)=c("draw","whether_independent","cluster")
      identifier[,1]=1
      identifier[,2]=1
      identifier[,3]=site.list[pick[1]]
      subsample[[l]][[t]][[i]]=cbind(subsample[[l]][[t]][[i]],identifier)
      
      distance=geo.dist[pick[1],]#the distance of all sites from this first draw - to be updated in each j loop
      
      # sites dropped from this first draw
      pool.drop=site.list[which(distance<=cutoff & distance>0)]# the list of sites viable for further draw - to be updated in each j loop
      pool.size.drop = length(pool.drop)
      
      if(pool.size.drop>0){
        data.drop=data[data$SiteID==pool.drop[1],]# a stub of data.drop to be potentially added on
        if (pool.size.drop>1){
          for(p in 1:(pool.size.drop-1)){
            data.drop=rbind(data.drop,data[data$SiteID==pool.drop[1+p],])
          }
        }else{
          data.drop=data.drop
        }
        
        identifier.drop = data.frame(matrix(NA,ncol=3,nrow=nrow(data.drop)));names(identifier.drop)=c("draw","whether_independent","cluster")
        identifier.drop[,1]=1
        identifier.drop[,2]=0
        identifier.drop[,3]=site.list[pick[1]]
        
        data.drop=cbind(data.drop,identifier.drop)
        
        # site picked and the ones dropped
        subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],data.drop)
      }else{
        subsample[[l]][[t]][[i]]=subsample[[l]][[t]][[i]]
      }
      
      # viable sites after the first draw
      pool=site.list[which(distance>cutoff)]# the list of sites viable for further draw - to be updated in each j loop
      pool.size=length(pool)# the number of sites viable for further draw - to be updated in each j loop
      pool.dist=geo.dist[which(distance>cutoff),which(distance>cutoff)]# the distance matrix only for sites kept after applying the cutoff - to be updated in each j loop
      
      # subsequent draws based on these viable sites - viable sites are to be updated after each draw
      for(j in 1:(site.number-1)){#the max number of samples that could be drawn into the subsample
        print(j)
        if (pool.size > 1) {
          # a new draw
          pick[1+j]=sample(x=1:pool.size,size=1)#generate "resample"-number of start sites (in the form of its numbering postion)
          subsample.temp=data[data$SiteID==pool[pick[1+j]],]
          identifier = data.frame(matrix(NA,ncol=3,nrow=nrow(subsample.temp)));names(identifier)=c("draw","whether_independent","cluster")
          identifier[,1]=1+j
          identifier[,2]=1
          identifier[,3]=pool[pick[1+j]]
          subsample.temp=cbind(subsample.temp,identifier)
          
          subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],subsample.temp)
          
          distance=pool.dist[pick[1+j],]
          
          # sites dropped from this previous draw
          pool.drop=pool[which(distance<=cutoff & distance>0)]# the list of sites viable for further draw - to be updated in each j loop
          pool.size.drop = length(pool.drop)
          
          if(pool.size.drop>0){
            data.drop=data[data$SiteID==pool.drop[1],]# a stub of data.drop to be potentially added on
            if (pool.size.drop>1){
              for(p in 1:(pool.size.drop-1)){
                data.drop=rbind(data.drop,data[data$SiteID==pool.drop[1+p],])
              }
            }else{
              data.drop=data.drop
            }
            identifier.drop = data.frame(matrix(NA,ncol=3,nrow=nrow(data.drop)));names(identifier.drop)=c("draw","whether_independent","cluster")
            identifier.drop[,1]=1+j
            identifier.drop[,2]=0
            identifier.drop[,3]=pool[pick[1+j]]
            
            data.drop=cbind(data.drop,identifier.drop)
            
            # site picked and the ones dropped
            subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],data.drop)
          }else{
            subsample[[l]][[t]][[i]]=subsample[[l]][[t]][[i]]
          }
          
          # viable sites after the previous draw
          pool=pool[which(distance>cutoff)]
          pool.size=length(pool)# the number of sites viable for further draw - to be updated in each j loop
          pool.dist=pool.dist[which(distance>cutoff),which(distance>cutoff)]# the distance matrix only for sites kept after applying the cutoff - to be updated in each j loop
        }else{
          if(pool.size==1){
            # the only possible new draw
            pick[1+j]=1
            
            subsample.temp=data[data$SiteID==pool[pick[1+j]],]
            identifier = data.frame(matrix(NA,ncol=3,nrow=nrow(subsample.temp)));names(identifier)=c("draw","whether_independent","cluster")
            identifier[,1]=1+j
            identifier[,2]=1
            identifier[,3]=pool[pick[1+j]]
            subsample.temp=cbind(subsample.temp,identifier)
            
            subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],subsample.temp)
            
            pool.size=0
          }else{
            subsample[[l]][[t]][[i]]=subsample[[l]][[t]][[i]]
          }
        }  
      }
      
      filename = paste(l,taxa[[t]],i, "subsample.csv")
      write.csv(subsample[[l]][[t]][[i]], filename)
      
      # double check if sites too close to each other have been elimiated
      independent = subsample[[l]][[t]][[i]][subsample[[l]][[t]][[i]]$whether_independent==1,]
      independent.number=length(unique(independent$SiteID))
      new.site.geo=as.matrix(independent[!duplicated(independent[,5:6]),5:6]) #the geo location of all sites, keeping only lat long
      new.geo.dist=spDists(x=new.site.geo,longlat=TRUE)#a matrix of distance between all point pairs, to be referenced to in checking distanc
      
      if(sum(new.geo.dist>cutoff)==(dim(new.geo.dist)[1]^2-dim(new.geo.dist)[1])){#a logical test whether all distances except for the diagonals are greater than cutoff)
        check.independent=1
      }else{
        check.independent=0
      }
      
      # double check if the produced data set is of the same size as the original data set
      if(length(unique(data$SiteID))==length(unique(subsample[[l]][[t]][[i]]$SiteID))){
        check.complete=1
      }else{
        check.complete=0
      }
      
      book.keeping.temporary=as.data.frame(matrix(NA,ncol=8,nrow=1))
      book.keeping.temporary[1]=m
      book.keeping.temporary[2]=taxa[t]
      book.keeping.temporary[3]=site.number
      book.keeping.temporary[4]=i
      book.keeping.temporary[5]=independent.number
      book.keeping.temporary[6]=length(unique(subsample[[l]][[t]][[i]]$SiteID))
      book.keeping.temporary[7]=check.independent
      book.keeping.temporary[8]=check.complete
      
      names(book.keeping.temporary)=c("Criterion","Taxon","OriginalSampleSize","Subsample","IndependentSampleSize","AllSampleSize_post","SatisfyMinimumDistance","SatisfyDataCompletion")
      
      book.keeping=rbind(book.keeping,book.keeping.temporary)
    }
  }
}

## double check quality of subsamples
book.keeping[book.keeping$IndependentSampleSize>book.keeping$OriginalSampleSize,]
# returns only the NA row, which is the stub/first row -> no problem
unique(book.keeping$SatisfyMinimumDistance)
# returns only 1 and NA -> all good
# write book.keeping data frame to .csv
write.csv(book.keeping, "bookkeeping.csv")

### For data on length of temporal occurrence span during breeding season (length of temporal occurrence span represented as "span" in data set)
stay.breeding.date  = read.csv('Breedingstay.csv', header=T)
stay.breeding.date$SiteID = as.factor(stay.breeding.date $SiteID)
stay.breeding.date$Year = as.numeric(stay.breeding.date $Year)
head(stay.breeding.date)

resample=1000
taxa.all = list(c("Cicada","Frog","Goose","Swallow"),
                c("Cicada","Cricket","Frog","Goose","Swallow"),
                c("Goose","Swallow","Toad"),
                c("Frog","Goose","Swallow"))
cutoff.all=list(c(400,1000,600,800),
                c(200,200,1000,600,1200),
                c(600,1200,200),
                c(1000,600,1200))
m.all=c(5,6,10,11)
subsample=list(NA)

book.keeping=as.data.frame(matrix(NA,ncol=8,nrow=1))
names(book.keeping)=c("Criterion","Taxon","OriginalSampleSize","Subsample","IndependentSampleSize","AllSampleSize_post","SatisfyMinimumDistance","SatisfyDataCompletion")

for (l in 1:4){# the large loop of criteria
  print(l)
  taxa = taxa.all[[l]]
  cutoff.suball = cutoff.all[[l]]
  m=m.all[l]
  taxa.number=length(taxa)
  subsample[[l]]=list(NA)
  
  for (t in 1:taxa.number){#looping through each taxa
    print(t)
    subsample[[l]][[t]]=list(NA)
    cutoff=cutoff.suball[[t]]
    
    # information about the original dataset that is not going to be updated in subsequent loops
    data=stay.breeding.date[stay.breeding.date$AnalysisGroup==taxa[[t]] & stay.breeding.date$m==m,]
    site.list=unique(data$SiteID)#the list of viable site list - to be updated in each j loop
    site.number=length(site.list) # the number of all sites - remains unchanged
    site.geo=as.matrix(data[!duplicated(data[,5:6]),5:6]) #the geo location of all sites, keeping only lat long
    geo.dist=spDists(x=site.geo,longlat=TRUE)#a matrix of distance between all point pairs, to be referenced to in checking distanc
    
    for (i in 1:resample){
      print(i)
      
      pick=rep(NA,site.number)#an empty vector to store all the picks
      
      # the first draw
      pick[1]=sample(x=1:site.number,size=1)#generate 1st "resample" draw - to be updated in each j loop
      subsample[[l]][[t]][[i]]=data[data$SiteID==site.list[pick[1]],]#put data from the first draw into sub - to be udpated in each j loop
      # add an identifier to note this is data from draws
      identifier = data.frame(matrix(NA,ncol=3,nrow=nrow(subsample[[l]][[t]][[i]])));names(identifier)=c("draw","whether_independent","cluster")
      identifier[,1]=1
      identifier[,2]=1
      identifier[,3]=site.list[pick[1]]
      subsample[[l]][[t]][[i]]=cbind(subsample[[l]][[t]][[i]],identifier)
      
      distance=geo.dist[pick[1],]#the distance of all sites from this first draw - to be updated in each j loop
      
      # sites dropped from this first draw
      pool.drop=site.list[which(distance<=cutoff & distance>0)]# the list of sites viable for further draw - to be updated in each j loop
      pool.size.drop = length(pool.drop)
      
      if(pool.size.drop>0){
        data.drop=data[data$SiteID==pool.drop[1],]# a stub of data.drop to be potentially added on
        if (pool.size.drop>1){
          for(p in 1:(pool.size.drop-1)){
            data.drop=rbind(data.drop,data[data$SiteID==pool.drop[1+p],])
          }
        }else{
          data.drop=data.drop
        }
        
        identifier.drop = data.frame(matrix(NA,ncol=3,nrow=nrow(data.drop)));names(identifier.drop)=c("draw","whether_independent","cluster")
        identifier.drop[,1]=1
        identifier.drop[,2]=0
        identifier.drop[,3]=site.list[pick[1]]
        
        data.drop=cbind(data.drop,identifier.drop)
        
        # site picked and the ones dropped
        subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],data.drop)
      }else{
        subsample[[l]][[t]][[i]]=subsample[[l]][[t]][[i]]
      }
      
      # viable sites after the first draw
      pool=site.list[which(distance>cutoff)]# the list of sites viable for further draw - to be updated in each j loop
      pool.size=length(pool)# the number of sites viable for further draw - to be updated in each j loop
      pool.dist=geo.dist[which(distance>cutoff),which(distance>cutoff)]# the distance matrix only for sites kept after applying the cutoff - to be updated in each j loop
      
      # subsequent draws based on these viable sites - viable sites are to be updated after each draw
      for(j in 1:(site.number-1)){#the max number of samples that could be drawn into the subsample
        print(j)
        if (pool.size > 1) {
          # a new draw
          pick[1+j]=sample(x=1:pool.size,size=1)#generate "resample"-number of start sites (in the form of its numbering postion)
          subsample.temp=data[data$SiteID==pool[pick[1+j]],]
          identifier = data.frame(matrix(NA,ncol=3,nrow=nrow(subsample.temp)));names(identifier)=c("draw","whether_independent","cluster")
          identifier[,1]=1+j
          identifier[,2]=1
          identifier[,3]=pool[pick[1+j]]
          subsample.temp=cbind(subsample.temp,identifier)
          
          subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],subsample.temp)
          
          distance=pool.dist[pick[1+j],]
          
          # sites dropped from this previous draw
          pool.drop=pool[which(distance<=cutoff & distance>0)]# the list of sites viable for further draw - to be updated in each j loop
          pool.size.drop = length(pool.drop)
          
          if(pool.size.drop>0){
            data.drop=data[data$SiteID==pool.drop[1],]# a stub of data.drop to be potentially added on
            if (pool.size.drop>1){
              for(p in 1:(pool.size.drop-1)){
                data.drop=rbind(data.drop,data[data$SiteID==pool.drop[1+p],])
              }
            }else{
              data.drop=data.drop
            }
            identifier.drop = data.frame(matrix(NA,ncol=3,nrow=nrow(data.drop)));names(identifier.drop)=c("draw","whether_independent","cluster")
            identifier.drop[,1]=1+j
            identifier.drop[,2]=0
            identifier.drop[,3]=pool[pick[1+j]]
            
            data.drop=cbind(data.drop,identifier.drop)
            
            # site picked and the ones dropped
            subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],data.drop)
          }else{
            subsample[[l]][[t]][[i]]=subsample[[l]][[t]][[i]]
          }
          
          # viable sites after the previous draw
          pool=pool[which(distance>cutoff)]
          pool.size=length(pool)# the number of sites viable for further draw - to be updated in each j loop
          pool.dist=pool.dist[which(distance>cutoff),which(distance>cutoff)]# the distance matrix only for sites kept after applying the cutoff - to be updated in each j loop
        }else{
          if(pool.size==1){
            # the only possible new draw
            pick[1+j]=1
            
            subsample.temp=data[data$SiteID==pool[pick[1+j]],]
            identifier = data.frame(matrix(NA,ncol=3,nrow=nrow(subsample.temp)));names(identifier)=c("draw","whether_independent","cluster")
            identifier[,1]=1+j
            identifier[,2]=1
            identifier[,3]=pool[pick[1+j]]
            subsample.temp=cbind(subsample.temp,identifier)
            
            subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],subsample.temp)
            
            pool.size=0
          }else{
            subsample[[l]][[t]][[i]]=subsample[[l]][[t]][[i]]
          }
        }  
      }
      
      filename = paste(l,taxa[[t]],i, "subsample.csv")
      write.csv(subsample[[l]][[t]][[i]], filename)
      
      # double check if sites too close to each other have been elimiated
      independent = subsample[[l]][[t]][[i]][subsample[[l]][[t]][[i]]$whether_independent==1,]
      independent.number=length(unique(independent$SiteID))
      new.site.geo=as.matrix(independent[!duplicated(independent[,5:6]),5:6]) #the geo location of all sites, keeping only lat long
      new.geo.dist=spDists(x=new.site.geo,longlat=TRUE)#a matrix of distance between all point pairs, to be referenced to in checking distanc
      
      if(sum(new.geo.dist>cutoff)==(dim(new.geo.dist)[1]^2-dim(new.geo.dist)[1])){#a logical test whether all distances except for the diagonals are greater than cutoff)
        check.independent=1
      }else{
        check.independent=0
      }
      
      # double check if the produced data set is of the same size as the original data set
      if(length(unique(data$SiteID))==length(unique(subsample[[l]][[t]][[i]]$SiteID))){
        check.complete=1
      }else{
        check.complete=0
      }
      
      book.keeping.temporary=as.data.frame(matrix(NA,ncol=8,nrow=1))
      book.keeping.temporary[1]=m
      book.keeping.temporary[2]=taxa[t]
      book.keeping.temporary[3]=site.number
      book.keeping.temporary[4]=i
      book.keeping.temporary[5]=independent.number
      book.keeping.temporary[6]=length(unique(subsample[[l]][[t]][[i]]$SiteID))
      book.keeping.temporary[7]=check.independent
      book.keeping.temporary[8]=check.complete
      
      names(book.keeping.temporary)=c("Criterion","Taxon","OriginalSampleSize","Subsample","IndependentSampleSize","AllSampleSize_post","SatisfyMinimumDistance","SatisfyDataCompletion")
      
      book.keeping=rbind(book.keeping,book.keeping.temporary)
    }
  }
}

## double check quality of subsamples
book.keeping[book.keeping$IndependentSampleSize>book.keeping$OriginalSampleSize,]
# returns only the NA row, which is the stub/first row -> no problem
unique(book.keeping$SatisfyMinimumDistance)
# returns only 1 and NA -> all good
# write book.keeping data frame to .csv
write.csv(book.keeping, "bookkeeping.csv")

### For temporal overlap data during breeding season
overlap.breeding.date  = read.csv('CommunityOverlap.csv', header=T)
overlap.breeding.date$SiteID = as.factor(overlap.breeding.date $SiteID)
overlap.breeding.date$Year = as.numeric(overlap.breeding.date $Year)
head(overlap.breeding.date)

resample=1000
taxa.all = list(c("Cicada_Cricket","Cicada_Cuckoo","Cicada_Swallow","Frog_Swallow","Frog_Toad","Goose_Swallow"),
                c("Cicada_Swallow","Cricket_Frog","Frog_Swallow","Frog_Toad","Goose_Swallow"),
                c("Cicada_Swallow","Frog_Swallow"))
cutoff.all=list(c(200,200,200,1000,200,400),
                c(200,200,800,200,400),
                c(200,200))
m.all=c(5,6,11)
subsample=list(NA)

book.keeping=as.data.frame(matrix(NA,ncol=8,nrow=1))
names(book.keeping)=c("Criterion","Taxon","OriginalSampleSize","Subsample","IndependentSampleSize","AllSampleSize_post","SatisfyMinimumDistance","SatisfyDataCompletion")

for (l in 1:3){# the large loop of criteria
  print(l)
  taxa = taxa.all[[l]]
  cutoff.suball = cutoff.all[[l]]
  m=m.all[l]
  taxa.number=length(taxa)
  subsample[[l]]=list(NA)
  
  for (t in 1:taxa.number){#looping through each taxa
    print(t)
    subsample[[l]][[t]]=list(NA)
    cutoff=cutoff.suball[[t]]
    
    # information about the original dataset that is not going to be updated in subsequent loops
    data=overlap.breeding.date[overlap.breeding.date$TaxonPair==taxa[[t]] & overlap.breeding.date$m==m,]
    site.list=unique(data$SiteID)#the list of viable site list - to be updated in each j loop
    site.number=length(site.list) # the number of all sites - remains unchanged
    site.geo=as.matrix(data[!duplicated(data[,4:5]),4:5]) #the geo location of all sites, keeping only lat long
    geo.dist=spDists(x=site.geo,longlat=TRUE)#a matrix of distance between all point pairs, to be referenced to in checking distanc
    
    for (i in 1:resample){
      print(i)
      
      pick=rep(NA,site.number)#an empty vector to store all the picks
      
      # the first draw
      pick[1]=sample(x=1:site.number,size=1)#generate 1st "resample" draw - to be updated in each j loop
      subsample[[l]][[t]][[i]]=data[data$SiteID==site.list[pick[1]],]#put data from the first draw into sub - to be udpated in each j loop
      # add an identifier to note this is data from draws
      identifier = data.frame(matrix(NA,ncol=3,nrow=nrow(subsample[[l]][[t]][[i]])));names(identifier)=c("draw","whether_independent","cluster")
      identifier[,1]=1
      identifier[,2]=1
      identifier[,3]=site.list[pick[1]]
      subsample[[l]][[t]][[i]]=cbind(subsample[[l]][[t]][[i]],identifier)
      
      distance=geo.dist[pick[1],]#the distance of all sites from this first draw - to be updated in each j loop
      
      # sites dropped from this first draw
      pool.drop=site.list[which(distance<=cutoff & distance>0)]# the list of sites viable for further draw - to be updated in each j loop
      pool.size.drop = length(pool.drop)
      
      if(pool.size.drop>0){
        data.drop=data[data$SiteID==pool.drop[1],]# a stub of data.drop to be potentially added on
        if (pool.size.drop>1){
          for(p in 1:(pool.size.drop-1)){
            data.drop=rbind(data.drop,data[data$SiteID==pool.drop[1+p],])
          }
        }else{
          data.drop=data.drop
        }
        
        identifier.drop = data.frame(matrix(NA,ncol=3,nrow=nrow(data.drop)));names(identifier.drop)=c("draw","whether_independent","cluster")
        identifier.drop[,1]=1
        identifier.drop[,2]=0
        identifier.drop[,3]=site.list[pick[1]]
        
        data.drop=cbind(data.drop,identifier.drop)
        
        # site picked and the ones dropped
        subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],data.drop)
      }else{
        subsample[[l]][[t]][[i]]=subsample[[l]][[t]][[i]]
      }
      
      # viable sites after the first draw
      pool=site.list[which(distance>cutoff)]# the list of sites viable for further draw - to be updated in each j loop
      pool.size=length(pool)# the number of sites viable for further draw - to be updated in each j loop
      pool.dist=geo.dist[which(distance>cutoff),which(distance>cutoff)]# the distance matrix only for sites kept after applying the cutoff - to be updated in each j loop
      
      # subsequent draws based on these viable sites - viable sites are to be updated after each draw
      for(j in 1:(site.number-1)){#the max number of samples that could be drawn into the subsample
        print(j)
        if (pool.size > 1) {
          # a new draw
          pick[1+j]=sample(x=1:pool.size,size=1)#generate "resample"-number of start sites (in the form of its numbering postion)
          subsample.temp=data[data$SiteID==pool[pick[1+j]],]
          identifier = data.frame(matrix(NA,ncol=3,nrow=nrow(subsample.temp)));names(identifier)=c("draw","whether_independent","cluster")
          identifier[,1]=1+j
          identifier[,2]=1
          identifier[,3]=pool[pick[1+j]]
          subsample.temp=cbind(subsample.temp,identifier)
          
          subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],subsample.temp)
          
          distance=pool.dist[pick[1+j],]
          
          # sites dropped from this previous draw
          pool.drop=pool[which(distance<=cutoff & distance>0)]# the list of sites viable for further draw - to be updated in each j loop
          pool.size.drop = length(pool.drop)
          
          if(pool.size.drop>0){
            data.drop=data[data$SiteID==pool.drop[1],]# a stub of data.drop to be potentially added on
            if (pool.size.drop>1){
              for(p in 1:(pool.size.drop-1)){
                data.drop=rbind(data.drop,data[data$SiteID==pool.drop[1+p],])
              }
            }else{
              data.drop=data.drop
            }
            identifier.drop = data.frame(matrix(NA,ncol=3,nrow=nrow(data.drop)));names(identifier.drop)=c("draw","whether_independent","cluster")
            identifier.drop[,1]=1+j
            identifier.drop[,2]=0
            identifier.drop[,3]=pool[pick[1+j]]
            
            data.drop=cbind(data.drop,identifier.drop)
            
            # site picked and the ones dropped
            subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],data.drop)
          }else{
            subsample[[l]][[t]][[i]]=subsample[[l]][[t]][[i]]
          }
          
          # viable sites after the previous draw
          pool=pool[which(distance>cutoff)]
          pool.size=length(pool)# the number of sites viable for further draw - to be updated in each j loop
          pool.dist=pool.dist[which(distance>cutoff),which(distance>cutoff)]# the distance matrix only for sites kept after applying the cutoff - to be updated in each j loop
        }else{
          if(pool.size==1){
            # the only possible new draw
            pick[1+j]=1
            
            subsample.temp=data[data$SiteID==pool[pick[1+j]],]
            identifier = data.frame(matrix(NA,ncol=3,nrow=nrow(subsample.temp)));names(identifier)=c("draw","whether_independent","cluster")
            identifier[,1]=1+j
            identifier[,2]=1
            identifier[,3]=pool[pick[1+j]]
            subsample.temp=cbind(subsample.temp,identifier)
            
            subsample[[l]][[t]][[i]]=rbind(subsample[[l]][[t]][[i]],subsample.temp)
            
            pool.size=0
          }else{
            subsample[[l]][[t]][[i]]=subsample[[l]][[t]][[i]]
          }
        }  
      }
      
      filename = paste(l,taxa[[t]],i, "subsample.csv")
      write.csv(subsample[[l]][[t]][[i]], filename)
      
      # double check if sites too close to each other have been elimiated
      independent = subsample[[l]][[t]][[i]][subsample[[l]][[t]][[i]]$whether_independent==1,]
      independent.number=length(unique(independent$SiteID))
      new.site.geo=as.matrix(independent[!duplicated(independent[,4:5]),4:5]) #the geo location of all sites, keeping only lat long
      new.geo.dist=spDists(x=new.site.geo,longlat=TRUE)#a matrix of distance between all point pairs, to be referenced to in checking distanc
      
      if(sum(new.geo.dist>cutoff)==(dim(new.geo.dist)[1]^2-dim(new.geo.dist)[1])){#a logical test whether all distances except for the diagonals are greater than cutoff)
        check.independent=1
      }else{
        check.independent=0
      }
      
      # double check if the produced data set is of the same size as the original data set
      if(length(unique(data$SiteID))==length(unique(subsample[[l]][[t]][[i]]$SiteID))){
        check.complete=1
      }else{
        check.complete=0
      }
      
      book.keeping.temporary=as.data.frame(matrix(NA,ncol=8,nrow=1))
      book.keeping.temporary[1]=m
      book.keeping.temporary[2]=taxa[t]
      book.keeping.temporary[3]=site.number
      book.keeping.temporary[4]=i
      book.keeping.temporary[5]=independent.number
      book.keeping.temporary[6]=length(unique(subsample[[l]][[t]][[i]]$SiteID))
      book.keeping.temporary[7]=check.independent
      book.keeping.temporary[8]=check.complete
      
      names(book.keeping.temporary)=c("Criterion","Taxon","OriginalSampleSize","Subsample","IndependentSampleSize","AllSampleSize_post","SatisfyMinimumDistance","SatisfyDataCompletion")
      
      book.keeping=rbind(book.keeping,book.keeping.temporary)
    }
  }
}

## double check quality of subsamples
book.keeping[book.keeping$IndependentSampleSize>book.keeping$OriginalSampleSize,]
# returns only the NA row, which is the stub/first row -> no problem
unique(book.keeping$SatisfyMinimumDistance)
# returns only 1 and NA -> all good
# write book.keeping data frame to .csv
write.csv(book.keeping, "bookkeeping.csv")
