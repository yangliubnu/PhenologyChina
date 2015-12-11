### 1st type of scenario: for both species, TOW position changes; TOW span does not change
##  Patterns of TOW change for each species under climate change
begin.a.change = rnorm(n=10000,mean=0,sd=1)
begin.b.change = rnorm(n=10000,mean=0,sd=1)
length.a.change = 0
length.b.change = 0

##  Updated patterns of TOW for each species
begin.a = begin.a.original+begin.a.change
begin.b = begin.b.original+begin.b.change
length.a = length.a.original+length.a.change
length.b = length.b.original+length.b.change
end.a = begin.a+length.a
end.b = begin.b+length.b

##  Updated patterns of temporal overlap between species a and b
overlap = rep(NA,10000)
for(i in 1:10000){
  overlap[i] = min(end.a[i], end.b[i])-max(begin.a[i], begin.b[i])
  if(overlap[i]<0){
    overlap[i]=0}
  }

##  Plot kernel density curve
plot(density(overlap),xlim=c(0,4),ylim=c(0,2),lwd=3,
     col="grey30") # Plot exported at 750*450

### 2nd type of scenario: TOW pisition changes for both species; TOW span changes by +1 for species a, remains unchanged for species b
##  Patterns of TOW change for each species under climate change
begin.a.change = rnorm(n=10000,mean=0,sd=1)
begin.b.change = rnorm(n=10000,mean=0,sd=1)
length.a.change = 1
length.b.change = 0

##  Updated patterns of TOW for each species
begin.a = begin.a.original+begin.a.change
begin.b = begin.b.original+begin.b.change
length.a = length.a.original+length.a.change
length.b = length.b.original+length.b.change
end.a = begin.a+length.a
end.b = begin.b+length.b

##  Updated patterns of temporal overlap between species a and b
overlap.1 = rep(NA,10000)
for(i in 1:10000){
  overlap.1[i] = min(end.a[i], end.b[i])-max(begin.a[i], begin.b[i])
  if(overlap.1[i]<0){
    overlap.1[i]=0}
}

##  Plotting of kernel density curve on hold - to be plotted together with the next two plots

### 3rd type of scenario: TOW pisition changes for both species; TOW span changes by +1 for species b, remains unchanged for species a
##  Patterns of TOW change for each species under climate change
begin.a.change = rnorm(n=10000,mean=0,sd=1)
begin.b.change = rnorm(n=10000,mean=0,sd=1)
length.a.change = 0
length.b.change = 1

##  Updated patterns of TOW for each species
begin.a = begin.a.original+begin.a.change
begin.b = begin.b.original+begin.b.change
length.a = length.a.original+length.a.change
length.b = length.b.original+length.b.change
end.a = begin.a+length.a
end.b = begin.b+length.b

##  Updated patterns of temporal overlap between species a and b
overlap.2 = rep(NA,10000)
for(i in 1:10000){
  overlap.2[i] = min(end.a[i], end.b[i])-max(begin.a[i], begin.b[i])
  if(overlap.2[i]<0){
    overlap.2[i]=0}
}

##  Plotting of kernel density curve on hold - to be plotted together with the next plot

### 4th type of scenario: TOW pisition changes for both species; TOW span changes by +1 for both species
##  Patterns of TOW change for each species under climate change
begin.a.change = rnorm(n=10000,mean=0,sd=1)
begin.b.change = rnorm(n=10000,mean=0,sd=1)
length.a.change = 1
length.b.change = 1

##  Updated patterns of TOW for each species
begin.a = begin.a.original+begin.a.change
begin.b = begin.b.original+begin.b.change
length.a = length.a.original+length.a.change
length.b = length.b.original+length.b.change
end.a = begin.a+length.a
end.b = begin.b+length.b

##  Updated patterns of temporal overlap between species a and b
overlap.3 = rep(NA,10000)
for(i in 1:10000){
  overlap.3[i] = min(end.a[i], end.b[i])-max(begin.a[i], begin.b[i])
  if(overlap.3[i]<0){
    overlap.3[i]=0}
}

##  Plot kernel density curve: for all three preceding plots (scenarios 2-4)
plot(density(overlap.1),xlim=c(0,4),ylim=c(0,2),lwd=3,
     col="coral")
lines(density(overlap.2, bw=0.1649), lwd=3,
      col="blue",)
lines(density(overlap.3, bw=0.1649),lwd=3,
      col="grey30")

### 5th type of scenario: TOW pisition changes for both species; TOW span changes by -1 for species a, remains unchanged for species b
##  Patterns of TOW change for each species under climate change
begin.a.change = rnorm(n=10000,mean=0,sd=1)
begin.b.change = rnorm(n=10000,mean=0,sd=1)
length.a.change = -1
length.b.change = 0

##  Updated patterns of TOW for each species
begin.a = begin.a.original+begin.a.change
begin.b = begin.b.original+begin.b.change
length.a = length.a.original+length.a.change
length.b = length.b.original+length.b.change
end.a = begin.a+length.a
end.b = begin.b+length.b

##  Updated patterns of temporal overlap between species a and b
overlap.1 = rep(NA,10000)
for(i in 1:10000){
  overlap.1[i] = min(end.a[i], end.b[i])-max(begin.a[i], begin.b[i])
  if(overlap.1[i]<0){
    overlap.1[i]=0}
}

##  Plotting of kernel density curve on hold - to be plotted together with the next two plots

### 6th type of scenario: TOW pisition changes for both species; TOW span changes by -1 for species b, remains unchanged for species a
##  Patterns of TOW change for each species under climate change
begin.a.change = rnorm(n=10000,mean=0,sd=1)
begin.b.change = rnorm(n=10000,mean=0,sd=1)
length.a.change = 0
length.b.change = -1

##  Updated patterns of TOW for each species
begin.a = begin.a.original+begin.a.change
begin.b = begin.b.original+begin.b.change
length.a = length.a.original+length.a.change
length.b = length.b.original+length.b.change
end.a = begin.a+length.a
end.b = begin.b+length.b

##  Updated patterns of temporal overlap between species a and b
overlap.2 = rep(NA,10000)
for(i in 1:10000){
  overlap.2[i] = min(end.a[i], end.b[i])-max(begin.a[i], begin.b[i])
  if(overlap.2[i]<0){
    overlap.2[i]=0}
}

##  Plotting of kernel density curve on hold - to be plotted together with the next plot

### 7th type of scenario: TOW pisition changes for both species; TOW span changes by -1 for both species
##  Patterns of TOW change for each species under climate change
begin.a.change = rnorm(n=10000,mean=0,sd=1)
begin.b.change = rnorm(n=10000,mean=0,sd=1)
length.a.change = -1
length.b.change = -1

##  Updated patterns of TOW for each species
begin.a = begin.a.original+begin.a.change
begin.b = begin.b.original+begin.b.change
length.a = length.a.original+length.a.change
length.b = length.b.original+length.b.change
end.a = begin.a+length.a
end.b = begin.b+length.b

##  Updated patterns of temporal overlap between species a and b
overlap.3 = rep(NA,10000)
for(i in 1:10000){
  overlap.3[i] = min(end.a[i], end.b[i])-max(begin.a[i], begin.b[i])
  if(overlap.3[i]<0){
    overlap.3[i]=0}
}

##  Plot kernel density curve: for all three preceding plots (scenarios 5-7)
plot(density(overlap.1),xlim=c(0,4),ylim=c(0,2),lwd=3,
     col="coral")
lines(density(overlap.2, bw=0.09882), lwd=3,
      col="blue",)
lines(density(overlap.3, bw=0.09882),lwd=3,
      col="grey30")
