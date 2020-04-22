---
title: "Biophysics Pset 1"
output:
  pdf_document: default
  html_notebook: default
  word_document: default
  html_document:
    df_print: paged
---
Maya Waarts Biophysics Pset 1

Problem 1a-c. 

```r
n<-1000
k<-list(0.001,0.005,0.01,0.05,0.1)
Means<-list(0,0,0,0)
Variances<-list(0,0,0,0)


for (m in 1:length(k)) {
N_in_k<-rep(0, 100)
for (i in 1:100) {
  list <- rep(0, n)
  for (j in 1:n) {  #generate list of random numbers
    rand<-runif(1)
    list[j]<-rand
    if (rand < k[m]) {  #count how many of the random numbers are in [0,k]
      N_in_k[i]<-N_in_k[i]+1
    }
  }
Means[m]<-mean(N_in_k)
Variances[m] <-var(N_in_k)
}

}

plot(main="Variance Dependence on Mean for Smaller k", Means, Variances)
```

![](biophysics-_files/figure-latex/unnamed-chunk-1-1.pdf)<!-- --> 
For very small values of k, the variance of N equals the mean of N. This is what one would expect from a poisson process. The process of randomly and independently selecting random numbers within an interval is itself not a poisson process but rather a discrete, binomial process for n events, each with an independent probability of success, p=k. A poisson distribution in contrast to a binomial distribution is continuous, but a binomial distribution approaches a poisson distribution as n-> infinity and p-> 0 such that np-> lambda where lambda is the success rate. N was chosen to be sufficiently large here to approximate the first condition (n-> infinity), and when p is close to 0, as it is in this set of k values (satisfying the second condition) the approximation is valid. A common rule of thumb is that the approximation is reasonable when np<10. In this example n = 100 and p ranges from 0.001 to 0.1 giving a maximum np of 10. This explains why the distribution appears to be explained by a poisson process for these values. 

1d. 

```r
n<-1000
k<-list(0.1,0.2,0.4,0.8,1)
Means<-list(0,0,0,0)
Variances<-list(0,0,0,0)

for (m in 1:length(k)) {
N_in_k<-rep(0, 100)
for (i in 1:100) {
  list <- rep(0, n)
  for (j in 1:n) {   #generate list of random numbers
    rand<-runif(1)
    list[j]<-rand
    if (rand < k[m]) {  #count how many of the random numbers are within [0,k]
      N_in_k[i]<-N_in_k[i]+1
    }
  }
Means[m]<-mean(N_in_k)
Variances[m] <-var(N_in_k)
}

}

#plot(type="l", xlab="Time", ylab="density", ylim=c(0,65), col ="red")
plot(main="Variance Dependence on Mean for Larger k", Means, Variances)
```

![](biophysics-_files/figure-latex/unnamed-chunk-2-1.pdf)<!-- --> 
For larger k's the results show a binomial distribution. This is because the k values are too large for the process to be accurately approximated as a poisson process. The probability of success is too far from 0 and so it can no longer be thought of as sampling from an infinate range. This is most clear for the k=1 case. A poisson process would predict a mean of n (the number of events), and a variance also of n. This cannot be the case, however, because the points were taken from a range between 0 and 1 so all the points must be less than one and the variance must be 0 as it is in the plot. 



Problem 2, L_run/L_c = 0.01

```r
x<-seq(0,6,.01)
time<-seq(0,499,1)  #define length of time to track the e. coli 
L<- 0.046  # define L_run to be 0.01 of L_c (L_c is 4.6 based on profile defined later)
totalTime<-500
runSize<-0 # initialize run size
n<-100 # set number of trials for the e. coli run 
percentTimeAtHighPeak<-rep(0,n) 
percentTimeAtLowPeak<-rep(0,n) 
startingPoint<- rep(0,n) #initialize starting point

conc<-function(X) {   # calculate the concentration profile at a given point
  y<- -(X-5)^4 - 6*(X-5)^3 - 8*(X-5)^2 + 7 # set profile to be a two - peaked curve
  return(y)
}

derivativeCalc<- function(p, step, Direction) { # calculate the concentration gradient 
 if (Direction==1){  # calculate gradient in the direction of motion, "1" = positive x direction
  g<- (conc(p+step) - conc(p))/step
 }
 else {
    g<- (conc(p-step) - conc(p))/step
  }
  return(g)
}

runSizeCalc<-function(G) { # calculate the run size proportional to the size of the gradient
  possibleGradients<- derivativeCalc(x,L,0)  
  maxGradient<-max(possibleGradients) #determine maximum possible gradient to use as normalizer
  minGradient<-0
  
  normalizedRunSize<- L+(G/maxGradient)*(2*L) #determines a run size between L and 3L
  return(normalizedRunSize)
  }

plot(main = "Concentration profile", xlab="Position", ylab="Concentration", x,conc(x), ylim=c(0,18), xlim=c(0, 6), type="l")

for (j in 1:n) {
#define e. coli starting position and initialize list of positions for each time point
yep<-rep(0,totalTime)
startingPoint[j]<-runif(1,1,5.6) # randomly select a starting point
yep[1]<-startingPoint[j]

for (i in 1:499) { #run 
  direction<-sample(0:1, 1) #randomly chose direction 
  gradient <-derivativeCalc(yep[i], L, direction)
  if (gradient>0) { #if gradient is possitive, choose run size based on magnitude of gradient
    runSize<-runSizeCalc(gradient)
  }
  else { #if gradient is negative, default to run size L_run
   runSize<-L
  }
  if (direction==0) { #move left
   yep[i+1]=yep[i]-runSize
  }
  else { #move right
    yep[i+1]=yep[i]+runSize
  }
}

timeAtHighPeak<-0
timeAtLowPeak <-0
for (i in 1:length(yep)) { #sum up time spent at each peak
  if (yep[i] > 1.5 && yep[i]<2.0) { #define range of high peak
    timeAtHighPeak<-timeAtHighPeak+1
  }
  else if (yep[i] > 4.8 && yep [i]< 5.3) { #define range of low peak
    timeAtLowPeak<-timeAtLowPeak+1
  }
}

percentTimeAtHighPeak[j] <- (timeAtHighPeak/totalTime)*100
percentTimeAtLowPeak[j] <- (timeAtLowPeak/totalTime)*100
}

#points(yep, conc(yep))
abline (v=1.5, col="red")
abline (v=2.0, col="red")
abline (v=4.8, col="red")
abline (v=5.3, col="red")
```

![](biophysics-_files/figure-latex/unnamed-chunk-3-1.pdf)<!-- --> 

```r
#plot(time, yep, ylim =c(1,5.6), type = 'l')

plot(main="Time spent at peak vs. starting position (small step size)", xlab="Starting Point", ylab="Percent Time", startingPoint, percentTimeAtHighPeak, col="red")
legend(4.5,30, legend =c("High Peak", "Low Peak"), col =c("red","blue"), pch=1)
points(startingPoint, percentTimeAtLowPeak, col="blue")
```

![](biophysics-_files/figure-latex/unnamed-chunk-3-2.pdf)<!-- --> 
In this simulation, L_c = 4.6 and L_run = 0.046. The ratio of L_run/L_c is therefore 0.01. Time at a peak is defined as being within 0.5 units of the maximum as illustrated by the vertical red lines on the concentration profile. In almost all cases, the e. coli will spend some significant (>20% total time) time at a peak, but the amount of time it spends at a peak depends on its starting position. If the e. coli starts at the trough (position ~3-4) it spends between 10% and 60% (with rare case <20%) of its time at a peak. If the e. coli starts in a trough it may find either of the two peaks, not necessarily only the one closer to it although it more commonly finds the one closer to it. If the e. coli starts near a peak on the other hand, it will spend a greater percentage of its' time at that peak (60-80%) and will almost never spend any significant time at the other peak even when the other peak is higher. There are, however, a few cases in which the e. coli starts near the low peak and ends up spending a greater percentage of time at the higher peak. There are far less cases where the reverse happens (ie the e. coli starts near the higher peak and spends a greater percentage of time at the lower peak). This is expected since an e. coli that starts near a high peak would have no advantage in finding a lower peak. 


Problem 2, L_run/L_c = 0.5

```r
x<-seq(0,6,.01)
time<-seq(0,499,1)  #define length of time to track the e. coli 
L<- 2.3  # define L_run to be 0.5 of L_c (L_c is 4.6 based on profile defined later)
totalTime<-500
runSize<-0 # initialize run size
n<-100 # set number of trials for the e. coli run 
percentTimeAtHighPeak<-rep(0,n) 
percentTimeAtLowPeak<-rep(0,n) 
startingPoint<- rep(0,n) #initialize starting point

conc<-function(X) {   # calculate the concentration profile at a given point
  y<- -(X-5)^4 - 6*(X-5)^3 - 8*(X-5)^2 + 7 # set profile to be a two - peaked curve
  return(y)
}

derivativeCalc<- function(p, step, Direction) { # calculate the concentration gradient 
 if (Direction==1){  # calculate gradient in the direction of motion, "1" = positive x direction
  g<- (conc(p+step) - conc(p))/step
 }
 else {
    g<- (conc(p-step) - conc(p))/step
  }
  return(g)
}

runSizeCalc<-function(G) { # calculate the run size proportional to the size of the gradient
  possibleGradients<- derivativeCalc(x,L,0)  
  maxGradient<-max(possibleGradients) #determine maximum possible gradient to use as normalizer
  minGradient<-0
  
  normalizedRunSize<- L+(G/maxGradient)*(2*L) #determines a run size between L and 3L
  return(normalizedRunSize)
  }

plot(main = "Concentration profile", xlab="Position", ylab="Concentration", x,conc(x), ylim=c(0,18), xlim=c(0, 6), type="l")

for (j in 1:n) {
#define e. coli starting position and initialize list of positions for each time point
yep<-rep(0,totalTime)
startingPoint[j]<-runif(1,1,5.6) # randomly select a starting point
yep[1]<-startingPoint[j]

for (i in 1:499) { #run 
  direction<-sample(0:1, 1) #randomly chose direction 
  gradient <-derivativeCalc(yep[i], L, direction)
  if (gradient>0) { #if gradient is possitive, choose run size based on magnitude of gradient
    runSize<-runSizeCalc(gradient)
  }
  else { #if gradient is negative, default to run size L_run
   runSize<-L
  }
  if (direction==0) { #move left
   yep[i+1]=yep[i]-runSize
  }
  else { #move right
    yep[i+1]=yep[i]+runSize
  }
}

timeAtHighPeak<-0
timeAtLowPeak <-0
for (i in 1:length(yep)) { #sum up time spent at each peak
  if (yep[i] > 1.5 && yep[i]<2.0) { #define range of high peak
    timeAtHighPeak<-timeAtHighPeak+1
  }
  else if (yep[i] > 4.8 && yep [i]< 5.3) { #define range of low peak
    timeAtLowPeak<-timeAtLowPeak+1
  }
}

percentTimeAtHighPeak[j] <- (timeAtHighPeak/totalTime)*100
percentTimeAtLowPeak[j] <- (timeAtLowPeak/totalTime)*100
}

#points(yep, conc(yep))
abline (v=1.5, col="red")
abline (v=2.0, col="red")
abline (v=4.8, col="red")
abline (v=5.3, col="red")
```

![](biophysics-_files/figure-latex/unnamed-chunk-4-1.pdf)<!-- --> 

```r
#plot(time, yep, ylim =c(1,5.6), type = 'l')

plot(main="Time spent at peak vs. starting position (large step size)", xlab="Starting Point", ylab="Percent Time", startingPoint, percentTimeAtHighPeak, col="red", ylim = c(0,100))
legend(4.5,30, legend =c("High Peak", "Low Peak"), col =c("red","blue"), pch=1)
points(startingPoint, percentTimeAtLowPeak, col="blue")
```

![](biophysics-_files/figure-latex/unnamed-chunk-4-2.pdf)<!-- --> 
In this simulation, L_c = 4.6 and L_run = 2.3. The ratio of L_run/L_c is therefore 0.5. Since this ratio is so large the motion of the e. coli is essentially random as even once it finds the peak, its' next step is long enough to bring it away already so  it is unable to stay there. This is especially true since the peaks in this concentration profile span about 2.3 distance units so a single step of minimum size is enough to get away from the peak. Therefore the e. coli spends no significant amount of time at either peak regardless of its' starting position. 


