require(lubridate)

maxTiter <<- 15

#beta constrained between 1 and Nneg/Ct
#Enforce Nneg-beta*Ct >0  ->  Nneg/Ct>beta

Nneg = 100
Npos = c(1,3,3,4,1,2,4,3,3,4,6,2,2,1,0)
beta = 2
Ct=10

#maxweights = sort(runif(10))
#maxweights = maxweights/sum(maxweights)

maxweights= c(1,1,2,3,4,6,4,2,1,.5)
maxweights = maxweights/sum(maxweights)

cutoffdecay = 6
lambdahi=.42
lambdalow=0.05



#maxweights starts at cutoffdecay weight goes to 15 weight
#it controls the possibilities of max titer by how many elements it has

onestep <- function(Npos,Nneg,beta, Ct, 
                    maxweights,cutoffdecay,
                    lambdahi,lambdalow, cutoffinfect=1){
    
    #First we decay
    
    highTiters = seq(15-length(maxweights)+1,15)   
    for (jj in highTiters){
        decay = rbinom(1,Npos[jj],lambdahi)
        Npos[jj] = Npos[jj]-decay
        Npos[jj-1] = Npos[jj-1] + decay
        
    }
    
    lowTiters = 1:(15-length(highTiters))
    for (jj in lowTiters){
        decay = rbinom(1,Npos[jj],lambdalow)
        Npos[jj] = Npos[jj]-decay
        if(jj-1 == 0) {
            Nneg = Nneg+decay
        } else { Npos[jj-1] = Npos[jj-1] + decay }
    }
    
    ######Calculate number of new infections
    
    #dataIs = beta*Ct
    #modelIs = sum(Npos[highTiters])
    newIs = beta*Ct
    
    if (newIs > 0){
    #Infect!#################
        Nneg = Nneg - newIs
        for (ii in seq_along(newIs)){
            max=sample(highTiters,1,prob=maxweights)
            Npos[max] = Npos[max]+1
        }
    }
        #if (Nneg < 0) message("N0:", Nneg, " Beta:", beta, " Ct:",Ct)
    return(list(Npos=Npos,Nneg=Nneg))
}


likeli <- function(Npos,Nneg,obspos,obsneg){
    ppos = Npos/sum(Npos,Nneg)
    pneg = Nneg/sum(Npos,Nneg)
    if (pneg < 0){
        return(-100)
    }
    const = factorial(sum(obspos,obsneg))/
        prod(factorial(obspos),factorial(obsneg))
    likeli= const*prod(ppos^obspos)*pneg^obsneg
    return(likeli)
}

loglikelifunc <- function(Npos,Nneg,obspos,obsneg){
    ppos = Npos/sum(Npos,Nneg)
    pneg = Nneg/sum(Npos,Nneg)
    if (pneg <= 0) {
        return(-10000)
    }
    logconst = log(factorial(sum(obspos,obsneg)))-
            sum(log(factorial(obspos)),log(factorial(obsneg)))
    ppos[!ppos]=1
    loglikeli = logconst+sum(log(ppos)*obspos,log(pneg)*obsneg,na.rm=T)
    return(loglikeli)
}



########################
#TODO

#function that samples betas
#acquire obspos and obsneg (which gives you Ct)



################################################################################
################################################################################
################################################################################


maxweights = sort(runif(10))
cutoffdecay = 6
lambdahi=.42
lambdalow=0.05

Nneg = 5000
#Npos = sample(10:50,15,replace = T)
Npos = rep(10,15)



obspos = sample(1:10,15,replace = T)
obsneg = 10   

Ct=0
beta=1

####beta is one per year
####Ct is one per week x per year
runmodel <- function(times,beta=rep(0,times),
                     Ct=matrix(1,nrow=52,ncol=times), ...){
    data <- data.frame(week=rep(1:52,each=15),rank=rep(1:15,52))
    
    
    for (tt in 1:times){
        simPos = matrix(NA, nrow = 15, ncol = 52)
        simPos[,1] = c(Npos)
        
        simNeg = rep(NA,52)
        simNeg[1]=Nneg

        for (ii in 2:52){
            temp=onestep(Npos,Nneg, beta[tt], Ct[ii,tt], maxweights, cutoffdecay=6,
                         lambdahi,lambdalow)
            simNeg[ii] = temp$Nneg
            simPos[,ii]= temp$Npos
        }
        
        num = NULL
        for (i in 1:52) num = c(num, simPos[,i])
        data[,(tt+2)]=num
    }
    if (times==1) {data$mean = data[,3]
                   data$density = data$mean/sum(data$mean)
    } else {data$mean = rowMeans(data[,3:(times+2)])
            data$density = data$mean/sum(data[1,3:12])}
    return(data)
}


data <- runmodel(30,Npos=rep(10,15))
####Average of 10 runs of 1 year with no new infections
ggplot(data, aes(week,rank))+
    geom_raster(aes(fill = density))+
    guides(fill=guide_legend(title=NULL))
 
###histogram of prevalence of middle values

####beta is 10 every year
####Ct is 1 per week x per year

####So 10 new infections per week, randomly distributed with
#####reasonable weights
data <- runmodel(10,beta=rep(1,10),Neg=10000)

ggplot(data, aes(week,rank))+
    geom_raster(aes(fill = density))+
    guides(fill=guide_legend(title=NULL))+ 
    scale_fill_gradient(low = "white",
                        high = "steelblue4")  



################################################################################






















Nneg = 100000
#Npos = sample(1:10,15,replace = T)
Npos = rep(2,15)
#maxweights = sort(runif(10))
#maxweights = maxweights/sum(maxweights)

maxweights= c(1,1,2,3,4,6,4,2,1,.5)
maxweights = maxweights/sum(maxweights)

cutoffdecay = 6
lambdahi=.42
lambdalow=0.05

obsneg=1
obspos=rep(1,15)

beta <- runif(1,1,30)
Ct <- 1

newState <- onestep(Npos,Nneg,beta, Ct, 
            maxweights,cutoffdecay,
            lambdahi,lambdalow)

Npos = newState$Npos
Nneg = newState$Nneg

logLike = loglikelifunc(Npos = newState$Npos, 
                        Nneg = newState$Nneg, 
                        obspos, obsneg)

logLike


################################################################################
################################################################################
################################################################################








startDate = "2008-08-25"
endDate = "2008-10-25"
birthday = "June 15"
birthdayweek = week("2015-06-15")

sim <- function(betas,startDate,endDate, byYr=T,
                betaby=c("month","season","hiandlo","days","weeks")){
    
    startDate=as.Date(startDate)
    endDate=as.Date(endDate)
    betaby = match.args(betaby)
    ####see if betas is compatible####
    
    num = endDate-startDate
    
    
    bt = data.frame(NA,ncols=year(endDate)-year(startDate),nrow=num)
    
    years = seq(year(startDate),year(endDate))
    betanum = 0
    
    
    
    
    
        
        
    }

    
    days =as.numeric( endDate-startDate )
    daily=len
    years = seq(year(startDate),year(endDate))
    weekly = week(startDate)
    
    return(len)

}



sim = function( betaby=c(T,F)){
    betaby=match.arg(betaby)
    print(betaby)
}


allbetas = rep(NA,len)
match(allbetas, betas)


################################################################################
################################################################################
################################################################################


sero.model<-function(y,sigma,sigmab,rho,rhob,delta,alpha,Ct,recruits){
    
    N0 = initial$neg+recruits-sigma*Ct
    N1 = #decay from N2 - decay to N0
    N2 = #
    N3 = 
    N4 =
    N5 = 
    N6 = #decay from N7 - decay to N5 + fraction of new infection
    N7 =
    N8 =
    N9 =
    N10 =
    N11 = 
    N12 = 
    N13 = 
    N14 = 
    N15 = #decay to N14 + fraction of new infection density
    
    return(c(N0,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15))
}


run.serology.model<-function(initialvals,
                             sigma1=rep(200,156),sigmab1=10,rho1=0.6,
                             rhob1=0.8,delta1=0.008,alpha1,btime=7,TT=156,Ct1,
                             recruits1,sero.model1=sero.model){
    
    
    for(i in 1:TT){
        y[,i+1]=sero.model1(*************)
    }
    return(y[,-1])
}
