###################################################################
# Multiple Frogs (15) and Toads (15) Simulation
###################################################################

#-------------------------------------- setup --------------------------------------------------

n.days <- 30     #no.of days in simulation
#sigma <- .5/30   # rate of zoospore exposure due to contact with zoospores in water

max.z.frog <- 410000					# max amt of zoospores frog can tolerate
max.z.toad <- 350000        # max amt of zoospores toad can tolerate

#sf <- 0.8     #shedding rate of zoospores for frogs
#st <- 0.35     #shedding rate of zoospores for toads
sf <- .5
st <- 0.073
  
  
df <- 0.3        # death rate of zoospores on frogs
dt <- 0.3        # death rate of zoospores on toads

rf <- 2        #growth rate of zoospores on host 
rt <- 1.6

prob.die <- .9  #probability leftover zoospores will die in pond
prob.nat.mort <- 0.005   # probability amphibian will die of natural mortality

P <- rep(NA, length = n.days)     #set up empty vector for pond zoospores
P[1] <- 1000      #set starting zoospore load in pond



#----- Frogs Matrix
Frogs <- matrix(NA,nrow=n.days+1,ncol=15)
rownames(Frogs) <- paste("t", 1:31,sep="")
colnames(Frogs) <- paste("F", 1:15,sep="")
Frogs[1,] <- 0
#Frogs[1,2] <- 1600

#----- Toads Matrix
Toads <- matrix(NA,nrow=n.days+1,ncol=15)
rownames(Toads) <- paste("t",1:31,sep="")
colnames(Toads) <- paste("T",1:15,sep="")
Toads[1,] <- 0
#Toads[1,2] <- 1600

#----- Zoospore Probabilities Matrix to be applied to multinomial matrix
z.probs <- matrix(NA, nrow = n.days +1, ncol = 3)       
rownames(z.probs) <- paste("t", 1:31,sep="")
colnames(z.probs) <- c("Else","F","T")

#------ Zoospore Multinomial Matrix to Determine where zoospores in P go
where.go.zoos <- matrix(NA, nrow = n.days +1, ncol = 3)
rownames(where.go.zoos) <- paste("t",1:31,sep="")
colnames(where.go.zoos) <- c("Else","F","T")

#zoos.to.stay <- rep(NA, length = n.days+1)
#zoos.to.move.F <- rep(NA, length = n.days+1)
#zoos.to.move.T <-  rep(NA, length = n.days+1)
#zoos.to.die <-  rep(NA, length = n.days+1)

# viable frogs and toads (create vectors of viable frogs and toads that aren't dead)
F <- rep(NA, length = n.days+1)
T <- rep(NA, length = n.days+1)
mean.int.f <- rep(NA, length = n.days+1)
mean.int.t <- rep(NA, length = n.days+1)

sigma <- rep(NA, length = n.days+1)









for(t in 2:(n.days+1)){


pos.F <- which(!is.na(Frogs[t-1,]))
F[t] <- length(which(!is.na(Frogs[t-1,])))

pos.T <- which(!is.na(Toads[t-1,]))
T[t] <- length(which(!is.na(Toads[t-1,])))

sigma[t] <- .5/((F[t] + T[t])*2)
# Second, calculate probabilities of zoospores to move in z.probs with pos.F and pos.T


z.probs[t,] <- c(1-(sigma[t]*F[t] + sigma[t]*T[t]), sigma[t]*F[t], sigma[t]*T[t])
   

# Third, calculate zoos.to.move to  F and T and zoos.to.stay and zoos.to.die in P

where.go.zoos <- rmultinom(1, size = P[t-1], z.probs[t,])

zoos.to.move.F <- where.go.zoos[2]
zoos.to.move.T <- where.go.zoos[3]

# of those that didn't infect a frog or toad, do they live/stay or die?
zoos.left <- where.go.zoos[1]
#zoos.die <- rbinom(1, zoos.left, prob.die) # rest die and are not kept track of
zoos.survive <- rbinom(1, zoos.left, 1- prob.die) 
# Fourth, figure out who.T and who.F 

who.F <- sample(pos.F, size = zoos.to.move.F, replace = TRUE)
who.T <- sample(pos.T, size = zoos.to.move.T, replace = TRUE)

  
   
   P[t] <- sum(sf * Frogs[t-1,],na.rm=TRUE) + sum(st * Toads[t-1,],na.rm=TRUE) + zoos.survive    ## OR  P[t] <- P[t-1] + sf * Frogs[t-1] + st * Toads[t-1] - zoos.to.die[t] - zoos.to.move.F[t] - zoos.to.move.T[t]
   mean.int.f[t] <- mean(Frogs[t-1,], na.rm = TRUE)
   mean.int.t[t] <- mean(Toads[t-1,], na.rm = TRUE)
   
    for(i in 1:15){

      Fnewzoos <- rf*Frogs[t-1,i] - sf* Frogs[t-1,i] - df* Frogs[t-1,i] + length(which(who.F==i))
    Frogs[t,i] <- ifelse(rbinom(1, 1, prob.nat.mort)==1, NA,
                    ifelse(  Fnewzoos > max.z.frog, NA,
                      ifelse(Fnewzoos < 0, 0, Fnewzoos)))
                            
    
    
      Tnewzoos <- rt*Toads[t-1,i] - st* Toads[t-1,i] - dt* Toads[t-1,i] + length(which(who.T==i))
    Toads[t,i] <- ifelse(rbinom(1, 1, prob.nat.mort)==1, NA,                  
                  ifelse(  Tnewzoos > max.z.toad, NA, 
                         ifelse(Tnewzoos < 0, 0, Tnewzoos)))
        
   
    }
}
    


P
Frogs
Toads
F
T

mean.int.t
mean.int.f
library(dplyr)
library(tidyr)

ToadToad <- as.data.frame(Toads)
Froggy <- as.data.frame(Frogs)
Pond <- as.data.frame(P)


Ponds <- Pond %>%
  rename(Zoospores= "P") %>%
  mutate(timestep = 1:31)%>%
  mutate(type = "Pond") 
  
  
FrogLong <- Froggy %>%
  gather("Amphib", "Zoospores") %>%
  mutate(timestep = rep(1:31, 15)) %>%
  mutate(type = "Frog")


Amphed <- ToadToad %>%
  gather("Amphib", "Zoospores") %>%
  mutate(timestep = rep(1:31, 15)) %>%
  mutate(type = "Toad") %>%
  bind_rows(FrogLong) %>%
  bind_rows(Ponds)

Toad2 <- ToadLong %>%
  group_by(timestep) %>%
  summarize(mean = mean(Zoospores, na.rm = TRUE), 
            sd = sd(Zoospores, na.rm = TRUE))


ggplot(Amphed, aes(x = timestep, y = Zoospores, col = type))+
  geom_point()+
  stat_smooth()+
  theme_classic()+
  theme(legend.position = c(.2,.8))


plot(1:31, P, xlab = "Days", ylab = "No. of Zoospores", col = "blue", type = "b", pch = 20, main = "15 Frogs + 15 Toads", ylim = )
lines(0:30, mean.int.f, col = "green",type = "b", pch = 20)
lines(0:30, mean.int.t, col = "orange",type = "b", pch = 20)
legend("topleft", legend = c("Pond", "Frog", "Toad"), col = c("blue", "green", "orange"), lty = 1)

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#  FROGS ONLY (30)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------- setup 

n.days <- 30     #no.of days in simulation
#sigma <- .5/30   # rate of zoospore exposure due to contact with zoospores in water

max.z.frog <- 410000					# max amt of zoospores frog can tolerate
max.z.toad <- 350000        # max amt of zoospores toad can tolerate

#sf <- 0.8     #shedding rate of zoospores for frogs
#st <- 0.35     #shedding rate of zoospores for toads
sf <- .5
st <- 0.073


df <- 0.3        # death rate of zoospores on frogs
dt <- 0.3        # death rate of zoospores on toads

rf <- 2        #growth rate of zoospores on host 
rt <- 1.6

prob.die <- .9  #probability leftover zoospores will die in pond
prob.nat.mort <- 0.005   # probability amphibian will die of natural mortality

P <- rep(NA, length = n.days)     #set up empty vector for pond zoospores
P[1] <- 1000      #set starting zoospore load in pond



#----- Frogs Matrix
Frogs <- matrix(NA,nrow=n.days+1,ncol=30)
rownames(Frogs) <- paste("t", 1:31,sep="")
colnames(Frogs) <- paste("F", 1:30,sep="")
Frogs[1,] <- 0
#Frogs[1,2] <- 1600

#----- Toads Matrix
Toads <- matrix(NA,nrow=n.days+1,ncol=30)
rownames(Toads) <- paste("t",1:31,sep="")
colnames(Toads) <- paste("T",1:30,sep="")
Toads[1,] <- NA
#Toads[1,2] <- 1600

#----- Zoospore Probabilities Matrix to be applied to multinomial matrix
z.probs <- matrix(NA, nrow = n.days +1, ncol = 3)       
rownames(z.probs) <- paste("t", 1:31,sep="")
colnames(z.probs) <- c("Else","F","T")

#------ Zoospore Multinomial Matrix to Determine where zoospores in P go
where.go.zoos <- matrix(NA, nrow = n.days +1, ncol = 3)
rownames(where.go.zoos) <- paste("t",1:31,sep="")
colnames(where.go.zoos) <- c("Else","F","T")

# viable frogs and toads (create vectors of viable frogs and toads that aren't dead)
F <- rep(NA, length = n.days+1)
T <- rep(NA, length = n.days+1)
mean.int.f <- rep(NA, length = n.days+1)
mean.int.t <- rep(NA, length = n.days+1)
sigma <- rep(NA, length = n.days+1)




for(t in 2:(n.days+1)){
  
  
  pos.F <- which(!is.na(Frogs[t-1,]))
  F[t] <- length(which(!is.na(Frogs[t-1,])))
  
  pos.T <- which(!is.na(Toads[t-1,]))
  T[t] <- length(which(!is.na(Toads[t-1,])))
  
  sigma[t] <- .5/((F[t] + T[t])*2)
  # Second, calculate probabilities of zoospores to move in z.probs with pos.F and pos.T
  
  
  z.probs[t,] <- c(1-(sigma[t]*F[t] + sigma[t]*T[t]), sigma[t]*F[t], sigma[t]*T[t])
  
  
  # Third, calculate zoos.to.move to  F and T and zoos.to.stay and zoos.to.die in P
  
  where.go.zoos <- rmultinom(1, size = P[t-1], z.probs[t,])
  
  zoos.to.move.F <- where.go.zoos[2]
  zoos.to.move.T <- where.go.zoos[3]
  
  # of those that didn't infect a frog or toad, do they live/stay or die?
  zoos.left <- where.go.zoos[1]
  #zoos.die <- rbinom(1, zoos.left, prob.die) # rest die and are not kept track of
  zoos.survive <- rbinom(1, zoos.left, 1- prob.die) 
  # Fourth, figure out who.T and who.F 
  
  who.F <- sample(pos.F, size = zoos.to.move.F, replace = TRUE)
  who.T <- sample(pos.T, size = zoos.to.move.T, replace = TRUE)
  
  
  
  P[t] <- sum(sf * Frogs[t-1,],na.rm=TRUE) + sum(st * Toads[t-1,],na.rm=TRUE) + zoos.survive    ## OR  P[t] <- P[t-1] + sf * Frogs[t-1] + st * Toads[t-1] - zoos.to.die[t] - zoos.to.move.F[t] - zoos.to.move.T[t]
  mean.int.f[t] <- mean(Frogs[t-1,], na.rm = TRUE)
  mean.int.t[t] <- mean(Toads[t-1,], na.rm = TRUE)
  
  for(i in 1:30){
    
    Fnewzoos <- rf*Frogs[t-1,i] - sf* Frogs[t-1,i] - df* Frogs[t-1,i] + length(which(who.F==i))
    Frogs[t,i] <- ifelse(rbinom(1, 1, prob.nat.mort)==1, NA,
                         ifelse(  Fnewzoos > max.z.frog, NA,
                                  ifelse(Fnewzoos < 0, 0, Fnewzoos)))
    
    
    
    Tnewzoos <- rt*Toads[t-1,i] - st* Toads[t-1,i] - dt* Toads[t-1,i] + length(which(who.T==i))
    Toads[t,i] <- ifelse(rbinom(1, 1, prob.nat.mort)==1, NA,                  
                         ifelse(  Tnewzoos > max.z.toad, NA, 
                                  ifelse(Tnewzoos < 0, 0, Tnewzoos)))
    
    
  }
}



P
Frogs
Toads
F
T

mean.int.t
mean.int.f




plot(1:31, P, xlab = "Days", ylab = "No. of Zoospores", col = "blue", type = "b", pch = 20, main = "Only 30 Frogs")
lines(0:30, mean.int.f, col = "green",type = "b", pch = 20)
lines(0:30, mean.int.t, col = "orange",type = "b", pch = 20)
legend("topleft", legend = c("Pond", "Frog"), col = c("blue", "green"), lty = 1)

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#  TOADS ONLY (30)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------- setup 

n.days <- 30     #no.of days in simulation
#sigma <- .5/30   # rate of zoospore exposure due to contact with zoospores in water

max.z.frog <- 410000					# max amt of zoospores frog can tolerate
max.z.toad <- 350000        # max amt of zoospores toad can tolerate

#sf <- 0.8     #shedding rate of zoospores for frogs
#st <- 0.35     #shedding rate of zoospores for toads
sf <- .5
st <- 0.073


df <- 0.3        # death rate of zoospores on frogs
dt <- 0.3        # death rate of zoospores on toads

rf <- 2        #growth rate of zoospores on host 
rt <- 1.6

prob.die <- .9  #probability leftover zoospores will die in pond
prob.nat.mort <- 0.005   # probability amphibian will die of natural mortality

P <- rep(NA, length = n.days)     #set up empty vector for pond zoospores
P[1] <- 1000      #set starting zoospore load in pond



#----- Frogs Matrix
Frogs <- matrix(NA,nrow=n.days+1,ncol=30)
rownames(Frogs) <- paste("t", 1:31,sep="")
colnames(Frogs) <- paste("F", 1:30,sep="")
Frogs[1,] <- NA
#Frogs[1,2] <- 1600

#----- Toads Matrix
Toads <- matrix(NA,nrow=n.days+1,ncol=30)
rownames(Toads) <- paste("t",1:31,sep="")
colnames(Toads) <- paste("T",1:30,sep="")
Toads[1,] <-0
#Toads[1,2] <- 1600

#----- Zoospore Probabilities Matrix to be applied to multinomial matrix
z.probs <- matrix(NA, nrow = n.days +1, ncol = 3)       
rownames(z.probs) <- paste("t", 1:31,sep="")
colnames(z.probs) <- c("Else","F","T")

#------ Zoospore Multinomial Matrix to Determine where zoospores in P go
where.go.zoos <- matrix(NA, nrow = n.days +1, ncol = 3)
rownames(where.go.zoos) <- paste("t",1:31,sep="")
colnames(where.go.zoos) <- c("Else","F","T")

# viable frogs and toads (create vectors of viable frogs and toads that aren't dead)
F <- rep(NA, length = n.days+1)
T <- rep(NA, length = n.days+1)
mean.int.f <- rep(NA, length = n.days+1)
mean.int.t <- rep(NA, length = n.days+1)
sigma <- rep(NA, length = n.days+1)




for(t in 2:(n.days+1)){
  
  
  pos.F <- which(!is.na(Frogs[t-1,]))
  F[t] <- length(which(!is.na(Frogs[t-1,])))
  
  pos.T <- which(!is.na(Toads[t-1,]))
  T[t] <- length(which(!is.na(Toads[t-1,])))
  
  sigma[t] <- .5/((F[t] + T[t])*2)
  # Second, calculate probabilities of zoospores to move in z.probs with pos.F and pos.T
  
  
  z.probs[t,] <- c(1-(sigma[t]*F[t] + sigma[t]*T[t]), sigma[t]*F[t], sigma[t]*T[t])
  
  
  # Third, calculate zoos.to.move to  F and T and zoos.to.stay and zoos.to.die in P
  
  where.go.zoos <- rmultinom(1, size = P[t-1], z.probs[t,])
  
  zoos.to.move.F <- where.go.zoos[2]
  zoos.to.move.T <- where.go.zoos[3]
  
  # of those that didn't infect a frog or toad, do they live/stay or die?
  zoos.left <- where.go.zoos[1]
  #zoos.die <- rbinom(1, zoos.left, prob.die) # rest die and are not kept track of
  zoos.survive <- rbinom(1, zoos.left, 1- prob.die) 
  # Fourth, figure out who.T and who.F 
  
  who.F <- sample(pos.F, size = zoos.to.move.F, replace = TRUE)
  who.T <- sample(pos.T, size = zoos.to.move.T, replace = TRUE)
  
  
  
  P[t] <- sum(sf * Frogs[t-1,],na.rm=TRUE) + sum(st * Toads[t-1,],na.rm=TRUE) + zoos.survive    ## OR  P[t] <- P[t-1] + sf * Frogs[t-1] + st * Toads[t-1] - zoos.to.die[t] - zoos.to.move.F[t] - zoos.to.move.T[t]
  mean.int.f[t] <- mean(Frogs[t-1,], na.rm = TRUE)
  mean.int.t[t] <- mean(Toads[t-1,], na.rm = TRUE)
  
  for(i in 1:30){
    
    Fnewzoos <- rf*Frogs[t-1,i] - sf* Frogs[t-1,i] - df* Frogs[t-1,i] + length(which(who.F==i))
    Frogs[t,i] <- ifelse(rbinom(1, 1, prob.nat.mort)==1, NA,
                         ifelse(  Fnewzoos > max.z.frog, NA,
                                  ifelse(Fnewzoos < 0, 0, Fnewzoos)))
    
    
    
    Tnewzoos <- rt*Toads[t-1,i] - st* Toads[t-1,i] - dt* Toads[t-1,i] + length(which(who.T==i))
    Toads[t,i] <- ifelse(rbinom(1, 1, prob.nat.mort)==1, NA,                  
                         ifelse(  Tnewzoos > max.z.toad, NA, 
                                  ifelse(Tnewzoos < 0, 0, Tnewzoos)))
    
    
  }
}



P
Frogs
Toads
F
T

mean.int.t
mean.int.f




plot(1:31, P, xlab = "Days", ylab = "No. of Zoospores", col = "blue", type = "b", pch = 20, main = "Only 30 Toads")
lines(0:30, mean.int.f, col = "green",type = "b", pch = 20)
lines(0:30, mean.int.t, col = "orange",type = "b", pch = 20)
legend("topleft", legend = c("Pond", "Toad"), col = c("blue", "orange"), lty = 1)


#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#  MOSTLY FROGS (25) and 5 TOADS
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------- setup 

n.days <- 30     #no.of days in simulation
#sigma <- .5/30   # rate of zoospore exposure due to contact with zoospores in water

max.z.frog <- 410000					# max amt of zoospores frog can tolerate
max.z.toad <- 350000        # max amt of zoospores toad can tolerate

#sf <- 0.8     #shedding rate of zoospores for frogs
#st <- 0.35     #shedding rate of zoospores for toads
sf <- .5
st <- 0.073


df <- 0.3        # death rate of zoospores on frogs
dt <- 0.3        # death rate of zoospores on toads

rf <- 2        #growth rate of zoospores on host 
rt <- 1.6

prob.die <- .9  #probability leftover zoospores will die in pond
prob.nat.mort <- 0.005   # probability amphibian will die of natural mortality

P <- rep(NA, length = n.days)     #set up empty vector for pond zoospores
P[1] <- 1000      #set starting zoospore load in pond



#----- Frogs Matrix
Frogs <- matrix(NA,nrow=n.days+1,ncol=30)
rownames(Frogs) <- paste("t", 1:31,sep="")
colnames(Frogs) <- paste("F", 1:30,sep="")
Frogs[1,] <- 0 
Frogs[1,1:5] <- NA

#----- Toads Matrix
Toads <- matrix(NA,nrow=n.days+1,ncol=30)
rownames(Toads) <- paste("t",1:31,sep="")
colnames(Toads) <- paste("T",1:30,sep="")
Toads[1,] <- NA
Toads[1,1:5] <- 0

#----- Zoospore Probabilities Matrix to be applied to multinomial matrix
z.probs <- matrix(NA, nrow = n.days +1, ncol = 3)       
rownames(z.probs) <- paste("t", 1:31,sep="")
colnames(z.probs) <- c("Else","F","T")

#------ Zoospore Multinomial Matrix to Determine where zoospores in P go
where.go.zoos <- matrix(NA, nrow = n.days +1, ncol = 3)
rownames(where.go.zoos) <- paste("t",1:31,sep="")
colnames(where.go.zoos) <- c("Else","F","T")

# viable frogs and toads (create vectors of viable frogs and toads that aren't dead)
F <- rep(NA, length = n.days+1)
T <- rep(NA, length = n.days+1)
mean.int.f <- rep(NA, length = n.days+1)
mean.int.t <- rep(NA, length = n.days+1)
sigma <- rep(NA, length = n.days+1)




for(t in 2:(n.days+1)){
  
  
  pos.F <- which(!is.na(Frogs[t-1,]))
  F[t] <- length(which(!is.na(Frogs[t-1,])))
  
  pos.T <- which(!is.na(Toads[t-1,]))
  T[t] <- length(which(!is.na(Toads[t-1,])))
  
  sigma[t] <- .5/((F[t] + T[t])*2)
  # Second, calculate probabilities of zoospores to move in z.probs with pos.F and pos.T
  
  
  z.probs[t,] <- c(1-(sigma[t]*F[t] + sigma[t]*T[t]), sigma[t]*F[t], sigma[t]*T[t])
  
  
  # Third, calculate zoos.to.move to  F and T and zoos.to.stay and zoos.to.die in P
  
  where.go.zoos <- rmultinom(1, size = P[t-1], z.probs[t,])
  
  zoos.to.move.F <- where.go.zoos[2]
  zoos.to.move.T <- where.go.zoos[3]
  
  # of those that didn't infect a frog or toad, do they live/stay or die?
  zoos.left <- where.go.zoos[1]
  #zoos.die <- rbinom(1, zoos.left, prob.die) # rest die and are not kept track of
  zoos.survive <- rbinom(1, zoos.left, 1- prob.die) 
  # Fourth, figure out who.T and who.F 
  
  who.F <- sample(pos.F, size = zoos.to.move.F, replace = TRUE)
  who.T <- sample(pos.T, size = zoos.to.move.T, replace = TRUE)
  
  
  
  P[t] <- sum(sf * Frogs[t-1,],na.rm=TRUE) + sum(st * Toads[t-1,],na.rm=TRUE) + zoos.survive    ## OR  P[t] <- P[t-1] + sf * Frogs[t-1] + st * Toads[t-1] - zoos.to.die[t] - zoos.to.move.F[t] - zoos.to.move.T[t]
  mean.int.f[t] <- mean(Frogs[t-1,], na.rm = TRUE)
  mean.int.t[t] <- mean(Toads[t-1,], na.rm = TRUE)
  
  for(i in 1:30){
    
    Fnewzoos <- rf*Frogs[t-1,i] - sf* Frogs[t-1,i] - df* Frogs[t-1,i] + length(which(who.F==i))
    Frogs[t,i] <- ifelse(rbinom(1, 1, prob.nat.mort)==1, NA,
                         ifelse(  Fnewzoos > max.z.frog, NA,
                                  ifelse(Fnewzoos < 0, 0, Fnewzoos)))
    
    
    
    Tnewzoos <- rt*Toads[t-1,i] - st* Toads[t-1,i] - dt* Toads[t-1,i] + length(which(who.T==i))
    Toads[t,i] <- ifelse(rbinom(1, 1, prob.nat.mort)==1, NA,                  
                         ifelse(  Tnewzoos > max.z.toad, NA, 
                                  ifelse(Tnewzoos < 0, 0, Tnewzoos)))
    
    
  }
}



P
Frogs
Toads
F
T



mean.int.t
mean.int.f


i = 1:30

plot(1:31, Toads[1:31,1], pch = 20, col = "orange")
points(1:31, Toads[1:31,2], pch = 20, col = "orange")
points(1:31, Toads[1:31,3], pch = 20, col = "orange")

plot(1:31, mean.int.t, type = "l", col = "red")
lines(Toads[1:31,1], col = "orange")

plot(1:31, Toads[,T1:T30])

matplot(output$level, output[,-1], t="l", lty=1)
library(ggplot2) ; qplot(level, value, colour=variable, geom="line", data=melt(output, id="level"))

t.df <- data.frame(Toads)

plot(t ~ T, t.df)


plot(1:31, Toads[1:31,1:30])
plot(NA, NA, ylim = c(0,7000), xlim = c(1,31))
points(length(Toads))
plot(1:30, Toads[t,i])


plot(1:31, P, xlab = "Days", ylab = "No. of Zoospores", col = "blue", type = "b", pch = 20, main = "25 Frogs + 5 Toads")
lines(0:30, mean.int.f, col = "green",type = "b", pch = 20)
lines(0:30, mean.int.t, col = "orange",type = "b", pch = 20)
legend("topleft", legend = c("Pond", "Frog" ,"Toad"), col = c("blue", "green", "orange"), lty = 1)

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#  MOSTLY TOADS (25) and 5 FROGS
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------- setup 

n.days <- 30     #no.of days in simulation
#sigma <- .5/30   # rate of zoospore exposure due to contact with zoospores in water

max.z.frog <- 410000					# max amt of zoospores frog can tolerate
max.z.toad <- 350000        # max amt of zoospores toad can tolerate

#sf <- 0.8     #shedding rate of zoospores for frogs
#st <- 0.35     #shedding rate of zoospores for toads
sf <- .5
st <- 0.073


df <- 0.3        # death rate of zoospores on frogs
dt <- 0.3        # death rate of zoospores on toads

rf <- 2        #growth rate of zoospores on host 
rt <- 1.6

prob.die <- .9  #probability leftover zoospores will die in pond
prob.nat.mort <- 0.005   # probability amphibian will die of natural mortality

P <- rep(NA, length = n.days)     #set up empty vector for pond zoospores
P[1] <- 1000      #set starting zoospore load in pond



#----- Frogs Matrix
Frogs <- matrix(NA,nrow=n.days+1,ncol=30)
rownames(Frogs) <- paste("t", 1:31,sep="")
colnames(Frogs) <- paste("F", 1:30,sep="")
Frogs[1,] <- NA 
Frogs[1,1:5] <- 0

#----- Toads Matrix
Toads <- matrix(NA,nrow=n.days+1,ncol=30)
rownames(Toads) <- paste("t",1:31,sep="")
colnames(Toads) <- paste("T",1:30,sep="")
Toads[1,] <- 0
Toads[1,1:5] <- NA

#----- Zoospore Probabilities Matrix to be applied to multinomial matrix
z.probs <- matrix(NA, nrow = n.days +1, ncol = 3)       
rownames(z.probs) <- paste("t", 1:31,sep="")
colnames(z.probs) <- c("Else","F","T")

#------ Zoospore Multinomial Matrix to Determine where zoospores in P go
where.go.zoos <- matrix(NA, nrow = n.days +1, ncol = 3)
rownames(where.go.zoos) <- paste("t",1:31,sep="")
colnames(where.go.zoos) <- c("Else","F","T")

# viable frogs and toads (create vectors of viable frogs and toads that aren't dead)
F <- rep(NA, length = n.days+1)
T <- rep(NA, length = n.days+1)
mean.int.f <- rep(NA, length = n.days+1)
mean.int.t <- rep(NA, length = n.days+1)
sigma <- rep(NA, length = n.days+1)




for(t in 2:(n.days+1)){
  
  
  pos.F <- which(!is.na(Frogs[t-1,]))
  F[t] <- length(which(!is.na(Frogs[t-1,])))
  
  pos.T <- which(!is.na(Toads[t-1,]))
  T[t] <- length(which(!is.na(Toads[t-1,])))
  
  sigma[t] <- .5/((F[t] + T[t])*2)
  # Second, calculate probabilities of zoospores to move in z.probs with pos.F and pos.T
  
  
  z.probs[t,] <- c(1-(sigma[t]*F[t] + sigma[t]*T[t]), sigma[t]*F[t], sigma[t]*T[t])
  
  
  # Third, calculate zoos.to.move to  F and T and zoos.to.stay and zoos.to.die in P
  
  where.go.zoos <- rmultinom(1, size = P[t-1], z.probs[t,])
  
  zoos.to.move.F <- where.go.zoos[2]
  zoos.to.move.T <- where.go.zoos[3]
  
  # of those that didn't infect a frog or toad, do they live/stay or die?
  zoos.left <- where.go.zoos[1]
  #zoos.die <- rbinom(1, zoos.left, prob.die) # rest die and are not kept track of
  zoos.survive <- rbinom(1, zoos.left, 1- prob.die) 
  # Fourth, figure out who.T and who.F 
  
  who.F <- sample(pos.F, size = zoos.to.move.F, replace = TRUE)
  who.T <- sample(pos.T, size = zoos.to.move.T, replace = TRUE)
  
  
  
  P[t] <- sum(sf * Frogs[t-1,],na.rm=TRUE) + sum(st * Toads[t-1,],na.rm=TRUE) + zoos.survive    ## OR  P[t] <- P[t-1] + sf * Frogs[t-1] + st * Toads[t-1] - zoos.to.die[t] - zoos.to.move.F[t] - zoos.to.move.T[t]
  mean.int.f[t] <- mean(Frogs[t-1,], na.rm = TRUE)
  mean.int.t[t] <- mean(Toads[t-1,], na.rm = TRUE)
  
  for(i in 1:30){
    
    Fnewzoos <- rf*Frogs[t-1,i] - sf* Frogs[t-1,i] - df* Frogs[t-1,i] + length(which(who.F==i))
    Frogs[t,i] <- ifelse(rbinom(1, 1, prob.nat.mort)==1, NA,
                         ifelse(  Fnewzoos > max.z.frog, NA,
                                  ifelse(Fnewzoos < 0, 0, Fnewzoos)))
    
    
    
    Tnewzoos <- rt*Toads[t-1,i] - st* Toads[t-1,i] - dt* Toads[t-1,i] + length(which(who.T==i))
    Toads[t,i] <- ifelse(rbinom(1, 1, prob.nat.mort)==1, NA,                  
                         ifelse(  Tnewzoos > max.z.toad, NA, 
                                  ifelse(Tnewzoos < 0, 0, Tnewzoos)))
    
    
  }
}



P
Frogs
Toads
F
T

mean.int.t
mean.int.f

rowMeans(Toads, na.rm = TRUE)



plot(1:31, P, xlab = "Days", ylab = "No. of Zoospores", col = "blue", type = "b", pch = 20, main = "5 Frogs + 25 Toads")
lines(0:30, mean.int.f, col = "green",type = "b", pch = 20)
lines(0:30, mean.int.t, col = "orange",type = "b", pch = 20)
legend("topleft", legend = c("Pond", "Frog" ,"Toad"), col = c("blue", "green", "orange"), lty = 1)

