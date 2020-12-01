library(eSIR)
library(devtools)
library(tidyverse)
library(dplyr)
library(rjags)
library(chron)
library(glue)
library(nimble)
#########set seed and load raw data
set.seed(20200601)
load('base.RData')
N_SA=57779622
#############################
############initialization of states
#############################
k=10####ascertained rate%
start_point=which(SA$Date==as.Date('2020/03/15'))
propA0 <- (1-0.01*k)/(0.01*k)
E0 <- (sum(cases_SA[start_point:(start_point+4)])) * (1 + propA0)
R0 <- SA[start_point,3]
D0 <- SA[start_point,2]
H0 <- SA[start_point-7,1]*0.5
if(R0+D0+H0>=SA[start_point,1]){
  I0=0
}else{
  I0 <- SA[start_point,1]-H0-R0-D0
}
A0 <- I0 * propA0
S0 <- N_SA - E0 - I0 - A0 - H0 - R0 - D0
yinit <- c(S = S0, E = E0, I = I0, R = R0, H = H0, A = A0, D = D0)
#############################
####format training data
#############################
time=seq.Date(from=min(SA$Date),to=as.Date("2020/07/31"),by='days')
Y=SA[start_point:length(time),1]
R=SA[start_point:length(time),3]
D=SA[start_point:length(time),2]
Y=Y/N_SA-(R+D)/N_SA
R=(R+D)/N_SA
eps = 1e-10
Y <- pmax(Y,eps)
R <- pmax(R,eps)
#############################
######mcmc setting 
#############################
begin_str="03/15/2020"
T_fin=200
nchain = 4
M=1e5
nburnin = 5e4
thn = 10
nadapt = 10000
len <- round(M/thn)*nchain #number of MCMC draws in total
T_prime <-length(Y)
begin <- chron(dates.=begin_str)
chron_ls <- chron(begin:(begin+T_fin))
end <- chron(begin:(begin+T_fin))[T_fin]
message(paste0("The follow-up is from ",begin," to ",
               end," and the last observed date is ", chron_ls [T_prime],".") )
#############################
####prior
#############################
beta0 = 0.26272
gamma0 = 0.0821
gamma0_sd = 0.1
R0_sd = 1
R0 = beta0/gamma0
beta0 <- R0*gamma0

lognorm.parm <- function(mu0,var0){
  var <- log(var0 / mu0^2 + 1)
  mu <- log(mu0) - var / 2
  list(mu = mu, var = var)
}
gamma_var <- gamma0_sd^2
lognorm_gamma_parm <- lognorm.parm(gamma0,gamma_var)
R0_var <- R0_sd^2
lognorm_R0_parm <- lognorm.parm(R0,R0_var)
pi <- c(rep(1,length(seq.Date(from=as.Date("2020/03/15"),to=as.Date("2020/03/26"),by='days'))),rep(0.75,T_fin))
#############################
################ MCMC ##########
#############################
########specify model
model1.string <- paste0("
  model{
     for(t in 2:(T_prime+1)){
       Km[t-1,1] <- -beta*pi[t-1]*theta[t-1,1]*theta[t-1,2]
       Km[t-1,9] <- gamma*theta[t-1,2]
       Km[t-1,5] <- -Km[t-1,1]-Km[t-1,9]
       Km[t-1,2] <- -beta*pi[t-1]*(theta[t-1,1]+0.5*Km[t-1,1])*(theta[t-1,2]+0.5*Km[t-1,5])
       Km[t-1,10] <- gamma*(theta[t-1,2]+0.5*Km[t-1,5])
       Km[t-1,6] <- -Km[t-1,2]-Km[t-1,10]
       Km[t-1,3] <- -beta*pi[t-1]*(theta[t-1,1]+0.5*Km[t-1,2])*(theta[t-1,2]+0.5*Km[t-1,6])
       Km[t-1,11] <- gamma*(theta[t-1,2]+0.5*Km[t-1,6])
       Km[t-1,7] <- -Km[t-1,3]-Km[t-1,11]
       Km[t-1,4] <- -beta*pi[t-1]*(theta[t-1,1]+Km[t-1,3])*(theta[t-1,2]+Km[t-1,7])
       Km[t-1,12] <- gamma*(theta[t-1,2]+Km[t-1,7])
       Km[t-1,8] <- -Km[t-1,4]-Km[t-1,12]
       alpha[t-1,1] <- theta[t-1,1]+(Km[t-1,1]+2*Km[t-1,2]+2*Km[t-1,3]+Km[t-1,4])/6
       alpha[t-1,2] <- theta[t-1,2]+(Km[t-1,5]+2*Km[t-1,6]+2*Km[t-1,7]+Km[t-1,8])/6
       alpha[t-1,3] <- theta[t-1,3]+(Km[t-1,9]+2*Km[t-1,10]+2*Km[t-1,11]+Km[t-1,12])/6
       theta[t,1:3] ~ ddirch(k*alpha[t-1,1:3])
       Y[t-1] ~ dbeta(lambdaY*theta[t,2],lambdaY*(1-theta[t,2]))
       R[t-1] ~ dbeta(lambdaR*theta[t,3],lambdaR*(1-theta[t,3]))
     }
    theta0[1:3]<-c(",1-Y[1]-R[1],",",Y[1],",", R[1],")
    theta[1,1:3] ~ ddirch(theta0[1:3])
    gamma ~  dlnorm(", lognorm_gamma_parm$mu, ",", 1 / lognorm_gamma_parm$var,")
    R0 ~ dlnorm(", lognorm_R0_parm$mu, ",", 1 / lognorm_R0_parm$var,")
    beta <- R0*gamma
    k ~  dgamma(2,0.0001)
    lambdaY ~ dgamma(2,0.0001)
    lambdaR ~ dgamma(2,0.0001)
  }
")

model.spec <- textConnection(model1.string)

posterior <- jags.model(
  model.spec,
  data = list(
    'Y' = Y,
    'R' = R,
    'T_prime' = T_prime,
    'pi' = pi
  ),
  n.chains = nchain,
  n.adapt = nadapt
)

#burn-in
update(posterior, nburnin)
#sampling
jags_sample <- jags.samples(
  posterior,
  c('theta','gamma','R0','beta','Y','lambdaY','lambdaR','k'),
  n.iter = M,
  thin = thn
)
###################################
#############organize the posterior of parameters
#############################
R0_p <- unlist(as.mcmc.list(jags_sample$R0))
gamma_p <- unlist(as.mcmc.list(jags_sample$gamma))
beta_p <- unlist(as.mcmc.list(jags_sample$beta))
lambdaY_p <- unlist(as.mcmc.list(jags_sample$lambdaY))
lambdaR_p <- unlist(as.mcmc.list(jags_sample$lambdaR))
k_p <- unlist(as.mcmc.list(jags_sample$k))
theta_p <- array(Reduce(rbind,as.mcmc.list(jags_sample$theta)),dim=c(len,T_prime+1,3))

############posterior means and 95% credible intervals
theta_p_mean <- apply(theta_p[,T_prime+1,],2,mean)
theta_p_ci <- as.vector(apply(theta_p[,T_prime+1,],2,quantile,c(0.025,0.5,0.975)))
R0_p_mean <- mean(R0_p)
R0_p_ci <- quantile(R0_p,c(0.025,0.5,0.975))

gamma_p_mean <- mean(gamma_p)
gamma_p_ci <- quantile(gamma_p,c(0.025,0.5,0.975))

beta_p_mean <- mean(beta_p)
beta_p_ci <- quantile(beta_p,c(0.025,0.5,0.975))

lambdaY_p_mean <- mean(lambdaY_p)
lambdaY_p_ci <- quantile(lambdaY_p,c(0.025,0.5,0.975))

lambdaR_p_mean <- mean(lambdaR_p)
lambdaR_p_ci <- quantile(lambdaR_p,c(0.025,0.5,0.975))

k_p_mean <- mean(k_p)
k_p_ci <- quantile(k_p,c(0.025,0.5,0.975))

###########################
#### Forecast ####
#############################
theta_pp <- array(0,dim=c(len,T_fin-T_prime,3))
Y_pp <- matrix(NA,nrow=len,ncol=T_fin-T_prime)
R_pp <- matrix(NA,nrow=len,ncol=T_fin-T_prime)
for(l in 1:len){
  thetalt1 <- theta_p[l,T_prime+1,1]
  thetalt2 <- theta_p[l,T_prime+1,2]
  thetalt3 <- theta_p[l,T_prime+1,3]
  betal <- c(beta_p)[l]
  gammal <- c(gamma_p)[l]
  kt <- c(k_p)[l]
  lambdaYl <- c(lambdaY_p)[l]
  lambdaRl <- c(lambdaR_p)[l]
  if(betal<0 |gammal<0 |thetalt1<0 |thetalt2<0 |thetalt3<0) next
  for(t in 1:(T_fin-T_prime)){
    Km <- NULL
    alpha_pp <- NULL
    Km[1] <- -betal*pi[t+T_prime]*thetalt1*thetalt2
    Km[9] <- gammal*thetalt2
    Km[5] <- -Km[1]-Km[9]
    
    Km[2] <- -betal*pi[t+T_prime]*(thetalt1+0.5*Km[1])*(thetalt2+0.5*Km[5])
    Km[10] <- gammal*(thetalt2+0.5*Km[5])
    Km[6] <- -Km[2]-Km[10]
    
    Km[3] <- -betal*pi[t+T_prime]*(thetalt1+0.5*Km[2])*(thetalt2+0.5*Km[6])
    Km[11] <- gammal*(thetalt2+0.5*Km[6])
    Km[7] <- -Km[3]-Km[11]
    
    Km[4] <- -betal*pi[t+T_prime]*(thetalt1+Km[3])*(thetalt2+Km[7])
    Km[12] <- gammal*(thetalt2+Km[7])
    Km[8] <- -Km[4]-Km[12]
    
    alpha_pp[1] <- thetalt1+(Km[1]+2*Km[2]+2*Km[3]+Km[4])/6
    alpha_pp[2] <- thetalt2+(Km[5]+2*Km[6]+2*Km[7]+Km[8])/6
    alpha_pp[3] <- thetalt3+(Km[9]+2*Km[10]+2*Km[11]+Km[12])/6
    
    thetalt_tmp <- rdirichlet(1,kt*c(alpha_pp))
    thetalt1<-theta_pp[l,t,1] <- thetalt_tmp[1]
    thetalt2<-theta_pp[l,t,2] <- thetalt_tmp[2]
    thetalt3<-theta_pp[l,t,3] <- thetalt_tmp[3]
    
    Y_pp[l,t] <- rbeta(1,lambdaYl*thetalt2,lambdaYl*(1-thetalt2))
    R_pp[l,t] <- rbeta(1,lambdaRl*thetalt3,lambdaRl*(1-thetalt3))
  }
}
#############################
#### save image ####
#############################
save.image('SIR.RData')


