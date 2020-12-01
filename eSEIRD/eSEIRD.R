library(devtools)
library(tidyverse)
library(dplyr)
library(rjags)
library(chron)
library(glue)
#####################
######################set seed and load raw data
####################
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
data_SA=cbind.data.frame(cases_SA,recovered_SA,death_SA)[c(start_point:length(time)),]

Y=data_SA$cases_SA
R=data_SA$recovered_SA
D=data_SA$death_SA
#############################
######mcmc setting 
#############################
M=1e5
thn=10
nburnin=5e4
nchain=4
nadapt=1e4
T_fin=200
begin_str="03/15/2020"
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
R0=3.2
R0_sd=1
lognorm.parm <- function(mu0,var0){
  var <- log(var0 / mu0^2 + 1)
  mu <- log(mu0) - var / 2
  list(mu = mu, var = var)
}
R0_var <- R0_sd^2
lognorm_R0_parm <- lognorm.parm(R0,R0_var)
pi <- c(rep(1,length(seq.Date(from=as.Date("2020/03/15"),to=as.Date("2020/03/26"),by='days'))),rep(0.75,T_fin))
n <- c(rep(0.0004*N_SA,length(seq.Date(from=as.Date("2020/03/15"),to=as.Date("2020/03/25"),by='days'))),rep(0,T_fin-length(seq.Date(from=as.Date("2020/03/15"),to=as.Date("2020/03/25"),by='days'))))

#############################
################ MCMC ##########
#############################
########specify model
model1.string <- paste0("
  model{
     for(t in 2:(T_prime+1)){
       Km[t-1,1] <- -Km[t-1,5]-Km[t-1,9]-Km[t-1,13]-Km[t-1,17]-Km[t-1,21]-Km[t-1,25]
       Km[t-1,2] <- -Km[t-1,6]-Km[t-1,10]-Km[t-1,14]-Km[t-1,18]-Km[t-1,22]-Km[t-1,26]
       Km[t-1,3] <- -Km[t-1,7]-Km[t-1,11]-Km[t-1,15]-Km[t-1,19]-Km[t-1,23]-Km[t-1,27]
       Km[t-1,4] <- -Km[t-1,8]-Km[t-1,12]-Km[t-1,16]-Km[t-1,20]-Km[t-1,24]-Km[t-1,28]
       
       Km[t-1,5] <- beta*pi[t-1]*theta[t-1,1]*(theta[t-1,4]+a*theta[t-1,3])/N-theta[t-1,2]/De-n[t-1]*theta[t-1,2]/N
       Km[t-1,6] <- beta*pi[t-1]*(theta[t-1,1]+0.5*Km[t-1,1])*(theta[t-1,4]+0.5*Km[t-1,13]+a*(theta[t-1,3]+0.5*Km[t-1,9]))/(N)-(theta[t-1,2]+0.5*Km[t-1,5])/De-n[t-1]*(theta[t-1,2]+0.5*Km[t-1,5])/N
       Km[t-1,7] <- beta*pi[t-1]*(theta[t-1,1]+0.5*Km[t-1,2])*(theta[t-1,4]+0.5*Km[t-1,14]+a*(theta[t-1,3]+0.5*Km[t-1,10]))/(N)-(theta[t-1,2]+0.5*Km[t-1,6])/De-n[t-1]*(theta[t-1,2]+0.5*Km[t-1,6])/(N)
       Km[t-1,8] <- beta*pi[t-1]*(theta[t-1,1]+Km[t-1,3])*(theta[t-1,4]+Km[t-1,15]+a*(theta[t-1,3]+Km[t-1,11]))/(N)-(theta[t-1,2]+Km[t-1,7])/De-n[t-1]*(theta[t-1,2]+Km[t-1,7])/N
      
       Km[t-1,9] <- (1-r)*theta[t-1,2]/De-((1-rate1)/2.9+rate1/2.9)*theta[t-1,3]-n[t-1]*theta[t-1,3]/N
       Km[t-1,10] <- (1-r)*(theta[t-1,2]+0.5*Km[t-1,5])/De-((1-rate1)/2.9+rate1/2.9)*(theta[t-1,3]+0.5*Km[t-1,9])-n[t-1]*(theta[t-1,3]+0.5*Km[t-1,9])/N
       Km[t-1,11] <- (1-r)*(theta[t-1,2]+0.5*Km[t-1,6])/De-((1-rate1)/2.9+rate1/2.9)*(theta[t-1,3]+0.5*Km[t-1,10])-n[t-1]*(theta[t-1,3]+0.5*Km[t-1,10])/N
       Km[t-1,12] <- (1-r)*(theta[t-1,2]+Km[t-1,7])/De-((1-rate1)/2.9+rate1/2.9)*(theta[t-1,3]+Km[t-1,11])-n[t-1]*(theta[t-1,3]+Km[t-1,11])/N
      
       Km[t-1,13] <- r*theta[t-1,2]/De-(((1-rate1)/2.9+rate1/2.9)+1/Dh)*theta[t-1,4]
       Km[t-1,14] <- r*(theta[t-1,2]+0.5*Km[t-1,5])/De-(((1-rate1)/2.9+rate1/2.9)+1/Dh)*(theta[t-1,4]+0.5*Km[t-1,13])
       Km[t-1,15] <- r*(theta[t-1,2]+0.5*Km[t-1,6])/De-(((1-rate1)/2.9+rate1/2.9)+1/Dh)*(theta[t-1,4]+0.5*Km[t-1,14])
       Km[t-1,16] <- r*(theta[t-1,2]+Km[t-1,7])/De-(((1-rate1)/2.9+rate1/2.9)+1/Dh)*(theta[t-1,4]+Km[t-1,15])
      
       Km[t-1,17] <- 1/Dh*theta[t-1,4]-(rate2/8.6+(1-rate2)/8.6)*theta[t-1,5]
       Km[t-1,18] <- 1/Dh*(theta[t-1,4]+0.5*Km[t-1,13])-(rate2/8.6+(1-rate2)/8.6)*(theta[t-1,5]+0.5*Km[t-1,17])
       Km[t-1,19] <- 1/Dh*(theta[t-1,4]+0.5*Km[t-1,14])-(rate2/8.6+(1-rate2)/8.6)*(theta[t-1,5]+0.5*Km[t-1,18])
       Km[t-1,20] <- 1/Dh*(theta[t-1,4]+Km[t-1,15])-(rate2/8.6+(1-rate2)/8.6)*(theta[t-1,5]+Km[t-1,19])
      
       Km[t-1,21] <- (1-rate1)/2.9*(theta[t-1,4]+theta[t-1,3])+(1-rate2)/8.6*theta[t-1,5]-n[t-1]*theta[t-1,6]/N
       Km[t-1,22] <- (1-rate1)/2.9*((theta[t-1,4]+0.5*Km[t-1,13])+(theta[t-1,3]+0.5*Km[t-1,9]))+(1-rate2)/8.6*(theta[t-1,5]+0.5*Km[t-1,17])-n[t-1]*(theta[t-1,6]+0.5*Km[t-1,21])/N
       Km[t-1,23] <- (1-rate1)/2.9*((theta[t-1,4]+0.5*Km[t-1,14])+(theta[t-1,3]+0.5*Km[t-1,10]))+(1-rate2)/8.6*(theta[t-1,5]+0.5*Km[t-1,18])-n[t-1]*(theta[t-1,6]+0.5*Km[t-1,22])/N
       Km[t-1,24] <- (1-rate1)/2.9*((theta[t-1,4]+Km[t-1,15])+(theta[t-1,3]+Km[t-1,11]))+(1-rate2)/8.6*(theta[t-1,5]+Km[t-1,19])-n[t-1]*(theta[t-1,6]+Km[t-1,23])/N
      
       Km[t-1,25] <- rate1/2.9*(theta[t-1,3]+theta[t-1,4])+rate2/8.6*theta[t-1,5]
       Km[t-1,26] <- rate1/2.9*((theta[t-1,3]+0.5*Km[t-1,9])+(theta[t-1,4]+0.5*Km[t-1,13]))+rate2/8.6*(theta[t-1,5]+0.5*Km[t-1,17])
       Km[t-1,27] <- rate1/2.9*((theta[t-1,3]+0.5*Km[t-1,10])+(theta[t-1,4]+0.5*Km[t-1,14]))+rate2/8.6*(theta[t-1,5]+0.5*Km[t-1,18])
       Km[t-1,28] <- rate1/2.9*((theta[t-1,3]+Km[t-1,11])+(theta[t-1,4]+Km[t-1,15]))+rate2/8.6*(theta[t-1,5]+Km[t-1,19])
      
       alpha[t-1,1] <- theta[t-1,1]+(Km[t-1,1]+2*Km[t-1,2]+2*Km[t-1,3]+Km[t-1,4])/6
       alpha[t-1,2] <- theta[t-1,2]+(Km[t-1,5]+2*Km[t-1,6]+2*Km[t-1,7]+Km[t-1,8])/6
       alpha[t-1,3] <- theta[t-1,3]+(Km[t-1,9]+2*Km[t-1,10]+2*Km[t-1,11]+Km[t-1,12])/6
       alpha[t-1,4] <- theta[t-1,4]+(Km[t-1,13]+2*Km[t-1,14]+2*Km[t-1,15]+Km[t-1,16])/6
       alpha[t-1,5] <- theta[t-1,5]+(Km[t-1,17]+2*Km[t-1,18]+2*Km[t-1,19]+Km[t-1,20])/6
       alpha[t-1,6] <- theta[t-1,6]+(Km[t-1,21]+2*Km[t-1,22]+2*Km[t-1,23]+Km[t-1,24])/6
       alpha[t-1,7] <- theta[t-1,7]+(Km[t-1,25]+2*Km[t-1,26]+2*Km[t-1,27]+Km[t-1,28])/6
       
      for(j in c(3,2,4,5,6,7)){ 
      theta[t,j] ~ dpois(alpha[t-1,j])
      }
      
      theta[t,1]= N-theta[t,7]-theta[t,2]-theta[t,3]-theta[t,4]-theta[t,5]-theta[t,6]
      
       Y[t-1] ~ dpois(r*theta[t-1,2]/De)
       R[t-1] ~ dpois((1-rate1)/2.9*theta[t-1,4]+(1-rate2)/8.6*theta[t-1,5])
       D[t-1] ~ dpois(rate1/2.9*theta[t-1,4]+rate2/8.6*theta[t-1,5])
     }
    theta0[1:7]<-c(",yinit[1],",",yinit[2],",", yinit[6],",", yinit[3],",",yinit[5],",", yinit[4],",", yinit[7],")
    for (i in c(3,2,4,5,6,7)){
      theta[1,i] ~ dpois(theta0[i])
    }
    theta[1,1]= N-theta0[3]-theta0[2]-theta0[7]-theta0[4]-theta0[5]-theta0[6]
    r ~ dbeta(10,90)
    R0 ~ dlnorm(", lognorm_R0_parm$mu, ",", 1 / lognorm_R0_parm$var,")
    rate1 ~ dbeta(0.0296,2.9304)
    rate2 ~ dbeta(0.44,1.76)
    beta <- R0/(a*(1-r)/De/((1/De)*((1-rate1)/2.9+rate1/2.9))+r/De/((1/De)*(((1-rate1)/2.9+rate1/2.9)+1/Dh)))
  }
")

model.spec <- textConnection(model1.string)

posterior1 <- jags.model(
  model.spec,
  data = list(
    'Y' = Y,
    'R' = R,
    'D' = D,
    'T_prime' = T_prime,
    'pi' = pi,
    'N'=N_SA,
    'n'=n,
    'a' = 0.55,
    'De' = 5.2,
    'Dh'=7),
  n.chains = nchain,
  n.adapt = nadapt
)
#burn-in
update(posterior1, nburnin) 
#sampling
jags_sample <- jags.samples(
  posterior1,
  c('theta','r','R0','rate1','rate2','beta','Y','R','D'),
  n.iter = M,
  thin = thn
)

###################################
#############organize the posterior of parameters
#############################
R0_p <- unlist(as.mcmc.list(jags_sample$R0))
r_p <- unlist(as.mcmc.list(jags_sample$r))
beta_p <- unlist(as.mcmc.list(jags_sample$beta))
rate1_p <- unlist(as.mcmc.list(jags_sample$rate1))
rate2_p <- unlist(as.mcmc.list(jags_sample$rate2))
theta_p <- array(Reduce(rbind,as.mcmc.list(jags_sample$theta)),dim=c(len,T_prime+1,7))

############posterior means and 95% credible intervals
theta_p_mean <- apply(theta_p[,T_prime+1,],2,mean)
theta_p_ci <- as.vector(apply(theta_p[,T_prime+1,],2,quantile,c(0.025,0.5,0.975)))
R0_p_mean <- mean(R0_p)
R0_p_ci <- quantile(R0_p,c(0.025,0.5,0.975))
r_p_mean <- mean(r_p)
r_p_ci <- quantile(r_p,c(0.025,0.5,0.975))
beta_p_mean <- mean(beta_p)
beta_p_ci <- quantile(beta_p,c(0.025,0.5,0.975))
rate1_p_mean <- mean(rate1_p)
rate1_p_ci <- quantile(rate1_p,c(0.025,0.5,0.975))
rate2_p_mean <- mean(rate2_p)
rate2_p_ci <- quantile(rate2_p,c(0.025,0.5,0.975))

###########################
#### Forecast ####
#############################
T_fin=200
theta_pp <- array(0,dim=c(len,T_fin,7))
Y_pp <- matrix(NA,nrow=len,ncol=T_fin)
R_pp <- matrix(NA,nrow=len,ncol=T_fin)
D_pp <- matrix(NA,nrow=len,ncol=T_fin)
N=N_SA
a = 0.55
De = 5.2
Dh=7

for(l in 1:len){
  thetalt1 <- theta_p[l,1,1]
  thetalt2 <- theta_p[l,1,2]
  thetalt3 <- theta_p[l,1,3]
  thetalt4 <- theta_p[l,1,4]
  thetalt5 <- theta_p[l,1,5]
  thetalt6 <- theta_p[l,1,6]
  thetalt7 <- theta_p[l,1,7]
  betal <- c(beta_p)[l]
  rl <- c(r_p)[l]
  rate1l <- c(rate1_p)[l]
  rate2l <- c(rate2_p)[l]
  if(betal<0 |rl<0 |rate1l<0 |rate2l<0 |thetalt1<0 |thetalt2<0 |thetalt3<0|thetalt4<0|thetalt5<0|thetalt6<0|thetalt7<0) next
  for(t in 1:(T_fin)){
    Km <- NULL
    alpha_pp <- NULL
    
    Km[5] <- betal*pi[t]*thetalt1*(thetalt4+a*thetalt3)/N-thetalt2/De-n[t]*thetalt2/N
    Km[13] <- rl*thetalt2/De-(((1-rate1l)/2.9+(rate1l)/2.9)+1/Dh)*thetalt4
    Km[17] <- 1/Dh*thetalt4-(rate2l/8.6+(1-rate2l)/8.6)*thetalt5
    Km[21] <- (1-rate1l)/2.9*(thetalt4+thetalt3)+(1-rate2l)/8.6*thetalt5-n[t]*thetalt6/N
    Km[25] <- (rate1l)/2.9*(thetalt3+thetalt4)+rate2l/8.6*thetalt5
    Km[9] <- (1-rl)*thetalt2/De-((1-rate1l)/2.9+(rate1l)/2.9)*thetalt3-n[t]*thetalt3/N
    Km[1] <- -Km[9]-Km[5]-Km[13]-Km[17]-Km[21]-Km[25]
    
    Km[10] <- (1-rl)*(thetalt2+0.5*Km[5])/De-((1-rate1l)/2.9+(rate1l)/2.9)*(thetalt3+0.5*Km[9])-n[t]*(thetalt3+0.5*Km[9])/N
    Km[6] <- betal*pi[t]*(thetalt1+0.5*Km[1])*(thetalt4+0.5*Km[13]+a*(thetalt3+0.5*Km[9]))/N-(thetalt2+0.5*Km[5])/De-n[t]*(thetalt2+0.5*Km[5])/N
    Km[14] <- rl*(thetalt2+0.5*Km[5])/De-(((1-rate1l)/2.9+(rate1l)/2.9)+1/Dh)*(thetalt4+0.5*Km[13])
    Km[18] <- 1/Dh*(thetalt4+0.5*Km[13])-(rate2l/8.6+(1-rate2l)/8.6)*(thetalt5+0.5*Km[17])
    Km[22] <- (1-rate1l)/2.9*((thetalt4+0.5*Km[13])+(thetalt3+0.5*Km[9]))+(1-rate2l)/8.6*(thetalt5+0.5*Km[17])-n[t]*(thetalt6+0.5*Km[21])/N
    Km[26] <- (rate1l)/2.9*((thetalt3+0.5*Km[9])+(thetalt4+0.5*Km[13]))+rate2l/8.6*(thetalt5+0.5*Km[17])
    Km[2] <- -Km[10]-Km[6]-Km[14]-Km[18]-Km[22]-Km[26]
    
    Km[11] <- (1-rl)*(thetalt2+0.5*Km[6])/De-((1-rate1l)/2.9+(rate1l)/2.9)*(thetalt3+0.5*Km[10])-n[t]*(thetalt3+0.5*Km[10])/N
    Km[7] <- betal*pi[t]*(thetalt1+0.5*Km[2])*(thetalt4+0.5*Km[14]+a*(thetalt3+0.5*Km[10]))/N-(thetalt2+0.5*Km[6])/De-n[t]*(thetalt2+0.5*Km[6])/N
    Km[15] <- rl*(thetalt2+0.5*Km[6])/De-(((1-rate1l)/2.9+(rate1l)/2.9)+1/Dh)*(thetalt4+0.5*Km[14])
    Km[19] <- 1/Dh*(thetalt4+0.5*Km[14])-(rate2l/8.6+(1-rate2l)/8.6)*(thetalt5+0.5*Km[18])
    Km[23] <- (1-rate1l)/2.9*((thetalt4+0.5*Km[14])+(thetalt3+0.5*Km[10]))+(1-rate2l)/8.6*(thetalt5+0.5*Km[18])-n[t]*(thetalt6+0.5*Km[22])/N
    Km[27] <- (rate1l)/2.9*((thetalt3+0.5*Km[10])+(thetalt4+0.5*Km[14]))+rate2l/8.6*(thetalt5+0.5*Km[18])
    Km[3] <- -Km[11]-Km[7]-Km[15]-Km[19]-Km[23]-Km[27]
    
    Km[12] <- (1-rl)*(thetalt2+Km[7])/De-((1-rate1l)/2.9+(rate1l)/2.9)*(thetalt3+Km[11])-n[t]*(thetalt3+Km[11])/N
    Km[8] <- betal*pi[t]*(thetalt1+Km[3])*(thetalt4+Km[15]+a*(thetalt3+Km[11]))/N-(thetalt2+Km[7])/De-n[t]*(thetalt2+Km[7])/N
    Km[16] <- rl*(thetalt2+Km[7])/De-(((1-rate1l)/2.9+(rate1l)/2.9)+1/Dh)*(thetalt4+Km[15])
    Km[20] <- 1/Dh*(thetalt4+Km[15])-(rate2l/8.6+(1-rate2l)/8.6)*(thetalt5+Km[19])
    Km[24] <- (1-rate1l)/2.9*((thetalt4+Km[15])+(thetalt3+Km[11]))+(1-rate2l)/8.6*(thetalt5+Km[19])-n[t]*(thetalt6+Km[23])/N
    Km[28] <- (rate1l)/2.9*((thetalt3+Km[11])+(thetalt4+Km[15]))+rate2l/8.6*(thetalt5+Km[19])
    Km[4] <- -Km[12]-Km[8]-Km[16]-Km[20]-Km[24]-Km[28]
    
    alpha_pp[1] <- thetalt1+(Km[1]+2*Km[2]+2*Km[3]+Km[4])/6
    alpha_pp[2] <- thetalt2+(Km[5]+2*Km[6]+2*Km[7]+Km[8])/6
    alpha_pp[3] <- thetalt3+(Km[9]+2*Km[10]+2*Km[11]+Km[12])/6
    alpha_pp[4] <- thetalt4+(Km[13]+2*Km[14]+2*Km[15]+Km[16])/6
    alpha_pp[5] <- thetalt5+(Km[17]+2*Km[18]+2*Km[19]+Km[20])/6
    alpha_pp[6] <- thetalt6+(Km[21]+2*Km[22]+2*Km[23]+Km[24])/6
    alpha_pp[7] <- thetalt7+(Km[25]+2*Km[26]+2*Km[27]+Km[28])/6
    thetalt_tmp=vector(length = 7)
    for (j in c(3,2,4,5,6,7)){
      thetalt_tmp[j] <- rpois(1,alpha_pp[j])
    }
    thetalt_tmp[1] <- N_SA-thetalt_tmp[3]-thetalt_tmp[2]-thetalt_tmp[7]-thetalt_tmp[4]-thetalt_tmp[5]-thetalt_tmp[6]
    thetalt1<-theta_pp[l,t,1] <- thetalt_tmp[1]
    thetalt2<-theta_pp[l,t,2] <- thetalt_tmp[2]
    thetalt3<-theta_pp[l,t,3] <- thetalt_tmp[3]
    thetalt4<-theta_pp[l,t,4] <- thetalt_tmp[4]
    thetalt5<-theta_pp[l,t,5] <- thetalt_tmp[5]
    thetalt6<-theta_pp[l,t,6] <- thetalt_tmp[6]
    thetalt7<-theta_pp[l,t,7] <- thetalt_tmp[7]
    
    Y_pp[l,t] <- rpois(1,rl*thetalt2/De)
    R_pp[l,t] <- rpois(1,(1-rate1l)/2.9*thetalt4+(1-rate2l)/8.6*thetalt5)
    D_pp[l,t] <- rpois(1,(rate1l)/2.9*thetalt4+rate2l/8.6*thetalt5)
  }
}
#############################
#### save image ####
#############################
save.image('eSEIRD.RData')

