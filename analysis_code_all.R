# Analysis Code - Data Analysis in 'Leveraging Multiple Statistical Methods for Inverse Prediction in Nuclear Forensics Applications' (LMSMIPNFA) v. 1.0

# Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.

# Code to reproduce the data analysis in Leveraging Multiple Statistical Methods for Inverse Prediction in Nuclear Forensics Applications by Lewis, Zhang, Anderson-Cook (2017)

sessionInfo()

# 
# Auxilary Functions to make predictions under each method ----
#

#Forward Modeling Classical Linear Regression----


classical_forward<-function(df_train, df_test,weights,n_inits=3, sd=0){
  
  prior_range<-rbind(rgHNO3,rgPu,rgTemp)
  #rbind(rep(max(abs(rgHNO3)),2)*c(-1,1),rep(max(abs(rgPu)),2)*c(-1,1),rep(max(abs(rgTemp)),2)*c(-1,1))
  
  full_fits<-apply(df_train[,y_names], 2, function(y) lm(y~HNO3+Pu+Temp+HNO3*Pu+HNO3*Temp+Pu*Temp+I(HNO3^2)+I(Pu^2)+I(Temp^2), data=df_train))
  weights_full<-sapply(full_fits, FUN=function(f) summary(f)$sigma)
  
  
  step_fits<-lapply(full_fits, function(l) step(l, trace=0))
  weights_step<-sapply(step_fits, FUN=function(f) summary(f)$sigma)
  ypred <- function(x, model_list){
    x <- data.frame(t(x))
    names(x)<-x_names
    numResp<-length(model_list)
    yp<-numeric(numResp)
    for(i in 1:numResp){
      yp[i] <- predict(model_list[[i]], newdata=x)
    }
    yp
  }
  
  func_ss <- function(x, ystar, model_list, weights=NULL){
    yp <- ypred(x, model_list)
    if(is.null(weights)){weights<-1}
    sum(((ystar-yp)/weights)^2)
  }
  
  par_init<-sapply(1:n_inits, function(a) df_test[,x_names]+ rnorm(3, sd=sd))
  par_init<-apply(par_init, 2, function(x){
    if(any(x < prior_range[,1] | any(x > prior_range[,2]))){
      as.numeric(df_test[,x_names])
    } else {
      as.numeric(x)
    }
  })
  
  x_pred_full <- apply(par_init, 2, function(init) optim(par=init, fn=func_ss, method="Nelder-Mead", ystar=df_test[,y_names], model_list=full_fits,weights=weights))#, lower=c(rgHNO3[1],rgPu[1], rgTemp[1]), upper=c(rgHNO3[2],rgPu[2], rgTemp[2])))
  vals<-sapply(x_pred_full , function(x) x$value)
  # plot(vals)
  # abline(h=vals[which(vals==min(vals))[1]])
  x_pred_full<-x_pred_full[[which(vals==min(vals))[1]]]
  
  x_pred_step <- apply(par_init, 2, function(init) optim(par=init, fn=func_ss, method="Nelder-Mead", ystar=df_test[,y_names], model_list=step_fits,weights=weights))#,lower=c(rgHNO3[1],rgPu[1], rgTemp[1]), upper=c(rgHNO3[2],rgPu[2], rgTemp[2])))
  vals2<-sapply(x_pred_step , function(x) x$value)
  x_pred_step<-x_pred_step[[which(vals2==min(vals2))[1]]]
  
  out<-list()
  out$x_pred_full<-x_pred_full$par
  out$full_conv<-x_pred_full$convergence
  out$x_pred_step<-x_pred_step$par
  out$step_conv<-x_pred_step$convergence
  out$x_star<-df_test[,x_names] #true x
  out
}


#uses constrOptim instead of optim
classical_forward2<-function(df_train, df_test, weights,n_inits=3, sd=0){
  
  prior_range<-rbind(rgHNO3,rgPu,rgTemp)
  
  full_fits<-apply(df_train[,y_names], 2, function(y) lm(y~HNO3+Pu+Temp+HNO3*Pu+HNO3*Temp+Pu*Temp+I(HNO3^2)+I(Pu^2)+I(Temp^2), data=df_train))
  
  step_fits<-lapply(full_fits, function(l) step(l, trace=0, direction='both'))
  
  ypred <- function(x, model_list){
    x <- data.frame(t(x))
    names(x)<-x_names
    numResp<-length(model_list)
    yp<-numeric(numResp)
    for(i in 1:numResp){
      yp[i] <- predict(model_list[[i]], newdata=x)
    }
    yp
  }
  
  func_ss <- function(x, ystar, model_list, weights=NULL){
    yp <- ypred(x, model_list)
    if(is.null(weights)){weights<-1}
    sum(((ystar-yp)/weights)^2)
  }
  
  par_init<-sapply(1:n_inits, function(a) df_test[,x_names]+ rnorm(3, sd=sd))
  par_init<-apply(par_init, 2, function(x){
    if(any(x < prior_range[,1] | any(x > prior_range[,2]))){
      as.numeric(df_test[,x_names])
    } else {
      as.numeric(x)
    }
  })
  
  ui<-matrix(c(1,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,1,-1), nrow=6, ncol=3)    
  ci<-c(rgHNO3[1],-rgHNO3[2],rgPu[1],-rgPu[2],rgTemp[1],-rgTemp[2])
  
  x_pred_full <- apply(par_init, 2, function(init) constrOptim(theta=init, f=func_ss, ystar=df_test[,y_names], model_list=full_fits,weights=weights,ui=ui, ci=ci, grad=NULL))
  vals<-sapply(x_pred_full , function(x) x$value)
  # plot(vals)
  # abline(h=vals[which(vals==min(vals))[1]])
  
  x_pred_full<-x_pred_full[[which(vals==min(vals))[1]]]
  
  x_pred_step <- apply(par_init, 2, function(init) constrOptim(theta=init, f=func_ss, ystar=df_test[,y_names], model_list=step_fits,weights=weights,ui=ui, ci=ci, grad=NULL))
  vals<-sapply(x_pred_step , function(x) x$value)
  x_pred_step<-x_pred_step[[which(vals==min(vals))[1]]]
  
  out<-list()
  out$x_pred_full<-x_pred_full$par
  out$full_conv<-x_pred_full$convergence
  out$x_pred_step<-x_pred_step$par
  out$step_conv<-x_pred_step$convergence
  out$x_star<-df_test[,x_names] #true x
  out
}





# Forward Modeling, Bayesian Linear modeling using the MAP estimator----
#note - we use plug in mle's for beta's and sigmas and only solve the optimizization for x

bayesian_map <- function(df_train, df_test, n_inits=3,sd=0) {
  
  full_fits<-apply(df_train[,y_names], 2, function(y) lm(y~HNO3+Pu+Temp+HNO3*Pu+HNO3*Temp+Pu*Temp+I(HNO3^2)+I(Pu^2)+I(Temp^2), data=df_train))
  #AIC
  step_fits<-lapply(full_fits, function(l) step(l, trace=0))
  
  prior_range<-rbind(rgHNO3,rgPu,rgTemp)
  ystar<-df_test[,y_names]
  fn_log_likelihood<-function(model_list, ystar,xstar){
    xstar<-data.frame(t(xstar))
    names(xstar)<-x_names
    # log_lik<-lapply(model_list, function(l) {
    #    mu_hat<-predict(l)
    #    sigma_hat<-summary(l)$sigma
    #    y_obs<-l$residuals+l$fitted.values
    #    sum(dnorm(y_obs, mu_hat, sd=sigma_hat, log=TRUE))
    #    })
    # log_lik<-sum(unlist(log_lik)) #part not needed since all values fixed here
    
    log_lik_str<-mapply(FUN=function(l,y) {
      mu_hat<-predict(l, newdata=xstar)
      sigma_hat<-summary(l)$sigma
      sum(dnorm(y, mu_hat, sd=sigma_hat, log=TRUE))
    }, model_list, ystar)
    log_lik_str<-sum(unlist(log_lik_str))
    log_lik_str #+log_lik
  }
  
  log_prior_xstar<-function(xstar){
    if(any(xstar < prior_range[,1] | any(xstar > prior_range[,2]))){
      -Inf
    } else {
      0
    }
  }
  
  neg_log_posterior<-function(xstar,model_list, ystar){
    -fn_log_likelihood(model_list, ystar,xstar)-sum(log_prior_xstar(xstar))
  }  
  
  par_init<-sapply(1:n_inits, function(a) df_test[,x_names]+ rnorm(3, sd=sd))
  par_init<-apply(par_init, 2, function(x){
    if(any(x < prior_range[,1] | any(x > prior_range[,2]))){
      as.numeric(df_test[,x_names])
    } else {
      as.numeric(x)
    }
  })
  
  #par_init<-c(0,0,0)
  map_full<-apply(par_init, 2, function(init) optim(init, fn=neg_log_posterior,model_list=full_fits, ystar=ystar, control=list(maxit=500), method='Nelder-Mead'))
  vals<-sapply(map_full, function(x) x$value)
  map_full<-map_full[[which(vals==min(vals))[1]]]
  
  map_step<-apply(par_init, 2, function(init) optim(init, fn=neg_log_posterior,model_list=step_fits, ystar=ystar, control=list(maxit=500), method='Nelder-Mead'))
  vals<-sapply(map_step, function(x) x$value)
  map_step<-map_step[[which(vals==min(vals))[1]]]
  
  out<-list()
  out$map_full<-map_full$par
  out$map_full_conv<-map_full$convergence
  out$map_step<-map_step$par
  out$map_step_conv<-map_step$convergence
  out$x_star<-df_test[,x_names] #true x
  out
}  

#mean((map_full[[12]]$par-df_test[,x_names])^2)

# Inverse Modeling, lm,  Principal Components Regression (PCR) ----

lm_inverse<-function(df_train,df_test){
  df_train<-df_train[complete.cases(df_train),]
  resp<-as.matrix(df_train[,x_names])
  fit<-apply(resp,2, function(y) lm(y~Mode+Aggwt+Length+Width+Thickness, #+
                                    # I(Mode^2)+
                                    # I(Aggwt^2)+
                                    # I(Length^2)+
                                    # I(Width^2)+
                                    # I(Thickness^2),
                                    data=df_train))
  x_pred<-sapply(fit, function(f) as.numeric(predict(f, newdata=df_test)))
  x_pred_step<-sapply(fit, function(f) {
    f2<-step(f,trace = FALSE)
    as.numeric(predict(f2, newdata=df_test))
  })
  out<-list()
  out$x_pred_lm<-x_pred
  out$x_pred_lm_step<-x_pred_step
  out$x_star<-df_test[,x_names]
  out
}





library(pls)
pcr_inverse<-function(df_train,df_test, propVar=.95, scale=TRUE){
  resp<-as.matrix(df_train[,x_names])
  fit<-pcr(resp~Mode+Aggwt+Length+Width+Thickness,#^2+
           # I(Mode^2)+
           # I(Aggwt^2)+
           # I(Length^2)+
           # I(Width^2)+
           # I(Thickness^2),
           data=df_train, scale=scale)
  ncomp<-min(which(cumsum(fit$Xvar)/fit$Xtotvar>=propVar))
  x_pred<-predict(fit, newdata=df_test,ncomp=ncomp)
  out<-list()
  out$x_pred_pcr<-x_pred
  out$x_star<-df_test[,x_names]
  out
}



# Inverse Modeling, Partial Least Squares Regression (PLSR) ----
plsr_inverse<-function(df_train,df_test, propVar=.95, scale=TRUE){
  resp<-as.matrix(df_train[,x_names])
  fit<-plsr(resp~Mode+Aggwt+Length+Width+Thickness,#^2+
            # I(Mode^2)+
            # I(Aggwt^2)+
            # I(Length^2)+
            # I(Width^2)+
            # I(Thickness^2),
            data=df_train, method='oscorespls', scale=TRUE)
  ncomp<-min(which(cumsum(fit$Xvar)/fit$Xtotvar>=propVar))
  x_pred<-predict(fit, newdata=df_test,ncomp=ncomp)
  out<-list()
  out$x_pred_plsr<-x_pred
  out$x_star<-df_test[,x_names]
  out
}

# end Auxilary functions ------







# 
# Analysis of PuO2 data from Burney's 1980's report ----
#


library(Hmisc)
library(effects)
library(stargazer)
library(MASS)
library(rstan)
library(ggplot2)
library(GGally)
library(bayesplot)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#source('auxilary_functions.R')



simulatedData<-0 #use simulated data or no?


#Processing Variables ----
#HNO_3, M
#Pu, g/L
#Temperature, C

# Observed Responses
# mode, um
# Agg. Index, wt%
# Length, um
# Width, um
# Thickness, um

#scale_data<-1 # center and scale the data - keep as 1

par("pch"=19)
x_names<-c("HNO3" ,"Pu" ,"Temp" )
y_names<-c("Mode", "Aggwt" , "Length", "Width","Thickness")

pu<-read.csv('PuO2.csv')
names_pu<-names(pu)
#View(pu)
# apply(pu, 2, sd, na.rm=TRUE)

names(pu)<-c('isoSample', 'Mode \n ($\\mu$m)', 'Agglomeration \n (\\% weight)', 'Length \n  ($\\mu$m)', 'Width \n ($\\mu$m)', 'Thickness \n ($\\mu$m)', 'HNO$_3$ \n (M)', "Pu \n (g/L)", 'Temperature \n (C)')


# l_table <- latex(pu[,-1], caption="Complete data set from \\cite{burney1984} of experimental data characterizing the effects of three processing variables (column names HN03,  Pu, Temp)  on particle size and morphology parameters (Mode, Aggwt, Length, Width, Thickness) of produced plutonium oxide. Several types of missing/censored data is prevalent, indicated by either missing cells or letter codes. We do not consider the missing data mechanism in this analysis and treat all missing/censored values as unobserved.", label="burneyData", 
#                  rowlabel='Sample \n \\#',where='!htbp')
# 
# latexVerbatim(l_table)
# print(dvi(l_table, width = 5, height = 6.5), , width = 5, height = 6.5)

pu<-t(apply(pu, 1, function(x) as.numeric(as.character(x))))
# colLabs<-c('Mode', 'Agg. Index', 'Length', 'Width', 'Thickness', 'HNO', 'Pu', 'Temp.')

# 
# #pdf('figs/pairsPlot.pdf')
# ggpairs(as.data.frame(pu[-10,-1]), diag=list(continuous="barDiag"), lower = list(continuous = "cor"), upper = list(continuous = "points"), columnLabels = colLabs)+theme_bw()
# #dev.off()

colnames(pu)<-names_pu
pu<-as.data.frame(pu)

pu<-pu[,-1]
pu<-pu[-which(pu$Length==205),] #length 205??
pu<-pu[-which(pu$Temp==70),]
pu<-pu[-which(pu$Mode==12),]
pu[complete.cases(pu),]

#scale all variable to [-1,1] 
pu_scale<-apply(pu,2, function(x){
  middle<-(max(x, na.rm = TRUE)+min(x, na.rm = TRUE))/2
  dif<-(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))/2
  (x-middle)/dif
})



# if(scale_data){
pu<-as.data.frame(pu_scale)
weights<-1
# }
apply(pu, 2,range, na.rm=TRUE)
comp_cases<-complete.cases(pu)
pu_test<-pu[comp_cases,]
dim(pu_test)
pu_train<-pu[!comp_cases,]
dim(pu_train)
test_cases<-c(1:nrow(pu_test))  
# test_cases <- 1:2 #for testing
# pairs(pu, cex=.5)
# pairs(pu[,x_names], cex=.5)
# pairs(pu[,y_names], cex=.5)

set.seed(124)
# Using Simulated data -----
if(simulatedData){
  
  sd<-1
  pu$Mode<-pu$HNO3-2*pu$HNO3^2+pu$Pu+pu$Temp+1.5*pu$Pu*pu$Temp+rnorm(nrow(pu),sd=sd)
  pu$Aggwt<-3*pu$HNO3-pu$Pu^2+2*pu$Temp^2+rnorm(nrow(pu),sd=sd)
  pu$Length<-pu$HNO3+pu$Pu-2*pu$Temp+rnorm(nrow(pu),sd=sd)
  pu$Width<-pu$HNO3-pu$Pu+2*pu$HNO3*pu$Pu+rnorm(nrow(pu),sd=sd)
  pu$Thickness<-pu$HNO3-2*pu$HNO3^2-pu$Pu-pu$Temp+rnorm(nrow(pu),sd=sd)
  
  # pu$Mode<-pu$HNO3-2*pu$HNO3+pu$Pu+pu$Temp+rnorm(nrow(pu),sd=sd)
  # pu$Aggwt<-3*pu$HNO3-pu$Pu+2*pu$Temp+rnorm(nrow(pu),sd=sd)
  # pu$Length<-pu$HNO3+pu$Pu-2*pu$Temp+rnorm(nrow(pu),sd=sd)
  # pu$Width<-pu$HNO3-pu$Pu+2*pu$HNO3*pu$Pu+rnorm(nrow(pu),sd=sd)
  # pu$Thickness<-pu$HNO3-2*pu$HNO3-pu$Pu-pu$Temp+rnorm(nrow(pu),sd=sd)
  
  pairs(pu)
  comp_cases<-complete.cases(pu)
  pu_test<-pu[comp_cases,]
  dim(pu_test)
  pu_train<-pu[!comp_cases,]
  dim(pu_train)
  test_cases<-1:10
}



#Forward Modeling Classical Linear Regression----
ranges<-as.data.frame(apply(pu,2, range, na.rm=TRUE))
rgHNO3<-ranges$HNO3+c(-.5,.5)
rgPu<-ranges$Pu+c(-.5,.5)
rgTemp<-ranges$Temp+c(-.5,.5)

apply(pu,2, sd, na.rm=TRUE)

if(simulatedData){
  n_inits<-10 #number of starting points for optimization routine
  sd_init<-.1 #pertubation of x* for a starting point
} else{
  n_inits<-10
  sd_init<-.1
}

cl_foward_preds<-list()
j<-1 #make work when test_cases isn't sequential
for(i in test_cases){
  df_train<-rbind(pu_train, pu_test[-i,])
  df_test<-pu_test[i,]
  cl_foward_preds[[j]]<-classical_forward(df_train, df_test,weights = weights,n_inits = n_inits, sd=sd_init)
  j<-j+1
  print(i)
}


conv<-t(sapply(cl_foward_preds, function(y) c(y$full_conv,y$step_conv)))
any(conv)
full_preds<-t(sapply(cl_foward_preds, function(y) y$x_pred_full))
full_preds<-as.data.frame(full_preds)
names(full_preds)<-x_names
step_preds<-t(sapply(cl_foward_preds, function(y) y$x_pred_step))
step_preds<-as.data.frame(step_preds)
names(step_preds)<-x_names
x_star<-pu_test[test_cases,x_names]

# par(mfrow=c(1,3))
# for(r in 1:3){
# plot(full_preds[,r], step_preds[,r], main=r, cex=.1)
# text(full_preds[,r], step_preds[,r], test_cases)
# abline(0,1)
# }
# 
# r<-3
# sqrt(mean((full_preds[,r]-x_star[,r])^2))
# sqrt(mean((step_preds[,r]-x_star[,r])^2))
# 

# # Bayesian MAP -----
# if(i %in% c(12,13,15,16)){
#   n_inits<-10
#   sd_init<-.5
# } else{
#   n_inits<-1
#   sd_init<-0
# }
map_foward_preds<-list()
j<-1
for(i in test_cases){
  df_train<-rbind(pu_train, pu_test[-i,])
  df_test<-pu_test[i,]
  map_foward_preds[[j]]<-bayesian_map(df_train, df_test, n_inits = n_inits, sd=sd_init)
  j<-j+1
  print(i)
}
names(map_foward_preds[[1]])
conv<-t(sapply(map_foward_preds, function(y) c(y$map_full_conv,y$map_step_conv)))
any(conv)

map_full_preds<-t(sapply(map_foward_preds, function(y) y$map_full))
map_full_preds<-as.data.frame(map_full_preds)
names(map_full_preds)<-x_names
map_step_preds<-t(sapply(map_foward_preds, function(y) y$map_step))
map_step_preds<-as.data.frame(map_step_preds)
names(map_step_preds)<-x_names


# lm inverse ----

inverse_fits<-list()
j<-1
for(i in test_cases){
  df_train<-rbind(pu_train, pu_test[-i,])
  df_test<-pu_test[i,]
  inverse_fits[[j]]<-lm_inverse(df_train,df_test)
  j<-j+1
  print(i)
}

inverse_preds<-t(sapply(inverse_fits, function(y) y$x_pred_lm))
inverse_preds<-as.data.frame(inverse_preds)
names(inverse_preds)<-x_names

inverse_preds_step<-t(sapply(inverse_fits, function(y) y$x_pred_lm_step))
inverse_preds_step<-as.data.frame(inverse_preds_step)
names(inverse_preds_step)<-x_names

# PCR ----

propVar<-0.95
pcr_fits<-list()
j<-1
for(i in test_cases){
  df_train<-rbind(pu_train, pu_test[-i,])
  df_test<-pu_test[i,]
  pcr_fits[[j]]<-pcr_inverse(df_train,df_test, propVar=propVar)
  j<-j+1
  print(i)
}



pcr_preds<-t(sapply(pcr_fits, function(y) y$x_pred_pcr))
pcr_preds<-as.data.frame(pcr_preds)
names(pcr_preds)<-x_names
# x_star<-pu_test[,x_names]


# PLSR ----

plsr_fits<-list()
j<-1
for(i in test_cases){
  df_train<-rbind(pu_train, pu_test[-i,])
  df_test<-pu_test[i,]
  plsr_fits[[j]]<-plsr_inverse(df_train,df_test, propVar=propVar)
  j<-j+1
  print(i)
}


plsr_preds<-t(sapply(plsr_fits, function(y) y$x_pred_plsr))
plsr_preds<-as.data.frame(plsr_preds)
names(plsr_preds)<-x_names

apply(pu,2, function(x) sum(!is.na(x)))


#save.image(paste0('BurneyDataAnalysis_simData', simulatedData, '.RData'))



#
#Forward Full Bayes Modeling----
#


#rm(list=ls())
library(Hmisc)
library(effects)
library(stargazer)
library(MASS)
library(rstan)
library(ggplot2)
library(GGally)
library(bayesplot)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# Stan Model ------
burney.stan <- '
data {
int <lower=0> N_Mode;
int <lower=0> N_Aggwt;
int <lower=0> N_Length;
int <lower=0> N_Width;
int <lower=0> N_Thickness;
vector[N_Mode] y_Mode;
matrix[N_Mode, 3] X_Mode;
vector[N_Aggwt] y_Aggwt;
matrix[N_Aggwt, 3] X_Aggwt;
vector[N_Length] y_Length;
matrix[N_Length, 3] X_Length;
vector[N_Width] y_Width;
matrix[N_Width, 3] X_Width;
vector[N_Thickness] y_Thickness;
matrix[N_Thickness, 3] X_Thickness;
vector[5] y_star; // Mode, Aggwt, Length, Width, Thickness
vector[2] rgHNO3; // prior limits on xstar
vector[2] rgPu;
vector[2] rgTemp;
}


parameters  {
vector[3] x_star;
vector[10] beta_Mode;
vector[10] beta_Aggwt;
vector[10] beta_Length;
vector[10] beta_Width;
vector[10] beta_Thickness;
real <lower=0> sigma_Mode; 
real <lower=0> sigma_Aggwt; 
real <lower=0> sigma_Length; 
real <lower=0> sigma_Width; 
real <lower=0> sigma_Thickness; 
}

transformed parameters{
real mu_Mode[N_Mode];
real mu_Aggwt[N_Aggwt];
real mu_Length[N_Length];
real mu_Width[N_Width];
real mu_Thickness[N_Thickness];
real mu_Mode_star;
real mu_Aggwt_star;
real mu_Length_star;
real mu_Width_star;
real mu_Thickness_star;


for(j in 1:N_Mode){
mu_Mode[j]=beta_Mode[1]+
beta_Mode[2]*X_Mode[j,1]+
beta_Mode[3]*X_Mode[j,2]+
beta_Mode[4]*X_Mode[j,3]+
beta_Mode[5]*X_Mode[j,1]^2+
beta_Mode[6]*X_Mode[j,2]^2+
beta_Mode[7]*X_Mode[j,3]^2+
beta_Mode[8]*X_Mode[j,1]*X_Mode[j,2]+
beta_Mode[9]*X_Mode[j,1]*X_Mode[j,3]+
beta_Mode[10]*X_Mode[j,2]*X_Mode[j,3];
}

mu_Mode_star=beta_Mode[1]+
beta_Mode[2]*x_star[1]+
beta_Mode[3]*x_star[2]+
beta_Mode[4]*x_star[3]+
beta_Mode[5]*x_star[1]^2+
beta_Mode[6]*x_star[2]^2+
beta_Mode[7]*x_star[3]^2+
beta_Mode[8]*x_star[1]*x_star[2]+
beta_Mode[9]*x_star[1]*x_star[3]+
beta_Mode[10]*x_star[2]*x_star[3];

for(j in 1:N_Aggwt){
mu_Aggwt[j]=beta_Aggwt[1]+
beta_Aggwt[2]*X_Aggwt[j,1]+
beta_Aggwt[3]*X_Aggwt[j,2]+
beta_Aggwt[4]*X_Aggwt[j,3]+
beta_Aggwt[5]*X_Aggwt[j,1]^2+
beta_Aggwt[6]*X_Aggwt[j,2]^2+
beta_Aggwt[7]*X_Aggwt[j,3]^2+
beta_Aggwt[8]*X_Aggwt[j,1]*X_Aggwt[j,2]+
beta_Aggwt[9]*X_Aggwt[j,1]*X_Aggwt[j,3]+
beta_Aggwt[10]*X_Aggwt[j,2]*X_Aggwt[j,3];
}

mu_Aggwt_star=beta_Aggwt[1]+
beta_Aggwt[2]*x_star[1]+
beta_Aggwt[3]*x_star[2]+
beta_Aggwt[4]*x_star[3]+
beta_Aggwt[5]*x_star[1]^2+
beta_Aggwt[6]*x_star[2]^2+
beta_Aggwt[7]*x_star[3]^2+
beta_Aggwt[8]*x_star[1]*x_star[2]+
beta_Aggwt[9]*x_star[1]*x_star[3]+
beta_Aggwt[10]*x_star[2]*x_star[3];





for(j in 1:N_Length){
mu_Length[j]=beta_Length[1]+
beta_Length[2]*X_Length[j,1]+
beta_Length[3]*X_Length[j,2]+
beta_Length[4]*X_Length[j,3]+
beta_Length[5]*X_Length[j,1]^2+
beta_Length[6]*X_Length[j,2]^2+
beta_Length[7]*X_Length[j,3]^2+
beta_Length[8]*X_Length[j,1]*X_Length[j,2]+
beta_Length[9]*X_Length[j,1]*X_Length[j,3]+
beta_Length[10]*X_Length[j,2]*X_Length[j,3];
}


mu_Length_star=beta_Length[1]+
beta_Length[2]*x_star[1]+
beta_Length[3]*x_star[2]+
beta_Length[4]*x_star[3]+
beta_Length[5]*x_star[1]^2+
beta_Length[6]*x_star[2]^2+
beta_Length[7]*x_star[3]^2+
beta_Length[8]*x_star[1]*x_star[2]+
beta_Length[9]*x_star[1]*x_star[3]+
beta_Length[10]*x_star[2]*x_star[3];


for(j in 1:N_Width){
mu_Width[j]=beta_Width[1]+
beta_Width[2]*X_Width[j,1]+
beta_Width[3]*X_Width[j,2]+
beta_Width[4]*X_Width[j,3]+
beta_Width[5]*X_Width[j,1]^2+
beta_Width[6]*X_Width[j,2]^2+
beta_Width[7]*X_Width[j,3]^2+
beta_Width[8]*X_Width[j,1]*X_Width[j,2]+
beta_Width[9]*X_Width[j,1]*X_Width[j,3]+
beta_Width[10]*X_Width[j,2]*X_Width[j,3];
}


mu_Width_star=beta_Width[1]+
beta_Width[2]*x_star[1]+
beta_Width[3]*x_star[2]+
beta_Width[4]*x_star[3]+
beta_Width[5]*x_star[1]^2+
beta_Width[6]*x_star[2]^2+
beta_Width[7]*x_star[3]^2+
beta_Width[8]*x_star[1]*x_star[2]+
beta_Width[9]*x_star[1]*x_star[3]+
beta_Width[10]*x_star[2]*x_star[3];


for(j in 1:N_Thickness){
mu_Thickness[j]=beta_Thickness[1]+
beta_Thickness[2]*X_Thickness[j,1]+
beta_Thickness[3]*X_Thickness[j,2]+
beta_Thickness[4]*X_Thickness[j,3]+
beta_Thickness[5]*X_Thickness[j,1]^2+
beta_Thickness[6]*X_Thickness[j,2]^2+
beta_Thickness[7]*X_Thickness[j,3]^2+
beta_Thickness[8]*X_Thickness[j,1]*X_Thickness[j,2]+
beta_Thickness[9]*X_Thickness[j,1]*X_Thickness[j,3]+
beta_Thickness[10]*X_Thickness[j,2]*X_Thickness[j,3];
}

mu_Thickness_star=beta_Thickness[1]+
beta_Thickness[2]*x_star[1]+
beta_Thickness[3]*x_star[2]+
beta_Thickness[4]*x_star[3]+
beta_Thickness[5]*x_star[1]^2+
beta_Thickness[6]*x_star[2]^2+
beta_Thickness[7]*x_star[3]^2+
beta_Thickness[8]*x_star[1]*x_star[2]+
beta_Thickness[9]*x_star[1]*x_star[3]+
beta_Thickness[10]*x_star[2]*x_star[3];
}

model{
beta_Mode~normal(0,10);
beta_Aggwt~normal(0,10);
beta_Length~normal(0,10);
beta_Width~normal(0,10);
beta_Thickness~normal(0,10);
sigma_Mode~normal(0,10);
sigma_Aggwt~normal(0,10);
sigma_Length~normal(0,10);
sigma_Width~normal(0,10);
sigma_Thickness~normal(0,10);

x_star[1]~uniform(rgHNO3[1],rgHNO3[2]);
x_star[2]~uniform(rgPu[1],rgPu[2]);
x_star[3]~uniform(rgTemp[1],rgTemp[2]);
// x_star~normal(0,1000);


y_Mode~normal(mu_Mode, sigma_Mode);
y_Aggwt~normal(mu_Aggwt, sigma_Aggwt);
y_Length~normal(mu_Length, sigma_Length);
y_Width~normal(mu_Width, sigma_Width);
y_Thickness~normal(mu_Thickness, sigma_Thickness);
y_star[1]~normal(mu_Mode_star, sigma_Mode);
y_star[2]~normal(mu_Aggwt_star, sigma_Aggwt);
y_star[3]~normal(mu_Length_star, sigma_Length);
y_star[4]~normal(mu_Width_star, sigma_Width);
y_star[5]~normal(mu_Thickness_star, sigma_Thickness);
}
'




##load workspace
#simulatedData<-0 #0=actual data, 1=simulated data.
#load(paste0('BurneyDataAnalysis_simData', simulatedData, '.RData'))
if(!simulatedData){
  test_cases<-c(4)
} else{
  test_cases<-c(1)
}
test_cases
pu_test[test_cases,x_names]
bayes_preds<-matrix(NA, nrow=length(test_cases), ncol=3)
fits_all<-list()
n_chains<-3
if(simulatedData){
  iter<-5e4
  #iter<-1e3
} else{
  iter<-1e5
  #iter<-1e3
}

j<-1
#prior range for uniforms on x_star
prior_range<-rbind(rgHNO3,rgPu,rgTemp)
prior_range
for(i in test_cases){
  df_train<-rbind(pu_train, pu_test[-i,])
  df_test<-pu_test[i,]
  for(nm in 1:length(y_names)){
    assign(y_names[nm],df_train[, c(y_names[nm], x_names)])
    comp<-complete.cases(get(y_names[nm]))
    assign(y_names[nm],get(y_names[nm])[comp,])
    assign(paste0('y_',  y_names[nm]),get(y_names[nm])[,y_names[nm]])
    assign(paste0('X_',  y_names[nm]),get(y_names[nm])[,x_names])
    assign(paste0('N_',  y_names[nm]),length(get(paste0('y_',  y_names[nm]))))
  }
  y_star<-as.numeric(df_test[,y_names])
  stan_data<-list(
    N_Mode=N_Mode,
    N_Aggwt=N_Aggwt,
    N_Length=N_Length,
    N_Width=N_Width,
    N_Thickness=N_Thickness,
    y_Mode=y_Mode,
    X_Mode=X_Mode,
    y_Aggwt=y_Aggwt,
    X_Aggwt=X_Aggwt,
    y_Length=y_Length,
    X_Length=X_Length,
    y_Width=y_Width,
    X_Width=X_Width,
    y_Thickness=y_Thickness,
    X_Thickness=X_Thickness,
    y_star=y_star,
    rgHNO3= prior_range[1,],
    rgPu= prior_range[2,],
    rgTemp= prior_range[3,]
  )
  
  init<-list()
  for(tt in 1:n_chains){
    init[[tt]]<-list('x_star'=full_preds[j,]+rnorm(3, sd=.1))
  }
  rbind(rgHNO3,rgPu,rgTemp)
  fit <- stan(model_code = burney.stan, data = stan_data, init=init,iter = iter,warmup=floor(.75*iter), chains=n_chains, control=list(adapt_delta=.95))
  fits_all[[j]]<-fit
  #x_star<-extract(fit, "x_star")
  #rstan::traceplot(fit, pars='x_star')
  bayes_preds[j,]<-summary(fit, pars = c("x_star"), probs = c(0.1, 0.9))$summary[,"mean"]
  j<-j+1
  print(i)
}

bayes_preds<-as.data.frame(bayes_preds)
names(bayes_preds)<-x_names

#save.image(paste0('BurneyDataAnalysisBayes_simData', simulatedData, '.RData'))



# end full bayes ---



#
# Summarize the results of the burney analysis
#

library(Hmisc)
library(effects)
library(stargazer)
library(MASS)
library(rstan)
library(ggplot2)
library(GGally)
library(bayesplot)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(tikzDevice)

#figDir<-file.path('~', 'Documents', 'PuSignatures','methods_for_inv_pred_paper', 'figs')
plotCors<-brewer.pal(9, 'Set1')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
par(pch=19)
#load workspace
# simulatedData<-0#0=actual data, 1=simulated data.
#load(paste0('BurneyDataAnalysis_simData', simulatedData, '.RData'))

#pairs plot
colLabs<-c('Mode', 'Agg. Index', 'Length', 'Width', 'Thickness', 'HNO$_3$', 'Pu', 'Temp.')
pairs_plot<-ggpairs(as.data.frame(pu), diag=list(continuous="barDiag"), lower = list(continuous = "cor"), upper = list(continuous = wrap("points", size=.75, alpha=.5)), columnLabels = colLabs)+theme_bw()+theme(axis.text.x=element_text(angle=60,hjust=1))
#tikz(file.path(figDir, paste0('pairsPlot_',simulatedData, '.tex')),width = 6.1, height = 6.1)
##pdf(paste0('figs/pairsPlot_', simulatedData, '.pdf'))
pairs_plot
#dev.off() 

#Linear Models
names(pu)
fit_mode<-lm(Mode~(HNO3+Pu+Temp)^2+I(HNO3^2)+I(Pu^2)+I(Temp^2), data=pu)
fit_agg<-lm(Aggwt~(HNO3+Pu+Temp)^2+I(HNO3^2)+I(Pu^2)+I(Temp^2), data=pu)
fit_length<-lm(Length~(HNO3+Pu+Temp)^2+I(HNO3^2)+I(Pu^2)+I(Temp^2), data=pu)
fit_width<-lm(Width~(HNO3+Pu+Temp)^2+I(HNO3^2)+I(Pu^2)+I(Temp^2), data=pu)
fit_thickness<-lm(Thickness~(HNO3+Pu+Temp)^2+I(HNO3^2)+I(Pu^2)+I(Temp^2), data=pu)
summary(fit_mode)
summary(fit_agg)
summary(fit_length)
summary(fit_width)
summary(fit_thickness)
ggcoef(fit_thickness)
#SUMMARIES-----

rmse<-function(x,y){
  sqrt(mean((y-x)^2))
}


rmse_df<-data.frame(HNO3=c(rmse(x_star$HNO3, full_preds$HNO3),
                           rmse(x_star$HNO3, step_preds$HNO3),
                           rmse(x_star$HNO3, map_full_preds$HNO3),
                           rmse(x_star$HNO3, map_step_preds$HNO3),
                           #rmse(x_star$HNO3, bayes_preds$HNO3),
                           #rmse(x_star$HNO3, inverse_preds$HNO3),
                           #rmse(x_star$HNO3, inverse_preds_step$HNO3),
                           rmse(x_star$HNO3, pcr_preds$HNO3),
                           rmse(x_star$HNO3, plsr_preds$HNO3),
                           rmse(x_star$HNO3, inverse_preds$HNO3),
                           rmse(x_star$HNO3, rep(mean(rgHNO3),nrow(x_star)))),
                    Pu=c(rmse(x_star$Pu, full_preds$Pu),
                         rmse(x_star$Pu, step_preds$Pu),
                         rmse(x_star$Pu, map_full_preds$Pu),
                         rmse(x_star$Pu, map_step_preds$Pu),
                         #rmse(x_star$Pu, bayes_preds$Pu),
                         #rmse(x_star$Pu, inverse_preds$Pu),
                         #rmse(x_star$Pu, inverse_preds_step$Pu),
                         rmse(x_star$Pu, pcr_preds$Pu),
                         rmse(x_star$Pu, plsr_preds$Pu),
                         rmse(x_star$Pu, inverse_preds$Pu),
                         rmse(x_star$Pu, rep(mean(rgPu),nrow(x_star)))),
                    Temp=c(rmse(x_star$Temp, full_preds$Temp),
                           rmse(x_star$Temp, step_preds$Temp),
                           rmse(x_star$Temp, map_full_preds$Temp),
                           rmse(x_star$Temp, map_step_preds$Temp),
                           #rmse(x_star$Temp, bayes_preds$Temp),
                           #rmse(x_star$Temp, inverse_preds$Temp),
                           #rmse(x_star$Temp, inverse_preds_step$Temp),
                           rmse(x_star$Temp, pcr_preds$Temp),
                           rmse(x_star$Temp, plsr_preds$Temp),
                           rmse(x_star$Temp, inverse_preds$Temp),
                           rmse(x_star$Temp, rep(mean(rgTemp),nrow(x_star)))))
rmse_df
rownames(rmse_df)<-c('Forward full', 
                     'Forward step', 
                     'MAP full', 
                     'MAP step', 
                     #'Bayes',
                     #'inverse full',
                     #'inverse step',
                     'PCR',
                     'PLSR',
                     'Main effects',
                     'Prior mean')
rmse_df


rmse_df_melt<-melt(as.matrix(rmse_df), value.name = 'RMSE')
names(rmse_df_melt)<-c('Model', 'Precipitation Variable','RMSE')

rmse_plot<-ggplot(rmse_df_melt,aes(x=Model, y=RMSE,group=`Precipitation Variable`, col=`Precipitation Variable`))+
  geom_line(size=1)+
  geom_point(size=4)+
  scale_colour_discrete(name="Precipitation Variable", labels=c('HNO$_3$', "Pu", 'Temperature'))+
  theme_bw()+theme(axis.text.x=element_text(angle=60,hjust=1), plot.title = element_text(hjust = 0.5))+ggtitle('RMSE of Predictions')
rmse_plot

#tikz(file.path(figDir, paste0('rmsePlot_',simulatedData, '.tex')),width = 7, height = 5)
##pdf(file.path(figDir, paste0('rmsePlot_',simulatedData, '.tex')))
rmse_plot
#dev.off() 





#Summary HNO3----

HNO3_preds<-data.frame(x_star=x_star$HNO3, 
                       forward_full=full_preds$HNO3,
                       forward_stepwise=step_preds$HNO3,
                       map_full=map_full_preds$HNO3,
                       map_stepwise=map_step_preds$HNO3,
                       #Bayes=bayes_preds$HNO3,
                       #inverse_full=inverse_preds$HNO3,
                       #inverse_step=inverse_preds_step$HNO3,
                       PCR=pcr_preds$HNO3, 
                       PLSR=plsr_preds$HNO3)

#Summary Pu----


Pu_preds<-data.frame(x_star=x_star$Pu, 
                     forward_full=full_preds$Pu,
                     forward_stepwise=step_preds$Pu,
                     map_full=map_full_preds$Pu,
                     map_stepwise=map_step_preds$Pu,
                     #Bayes=bayes_preds$Pu,
                     #inverse_full=inverse_preds$Pu,
                     #inverse_step=inverse_preds_step$Pu,
                     PCR=pcr_preds$Pu, 
                     PLSR=plsr_preds$Pu)


#Summary Temp----

Temp_preds<-data.frame(x_star=x_star$Temp, 
                       forward_full=full_preds$Temp,
                       forward_stepwise=step_preds$Temp,
                       map_full=map_full_preds$Temp,
                       map_stepwise=map_step_preds$Temp,
                       #Bayes=bayes_preds$Temp,
                       #inverse_full=inverse_preds$Temp,
                       #inverse_step=inverse_preds_step$Temp,
                       PCR=pcr_preds$Temp, 
                       PLSR=plsr_preds$Temp)


all.equal(pu_test[test_cases,'HNO3'], HNO3_preds$x_star, check.attributes=TRUE)
all.equal(pu_test[test_cases,'Pu'], Pu_preds$x_star, check.attributes=TRUE)
all.equal(pu_test[test_cases,'Temp'], Temp_preds$x_star, check.attributes=TRUE)
all.equal(as.matrix(pu_test[test_cases,x_names]), cbind(HNO3_preds$x_star,Pu_preds$x_star,Temp_preds$x_star), check.attributes=FALSE)

x_star<-pu_test[test_cases,x_names]





my_line <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1,..., col=2, lty=2, lwd=2)
}
par(pch=19, cex=.5)
pairs(HNO3_preds, lower.panel = my_line, upper.panel = my_line, main=expression(HNO[3]), cex=.5)
pairs(Pu_preds, lower.panel = my_line, upper.panel = my_line, main=expression(Pu), cex=.5)
pairs(Temp_preds, lower.panel = my_line, upper.panel = my_line, main=expression(Tempurature), cex=.5)


#------
#tikz(file.path(figDir, paste0('comparePlot_',simulatedData, '.tex')),width = 6, height = 5)
##pdf(paste0('figs/comparePlot_', simulatedData, '.pdf'))
par(mfrow=c(1,1))
plot(HNO3_preds$forward_full, HNO3_preds$map_full, type='n', ylab='MAP full', xlab='Forward full', main='Comparison of Classical and MAP Estimates of HNO$_3$')
text(HNO3_preds$forward_full, HNO3_preds$map_full,label=1:nrow(Temp_preds),col='blue')
abline(0,1, col=2, lwd=2, lty=2)
grid(col=1)
#dev.off()


plot(HNO3_preds$forward_stepwise, HNO3_preds$map_stepwise, type='n')
text(HNO3_preds$forward_stepwise, HNO3_preds$map_stepwise,label=1:nrow(Temp_preds),col='blue')
abline(0,1)

plot(Pu_preds$forward_full, Pu_preds$map_full, type='n')
text(Pu_preds$forward_full, Pu_preds$map_full,label=1:nrow(Temp_preds),col='blue')
abline(0,1)

plot(Pu_preds$forward_stepwise, Pu_preds$map_stepwise, type='n')
text(Pu_preds$forward_stepwise, Pu_preds$map_stepwise,label=1:nrow(Temp_preds),col='blue')
abline(0,1)


plot(Temp_preds$forward_full, Temp_preds$map_full, type='n')
text(Temp_preds$forward_full, Temp_preds$map_full,label=1:nrow(Temp_preds),col='blue')
abline(0,1)

plot(Temp_preds$forward_stepwise, Temp_preds$map_stepwise, type='n')
text(Temp_preds$forward_stepwise, Temp_preds$map_stepwise,label=1:nrow(Temp_preds),col='blue')
abline(0,1)


#Bayesian Summary -----

#load(paste0('BurneyDataAnalysisBayes_simData', simulatedData, '.RData'))

# # Check convergence of Bayesian fits -----
par(mfrow=c(1,3))
n_mcmc<-dim(fit)[1]
# for(j in 1:length(test_cases)){
# # j<-2
j<-1
print(range(rstan::summary(fits_all[[j]], probs = c(0.1, 0.9))$summary[,"Rhat"]))
draws <- as.array(fits_all[[j]], pars="x_star")
par(mfrow=c(3,1))
for(k in 1:3){
  plot(draws[,1,k], type='l',col=plotCors[4])
  lines(draws[,2,k], type='l', col=plotCors[5])
  lines(draws[,3,k], type='l', col=plotCors[6])
}

# mcmc_trace_highlight(draws,  facet_args = list(ncol = 1, strip.position = "left"), highlight=2) 

#pdf(paste0('figs/traceplot_', simulatedData, '.pdf'))
mcmc_trace(draws,  facet_args = list(ncol = 1, strip.position = "left")) 
#dev.off()

summary(rhat(fits_all[[j]],regex_pars='mu'))
rhats<-rhat(fits_all[[j]],regex_pars='mu')
rhats<-melt(rhats)
ggRhat<-ggplot(rhats, aes(x=value))+geom_histogram(col='black', fill='white')+
  labs(title='Potential Scale Reduction Statistics')+theme_bw()
ggRhat

##tikz(paste0('figs/Rhat_', simulatedData, '.tex'))
##pdf(paste0('figs/Rhat_', simulatedData, '.pdf'))
ggRhat
##dev.off()


plot((draws[,1,3]), type='l')
abline(h=x_star[test_cases[j],1], col=2)
plot(draws[,1,2], type='l')
abline(h=x_star[j,2], col=2)
plot(draws[,1,3], type='l', main=j)
abline(h=x_star[j,3], col=2)


# #Check the classical beta hats are close to the posterior means
postBetas<-lsBetas<-array(NA, c(10, length(test_cases), length(y_names)))
for(i in 1:length(test_cases)){
  df_train<-rbind(pu_train, pu_test[-test_cases[i],])
  df_test<-pu_test[test_cases[i],]
  for(j in 1:length(y_names)){
    postBetas[,i,j]<-rstan::summary(fits_all[[i]],pars=paste0("beta_", y_names[j]), probs = c(0.1, 0.9))$summary[,'mean']
    full_fit<-lm(df_train[,y_names[j]]~HNO3+Pu+Temp+HNO3*Pu+HNO3*Temp+Pu*Temp+I(HNO3^2)+I(Pu^2)+I(Temp^2), data=df_train)
    lsBetas[,i,j]<-coef(full_fit)
  }
}
par(mfrow=c(1,1))
plot(postBetas,lsBetas, pch=19)
abline(0,1, col=2, lwd=2)
cor(postBetas,lsBetas)
summary(lm(postBetas~lsBetas))


for(j in 1:length(test_cases)){
  print(range(rstan::summary(fits_all[[j]],pars='x_star', probs = c(0.1, 0.9))$summary[,"Rhat"]))
}
draws <- rstan::extract(fits_all[[j]], pars="x_star")$x_star




ggdf<-data.frame(draws[,1],draws[,2],draws[,3])
names(ggdf)<-x_names
estimates<-data.frame(rbind(full_preds[test_cases[j],],
                            map_full_preds[test_cases[j],],
                            step_preds[test_cases[j],],
                            map_step_preds[test_cases[j],],
                            pcr_preds[test_cases[j],],
                            plsr_preds[test_cases[j],],
                            inverse_preds[test_cases[j],]),
                      Model=c('forward full', 'MAP full', 'forward step', 'MAP step', 'PCR', 'PLSR','Main Effects'))

# ggdf2<-melt(ggdf)
# 
# gg1<-ggplot(data=ggdf2, aes(x=value), facets~variable) +
#   geom_density()
#   labs(title='Marginal posterior of HNO3')+
#   geom_density()+theme_bw()+theme(legend.position='left')+
#   geom_point(data=estimates, aes(x=HNO3,y=0, col=Model), size=3, pch=19)+
#   geom_vline(data=priors, aes(xintercept=limits), lty=2) 

estimates$Model<-factor(estimates$Model, levels=estimates$Model)
priors<-data.frame(limits=c(-1.5,1.5))
gg1<-ggplot(ggdf, aes(x=HNO3)) +
  geom_density()+theme_bw()+theme(legend.position='none')+
  geom_point(data=estimates, aes(x=HNO3,y=0.05, col=Model),shape=1:7, size=5, pch=16)+
  geom_vline(data=priors, aes(xintercept=limits), lty=2) +
  labs(title=expression(paste('Marginal posterior of HNO$_3$')))+xlab(expression(paste(' HNO$_3$ ')))

gg2<-ggplot(ggdf, aes(x=Pu)) +
  labs(title='Marginal posterior of Pu')+
  geom_density()+theme_bw()+theme(legend.position='none')+
  geom_point(data=estimates, aes(x=Pu,y=0.05, col=Model),shape=1:7, size=5, pch=16)+
  geom_vline(data=priors, aes(xintercept=limits), lty=2)
gg3<-ggplot(ggdf, aes(x=Temp)) +
  labs(title='Marginal posterior of Temperature')+ xlab('Temperature')+
  geom_density()+theme_bw()+theme(legend.position='none')+
  geom_point(data=estimates, aes(x=Temp,y=0.05, col=Model),shape=1:7, size=5, pch=16)+
  geom_vline(data=priors, aes(xintercept=limits), lty=2)


#density with legend
gg4<-ggplot(ggdf, aes(x=HNO3)) +
  geom_density()+theme_bw()+theme(legend.position='left')+
  geom_point(data=estimates, aes(x=HNO3,y=0.05, col=Model),shape=1:7, size=5, pch=16)
#Extract Legend 
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 
legend <- g_legend(gg4) 

plot_grid(gg1, gg2,gg3,legend, align='h')
#dev.off()
#tikz(file.path(figDir, paste0('marginalPosterior_',simulatedData, '.tex')),width = 6, height = 5)
##pdf(paste0('figs/marginalPosterior_', simulatedData, '.pdf'))
plot_grid(gg1, gg2,gg3,legend, align='h')
#dev.off() 

# end Bayesian Summary -----













