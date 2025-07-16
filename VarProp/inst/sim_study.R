library(fields)
library(mrds)
library(dsm)
library(ggplot2)

#function to generate two correlated, spatially autocorrelated
#covariates using a co-regionalization model
gen_2_covs <- function(n_x,n_y,B1=0.5) {
  #First, generate two independent spatial processes
  Locs = expand.grid(y=c(1:n_y),x=c(1:n_x))-0.5  #locations of gridded centroids
  Cov_exp = fields::Exp.cov(Locs,Locs,aRange=5)
  L = chol(Cov_exp)
  n=n_x*n_y
  Mu1 = t(L) %*% rnorm(n)
  Mu2 = t(L) %*% rnorm(n)
  
  Covs= matrix(Mu1,n,2)
  Covs[,2]=B1*Mu1 + (1-abs(B1))*Mu2
  Covs 
}

lognormal_CI <- function(estimate,CV=NULL,SE=NULL){
  if(is.null(SE)==0){
     varlog = log(1+SE^2/estimate^2)
  }
  if(is.null(CV)==0){
    varlog = log(1+CV^2)
  }
  C = exp(1.959964*sqrt(varlog))
  return(c(estimate/C,estimate*C))
}

rqpois <- function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta-1))
}

#400 grid cells, each 2 x 2 (so half-width can be 1)
n_sim=1000
n_boot = 1000
n_x = 25
n_y = 25
n_s = n_x*n_y
set.seed(12345)
B1_vec = rnorm(n_sim,0,0.5)
Which_out = which(abs(B1_vec)>1)
B1_vec[Which_out]=0
Cor_vec = rep(NA,n_sim)
Cov_array =array(NA,dim=c(n_sim,n_s,2))
N_s = matrix(0,n_sim,n_s)
XY =  expand.grid(y=c(1:n_y),x=c(1:n_x))-0.5  #locations of gridded centroids
n_sample = 81
Counts = rep(0,n_sample)
det_intercept = -0.5
det_beta = 0.5
preddata = data.frame(XY)
preddata$area = 4
pred_df=preddata
pred_df$off.set=log(preddata$area)

Sampled_xy=expand.grid(y=c(1,4,7,10,13,16,19,22,25),x=c(1,4,7,10,13,16,19,22,25))
Sampled = (Sampled_xy$x-1)*25+Sampled_xy$y  #vector index

Ests = Bias = matrix(NA,5,n_sim) #k5, k8, k8-wt, k8-HT, varprop
EDF = matrix(NA,2,n_sim)
Est_array = array(NA,dim=c(n_sim,n_s,5))
Varprop_diagnostics = matrix(NA,3,n_sim)
P_mat = matrix(0,n_sim,n_sample)
CV = Coverage = matrix(NA,7,n_sim) #k5, k8, k8-wt, k8-HT, varprop, k8 boot, ht boot
Zeros_boot = rep(0,n_sample)
Ones_boot = rep(1,n_sample)
N_boot = N_boot_ht = rep(NA,n_boot)

for(isim in 1:n_sim){
  #Sampled = sample(c(1:n_s),n_sample)  #vector index
  #Sampled_xy=XY[Sampled,]
  Which_no_sample = c(1:n_s)
  Which_no_sample = Which_no_sample[-Sampled]
  
  Cov_array[isim,,] = gen_2_covs(n_x,n_y,B1_vec[isim])
  #preddata$covdens=Cov_array[isim,,1]  #not needed if only spatial model
  Cor_vec[isim]=cor(Cov_array[isim,,1],Cov_array[isim,,2])
  N_s[isim,]=rqpois(n_s,exp(2+1*Cov_array[isim,,1]),1.2)
  N_true = sum(N_s[isim,])
  N_true_uncovered = sum(N_s[isim,Which_no_sample])
  N_covered = N_s[isim,Sampled]
  n_covered = sum(N_covered)
  Dists = runif(n_covered)
  Zeros = rep(0,n_covered)
  Dist_covs = rep(Cov_array[isim,Sampled,2],N_covered)
  Sigma = exp(det_intercept+det_beta*Dist_covs)
  P=dnorm(Dists,Zeros,Sigma)/dnorm(Zeros,Zeros,Sigma)
  P_mat[isim,]=P[Sampled]
  Observed = 1*(runif(n_covered)<P)
  Which_obs = which(Observed==1)
  Dist_obs = Dists[Which_obs]  #for detection function fitting
  Cov_obs = Dist_covs[Which_obs]
  n_obs=length(Which_obs)
  dist_data <- data.frame("object"=c(1:n_obs),"observer"=rep(1,n_obs),
                          "detected"=rep(1,n_obs),"distance"=Dist_obs,
                          "covdet"=Cov_obs,"size"=1)
  sim_ddf <- mrds::ddf(dsmodel=~mcds(key="hn",formula=~covdet),
                       meta.data=list(width=1),
                       data=dist_data)
  Beta_hat = sim_ddf$par
  VC_hat = solve(sim_ddf$hessian)
  
  #fit base DSM
  obsdata = dist_data
  SegID_all_animals = rep(c(1:n_sample),N_covered)
  obsdata$Sample.Label = SegID_all_animals[Which_obs]
  obsdata$size=1
  segdata = data.frame(x=Sampled_xy$x,y=Sampled_xy$y,
                       Sample.Label=c(1:n_sample),
                       covdens=Cov_array[isim,Sampled,1],
                       covdet=Cov_array[isim,Sampled,2],
                       Effort=2)
  pred_det <- predict(sim_ddf,newdata=segdata)$fitted
  wts_pinv <- 1/(pred_det)
  wts_pinv = wts_pinv/mean(wts_pinv)  #make weights sum to n
  
  # dsm_flat <- dsm(count~1, ddf.obj=sim_ddf, 
  #                 segment.data=segdata, observation.data=obsdata, 
  #                 method="REML")
  # dsm_flat_pred <- predict(dsm_flat,preddata,preddata$area)
  dsm_k5 <- dsm(count~te(x,y,k=5), ddf.obj=sim_ddf, 
                  segment.data=segdata, observation.data=obsdata, 
                  method="REML")
  EDF[1,isim]= summary(dsm_k5)$edf
  dsm_k8 <- dsm(count~te(x,y,k=8), ddf.obj=sim_ddf, 
                  segment.data=segdata, observation.data=obsdata, 
                  method="REML")
  EDF[2,isim]= summary(dsm_k8)$edf
  dsm_k8_wt <- dsm(count~te(x,y,k=8), ddf.obj=sim_ddf, 
                segment.data=segdata, observation.data=obsdata, 
                method="REML",weights=wts_pinv)
  dsm_k8_ht <- dsm(abundance.est~te(x,y,k=8), ddf.obj=sim_ddf, 
                   segment.data=segdata, observation.data=obsdata, 
                   method="REML")
  
  #predicted abundance / bias
  dsm_k5_pred <- predict(dsm_k5, preddata, preddata$area)
  Ests[1,isim] = sum(dsm_k5_pred[Which_no_sample])
  Est_array[isim,,1]=dsm_k5_pred
  dsm_k8_pred <- predict(dsm_k8, preddata, preddata$area)
  Ests[2,isim] = sum(dsm_k8_pred[Which_no_sample])
  Est_array[isim,,2]=dsm_k8_pred
  dsm_k8_wt_pred <- predict(dsm_k8_wt, preddata, preddata$area)
  Ests[3,isim] = sum(dsm_k8_wt_pred[Which_no_sample])
  Est_array[isim,,3]=dsm_k8_wt_pred
  dsm_k8_ht_pred <- predict(dsm_k8_ht, preddata, preddata$area)
  Ests[4,isim]=sum(dsm_k8_ht_pred[Which_no_sample])
  Est_array[isim,,4]=dsm_k8_ht_pred
  for(imod in 1:4)Bias[imod,isim]=(Ests[imod,isim]-N_true_uncovered)/N_true_uncovered

  # cv, coverage: independent method
  preddata.varprop <- split(preddata[Which_no_sample,], 1:length(Which_no_sample))
  var1_k5 <- dsm_var_gam(dsm_k5, pred.data=preddata.varprop,
                         off.set=preddata$area[Which_no_sample])
  CV[1,isim]=summary(var1_k5)$cv
  cur_CI = lognormal_CI(Ests[1,isim],CV=CV[1,isim])
  Coverage[1,isim]=(N_true_uncovered>cur_CI[1] & N_true_uncovered<cur_CI[2])

  var1_k8 <- dsm_var_gam(dsm_k8, pred.data=preddata.varprop,
                                         off.set=preddata$area[Which_no_sample])
  CV[2,isim]=summary(var1_k8)$cv
  cur_CI = lognormal_CI(Ests[2,isim],CV=CV[2,isim])
  Coverage[2,isim]=(N_true_uncovered>cur_CI[1] & N_true_uncovered<cur_CI[2])
  
  var1_k8_wt <- dsm_var_gam(dsm_k8_wt, pred.data=preddata.varprop,
                         off.set=preddata$area[Which_no_sample])
  CV[3,isim]=summary(var1_k8_wt)$cv
  cur_CI = lognormal_CI(Ests[3,isim],CV=CV[3,isim])
  Coverage[3,isim]=(N_true_uncovered>cur_CI[1] & N_true_uncovered<cur_CI[2])
  
  var1_k8_ht <- dsm_var_gam(dsm_k8_ht, pred.data=preddata.varprop,
                            off.set=preddata$area[Which_no_sample])
  CV[4,isim]=summary(var1_k8_ht)$cv
  cur_CI = lognormal_CI(Ests[4,isim],CV=CV[4,isim])
  Coverage[4,isim]=(N_true_uncovered>cur_CI[1] & N_true_uncovered<cur_CI[2])
  
  
  # varprop method
  preddata$off.set=preddata$area
  var2_k8 <- dsm_varprop(dsm_k8,newdata=preddata[Which_no_sample,])
  Ests[5,isim] = sum(var2_k8$pred)
  CV[5,isim]=sqrt(var2_k8$var)/Ests[5,isim]
  cur_CI = lognormal_CI(Ests[5,isim],CV=CV[5,isim])
  Coverage[5,isim]=(N_true_uncovered>cur_CI[1] & N_true_uncovered<cur_CI[2])
  Bias[5,isim]=(Ests[5,isim]-N_true_uncovered)/N_true_uncovered
  Est_array[isim,,5]=dsm_varprop(dsm_k8,newdata=preddata)$pred
  Varprop_diagnostics[,isim]=(summary(var2_k8)$varprop_diagnostic)[[1]][,4]-
    (summary(var2_k8)$varprop_diagnostic)[[1]][,2]
  
  # #varprop using original fit for point estimate
  # CV[3,isim]=sqrt(var2_k8$var)/Ests[1,isim]
  # cur_CI = lognormal_CI(Ests[1,isim],CV=CV[3,isim])
  # Coverage[3,isim]=(N_true_uncovered>cur_CI[1] & N_true_uncovered<cur_CI[2])
  
  #constructive bootstrap
  gam_cv = summary(var1_k8)$gam.cv
  Count_df = data.frame("count"=tabulate(obsdata$Sample.Label,nbins=n_sample),
            "x"=Sampled_xy$x,"y"=Sampled_xy$y)
  # Sigma_hat = exp(Beta_hat[1]+Beta_hat[2]*Cov_array[isim,Sampled,2])
  # P_hat = (pnorm(Ones_boot,Zeros_boot,Sigma_hat)-0.5)/
  #   dnorm(Zeros_boot,Zeros_boot,Sigma_hat)
  # Count_df$off.set = log(P_hat)+log(4)
  # my_gam = gam(count~offset(off.set)+s(x,y),family="quasipoisson",
  #              data=Count_df,method="REML")
  # gam_pred_lp = predict(my_gam,newdata=pred_df[Which_no_sample,],type="lpmatrix")
  # Vbeta <- vcov(my_gam)
  # Vp <- gam_pred_lp %*% (Vbeta %*% t(gam_pred_lp))
  # Cur_est = predict(my_gam,newdata=pred_df[Which_no_sample,],type="response")
  # var_est = Cur_est %*% Vp %*% Cur_est #vector of derivatives for exp transformation is just the estimates
  # cv_est = sqrt(var_est)/sum(Cur_est)
  for(iboot in 1:n_boot){
    Beta_boot = rmvnorm(1,Beta_hat,VC_hat)
    Sigma_boot = exp(Beta_boot[1]+Beta_boot[2]*Cov_array[isim,Sampled,2])
    P_boot = (pnorm(Ones_boot,Zeros_boot,Sigma_boot)-0.5)/
      dnorm(Zeros_boot,Zeros_boot,Sigma_boot)
    if(any(Sigma_boot>1000000000))P_boot[which(Sigma_boot>1000000000)]=1.0 #numerical issues here
    
    # count model
    Count_df$off.set = log(P_boot)+log(4) # p * area on real scale
    gam_boot = gam(count~offset(off.set)+te(x,y,k=8),family="quasipoisson",
                   data=Count_df,method="REML")
    gam_pred = predict(gam_boot,newdata=pred_df[Which_no_sample,],type="response")
    N_boot[iboot]=sum(gam_pred)
    
    #H-T model
    Count_df$ht = Count_df$count/P_boot
    Count_df$off.set = log(4)
    gam_boot = gam(ht~offset(off.set)+te(x,y,k=8),family="quasipoisson",
                   data=Count_df,method="REML")
    gam_pred = predict(gam_boot,newdata=pred_df[Which_no_sample,],type="response")
    N_boot_ht[iboot]=sum(gam_pred)
  }
  var_Boot = var(N_boot)+(summary(var1_k8)$gam.cv*Ests[2,isim])^2
  CV[6,isim] = sqrt(var_Boot)/Ests[2,isim]
  cur_CI = lognormal_CI(Ests[2,isim],CV=CV[6,isim])
  Coverage[6,isim]=(N_true_uncovered>cur_CI[1] & N_true_uncovered<cur_CI[2])
  
  var_Boot = var(N_boot_ht)+(summary(var1_k8_ht)$gam.cv*Ests[4,isim])^2
  CV[7,isim] = sqrt(var_Boot)/Ests[4,isim]
  cur_CI = lognormal_CI(Ests[4,isim],CV=CV[7,isim])
  Coverage[7,isim]=(N_true_uncovered>cur_CI[1] & N_true_uncovered<cur_CI[2])
}

save.image('sim_results_psp.RData')

# Locs = expand.grid(y=c(1:n_y),x=c(1:n_x))-0.5  #locations of gridded centroids
# Plot_df = Locs
# Plot_df$Mu1 = N_s[isim,]
# Plot_df$Ests_sp = dsm_base_pred
# Plot_df$Ests_sp_k8 = dsm_base8_pred
# Plot_df$Ests_flat = dsm_flat_pred
# Plot_df$RelBias_sp = (dsm_base_pred-N_s[isim,])/N_s[isim,]
# Plot_df$RelBias_sp_k8 = (dsm_base8_pred-N_s[isim,])/N_s[isim,]
# Plot_df$RelBias_flat = (dsm_flat_pred-N_s[isim,])/N_s[isim,]
# P_df = data.frame(covdet=Cov_array[isim,,2],size=1)
# Plot_df$P = predict(sim_ddf,newdata=P_df)[[1]]
# 
# ggplot(Plot_df)+geom_tile(aes(x=x,y=y,fill=P))+scale_fill_viridis()
# ggplot(Plot_df)+geom_tile(aes(x=x,y=y,fill=Ests_sp))+scale_fill_viridis_c()
# ggplot(Plot_df)+geom_tile(aes(x=x,y=y,fill=Ests_sp_k8))+scale_fill_viridis_c()
# ggplot(Plot_df)+geom_tile(aes(x=x,y=y,fill=Mu1))+scale_fill_viridis_c()
# ggplot(Plot_df)+geom_tile(aes(x=x,y=y,fill=RelBias_sp))
# 

Emp_cor <- rep(NA,1000)
for(isim in 1:1000)Emp_cor[isim]=cor(Cov_array[isim,,2],Est_array[isim,,2])

Plot_df = data.frame("Emp_cor"=Emp_cor,"True_cor"=Cor_vec)
ggplot(Plot_df)+geom_point(aes(x=True_cor,y=Emp_cor))+xlab(expression(rho(bold(x)[hab],bold(x)[p])))+
  ylab(expression(hat(rho)(hat(bold(N))[s],bold(x)[p])))
png(filename="Cor_scatter.png",width=5,height=5,units='in',res=600)
ggplot(Plot_df)+geom_point(aes(x=True_cor,y=Emp_cor))+xlab(expression(rho(bold(x)[hab],bold(x)[p])))+
  ylab(expression(hat(rho)(hat(bold(N))[s],bold(x)[p])))
dev.off()



Cor_pred = c(-100:100)/100
n_pred = length(Cor_pred)
Cov_mod_list = CV_mod_list= vector("list",4)
Bias_mod_list =  vector("list",5)
new_data = data.frame(Cor=Cor_pred)
for(i in 1:5){
  gam_df = data.frame("Cor"=Cor_vec) #,"Coverage"=Coverage[i,],"CV"=CV[i,])
  gam_df$Bias = Bias[i,]
  gam_bias = gam(Bias~s(Cor,k=4),data=gam_df)
  Bias_mod_list[[i]]=predict(gam_bias,newdata=new_data,type="response")
}
for(i in 1:7){
  gam_df = data.frame("Cor"=Cor_vec,"Coverage"=Coverage[i,],"CV"=CV[i,])
  gam_cov = gam(Coverage~s(Cor,k=4),data=gam_df,family="binomial")
  gam_cv = gam(CV~s(Cor,k=4),data=gam_df)
  Cov_mod_list[[i]]=predict(gam_cov,newdata=new_data,type="response")
  CV_mod_list[[i]]=predict(gam_cv,newdata=new_data,type="response")
}

Plot_df = data.frame("Correlation"=rep(Cor_pred,7),
                     "Coverage"=unlist(Cov_mod_list),
                     "Procedure"=c(rep("k=5,delta",n_pred),rep("delta",n_pred),
                                   rep("wtd,delta",n_pred),rep("HT,delta",n_pred),
                                   rep("varprop",n_pred),rep("boot",n_pred),
                                   rep("HT,boot",n_pred)))
ggplot(Plot_df)+geom_line(aes(x=Correlation,y=Coverage,color=Procedure),size=1.5)+
  scale_color_viridis(discrete=TRUE)
png(filename="Sim_coverage.png",width=5,height=5,units='in',res=600)
ggplot(Plot_df)+geom_line(aes(x=Correlation,y=Coverage,color=Procedure),size=1.5)+
  scale_color_viridis(discrete=TRUE)
dev.off()

Plot_df = data.frame("Correlation"=rep(Cor_pred,7),
                     "CV"=unlist(CV_mod_list),
                     "Procedure"=c(rep("k=5,delta",n_pred),rep("delta",n_pred),
                                   rep("wtd,delta",n_pred),rep("HT,delta",n_pred),
                                   rep("varprop",n_pred),rep("boot",n_pred),
                                   rep("HT,boot",n_pred)))
ggplot(Plot_df)+geom_line(aes(x=Correlation,y=CV,color=Procedure),size=1.5)+
  scale_color_viridis(discrete=TRUE)

png(filename="Sim_CV.png",width=5,height=5,units='in',res=600)
ggplot(Plot_df)+geom_line(aes(x=Correlation,y=CV,color=Procedure),size=1.5)+
  scale_color_viridis(discrete=TRUE)
dev.off()

Plot_df = data.frame("Correlation"=rep(Cor_pred,5),
                     "Bias"=unlist(Bias_mod_list),
                     "Procedure"=c(rep("k=5",n_pred),rep("base",n_pred),
                                                rep("wtd",n_pred),rep("HT",n_pred),
                                                rep("varprop",n_pred)))
ggplot(Plot_df)+geom_line(aes(x=Correlation,y=Bias,color=Procedure),linewidth=1.5)+
  scale_color_viridis(discrete=TRUE)
png(filename="Sim_bias.png",width=5,height=5,units='in',res=600)
ggplot(Plot_df)+geom_line(aes(x=Correlation,y=Bias,color=Procedure),linewidth=1.5)+
  scale_color_viridis(discrete=TRUE)
dev.off()



#check out correlation between estimated surface and estimated detection
