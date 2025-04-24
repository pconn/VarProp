#plots to visualize simulation study setup

gen_2_covs <- function(n_x,n_y,B1=0.5) {
  library(fields)
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

rqpois <- function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta-1))
}

#400 grid cells, each 2 x 2 (so half-width can be 1)
library(mrds)
library(dsm)
n_x = 25
n_y = 25
n_s = n_x*n_y
set.seed(12345)

XY =  expand.grid(y=c(1:n_y),x=c(1:n_x))-0.5  #locations of gridded centroids
n_sample = 81
det_intercept = -0.5
det_beta = 0.5
preddata = data.frame(XY)
preddata$area = 4
pred_df=preddata
pred_df$off.set=log(preddata$area)

Sampled_xy=expand.grid(y=c(1,4,7,10,13,16,19,22,25),x=c(1,4,7,10,13,16,19,22,25))
Sampled = (Sampled_xy$x-1)*25+Sampled_xy$y  #vector index

Which_no_sample = c(1:n_s)
Which_no_sample = Which_no_sample[-Sampled]

Covs = gen_2_covs(n_x,n_y,0.6)
#preddata$covdens=Cov_array[isim,,1]  #not needed if only spatial model
cor(Covs[,1],Covs[,2])

N_s=rqpois(n_s,exp(2+1*Covs[,1]),1.2)



# residual mapsplot_N_map_xy<-function(N,XY,leg.title="Abundance"){
#devtools::install_github("kwstat/pals")
plot_N_map_xy<-function(N,XY,leg.title="Abundance",my_color="blue"){
  require(ggplot2)
  library(scales)
  library(viridis)
  #library(pals)
  Abundance=N
  Cur.df=cbind(data.frame(x=XY[,1],y=XY[,2],Abundance))
  colnames(Cur.df)=c("x","y","Abundance")
  tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
  if(my_color=="red")p1=ggplot(Cur.df)+aes(x,y,fill=Abundance)+geom_raster()+tmp.theme+scale_fill_gradient(low="white",high=muted("red"),name=leg.title)
  if(my_color=="blue")p1=ggplot(Cur.df)+aes(x,y,fill=Abundance)+geom_raster()+tmp.theme+scale_fill_gradient(low="white",high=muted("blue"),name=leg.title)
  if(my_color=="viridis")p1=ggplot(Cur.df)+aes(x,y,fill=Abundance)+geom_raster()+tmp.theme+scale_fill_viridis(name=leg.title)
      p1
}

plot_xp = plot_N_map_xy(N=Covs[,2],XY=XY,leg.title="Cov-p",my_color="blue")
plot_xh = plot_N_map_xy(N=Covs[,1],XY=XY,leg.title="Cov-N",my_color="red")
plot_N = plot_N_map_xy(N=N_s,XY=XY,leg.title="N-true",my_color="viridis")+scale_fill_viridis(limits=c(0,176),name="N-true")

#detection probability 
Dists = c(0:100)/100
Zeros = rep(0,101)
Pplot_covs = c(min(Covs[,2]),mean(Covs[,2]),max(Covs[,2]))
Sigma = cbind(exp(det_intercept+det_beta*Pplot_covs[1]),exp(det_intercept+det_beta*Pplot_covs[2]),exp(det_intercept+det_beta*Pplot_covs[3]))
P=dnorm(Dists,Zeros,Sigma)/dnorm(Zeros,Zeros,Sigma)
Pplot_df = data.frame("Distance"=rep(Dists,3),
                      "p"=c(dnorm(Dists,Zeros,Sigma[1])/dnorm(Zeros,Zeros,Sigma[1]),
                            dnorm(Dists,Zeros,Sigma[2])/dnorm(Zeros,Zeros,Sigma[2]),
                            dnorm(Dists,Zeros,Sigma[3])/dnorm(Zeros,Zeros,Sigma[3])),
                      "Type"=rep(c("low","medium","high"),each=101)
                      )
Pplot_df$Type = factor(Pplot_df$Type,levels=c("high","medium","low"))
Pplot = ggplot(Pplot_df)+geom_line(aes(x=Distance,y=p,color=Type),size=1.2)


#Counts plot
N_covered = N_s[Sampled]
Cell_ID = rep(Sampled,N_covered)
n_covered = sum(N_covered)
Dists = runif(n_covered)
Zeros = rep(0,n_covered)
Dist_covs = rep(Covs[Sampled,2],N_covered)
Sigma = exp(det_intercept+det_beta*Dist_covs)
P=dnorm(Dists,Zeros,Sigma)/dnorm(Zeros,Zeros,Sigma)
P_samp=P[Sampled]
Observed = 1*(runif(n_covered)<P)
Which_obs = which(Observed==1)
Dist_obs = Dists[Which_obs]  #for detection function fitting
Cov_obs = Dist_covs[Which_obs]
n_obs=length(Which_obs)
Counts = tabulate(Cell_ID[Which_obs])
Counts[Which_no_sample]=NA

Count_plot=plot_N_map_xy(N=Counts,XY=XY,leg.title="Counts",my_color="viridis")
 
#Estimated N plot
dist_data <- data.frame("object"=c(1:n_obs),"observer"=rep(1,n_obs),
                        "detected"=rep(1,n_obs),"distance"=Dist_obs,
                        "covdet"=Cov_obs,"size"=1)
sim_ddf <- mrds::ddf(dsmodel=~mcds(key="hn",formula=~covdet),
                     meta.data=list(width=1),
                     data=dist_data)
obsdata = dist_data
SegID_all_animals = rep(c(1:n_sample),N_covered)
obsdata$Sample.Label = SegID_all_animals[Which_obs]
obsdata$size=1
segdata = data.frame(x=Sampled_xy$x,y=Sampled_xy$y,
                     Sample.Label=c(1:n_sample),
                     covdens=Covs[Sampled,1],
                     covdet=Covs[Sampled,2],
                     Effort=2)
pred_det <- predict(sim_ddf,newdata=segdata)$fitted

dsm_k8 <- dsm(count~te(x,y,k=8), ddf.obj=sim_ddf, 
              segment.data=segdata, observation.data=obsdata, 
              method="REML")
dsm_k8_pred <- predict(dsm_k8, preddata, preddata$area)
Nest_plot=plot_N_map_xy(N=dsm_k8_pred,XY=XY,leg.title="N-est",my_color="viridis")
Nest_plot = Nest_plot+scale_fill_viridis(limits=c(0,176),name="N-est")


library(gridExtra)
plot_xp <- plot_xp + annotate("text", x = 1, y = 27, label = "A.")
plot_xh <- plot_xh + annotate("text", x = 1, y = 27, label = "B.")
plot_N <- plot_N + annotate("text", x = 1, y = 27, label = "C.")
Pplot <- Pplot + annotate("text", x = 0, y = 1.1, label = "D.")
Count_plot <- Count_plot + annotate("text", x = 1, y = 27, label = "E.")
Nest_plot <- Nest_plot + annotate("text", x = 1, y = 27, label = "F.")

p1 = ggplotGrob(plot_xp)
p2 = ggplotGrob(plot_xh)
p3=  ggplotGrob(plot_N)
p4=ggplotGrob(Pplot)
p5=ggplotGrob(Count_plot)
p6=ggplotGrob(Nest_plot)
png("DSM_sim_design.png",width=6,height=6,units="in",res=600)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
#cowplot::plot_grid(Bearded_plot,Ringed_plot,ncol=1,plot.margin=margin(0,10,0,0))
dev.off()


