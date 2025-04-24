#Script xy_tw_DSM_2022.R...Megan C. Ferguson...27 January 2025
#adapted by P. Conn for preferential sampling/variance propagation study April 2025

  #Notes
  #
  # 1. This script builds the density surface model with the bivariate isotropic 
  #    thin plate regression spline for Eastern Bering Sea belugas in 2022 
  #    from Ferguson et al. (2025). This DSM has the following characteristics:
  #     a. bivariate and isotropic smooths of x and y, comprising thin-plate  
  #        regression smoothing splines with shrinkage
  #     b. Tweedie pdf for counts
  #     c. log link
  # 
  # 2. Required input files:
  #    a. Data/FergusonEtal_20250125_EBS_Beluga_DSM_data.Rdata
  #    b. Data/FergusonEtal_20250125_EBS_Beluga_Nht_data.Rdata
  #    c. Data/NSDL_Hex_4strata.Rdata  (hexagon definitions for plotting)

    library(mgcv)
    library(dsm)
    library(tweedie)

      set.k <- 200
      
    #Input necessary objects
      load("./varprop/Data/FergusonEtal_20250125_EBS_Beluga_DSM_data.Rdata")
      load("./varprop/Data/FergusonEtal_20250125_EBS_Beluga_Nht_data.Rdata")


      #extra detection probability stuff for Horvitz-Thompson
      p_avail = 0.5
      p_line = 0.753
      cv_p_line = 0.015  
      
    #Create mgcv model
      
      b <- gam(formula = seg.ind ~ s(x, y, bs="ts", k=set.k) +
                                             offset(log(a.p)),
                                      family=tw(link="log"), 
                                      method="REML",
                                      data=gam.dat22)   
      Lp <- predict(b, newdata=predgrid.strat, type="lpmatrix")
      pred_mgcv <- predgrid.strat$a.p*exp(Lp%*%coef(b)) 
      
      
      #model for beaufort sea state
      gam.dat22$iBeaufMin1 = gam.dat22$iBeauf-1
      bss_mod <- gam(formula = iBeauf ~ s(x, y, bs="ts", k=set.k),
                     method="REML",
                     data=gam.dat22)  
      pred_bss <- predict(bss_mod,newdata=predgrid.strat,type="response")
      cor(pred_bss,pred_mgcv)
      
      #model for turbidity
      gam.dat22$iTurb = (gam.dat22$Turb=="yes")*1
      turb_mod <- gam(formula = iTurb ~ s(x, y, bs="ts", k=set.k),
                     family="binomial",
                     method="REML",
                     data=gam.dat22)  
      pred_turb <- predict(turb_mod,newdata=predgrid.strat,type="response")
      cor(pred_turb,pred_mgcv)  
      
      #fit DSM with H-T responses
      gam.dat22$ht = gam.dat22$seg.ind/gam.dat22$ddf.pred
      predgrid.strat$area=predgrid.strat$a.p
      ht_mod <- gam(formula = ht ~ s(x, y, bs="ts", k=set.k) +
                 offset(log(area)),
               family=tw(link="log"), 
               method="REML",
               data=gam.dat22)   
      Lp <- predict(ht_mod, newdata=predgrid.strat, type="lpmatrix")
      pred_ht <- predgrid.strat$a.p*exp(Lp%*%coef(ht_mod)) 
      
      #plots
      library(ggplot2)
      library(viridis)
      library(sf)
      load('./varprop/data/NSDL_Hex_4strata.Rdata')
      Grid_sf = st_as_sf(strat.buff.hexPols)
      Grid_sf$N1=pred_mgcv
      Grid_sf$Nht=pred_ht
      Grid_sf$Turb = pred_turb
      Grid_sf$BSS = pred_bss
      Grid_sf$Ndiff = Grid_sf$N1-Grid_sf$Nht
      
      plot_N1 <- ggplot(Grid_sf)+geom_sf(aes(fill=N1))+scale_fill_viridis(name=expression(hat(N)[s]))+
        annotate("text", -Inf, Inf, label = "A.", hjust = 0, vjust = 1)+xlab("Latitude")+ylab("Longitude")+
        scale_x_continuous(breaks=c(-166,-164,-162))
      plot_N2 <- ggplot(Grid_sf)+geom_sf(aes(fill=Nht))+scale_fill_viridis()
      plot_Ndiff <- ggplot(Grid_sf)+geom_sf(aes(fill=Ndiff))+scale_fill_viridis(name=expression(hat(Delta)(N[s])))+
        annotate("text", -Inf, Inf, label = "D.", hjust = 0, vjust = 1)+xlab("Latitude")+ylab("Longitude")+
        scale_x_continuous(breaks=c(-166,-164,-162))
      plot_turb <- ggplot(Grid_sf)+geom_sf(aes(fill=Turb))+scale_fill_viridis()+
        annotate("text", -Inf, Inf, label = "C.", hjust = 0, vjust = 1)+xlab("Latitude")+ylab("Longitude")+
          scale_x_continuous(breaks=c(-166,-164,-162))
      plot_BSS <- ggplot(Grid_sf)+geom_sf(aes(fill=BSS))+scale_fill_viridis()+
        annotate("text", -Inf, Inf, label = "B.", hjust = 0, vjust = 1)+xlab("Latitude")+ylab("Longitude")+
        scale_x_continuous(breaks=c(-166,-164,-162))
      
      p1 = ggplotGrob(plot_N1)
      p2 = ggplotGrob(plot_BSS)
      p3=  ggplotGrob(plot_turb)
      p4=ggplotGrob(plot_Ndiff)
      grid.arrange(p1,p2,p3,p4,ncol=2)

      png("beluga_ests.png",width=6,height=6,units="in",res=600)
      cowplot::plot_grid(plot_N1,plot_BSS,plot_turb,plot_Ndiff,ncol=2,panel_spacing = unit(0.05, "lines"))
      dev.off()
      
              
          
