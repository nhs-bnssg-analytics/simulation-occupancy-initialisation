rm(list=ls())

require(tidyverse)
require(pracma)
require(Runuran)
require(truncdist)
require(scales)
require(parallel)


setwd("S:/Finance/Shared Area/BNSSG - BI/8 Modelling and Analytics/working/rw/online simulation initialisation/")


################################################################################
################################################################################
# SECTION 3: RESULTS WITH TIME-HOMOGENEOUS ARRIVAL RATE

# formulae
gen_rand_timehomog<-function(sf,n) {
  ftr<-function(x) sf(x)/integrate(sf,lower=0,upper=Inf)$value
  gen<-pinv.new(pdf=ftr,lb=0,ub=Inf)
  ur(gen,n)
}

# verification through simulation
sim_homog<-function(sf,n_entries) {
  wup<-uq(pinv.new(cdf=function(x) 1-sf(x),lb=0,ub=Inf),0.999)
  sim<-data.frame(entry=seq(0,wup+100*wup,length.out=n_entries)) %>%
    mutate(end=entry+ur(pinv.new(cdf=function(x) 1-sf(x),lb=0,ub=Inf),nrow(.)))
  bind_rows(lapply(seq(wup,max(sim$entry),length.out=100), function(x) {
    sim %>%
      filter(entry<x & end>x) %>%
      mutate(partial_los_from_entry=x-entry) %>%
      mutate(partial_los_to_end=end-x) %>%
      select(-c(entry,end))
  }))
}

eqn_sim_compare<-function(sf,sim_n_entries,plot_title) {
  sim<-sim_homog(sf=sf,n_entries=sim_n_entries)
  eqn<-gen_rand_timehomog(sf=sf,n=nrow(sim))
  data.frame(expo_sim=sim,expo_eqn=eqn) %>%
    pivot_longer(cols=everything(),names_to="type",values_to="value") %>%
    ggplot(aes(x=value,col=type)) +
    stat_ecdf() +
    labs(title=paste0(plot_title)) +
    xlab("Time") +
    ylab("Probability distribution") +
    theme_bw() +
    theme(legend.position="bottom")
}


########################################
# SECTION 3.2

plt_3_2_dat<-bind_rows(
  sim_homog(sf=function(x) exp(-2*x),n_entries=1000000) %>%
    mutate(eqn=gen_rand_timehomog(sf=function(x) exp(-2*x),n=nrow(.))) %>%
    pivot_longer(cols=everything(),names_to="type",values_to="value") %>%
    mutate(distribution="exponential(2)"),
  sim_homog(sf=function(x) 0.5-0.5*erf((log(x)-((log(0.5)-0.5)))/(1*sqrt(2))),n_entries=1000000) %>%
    mutate(eqn=gen_rand_timehomog(sf=function(x) 0.5-0.5*erf((log(x)-((log(0.5)-0.5)))/(1*sqrt(2))),n=nrow(.))) %>%
    pivot_longer(cols=everything(),names_to="type",values_to="value") %>%
    mutate(distribution="lognormal(-1.19,1)"),
  sim_homog(sf=function(x) 0.5-0.5*erf((log(x)-((log(0.5)-4.5)))/(1.5*sqrt(2))),n_entries=1000000) %>%
    mutate(eqn=gen_rand_timehomog(sf=function(x) 0.5-0.5*erf((log(x)-((log(0.5)-4.5)))/(1.5*sqrt(2))),n=nrow(.))) %>%
    pivot_longer(cols=everything(),names_to="type",values_to="value") %>%
    mutate(distribution="lognormal(-1.47,1.25)")
)%>%
  mutate(type=recode(type,"eqn"="Formulae",
                     "partial_los_from_entry"="Simulation (r)",
                     "partial_los_to_end"="Simulation (p)")) %>%
  mutate(type=factor(type,levels=c("Formulae","Simulation (r)","Simulation (p)")))


png("fig2.png",height=3,width=5,unit="in",res=600)
print(plt_3_2_dat %>%
  ggplot(aes(x=value,colour=type)) +
    stat_ecdf(data = plt_3_2_dat %>% filter(type=="Formulae"), linetype=1) +
    stat_ecdf(data = plt_3_2_dat %>% filter(type=="Simulation (p)"),linetype=2) +
    stat_ecdf(data = plt_3_2_dat %>% filter(type=="Simulation (r)"), linetype=3) +
    facet_wrap(~distribution,nrow=1,scales="free_x") +
    xlab("Partial service time") +
    coord_cartesian(xlim=c(0,10)) +
    scale_x_continuous(breaks=c(0,2,4,6,8,10)) +
    ylab("Probability distribution") +
    theme_bw() +
    theme(legend.position="bottom",
          legend.title=element_blank()) +
    guides(size = "none")
  
)
dev.off()


########################################
# SECTION 3.3: impact on main simulation results

# exponential

sim_res_exp<-function(init_method) {
  sqraw<-seq(0,10,0.1)
  sq<-round(sqraw,1)
  cl<-makeCluster(detectCores()-1)
  clusterExport(cl=cl,varlist="gen_rand_timehomog")
  sim_res_tmp<-do.call("cbind",parLapply(cl,1:5000,function(rep) {
    require(Runuran)
    require(pracma)
    require(tidyverse)
    require(truncdist)
    set.seed(rep+ifelse(init_method=="formulae",1000000,ifelse(init_method=="method1",2000000,ifelse(init_method=="method2",3000000,4000000))))
    if(init_method=="formulae") dfx<-data.frame(arrs=-gen_rand_timehomog(sf=function(x) exp(-2*x),n=0.5*100))
    if(init_method=="method1") dfx<-data.frame(arrs=rep(-0.01,0.5*100))
    if(init_method=="method2") dfx<-data.frame(arrs=-rexp(0.5*100,rate=2))
    if(init_method=="method3") {
      tmp_mean<-1/2
      tmp_var<-1/(2^2)
      dfx<-data.frame(arrs=rep(-0.5*(tmp_mean+tmp_var/tmp_mean),0.5*100))
    }
    log_old_sims<-dfx %>%
      rowwise() %>%
      mutate(end=rtrunc(1,"exp",a=-arrs,b=Inf,rate=2)+arrs)
    log_new_sims<-data.frame(arrs=round(seq(0.005,10.005,0.01),3)) %>%
      mutate(end=arrs+rexp(nrow(.),rate=2))
    logx_sims<-rbind(log_old_sims,log_new_sims)
    res_sims<-data.frame(xname=sapply(sq,function(x) length(which(logx_sims$arrs<x & logx_sims$end>x))))
    #plot(res_sims$xname,type="l")
    names(res_sims)<-paste0("rep_",rep)
    return(res_sims)
  }))
  stopCluster(cl)
  sim_res_tmp %>%
    mutate(time=sq) %>%
    pivot_longer(cols=-time,names_to="rep",values_to="occupancy",names_prefix="rep_") %>%
    group_by(time) %>%
    summarise(mean_occ=mean(occupancy),occ_l80=quantile(occupancy,0.1),occ_u80=quantile(occupancy,0.9),
              occ_l95=quantile(occupancy,0.025),occ_u95=quantile(occupancy,0.975))
}

sim_res_exp_formulae<-sim_res_exp("formulae") %>% mutate(init_method="Formulae")
sim_res_exp_method1<-sim_res_exp("method1") %>% mutate(init_method="Method 1")
sim_res_exp_method2<-sim_res_exp("method2") %>% mutate(init_method="Method 2")
sim_res_exp_method3<-sim_res_exp("method3") %>% mutate(init_method="Method 3")
sim_res_exp_all<-bind_rows(sim_res_exp_formulae,sim_res_exp_method1,sim_res_exp_method2,sim_res_exp_method3)

sim_res_exp_all %>%
  ggplot(aes(x=time)) +
  geom_ribbon(aes(ymin=occ_l95,ymax=occ_u95),fill="lightblue3",alpha=0.5) +
  geom_ribbon(aes(ymin=occ_l80,ymax=occ_u80),fill="blue",alpha=0.5) +
  geom_line(aes(y=mean_occ)) +
  facet_wrap(~init_method,nrow=1) +
  theme_bw()

# lognormal (sigma = 1)

sim_res_lnorm<-function(init_method) {
  sqraw<-seq(0,10,0.1)
  sq<-round(sqraw,1)
  cl<-makeCluster(detectCores()-1)
  clusterExport(cl=cl,varlist="gen_rand_timehomog")
  sim_res_tmp<-do.call("cbind",parLapply(cl,1:5000,function(rep) {
    require(Runuran)
    require(pracma)
    require(tidyverse)
    require(truncdist)
    set.seed(rep+ifelse(init_method=="formulae",1000000,ifelse(init_method=="method1",2000000,ifelse(init_method=="method2",3000000,4000000))))
    if(init_method=="formulae") dfx<-data.frame(arrs=-gen_rand_timehomog(sf=function(x) 0.5-0.5*erf((log(x)-((log(0.5))-0.5))/(1*sqrt(2))),
                                                                        n=0.5*100))
    if(init_method=="method1") dfx<-data.frame(arrs=rep(-0.1,0.5*100))
    if(init_method=="method2") dfx<-data.frame(arrs=-rlnorm(0.5*100,meanlog=(log(0.5))-0.5,sdlog=1))
    if(init_method=="method3") {
      tmp_mean<-0.5
      tmp_var<-(exp(1^2)-1)*exp(2*((log(0.5))-0.5)+1^2)
      dfx<-data.frame(arrs=rep(-0.5*(tmp_mean+tmp_var/tmp_mean),0.5*100))
    }
    log_old_sims<-dfx %>%
      rowwise() %>%
      mutate(end=rtrunc(1,"lnorm",a=-arrs,b=Inf,meanlog=(log(0.5))-0.5,sdlog=1)+arrs)
    log_new_sims<-data.frame(arrs=round(seq(0.005,10.005,0.01),3)) %>%
      mutate(end=arrs+rlnorm(nrow(.),meanlog=(log(0.5))-0.5,sdlog=1))
    logx_sims<-rbind(log_old_sims,log_new_sims)
    res_sims<-data.frame(xname=sapply(sq,function(x) length(which(logx_sims$arrs<x & logx_sims$end>x))))
    names(res_sims)<-paste0("rep_",rep)
    return(res_sims)
  }))
  stopCluster(cl)
  sim_res_tmp %>%
    mutate(time=sq) %>%
    pivot_longer(cols=-time,names_to="rep",values_to="occupancy",names_prefix="rep_") %>%
    group_by(time) %>%
    summarise(mean_occ=mean(occupancy),occ_l80=quantile(occupancy,0.1),occ_u80=quantile(occupancy,0.9),
              occ_l95=quantile(occupancy,0.025),occ_u95=quantile(occupancy,0.975))
}

sim_res_lnorm_formulae<-sim_res_lnorm("formulae") %>% mutate(init_method="Formulae")
sim_res_lnorm_method1<-sim_res_lnorm("method1") %>% mutate(init_method="Method 1")
sim_res_lnorm_method2<-sim_res_lnorm("method2") %>% mutate(init_method="Method 2")
sim_res_lnorm_method3<-sim_res_lnorm("method3") %>% mutate(init_method="Method 3")
sim_res_lnorm_all<-bind_rows(sim_res_lnorm_formulae,sim_res_lnorm_method1,sim_res_lnorm_method2,sim_res_lnorm_method3)

sim_res_lnorm_all %>%
  ggplot(aes(x=time)) +
  geom_ribbon(aes(ymin=occ_l95,ymax=occ_u95),fill="lightblue3",alpha=0.5) +
  geom_ribbon(aes(ymin=occ_l80,ymax=occ_u80),fill="blue",alpha=0.5) +
  geom_line(aes(y=mean_occ)) +
  facet_wrap(~init_method,nrow=1) +
  theme_bw()

# lognormal (sigma = 1.25)

sim_res_lnorm2<-function(init_method) {
  sqraw<-seq(0,10,0.1)
  sq<-round(sqraw,1)
  cl<-makeCluster(detectCores()-1)
  clusterExport(cl=cl,varlist="gen_rand_timehomog")
  sim_res_tmp<-do.call("cbind",parLapply(cl,1:5000,function(rep) {
    require(Runuran)
    require(pracma)
    require(tidyverse)
    require(truncdist)
    set.seed(rep+ifelse(init_method=="formulae",1000000,ifelse(init_method=="method1",2000000,ifelse(init_method=="method2",3000000,4000000))))
    if(init_method=="formulae") dfx<-data.frame(arrs=-gen_rand_timehomog(sf=function(x) 0.5-0.5*erf((log(x)-((log(0.5))-0.78125))/(1.25*sqrt(2))),
                                                                         n=0.5*100))
    if(init_method=="method1") dfx<-data.frame(arrs=rep(-0.1,0.5*100))
    if(init_method=="method2") dfx<-data.frame(arrs=-rlnorm(0.5*100,meanlog=(log(0.5))-0.78125,sdlog=1.25))
    if(init_method=="method3") {
      tmp_mean<-0.5
      tmp_var<-(exp(1.25^2)-1)*exp(2*((log(0.5))-0.78125)+1.25^2)
      dfx<-data.frame(arrs=rep(-0.5*(tmp_mean+tmp_var/tmp_mean),0.5*100))
    }
    log_old_sims<-dfx %>%
      rowwise() %>%
      mutate(end=rtrunc(1,"lnorm",a=-arrs,b=Inf,meanlog=(log(0.5))-0.78125,sdlog=1.25)+arrs)
    log_new_sims<-data.frame(arrs=round(seq(0.005,10.005,0.01),3)) %>%
      mutate(end=arrs+rlnorm(nrow(.),meanlog=(log(0.5))-0.78125,sdlog=1.25))
    logx_sims<-rbind(log_old_sims,log_new_sims)
    res_sims<-data.frame(xname=sapply(sq,function(x) length(which(logx_sims$arrs<x & logx_sims$end>x))))
    names(res_sims)<-paste0("rep_",rep)
    return(res_sims)
  }))
  stopCluster(cl)
  sim_res_tmp %>%
    mutate(time=sq) %>%
    pivot_longer(cols=-time,names_to="rep",values_to="occupancy",names_prefix="rep_") %>%
    group_by(time) %>%
    summarise(mean_occ=mean(occupancy),occ_l80=quantile(occupancy,0.1),occ_u80=quantile(occupancy,0.9),
              occ_l95=quantile(occupancy,0.025),occ_u95=quantile(occupancy,0.975))
}

sim_res_lnorm2_formulae<-sim_res_lnorm2("formulae") %>% mutate(init_method="Formulae")
sim_res_lnorm2_method1<-sim_res_lnorm2("method1") %>% mutate(init_method="Method 1")
sim_res_lnorm2_method2<-sim_res_lnorm2("method2") %>% mutate(init_method="Method 2")
sim_res_lnorm2_method3<-sim_res_lnorm2("method3") %>% mutate(init_method="Method 3")
sim_res_lnorm2_all<-bind_rows(sim_res_lnorm2_formulae,sim_res_lnorm2_method1,sim_res_lnorm2_method2,sim_res_lnorm2_method3)

sim_res_lnorm2_all %>%
  ggplot(aes(x=time)) +
  geom_ribbon(aes(ymin=occ_l95,ymax=occ_u95),fill="lightblue3",alpha=0.5) +
  geom_ribbon(aes(ymin=occ_l80,ymax=occ_u80),fill="blue",alpha=0.5) +
  geom_line(aes(y=mean_occ)) +
  facet_wrap(~init_method,nrow=1) +
  theme_bw()


### combined plot

sim_res_exp_lnorm_lnorm2<-bind_rows(sim_res_exp_all %>% mutate(dist="exponential(2)"),
                             sim_res_lnorm_all %>% mutate(dist="lognormal(-1.19,1)"),
                             sim_res_lnorm2_all %>% mutate(dist="lognormal(-1.47,1.25)"))
saveRDS(sim_res_exp_lnorm_lnorm2,file="res_3_3.R")

#sim_res_exp_lnorm_lnorm2<-readRDS("res_3_3.R")

png("fig3.png",units="in",height=5,width=6,res=600)
print(
  sim_res_exp_lnorm_lnorm2 %>%
    ggplot(aes(x=time)) +
    geom_ribbon(aes(ymin=occ_l95,ymax=occ_u95),fill="skyblue",alpha=0.3) +
    geom_ribbon(aes(ymin=occ_l80,ymax=occ_u80),fill="skyblue4",alpha=0.3) +
    geom_line(aes(y=mean_occ)) +
    facet_grid(dist~init_method) +
    scale_x_continuous(breaks=c(0,2,4,6,8,10)) +
    xlab("Calendar time") +
    ylab("Occupancy") +
    theme_bw()
)
dev.off()

# errors for method 1 vs formulae
drops_method1<-sim_res_exp_lnorm_lnorm2 %>%
  filter(init_method=="Method 1") %>%
  filter(time>0 & time<1) %>%
  group_by(init_method,dist) %>%
  summarise(min=min(mean_occ),max=max(mean_occ))
drops_method1_ln1<-drops_method1 %>% filter(dist=="lognormal(-1.19,1)") %>% .$min
drops_method1_ln2<-drops_method1 %>% filter(dist=="lognormal(-1.47,1.25)") %>% .$min
50-drops_method1_ln1
(50-drops_method1_ln1)/50
50-drops_method1_ln2
(50-drops_method1_ln2)/50

# errors for method 2 vs formulae
drops_method2<-sim_res_exp_lnorm_lnorm2 %>%
  filter(init_method=="Method 2") %>%
  filter(time>0 & time<1) %>%
  group_by(init_method,dist) %>%
  summarise(min=min(mean_occ),max=max(mean_occ))
drops_method2_ln1<-drops_method2 %>% filter(dist=="lognormal(-1.19,1)") %>% .$max
drops_method2_ln2<-drops_method2 %>% filter(dist=="lognormal(-1.47,1.25)") %>% .$max
50-drops_method2_ln1
(50-drops_method2_ln1)/50
50-drops_method2_ln2
(50-drops_method2_ln2)/50

# errors for method 3 vs formulae
drops_method3<-sim_res_exp_lnorm_lnorm2 %>%
  filter(init_method=="Method 3") %>%
  filter(time>0.5 & time<1.5) %>%
  group_by(init_method,dist) %>%
  summarise(min=min(mean_occ),max=max(mean_occ))
drops_method3_ln1<-drops_method3 %>% filter(dist=="lognormal(-1.19,1)") %>% .$max
drops_method3_ln2<-drops_method3 %>% filter(dist=="lognormal(-1.47,1.25)") %>% .$max
50-drops_method3_ln1
(50-drops_method3_ln1)/50
50-drops_method3_ln2
(50-drops_method3_ln2)/50


################################################################################
################################################################################
# SECTION 4: RESULTS WITH TIME-INHOMOGENEOUS ARRIVAL RATE

########################################
# SECTION 4.1

png("fig4.png",units="in",height=2,width=2.5,res=600)
print(
  data.frame(x=seq(-10,0,0.01),y=200/2*(1+sin(-(seq(-10,0,0.01))*pi*2))) %>%
    ggplot(aes(x=x,y=y)) +
    geom_line() +
    xlab("Calendar time") +
    ylab("Arrivals") +
    scale_x_continuous(breaks=-c(10,8,6,4,2,0)) +
    theme_bw()
)
dev.off()

########################################
# SECTION 4.2

#### exponential

# simulation
tsq<-round(seq(-10,0,0.01),2)
tdf1<-data.frame(time=tsq,nstarts=round(100*200/2*(1+sin(-2*pi*tsq))))
inhomog_sim_expo<-bind_rows(lapply(1:nrow(tdf1),function(x) {
  tmpt<-tdf1$time[x]
  tmpn<-tdf1$nstarts[x]
  data.frame(start=rep(tmpt,tmpn),end=tmpt+rexp(tmpn,rate=2))
})) %>%
  filter(start<0 & end>0) %>%
  mutate(partial=0-start) %>%
  mutate(distn="exponential(2)",type="Simulation") %>%
  select(partial,distn,type)

# formulae
gen<-pinv.new(pdf=function(x) 1.55156*(1+sin(2*pi*x))*exp(-2*x),lb=0,ub=Inf)
inhomog_eqn_expo<-data.frame(partial=ur(gen,nrow(inhomog_sim_expo))) %>%
  mutate(distn="exponential(2)",type="Formulae")

inhomog_expo<-bind_rows(inhomog_sim_expo,inhomog_eqn_expo)

#### lognormal (sigma = 1)

# simulation
tsq<-round(seq(-10,0,0.01),2)
tdf1<-data.frame(time=tsq,nstarts=round(100*200/2*(1+sin(-2*pi*tsq))))
inhomog_sim_lnorm<-bind_rows(lapply(1:nrow(tdf1),function(x) {
  tmpt<-tdf1$time[x]
  tmpn<-tdf1$nstarts[x]
  data.frame(start=rep(tmpt,tmpn),end=tmpt+rlnorm(tmpn,meanlog=log(0.5)-0.5,sdlog=1))
})) %>%
  filter(start<0 & end>0) %>%
  mutate(partial=0-start) %>%
  mutate(distn="lognormal(-1.19,1)",type="Simulation") %>%
  select(partial,distn,type)

# formulae
ftr<-function(x) (1+sin(2*pi*x))*(0.5-0.5*erf((0.5+log(2*x))/sqrt(2)))/
                  integrate(function(y) (1+sin(2*pi*y))*(0.5-0.5*erf((0.5+log(2*y))/sqrt(2))),lower=0,upper=Inf)$value
gen<-pinv.new(pdf=ftr,lb=0.000001,ub=10)
inhomog_eqn_lnorm<-data.frame(partial=ur(gen,nrow(inhomog_sim_lnorm))) %>%
  mutate(distn="lognormal(-1.19,1)",type="Formulae")

inhomog_lnorm<-bind_rows(inhomog_sim_lnorm,inhomog_eqn_lnorm)

#### lognormal (sigma = 1.25)

# simulation
tsq<-round(seq(-10,0,0.01),2)
tdf1<-data.frame(time=tsq,nstarts=round(100*200/2*(1+sin(-2*pi*tsq))))
inhomog_sim_lnorm2<-bind_rows(lapply(1:nrow(tdf1),function(x) {
  tmpt<-tdf1$time[x]
  tmpn<-tdf1$nstarts[x]
  data.frame(start=rep(tmpt,tmpn),end=tmpt+rlnorm(tmpn,meanlog=log(0.5)-0.78125,sdlog=1.25))
})) %>%
  filter(start<0 & end>0) %>%
  mutate(partial=0-start) %>%
  mutate(distn="lognormal(-1.47,1.25)",type="Simulation") %>%
  select(partial,distn,type)

# formulae
ftr<-function(x) (1+sin(2*pi*x))*(0.5-0.5*erf((0.78125+log(2*x))/(1.25*sqrt(2))))/
  integrate(function(y) (1+sin(2*pi*y))*(0.5-0.5*erf((0.78125+log(2*y))/(1.25*sqrt(2)))),lower=0,upper=Inf)$value
gen<-pinv.new(pdf=ftr,lb=0.000001,ub=10)
inhomog_eqn_lnorm2<-data.frame(partial=ur(gen,nrow(inhomog_sim_lnorm2))) %>%
  mutate(distn="lognormal(-1.47,1.25)",type="Formulae")

inhomog_lnorm2<-bind_rows(inhomog_sim_lnorm2,inhomog_eqn_lnorm2)


#### all

plt_4_2_dat<-bind_rows(inhomog_expo,inhomog_lnorm,inhomog_lnorm2)

png("fig5.png",height=3,width=5,unit="in",res=600)
print(
  ggplot(plt_4_2_dat,aes(x=partial,colour=type)) +
    stat_ecdf(data = plt_4_2_dat %>% filter(type=="Formulae"), linetype=1) +
    stat_ecdf(data = plt_4_2_dat %>% filter(type=="Simulation"), linetype=2) +
    facet_wrap(~distn,nrow=1,scales="free_x") +
    xlab("Partial service time") +
    coord_cartesian(xlim=c(0,10)) +
    scale_x_continuous(breaks=c(0,2,4,6,8,10)) +
    ylab("Probability distribution") +
    theme_bw() +
    theme(legend.position="bottom",
          legend.title=element_blank())
)
dev.off()


########################################
# SECTION 4.3

# exponential

sim_res_exp<-function(init_method) {
  sqraw<-seq(0,10,0.1)
  sq<-round(sqraw,1)
  cl<-makeCluster(detectCores()-1)
  sim_res_tmp<-do.call("cbind",parLapply(cl,1:5000,function(rep) {
    require(Runuran)
    require(pracma)
    require(tidyverse)
    require(truncdist)
    set.seed(rep+ifelse(init_method=="formulae",1000000,ifelse(init_method=="method1",2000000,ifelse(init_method=="method2",3000000,4000000))))
    if(init_method=="formulae") {
      gen<-pinv.new(pdf=function(x) 1.55156*(1+sin(2*pi*x))*exp(-2*x),lb=0,ub=Inf)
      dfx<-data.frame(arrs=-ur(gen,64))
    }
    if(init_method=="method1") dfx<-data.frame(arrs=rep(-0.01,64))
    if(init_method=="method2") dfx<-data.frame(arrs=-rexp(64,rate=2))
    if(init_method=="method3") {
      tmp_mean<-1/2
      tmp_var<-1/(2^2)
      dfx<-data.frame(arrs=rep(-0.5*(tmp_mean+tmp_var/tmp_mean),64))
    }
    log_old_sims<-dfx %>%
      rowwise() %>%
      mutate(end=rtrunc(1,"exp",a=-arrs,b=Inf,rate=2)+arrs)
    tdf<-data.frame(time=round(seq(0,10,0.01),3),narrs=200/2*(1-sin(2*pi*round(seq(0,10,0.01),3)))/100) %>%
      rowwise() %>%
      mutate(narrs_rounded=floor(narrs)+ifelse(runif(1)<(narrs-floor(narrs)),1,0)) %>%
      select(-narrs)
    log_new_sims<-tdf %>%
      uncount(narrs_rounded) %>%
      rename(arrs=time) %>%
      mutate(end=arrs+rexp(nrow(.),rate=2))
    #log_new_sims<-data.frame(arrs=round(seq(0.005,10.005,0.01),3)) %>%
    #  mutate(end=arrs+rexp(nrow(.),rate=2))
    logx_sims<-rbind(log_old_sims,log_new_sims)
    res_sims<-data.frame(xname=sapply(sq,function(x) length(which(logx_sims$arrs<x & logx_sims$end>x))))
    #plot(res_sims$xname,type="l")
    names(res_sims)<-paste0("rep_",rep)
    return(res_sims)
  }))
  stopCluster(cl)
  sim_res_tmp %>%
    mutate(time=sq) %>%
    pivot_longer(cols=-time,names_to="rep",values_to="occupancy",names_prefix="rep_") %>%
    group_by(time) %>%
    summarise(mean_occ=mean(occupancy),occ_l80=quantile(occupancy,0.1),occ_u80=quantile(occupancy,0.9),
              occ_l95=quantile(occupancy,0.025),occ_u95=quantile(occupancy,0.975))
}

sim_res_exp_formulae<-sim_res_exp("formulae") %>% mutate(init_method="Formulae")
sim_res_exp_method1<-sim_res_exp("method1") %>% mutate(init_method="Method 1")
sim_res_exp_method2<-sim_res_exp("method2") %>% mutate(init_method="Method 2")
sim_res_exp_method3<-sim_res_exp("method3") %>% mutate(init_method="Method 3")
sim_res_exp_all<-bind_rows(sim_res_exp_formulae,sim_res_exp_method1,sim_res_exp_method2,sim_res_exp_method3)

sim_res_exp_all %>%
  ggplot(aes(x=time)) +
  geom_ribbon(aes(ymin=occ_l95,ymax=occ_u95),fill="lightblue3",alpha=0.5) +
  geom_ribbon(aes(ymin=occ_l80,ymax=occ_u80),fill="blue",alpha=0.5) +
  geom_line(aes(y=mean_occ)) +
  facet_wrap(~init_method,nrow=1) +
  theme_bw() 

# lognormal (sigma = 1)

sim_res_lnorm<-function(init_method) {
  sqraw<-seq(0,10,0.1)
  sq<-round(sqraw,1)
  cl<-makeCluster(detectCores()-1)
  sim_res_tmp<-do.call("cbind",parLapply(cl,1:5000,function(rep) {
    require(Runuran)
    require(pracma)
    require(tidyverse)
    require(truncdist)
    set.seed(rep+ifelse(init_method=="formulae",1000000,ifelse(init_method=="method1",2000000,ifelse(init_method=="method2",3000000,4000000))))
    if(init_method=="formulae") {
      ftr<-function(x) (1+sin(2*pi*x))*(0.5-0.5*erf((0.5+log(2*x))/sqrt(2)))/
        integrate(function(y) (1+sin(2*pi*y))*(0.5-0.5*erf((0.5+log(2*y))/sqrt(2))),lower=0,upper=Inf)$value
      gen<-pinv.new(pdf=ftr,lb=0.000001,ub=10)
      dfx<-data.frame(arrs=-ur(gen,64))
    }
    if(init_method=="method1") dfx<-data.frame(arrs=rep(-0.01,64))
    if(init_method=="method2") dfx<-data.frame(arrs=-rlnorm(64,meanlog=(log(0.5))-0.5,sdlog=1))
    if(init_method=="method3") {
      tmp_mean<-0.5
      tmp_var<-(exp(1^2)-1)*exp(2*((log(0.5))-0.5)+1^2)
      dfx<-data.frame(arrs=rep(-0.5*(tmp_mean+tmp_var/tmp_mean),64))
    }
    log_old_sims<-dfx %>%
      rowwise() %>%
      mutate(end=rtrunc(1,"lnorm",a=-arrs,b=Inf,meanlog=(log(0.5))-0.5,sdlog=1)+arrs)
    tdf<-data.frame(time=round(seq(0,10,0.01),3),narrs=200/2*(1-sin(2*pi*round(seq(0,10,0.01),3)))/100) %>%
      rowwise() %>%
      mutate(narrs_rounded=floor(narrs)+ifelse(runif(1)<(narrs-floor(narrs)),1,0)) %>%
      select(-narrs)
    log_new_sims<-tdf %>%
      uncount(narrs_rounded) %>%
      rename(arrs=time) %>%
      mutate(end=arrs+rlnorm(nrow(.),meanlog=(log(0.5)-0.5),sdlog=1))
    #log_new_sims<-data.frame(arrs=round(seq(0.005,10.005,0.01),3)) %>%
    #  mutate(end=arrs+rexp(nrow(.),rate=2))
    logx_sims<-rbind(log_old_sims,log_new_sims)
    res_sims<-data.frame(xname=sapply(sq,function(x) length(which(logx_sims$arrs<x & logx_sims$end>x))))
    #plot(res_sims$xname,type="l")
    names(res_sims)<-paste0("rep_",rep)
    return(res_sims)
  }))
  stopCluster(cl)
  sim_res_tmp %>%
    mutate(time=sq) %>%
    pivot_longer(cols=-time,names_to="rep",values_to="occupancy",names_prefix="rep_") %>%
    group_by(time) %>%
    summarise(mean_occ=mean(occupancy),occ_l80=quantile(occupancy,0.1),occ_u80=quantile(occupancy,0.9),
              occ_l95=quantile(occupancy,0.025),occ_u95=quantile(occupancy,0.975))
}

sim_res_lnorm_formulae<-sim_res_lnorm("formulae") %>% mutate(init_method="Formulae")
sim_res_lnorm_method1<-sim_res_lnorm("method1") %>% mutate(init_method="Method 1")
sim_res_lnorm_method2<-sim_res_lnorm("method2") %>% mutate(init_method="Method 2")
sim_res_lnorm_method3<-sim_res_lnorm("method3") %>% mutate(init_method="Method 3")
sim_res_lnorm_all<-bind_rows(sim_res_lnorm_formulae,sim_res_lnorm_method1,sim_res_lnorm_method2,sim_res_lnorm_method3)

sim_res_lnorm_all %>%
  ggplot(aes(x=time)) +
  geom_ribbon(aes(ymin=occ_l95,ymax=occ_u95),fill="lightblue3",alpha=0.5) +
  geom_ribbon(aes(ymin=occ_l80,ymax=occ_u80),fill="blue",alpha=0.5) +
  geom_line(aes(y=mean_occ)) +
  facet_wrap(~init_method,nrow=1) +
  theme_bw() 

# lognormal (sigma = 1.25)

sim_res_lnorm2<-function(init_method) {
  sqraw<-seq(0,10,0.1)
  sq<-round(sqraw,1)
  cl<-makeCluster(detectCores()-1)
  sim_res_tmp<-do.call("cbind",parLapply(cl,1:5000,function(rep) {
    require(Runuran)
    require(pracma)
    require(tidyverse)
    require(truncdist)
    set.seed(rep+ifelse(init_method=="formulae",1000000,ifelse(init_method=="method1",2000000,ifelse(init_method=="method2",3000000,4000000))))
    if(init_method=="formulae") {
      ftr<-function(x) (1+sin(2*pi*x))*(0.5-0.5*erf((0.78125+log(2*x))/(1.25*sqrt(2))))/
        integrate(function(y) (1+sin(2*pi*y))*(0.5-0.5*erf((0.78125+log(2*y))/(1.25*sqrt(2)))),lower=0,upper=Inf)$value
      gen<-pinv.new(pdf=ftr,lb=0.000001,ub=10)
      dfx<-data.frame(arrs=-ur(gen,64))
    }
    if(init_method=="method1") dfx<-data.frame(arrs=rep(-0.01,64))
    if(init_method=="method2") dfx<-data.frame(arrs=-rlnorm(64,meanlog=(log(0.5))-0.78125,sdlog=1.25))
    if(init_method=="method3") {
      tmp_mean<-0.5
      tmp_var<-(exp(1.25^2)-1)*exp(2*((log(0.5))-0.78125)+1.25^2)
      dfx<-data.frame(arrs=rep(-0.5*(tmp_mean+tmp_var/tmp_mean),64))
    }
    log_old_sims<-dfx %>%
      rowwise() %>%
      mutate(end=rtrunc(1,"lnorm",a=-arrs,b=Inf,meanlog=(log(0.5))-0.78125,sdlog=1.25)+arrs)
    tdf<-data.frame(time=round(seq(0,10,0.01),3),narrs=200/2*(1-sin(2*pi*round(seq(0,10,0.01),3)))/100) %>%
      rowwise() %>%
      mutate(narrs_rounded=floor(narrs)+ifelse(runif(1)<(narrs-floor(narrs)),1,0)) %>%
      select(-narrs)
    log_new_sims<-tdf %>%
      uncount(narrs_rounded) %>%
      rename(arrs=time) %>%
      mutate(end=arrs+rlnorm(nrow(.),meanlog=(log(0.5)-0.78125),sdlog=1.25))
    #log_new_sims<-data.frame(arrs=round(seq(0.005,10.005,0.01),3)) %>%
    #  mutate(end=arrs+rexp(nrow(.),rate=2))
    logx_sims<-rbind(log_old_sims,log_new_sims)
    res_sims<-data.frame(xname=sapply(sq,function(x) length(which(logx_sims$arrs<x & logx_sims$end>x))))
    #plot(res_sims$xname,type="l")
    names(res_sims)<-paste0("rep_",rep)
    return(res_sims)
  }))
  stopCluster(cl)
  sim_res_tmp %>%
    mutate(time=sq) %>%
    pivot_longer(cols=-time,names_to="rep",values_to="occupancy",names_prefix="rep_") %>%
    group_by(time) %>%
    summarise(mean_occ=mean(occupancy),occ_l80=quantile(occupancy,0.1),occ_u80=quantile(occupancy,0.9),
              occ_l95=quantile(occupancy,0.025),occ_u95=quantile(occupancy,0.975))
}

sim_res_lnorm2_formulae<-sim_res_lnorm2("formulae") %>% mutate(init_method="Formulae")
sim_res_lnorm2_method1<-sim_res_lnorm2("method1") %>% mutate(init_method="Method 1")
sim_res_lnorm2_method2<-sim_res_lnorm2("method2") %>% mutate(init_method="Method 2")
sim_res_lnorm2_method3<-sim_res_lnorm2("method3") %>% mutate(init_method="Method 3")
sim_res_lnorm2_all<-bind_rows(sim_res_lnorm2_formulae,sim_res_lnorm2_method1,sim_res_lnorm2_method2,sim_res_lnorm2_method3)

sim_res_lnorm2_all %>%
  ggplot(aes(x=time)) +
  geom_ribbon(aes(ymin=occ_l95,ymax=occ_u95),fill="lightblue3",alpha=0.5) +
  geom_ribbon(aes(ymin=occ_l80,ymax=occ_u80),fill="blue",alpha=0.5) +
  geom_line(aes(y=mean_occ)) +
  facet_wrap(~init_method,nrow=1) +
  theme_bw() 


# combined plot

sim_res_exp_lnorm_lnorm2<-bind_rows(sim_res_exp_all %>% mutate(dist="exponential(2)"),
                             sim_res_lnorm_all %>% mutate(dist="lognormal(-1.19,1)"),
                             sim_res_lnorm2_all %>% mutate(dist="lognormal(-1.47,1.25)"))
saveRDS(sim_res_exp_lnorm_lnorm2,file="res_4_3.R")

#sim_res_exp_lnorm_lnorm2<-readRDS("res_4_3.R")

png("fig6.png",units="in",height=5,width=6,res=600)
print(
  sim_res_exp_lnorm_lnorm2 %>%
    #group_by(init_method,dist) %>%
    #mutate(mean_occ_sm=c(NA,NA,NA,NA,NA,NA,NA,rollapply(zoo(mean_occ),15,function(x) c(1,2,3,4,5,6,7,8,7,6,5,4,3,2,1) %*% x/64),NA,NA,NA,NA,NA,NA,NA)) %>%
    ggplot(aes(x=time)) +
    geom_ribbon(aes(ymin=occ_l95,ymax=occ_u95),fill="skyblue",alpha=0.3) +
    geom_ribbon(aes(ymin=occ_l80,ymax=occ_u80),fill="skyblue4",alpha=0.3) +
    geom_line(aes(y=mean_occ)) +
    #geom_line(aes(y=mean_occ_sm),col="red") +
    facet_grid(dist~init_method) +
    scale_x_continuous(breaks=c(0,2,4,6,8,10)) +
    xlab("Calendar time") +
    ylab("Occupancy") +
    theme_bw()
)
dev.off()

# errors for method 1 vs formulae
drops_method1<-sim_res_exp_lnorm_lnorm2 %>%
  filter(init_method=="Method 1") %>%
  filter(time>0 & time<1) %>%
  group_by(init_method,dist) %>%
  summarise(min=min(mean_occ),max=max(mean_occ))
drops_method1_exp<-drops_method1 %>% filter(dist=="exponential(2)") %>% .$min
drops_method1_ln1<-drops_method1 %>% filter(dist=="lognormal(-1.19,1)") %>% .$min
drops_method1_ln2<-drops_method1 %>% filter(dist=="lognormal(-1.47,1.25)") %>% .$min
drops_method1_exp-drops_method1_ln1
(drops_method1_exp-drops_method1_ln1)/drops_method1_exp
drops_method1_exp-drops_method1_ln2
(drops_method1_exp-drops_method1_ln2)/drops_method1_exp

# errors for method 3 vs formulae
drops_method3<-sim_res_exp_lnorm_lnorm2 %>%
  filter(init_method=="Method 3") %>%
  filter(time>0.5 & time<1.5) %>%
  group_by(init_method,dist) %>%
  summarise(min=min(mean_occ),max=max(mean_occ))
drops_method3_exp<-drops_method3 %>% filter(dist=="exponential(2)") %>% .$max
drops_method3_ln1<-drops_method3 %>% filter(dist=="lognormal(-1.19,1)") %>% .$max
drops_method3_ln2<-drops_method3 %>% filter(dist=="lognormal(-1.47,1.25)") %>% .$max
drops_method3_exp-drops_method3_ln1
(drops_method3_exp-drops_method3_ln1)/drops_method3_exp
drops_method3_exp-drops_method3_ln2
(drops_method3_exp-drops_method3_ln2)/drops_method3_exp


########################################
# SECTION 5.1

SF <- function(x) 0.5-0.5*erf((log(x)-3)/(1.5*sqrt(2)))

samples<-function(SF, N) {
  eqn2 <- function(x) {SF(x) / integrate(SF, lower = 0, upper = Inf)$value}
  gen <- pinv.new(pdf = eqn2, lb = 0, ub = Inf)
  ur(gen, N)
}

round(samples(SF, 10), 2)










