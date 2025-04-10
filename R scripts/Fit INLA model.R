library(sp)
library(INLA)
###****************MAKE SEPARATE DATA SETS FOR EACH VILLAGE*******************************

bana_dat<-psc_dat1[which(psc_dat1[,"village"]=="Bana village"),]
bama_dat<-psc_dat1[which(psc_dat1[,"village"]=="Bana market"),]
pala_dat<-psc_dat1[which(psc_dat1[,"village"]=="Pala"),]
sour_dat<-psc_dat1[which(psc_dat1[,"village"]=="Souroukoudingan"),]

###****************SET UP INLA MODEL****************
#Data processing to remove NAs and extract geographic coordinates
bana_dat_inla <- na.omit(bana_dat[,c("count.f","month","longitude","latitude","site.id","year","village")])
bana_dat_inla$month <-as.numeric(bana_dat_inla$month)
bana_dat_inla$village <- as.character(bana_dat_inla$village)
long_lat_bana<-cbind(bana_dat_inla[,"longitude"],bana_dat_inla[,"latitude"])

pala_dat_inla <- na.omit(pala_dat[,c("count.f","month","longitude","latitude","site.id","year","village")])
pala_dat_inla$month <-as.numeric(pala_dat_inla$month)
pala_dat_inla$village <- as.character(pala_dat_inla$village)
long_lat_pala<-cbind(pala_dat_inla[,"longitude"],pala_dat_inla[,"latitude"])

bama_dat_inla <- na.omit(bama_dat[,c("count.f","month","longitude","latitude","site.id","year","village")])
bama_dat_inla$month <-as.numeric(bama_dat_inla$month)
bama_dat_inla$village <- as.character(bama_dat_inla$village)
long_lat_bama<-cbind(bama_dat_inla[,"longitude"],bama_dat_inla[,"latitude"])

sour_dat_inla <- na.omit(sour_dat[,c("count.f","month","longitude","latitude","site.id","year","village")])
sour_dat_inla$month <-as.numeric(sour_dat_inla$month)
sour_dat_inla$village <- as.character(sour_dat_inla$village)
long_lat_sour<-cbind(sour_dat_inla[,"longitude"],sour_dat_inla[,"latitude"])

long_lat_BMBS<-rbind(long_lat_bana,long_lat_pala,long_lat_bama,long_lat_sour)

#Define the INLA mesh 
m1.cutoff<-0.0004
m1.min.angle<-c(26,25)
m1.max.edge<-c(0.004,0.05)
tmesh.st<-1
tmesh.end<-12
tmesh.by<-1


bdy_bana<-inla.nonconvex.hull(rbind(long_lat_bana,long_lat_bama),convex=0.01)
#Combine bana and bama as they are close together
bdy_pala<-inla.nonconvex.hull(long_lat_pala,convex=0.01)
bdy_sour<-inla.nonconvex.hull(long_lat_sour,convex=0.01)

sp_bana=Polygon(bdy_bana$loc[,1:2],hole=FALSE)
sp_pala=Polygon(bdy_pala$loc[,1:2],hole=FALSE)
sp_sour=Polygon(bdy_sour$loc[,1:2],hole=FALSE)
sp_all <- SpatialPolygons(list(Polygons(list(sp_bana,sp_pala,sp_sour), '0')))
sp_inla_ms<-as.inla.mesh.segment(sp_all)
BMBS_mesh = inla.mesh.2d(rbind(long_lat_bana,long_lat_bama,long_lat_pala,long_lat_sour),boundary=sp_inla_ms,
                         cutoff=m1.cutoff,
                         min.angle=m1.min.angle,
                         max.edge=m1.max.edge,
                         max.n.strict = 10000,max.n=10000)
#Plot the mesh
dev.new(noRStudioGD = T)
plot(BMBS_mesh)

#Create the SPDE model using penalized complexity priors
BMBS_spde<- inla.spde2.pcmatern(mesh=BMBS_mesh, alpha = 2,
                                prior.range = c(0.001,0.1), 
                                prior.sigma = c(2,0.1))  

#Make the time mesh
mesh1d=inla.mesh.1d(seq(tmesh.st,tmesh.end,by=tmesh.by),interval=c(tmesh.st,tmesh.end),degree=2, boundary='free')

#Format the data for input into the INLA model
BMBS_dat_mod<-data.frame(y=c(bana_dat_inla$count.f,pala_dat_inla$count.f,bama_dat_inla$count.f,sour_dat_inla$count.f), 
                         w=rep(1,nrow(bana_dat_inla)+nrow(pala_dat_inla)+nrow(bama_dat_inla)+nrow(sour_dat_inla)), 
                         month = c(bana_dat_inla$month,pala_dat_inla$month,bama_dat_inla$month,sour_dat_inla$month), xcoo=long_lat_BMBS[,1],
                         ycoo=long_lat_BMBS[,2],site=c(bana_dat_inla$site.id,pala_dat_inla$site.id,bama_dat_inla$site.id,sour_dat_inla$site.id),
                         village=c(bana_dat_inla$village,pala_dat_inla$village,bama_dat_inla$village,sour_dat_inla$village))

#Make an index for the spatial Gaussian Markov random effect parameters
BMBS_iset1 <- inla.spde.make.index('BMBS.field', n.spde=BMBS_mesh$n, n.group=mesh1d$m)

#Make the projector matrix
BMBS_A1 <- inla.spde.make.A(BMBS_mesh,loc=cbind(BMBS_dat_mod$xcoo, BMBS_dat_mod$ycoo),group=c(BMBS_dat_mod$month),group.mesh=mesh1d)

#Combine everything into an INLA stack object
BMBS_stk1 <- inla.stack(tag='BMBS.data', data=list(y=BMBS_dat_mod$y), 
                        A=list(BMBS_A1,1,1,1),effects=list(BMBS_iset1, b0.1=BMBS_dat_mod$w,
                                                           site=BMBS_dat_mod$site, 
                                                           village=BMBS_dat_mod$village)) 

#Set prior for the iid random effect
eprec_iid<- list(prec=list(prior='pc.prec', param=c(0.09, 0.0001)))
eprec_iid2 <- list(prec=list(prior='pc.prec', param=c(2, 0.1)))

#Specify the model formula
poiss_form <- y ~ 0 + f(BMBS.field, model=BMBS_spde, group=BMBS.field.group,control.group=list(model='ar1')) +
  f(site, model='iid',hyper = eprec_iid) + f(village, model='iid',hyper = eprec_iid2)

###****************RUN THE INLA MODEL****************
poiss_PSC.res <- inla(poiss_form, family=c('Poisson'), data=inla.stack.data(BMBS_stk1), control.predictor=list(compute=TRUE,A=inla.stack.A(BMBS_stk1)),control.compute=list(config=TRUE,cpo=TRUE),control.inla=list(int.strategy='eb'))
