qmap("minnesota", zoom=6)















setwd("~/Documents/Projects/PolycultureBioenergy")
rm(list=ls()) #Cleans everything
source.with.encoding('~/Documents/Projects/HelperFunctions.R', encoding='UTF-8')
#This is modified from V5 following meeting with Adam on July 25, 2014
#One main difference is that many models do not include the data from 
#the first "production year" or year 2 of the planting because N was
#not applied during this year. We can also justify it by saying that the 
#plantings were still maturing.
library(nlme);library(car);library(faraway); library(ggplot2); library(lattice); library(reshape)
library(multcomp); library(vegan); library(lme4); library(AICcmodavg);
library(arm); library(coefplot); library(agricolae)
pdat<-read.table(file="PolyDat.txt", header=T)
pdat$Location<-as.factor(pdat$Location); pdat$fYear<-as.factor(pdat$Year);pdat$Plot<-as.factor(pdat$Plot);
pdat$Trt<-as.factor(pdat$Trt); pdat$Nfert<-as.factor(pdat$Nfert)
TrtNames<-c("1 spp: Switchgrass", "1 spp: Big Bluestem", "1 spp: Indiangrass", "1 spp: Canada Wild Rye",
            "4 spp: Grasses", "4 spp: Legumes", "4 spp: Forbs", "8 spp: Grass + Legumes",
            "8 spp: Grass + Forbs", "8 spp: Legumes + Forbs", "12 spp: Mix", "24 spp: Mix")
Trt_names<-list('1'="1 spp: Switchgrass", '2'="1 spp: Big Bluestem", '3'="1 spp: Indiangrass", '4'="1 spp: Canada Wild Rye",
                '5'="4 spp: Grasses", '6'="4 spp: Legumes", '7'="4 spp: Forbs", '8'="8 spp: Grass + Legumes",
                '9'="8 spp: Grass + Forbs", '10'="8 spp: Legumes + Forbs", '11'="12 spp: Mix", '12'="24 spp: Mix")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","gray30")
#Converting TotalBiomass from g / 0.5m2 to g / m2
pdat$TotalBiomass<-pdat$TotalBiomass*2
#This indicates the number of sown species for each treatment
sns<-c(1,1,1,1,4,4,4,8,8,8,12,24)
pdat$SR_P<-sns[as.numeric(pdat$Trt)]
#Calculates the number of non-NA and
#non-zeros in each row, be careful with indexing... only want to look through the species,
#so start with column 13 (Amorpha, first listed species) and go through 43 (zizia)
pdat$SR_A<-apply(pdat[,13:43],1, function(x) length(which(x>0)))
#Also, total cover does not always add up to weed cover plus all individual
#species, which is usually added in the field so this is understandable
#TotalCover2 is a sum of all the species and weed cover in the data set
#per row
pdat$TotalCover2<-rowSums(pdat[,c(8,13:43)], na.rm=T)
#Species with cover value of NA were not planted sown, but if the value is 0, it was
#sown but did not show up in the sample area. 
#Replicating pdat because we have to change NA to 0 for diversity
test2<-pdat  
#Make sure that you have 47 columns and that the last one is "TotalCover2"
colnames(test2)
#Now we need to divide the cover of each C4 plant by the total cover, and sum them together. 
#And this will be done for all other functional groups
pdat$C4_cover<-rowSums(cbind(test2[,14]/test2[,47],test2[,18]/test2[,47],test2[,33]/test2[,47],test2[,38]/test2[,47],test2[,41]/test2[,47]),na.rm=T)
pdat$C3_cover<-rowSums(cbind(test2[,22]/test2[,47],test2[,23]/test2[,47],test2[,24]/test2[,47]),na.rm=T)
pdat$Forb_cover<-rowSums(cbind(test2[,15]/test2[,47],test2[,21]/test2[,47],test2[,25]/test2[,47],test2[,27]/test2[,47],test2[,30]/test2[,47],test2[,32]/test2[,47],test2[,36]/test2[,47],test2[,37]/test2[,47],test2[,39]/test2[,47],test2[,40]/test2[,47],test2[,43]/test2[,47]),na.rm=T)
pdat$Legume_cover<-rowSums(cbind(test2[,13]/test2[,47],test2[,16]/test2[,47],test2[,17]/test2[,47],test2[,19]/test2[,47],test2[,20]/test2[,47],test2[,26]/test2[,47],test2[,28]/test2[,47],test2[,29]/test2[,47],test2[,31]/test2[,47],test2[,34]/test2[,47],test2[,35]/test2[,47],test2[,42]/test2[,47]),na.rm=T)
#We have to change NA for absent species to 0 to calculate diversity,
test2[is.na(test2)] <- 0
pdat$Div<-NULL
pdat$Div<-diversity(test2[,13:43])
#Adding the number of years since establishment. Note that this much be done for 
#southern and northern sites seperately
pdat$YrEst<-NA
pdat[as.character(pdat$Location)=="1",]$YrEst<-pdat[as.character(pdat$Location)=="1",]$Year-2006
pdat[as.character(pdat$Location)=="2",]$YrEst<-pdat[as.character(pdat$Location)=="2",]$Year-2007
pdat$Block<-as.factor(substr(as.numeric(as.character(pdat$Plot)), 2, 2))
#New to V5
pdat2<-droplevels(subset(pdat, pdat$YrEst>1))
#We should log transform the response. Need to clean data for that.
pdat2<-pdat[is.finite(log10(pdat$TotalBiomass)),]
pdat2$logbiomass<-log10(pdat2$TotalBiomass)
###################################################Data are ready!!!!!!!!!!!!!!!
##Reporting mean and SE of yield for certain treatments across sites and years
summarySE(subset(pdat, pdat$Trt=="1"|pdat$Trt=="2"|pdat$Trt=="5"|pdat$Trt=="8"|pdat$Trt=="11"|pdat$Trt=="12"),
                   measurevar="TotalBiomass", groupvars=c("Trt","Nfert"), na.rm=T)

tdat1<-subset(pdat, pdat$Trt=="1"|pdat$Trt=="2"|pdat$Trt=="5"|pdat$Trt=="8"|pdat$Trt=="11"|pdat$Trt=="12")
tdat2<-subset(tdat1, tdat1$Site=="R"&tdat1$Nfert=="0")
tdatmod<-lm(TotalBiomass~Trt, tdat2)
summary(tdatmod)
LSD.test(tdatmod, "Trt")
###
iter<-list(maxIter=500,msMaxIter=500, niterEM=500)

allmod<-lme(log10(TotalBiomass)~Trt*Nfert, random=~1|Site/Year, method="REML", data=pdat2, na.action=na.omit)
allmod2<-lme(log10(TotalBiomass)~Trt*Nfert*YrEst, random=~1+Nfert|Site/Year, method="REML", data=pdat2, na.action=na.omit)

#Effect of species treatments across N and time, accounting for
#variability in space
trtmod<-lme(log10(TotalBiomass)~Trt, random=~1|Site/Year, method="REML",
            data=subset(pdat2, pdat2$Nfert=="0"), na.action=na.omit)

#simout<-simulate(trtmodN, nsim=100)

simmod<-function(mod, nsim, dat) {
  
  parsout<-matrix(nrow=nsim, ncol=length(fixed.effects(mod)))
  colnames(parsout)<-names(fixed.effects(mod))
  
  for(i in 1:nsim) {
    modtmp<-update(mod, data=dat[sample(1:nrow(dat),nrow(dat), rep=T),])
    parsout[i,]<-fixed.effects(modtmp)
    
    print(i/nsim)
  }
  
  return(parsout)
}

parouttest<-simmod(trtmod, nsim=100, dat=subset(pdat2, pdat2$Nfert=="0"))
head(parouttest)
quantile(parouttest[,1]+parouttest[,2], c(0.025, 0.975))
quantile(parouttest[,1]+parouttest[,3], c(0.025, 0.975))

prs<-cbind(parouttest[,1], parouttest[,-1]+parouttest[,1])

plot(c(1,12), c(2.25, 3), type="n", xlab="trt", ylab="eff")
segments(1:12, apply(prs, 2, function(x) quantile(x, 0.025))
         , 1:12, apply(prs, 2, function(x) quantile(x, 0.975)))

#Attempting this same process above, only using the interaction
#model
allmod<-lme(log10(TotalBiomass)~Trt*Nfert, random=~1+Nfert|Site/Year, method="REML",
            data=pdat2, na.action=na.omit, control=lmeControl(opt='optim'))
poutall<-simmod(allmod, nsim=100, dat=pdat2)
head(poutall)

est1<-data.frame(poutall[,2:12]+poutall[,1])
est2<-data.frame(poutall[,13]+poutall[,14:24])
est3<-est1+est2

prsall<-cbind(poutall[,1], est1, poutall[,1]+poutall[,13], est3)
TrtNamesShort<-c("Switchgrass", "Big bluestem", "Indiangrass", "Canada wild rye",
                 "Grasses", "Legumes", "Forbs", "Grass/Legumes",
                 "Grass/Forbs", "Legumes/Forbs", "12 species mix", "High-diversity")

colnames(prsall)[1:12]<-TrtNamesShort
colnames(prsall)[13:24]<-paste(TrtNamesShort,rep("+N",12),sep="")
prstrans<-10^prsall
prstranslg<-melt(prstrans,variable_name="Treatment")
colnames(prstranslg)<-c("ID","OGname","Biomass")
prstranslg$Trt<-as.factor(rep(rep(TrtNamesShort, each=100), 2))
prstranslg$Nfert<-as.factor(c(rep("0", 1200),rep("1",1200)))
prssumx<-aggregate(Biomass~Nfert*Trt, data=prstranslg, mean)
prssumx2<-cbind(prssumx, aggregate(Biomass~Nfert*Trt, data=prstranslg, function(x) quantile(x, 0.025))[3],aggregate(Biomass~Nfert*Trt, data=prstranslg, function(x) quantile(x, 0.975))[3])
colnames(prssumx2)[4:5]<-c("Low","High")
#prssumx2$Nfert<-revalue(factor(prssumx2$Nfert), c("0"="0","1"="67"))
prssumx2$Trt<-factor(prssumx2$Trt, levels=c("Switchgrass", "Big bluestem", "Indiangrass", "Canada wild rye",
                                          "Grasses", "Legumes", "Forbs", "Grass/Legumes",
                                          "Grass/Forbs", "Legumes/Forbs", "12 species mix", "High-diversity"))
prssumx2$sig<-c("B", "b", "CD", "c", "F", "d", "E", "c", "DE", "bc", "A", "b", "B", "a", "B", "bc", "E", "d", "EF", "d", "D", "bc", "A", "a")

ggplot(prssumx2, aes(x=Trt, y=Biomass*0.01, fill=Nfert))+
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=Low*0.01, ymax=High*0.01),
                width=.2,position=position_dodge(.9))+
  geom_text(aes(label=sig, y=(High*0.01+.4)), size=4, position = position_dodge(width=1), show_guide=FALSE)+
  xlab(expression("Treatment " ))+
  ylab(expression("Biomass Yield " ~ (Mg ~ ha^{-1})))+
  scale_y_continuous(limits=c(0,8))+
  scale_fill_manual(values=cbPalette,
    labels=c(expression(paste("0 kg N ", ha^-1)), expression(paste("67 kg N ", ha^-1))))+
  guides(fill=guide_legend(title.theme=element_blank(), keywidth=.75, keyheight=.5))+
  theme(plot.title=element_text(size=12,face='bold', hjust=.02),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        legend.title=element_blank(),
        legend.position=c(.8,.9),
        legend.text=element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.x=element_text(size=12, angle=45, vjust=.5, color='black'),
        axis.text.y=element_text(size=12, color='black'),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank())

#Now we include year
###########
#Change in biomass yield through time. Not all linear trends.
iter2<-list(maxIter=10000,msMaxIter=10000, niterEM=10000)
modouttab<-matrix(data=NA, nrow=24, ncol=20)
modouttab<-data.frame(modouttab)
modlist<-vector("list",12)
for(j in c(1:12)){
  tt1<-droplevels(subset(pdat2, as.numeric(pdat2$Trt)==j))
  mmod4<-lme(logbiomass~Nfert*(YrEst+I(YrEst^2)+I(YrEst^3)), random=~1|Site/Plot,
             data=tt1, na.action=na.omit, control=iter2, method="ML")
fmod<-stepAIC(mmod4, direction="backward")
fmod<-update(fmod, .~., method="REML")
modlist[[j]]<-fmod
modouttab[(j*2)-1,1]<-j
modouttab[(j*2),1]<-j
modouttab[(j*2)-1,2:(length(fixef(fmod))+1)]<-names(fixef(fmod))
modouttab[(j*2),2:(length(fixef(fmod))+1)]<-round(fixef(fmod),3)
}

newdf<-data.frame(Nfert=factor(rep(c(0,1), each=7)), YrEst=1:7)
iter<-100
listres<-vector("list",12)

for(j in c(1:12)){
tt1<-droplevels(subset(pdat2, as.numeric(pdat2$Trt)==j))
datout<-matrix(ncol=nrow(newdf), nrow=iter, data=NA)
for(i in 1:iter) {
  subs<-sample(1:nrow(tt1), nrow(tt1), rep=T)
  newmod<-update(modlist[[j]], data=tt1[subs,])
  datout[i,]<-predict(newmod, newdata=newdf, level=0)
  print(i/iter)
}
listres[[j]]<-datout
}
allout<-data.frame(rbind(listres[[1]],listres[[2]],listres[[3]],listres[[4]],listres[[5]],
                         listres[[6]],listres[[7]],listres[[8]],listres[[9]],listres[[10]],
                         listres[[11]],listres[[12]]))
allout$Trt<-rep(1:12, each=100)
allout2<-melt(allout, id.vars="Trt" )
allout2$Nfert<-rep(c("0","1"), each=8400)
allout2$Year<-rep(1:7, each=1200)
sumallout<-aggregate(10^value~Trt*Nfert*Year, data=allout2, mean)
sumallout$UCI<-aggregate(10^value~Trt*Nfert*Year, data=allout2, function(x) quantile(x, 0.975))[,4]
sumallout$LCI<-aggregate(10^value~Trt*Nfert*Year, data=allout2, function(x) quantile(x, 0.025))[,4]
colnames(sumallout)[4]<-"value"
sumallout$Trt<-as.factor(sumallout$Trt)
sumallout$Trt2<-sumallout$Trt
levels(sumallout$Trt2)[levels(sumallout$Trt)==1:12]<-  c("Switchgrass", "Big bluestem", "Indiangrass", "Canada wild rye",
                                                         "Grasses", "Legumes", "Forbs", "Grass/Legumes",
                                                         "Grass/Forbs", "Legumes/Forbs", "12 species mix", "High-diversity")
#The data from one simulation have been saved below. It is stored in
#the Figures folder, nested within the AgronomyPaper folder.
#write.csv(sumallout,"Data4TrtTimePlot.csv")
#plotting 5.6 X 4.5 inches
ggplot(sumallout, aes(x=Year, y=value*0.01, color=Nfert, shape=Nfert, linetype=Nfert))+
  facet_wrap(~Trt2)+
  geom_line()+
  geom_point()+
  geom_pointrange(aes(ymax = UCI*0.01, ymin=LCI*0.01, color=Nfert), linetype=1)+
  xlab(expression("Year since establishment " ))+
  ylab(expression("Biomass Yield " ~ (Mg ~ ha^{-1})))+
  scale_y_continuous(limits=c(0,10))+
  scale_x_continuous(limits=c(1,7), breaks=c(1,3,5,7))+
  scale_color_manual(values=cbPalette,
                    labels=c(expression(paste("0 kg N ", ha^-1)), expression(paste("67 kg N ", ha^-1))))+
  scale_shape_manual(values=c(15,16),
                     labels=c(expression(paste("0 kg N ", ha^-1)), expression(paste("67 kg N ", ha^-1))))+
  scale_linetype_manual(values=c(1,2),
                        labels=c(expression(paste("0 kg N ", ha^-1)), expression(paste("67 kg N ", ha^-1))))+
  guides(color=guide_legend(title.theme=element_blank(), keywidth=.75, keyheight=.5))+
  theme(strip.text.y=element_text(size=10),
        strip.text.x=element_text(size=10),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA),
        legend.position=c(.88,.93),
        legend.key=element_rect(fill='white', colour='white'),
        legend.key.size=unit(.7, 'cm'),
        legend.text=element_text(size=9),
        axis.line = element_line(color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank())

#########################Percent Cover analysis
setwd("/Users/junge037/Documents/Projects/PolycultureBioenergy/AgronomyPaper/Figures/speciescomposition")
#tt<-pdat[,c(1:44,53)]
#tt<-tt[,c(1:6,45,7:44)]
pdatx<-pdat[,c(1:6,53,13:43)]
#Reporting average percent cover of species within treatments through time (years 1, 3 and 7)
realPCyr<-vector("list",3)
realPCtrt<-vector("list",8)
yrss<-c(1,3,7)
for(k in yrss){
  for(j in c(5:12)){
    tt1<-droplevels(subset(pdatx, pdatx$YrEst==k&as.numeric(pdatx$Trt)==j))
    tt2<-Filter(function(x)!all(is.na(x)),tt1)
    #tt2<-tt1[,-which(is.na(tt1[1,]))]
    tt3<-melt(tt2, id.vars=c(1:7), measure.vars=grep("_",colnames(tt2)),
              direction="long")
    tt4<-summarySE(tt3, measurevar="value",groupvars=c("Nfert", "variable"), na.rm=T)
    tt4$Trt<-(rep(tt3[1,5], length(tt4[,1])))
    tt4$Year<-(rep(tt3[1,7], length(tt4[,1])))
    #write.csv(tt4, file=paste("Revised-table for year",tt3[1,7],"and trt",tt3[1,5],".csv"))
    realPCtrt[[j-4]]<-tt4
  }
  realPCyr[[which(yrss==k)]]<-realPCtrt
}


#Looking at change in shannon diversity through time
cdiv_e<-subset(pdat[,c(1:12,45:53)], pdat$YrEst==1&as.numeric(pdat$Trt)>4)
cdiv_l<-subset(pdat[,c(1:12,45:53)], pdat$Year==2013&as.numeric(pdat$Trt)>4)
cdiv_d<-cdiv_l[,c(1:6,21)]
cdiv_d$divchange<-cdiv_l$Div-cdiv_e$Div
cdiv_s<-summarySE(cdiv_d, measurevar="divchange", groupvars=c("Trt", "Nfert"), na.rm=T)
div_modl<-lme(Div~Trt-1, random=~1|Location/Plot, data=cdiv_l)
div_mod<-lme(divchange~Trt/Nfert-1, random=~1|Plot, data=cdiv_d)
summary(div_mod)
divall<-simmod(div_mod, nsim=10000, dat=cdiv_d)
divall2<-simmod(div_modl, nsim=100, dat=cdiv_l)
head(divall)
divall2<-cbind(divall[,1:8], divall[,1:8]+divall[,9:16])
colnames(divall2)<-c(TrtNamesShort[5:12],paste(TrtNamesShort[5:12],"+N", sep=''))
divall3<-melt(divall2, variable_name="Treatment")
colnames(divall3)<-c("ID","OGname","Diversity")
divall3$Trt<-as.factor(rep(rep(TrtNamesShort[5:12], each=10000), 2))
divall3$Nfert<-as.factor(c(rep("0", 80000),rep("1",80000)))
divall4<-aggregate(Diversity~Nfert*Trt, data=divall3, mean)
divall5<-cbind(divall4, aggregate(Diversity~Nfert*Trt, data=divall3, function(x) quantile(x, 0.025))[3],aggregate(Diversity~Nfert*Trt, data=divall3, function(x) quantile(x, 0.975))[3])
colnames(divall5)[4:5]<-c("Low","High")
divall5$Trt2<-divall5$Trt
divall5$Trt2<-factor(divall5$Trt, levels=c("Grasses", "Legumes", "Forbs", "Grass/Legumes",
                                           "Grass/Forbs", "Legumes/Forbs", "12 species mix", "High-diversity"))
#Now plotting
ggplot(divall5, aes(x=Trt2, y=Diversity, fill=Nfert))+
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=Low, ymax=High),
                width=.2,position=position_dodge(.9))+
  xlab(expression("Treatment " ))+
  ylab("Change in Shannon diversity")+
  scale_fill_manual(values=cbPalette,
                    labels=c(expression(paste("0 kg N ", ha^-1)), expression(paste("67 kg N ", ha^-1))))+
  guides(fill=guide_legend(title.theme=element_blank(), keywidth=.75, keyheight=.5))+
  theme(plot.title=element_text(size=12,face='bold', hjust=.02),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        legend.title=element_blank(),
        legend.position=c(.8,.9),
        legend.text=element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.x=element_text(size=12, angle=45, vjust=.5, color='black'),
        axis.text.y=element_text(size=12, color='black'),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank())

#Model of change in relative abundance for each species within each treatment through time. 
hist(sqrt(subset(pdat, pdat$Andropogon_gerardii>2)$Andropogon_gerardii*0.01))
andro_5mod<-lm(sqrt(Andropogon_gerardii*0.01)~YrEst*Nfert, data=subset(pdat, pdat$Trt==5))
summary(lme(sqrt(Andropogon_gerardii*0.01)~YrEst*Nfert, random=~1|Site/Year, data=subset(pdat, pdat$Trt==5)))

zz1<-subset(pdat, pdat$Trt==5)
zz1$tab<-zz1[14]
summary(lme(sqrt(zz1[14]*0.01)~zz1[53]*zz1[6], random=~1|zz1[2]/zz1[3], data=zz1))
summary(lme(sqrt(colnames(zz1[14])*0.01)~YrEst*Nfert, random=~1|Site/Year, data=zz1))

plot(andro_5mod)
summary(andro_5mod)
xyplot(sqrt(Andropogon_gerardii)~YrEst|Nfert, group=Site, data=subset(pdat, pdat$Trt==5))

#Automate a process to calculate slope for all species*trt interactions through time
setwd("/Users/junge037/Documents/Projects/PolycultureBioenergy/AgronomyPaper/Figures/speciescompositionSLOPES")

pdatx<-pdat[,c(1:6,53,13:43)]

  for(j in c(5:12)){
    tt1<-droplevels(subset(pdatx, as.numeric(pdatx$Trt)==j))
    tt2<-Filter(function(x)!all(is.na(x)),tt1)
    numsp<-length(grep("_",colnames(tt2)))
    res<-data.frame("Trt"=rep(NA,numsp),"Species"=rep(NA,numsp),
                    "slope"=rep(NA,numsp), "lower"=rep(NA,numsp),
                    "upper"=rep(NA,numsp), "slope+N"=rep(NA,numsp),
                    "lower+N"=rep(NA,numsp), "upper+N"=rep(NA,numsp))
    
    for(k in 1:numsp){
    tt2$tab<-tt2[,7+k]
      #mod<-lm(sqrt(tab*0.01)~Nfert/YrEst-1, data=tt2)
      mod2<-lme(sqrt(tab)~Nfert/YrEst-1, random=~1|Site/Year,
                data=tt2, na.action=na.omit, control=iter2)
    ptt<-simmod(mod2, nsim=10, dat=tt2)  
    res[k,1]<-tt2[1,5]
    res[k,2]<-colnames(tt2[7+k])
    res[k,3]<-mean(ptt[,3], na.rm=T)
    res[k,4]<-quantile(ptt[,3], c(0.025, 0.975))[1]
    res[k,5]<-quantile(ptt[,3], c(0.025, 0.975))[2]
    res[k,6]<-mean(ptt[,4], na.rm=T)
    res[k,7]<-quantile(ptt[,4], c(0.025, 0.975))[1]
    res[k,8]<-quantile(ptt[,4], c(0.025, 0.975))[2]
    write.csv(res, file=paste("Update3-slopes for trt",tt2[1,5],".csv"))
    }
  }


#Selecting best fit model
setwd("/Users/junge037/Documents/Projects/PolycultureBioenergy/AgronomyPaper/Figures/speciescompchange2")
iter2<-list(maxIter=10000,msMaxIter=10000, niterEM=10000)

pdatx<-pdat[,c(1:6,53,13:43)]
outdat<-data.frame(Biomass=rep(NA, 14), Nfert=rep(c("0","1"),7), YrEst=rep(1:7, 2))
Trtmods<-vector("list",12)
Trtresponse<-vector("list",12)
iter<-10
newdat<-data.frame(Nfert=rep(c("0","1"),7), YrEst=rep(1:7, 2))
for(j in c(5:12)){
  tt1<-droplevels(subset(pdatx, as.numeric(pdatx$Trt)==j))
  tt2<-Filter(function(x)!all(is.na(x)),tt1)
  numsp<-length(grep("_",colnames(tt2)))
  splist<-vector("list", numsp)
  datout<-vector("list", numsp)
#need a way to exclude one species from Trt 12 (high-diversity),
#Liatris (k = 18) always produces an error and shuts everything down.
#Need a way to bypass this species without excluding it from beging so as to 
#preserve the order of species i.e. when k = 31, that's Zizia
  for(k in 1:numsp){
    tt2$tab<-tt2[,7+k]
    mmod4<-lme(sqrt(tab)~Nfert*(YrEst+I(YrEst^2)+I(YrEst^3)), random=~1|Site/Plot,
               data=tt2, na.action=na.omit, control=iter2, method="ML")
    fmod<-stepAIC(mmod4, direction="backward")
    fmod<-update(fmod, .~., method="REML")
    splist[[k]]<-fmod
    rdat<-matrix(ncol=nrow(newdat), nrow=iter, data=NA)
    sst<-data.frame(matrix(nrow=3, ncol=14))
    if(length(fixef(fmod))>1){
    for(i in 1:iter) {
      subs<-sample(1:nrow(tt2), nrow(tt2), rep=T)
      newmod<-update(fmod, data=tt2[subs,])
      rdat[i,]<-predict(newmod, newdata=newdf, level=0)
      print(i/iter)
    }
    for(hh in 1:14){
    #rdat[rdat<0]<-0
    sst[1,hh]<-mean(rdat[,hh], na.rm=T)^2
    sst[2,hh]<-quantile(rdat[,hh], 0.975, na.rm=T)^2
    sst[3,hh]<-quantile(rdat[,hh], 0.025, na.rm=T)^2
    }
    datout[[k]]<-sst
  }else{
    datout[[k]]<-paste("No effect")
      }
    }
  Trtmods[[j]]<-splist
  Trtresponse[[j]]<-datout
}

### Function to plot
plotspecies<-function(trt){
spresp<-merge_all(Trtresponse[[trt]])
spresp<-spresp[order(as.numeric(rownames(spresp))),]
tt1<-droplevels(subset(pdatx, as.numeric(pdatx$Trt)==trt))
tt2<-Filter(function(x)!all(is.na(x)),tt1)
specs<-colnames(tt2)[grep("_",colnames(tt2))]
#spresp$species<-rep(colnames(tt2)[grep("_",colnames(tt2))], each=3)
spresp2<-data.frame(t(spresp))
cnam<-rep(c("Cover","UCI","LCI"), length(specs))
cnam2<-rep(specs, each=3)
colnames(spresp2)<-paste(cnam,cnam2, sep='-')
spresp2$Nfert<-rep(c("0","1"), each=7)
spresp2$YrEst<-rep(c(1:7),2)
spresp3<-melt(spresp2, id.vars=c("Nfert", "YrEst"))
spresp3$response<-as.character(lapply(strsplit(as.character(spresp3$variable), split="-"), "[", 1))
spresp3$species<-as.character(lapply(strsplit(as.character(spresp3$variable), split="-"), "[", 2))
spresp4<-cast(spresp3[,-3], Nfert + YrEst + species~response)

ggplot(spresp4, aes(x=YrEst, y=Cover, color=species))+
  facet_wrap(~Nfert)+
  geom_line()+
  geom_point()+
  geom_pointrange(aes(ymax = UCI, ymin=LCI, color=species))+
  xlab(expression("Year since establishment " ))+
  ylab("Relative abundance (% cover)")+
  #scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(1,7))+
  #scale_color_manual(values=cbPalette)+                   
  guides(color=guide_legend(title.theme=element_blank(), keywidth=.75, keyheight=.5))+
  theme(plot.title=element_text(size=12,face='bold', hjust=.02),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        legend.key=element_rect(fill='white', colour='white'),
        legend.key.size=unit(.8, 'cm'),
        legend.text=element_text(size=10),
        axis.line = element_line(color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank())
}
 
plotspecies(12)
##Plotting spec. comp for grasses and grass/legume trts, and legumes
##Simply swap out the trt number in Trtresponse[[]] and two lines below
spresp<-merge_all(c(Trtresponse[[6]]))
spresp<-spresp[order(as.numeric(rownames(spresp))),]
tt1<-droplevels(subset(pdatx, as.numeric(pdatx$Trt)==6))
tt2<-Filter(function(x)!all(is.na(x)),tt1)
specs<-colnames(tt2)[grep("_",colnames(tt2))]
#spresp$species<-rep(colnames(tt2)[grep("_",colnames(tt2))], each=3)
spresp2<-data.frame(t(spresp))
cnam<-rep(c("Cover","UCI","LCI"), length(specs))
cnam2<-rep(specs, each=3)
colnames(spresp2)<-paste(cnam,cnam2, sep='-')
spresp2$Nfert<-rep(c("0","1"), each=7)
spresp2$YrEst<-rep(c(1:7),2)
spresp3<-melt(spresp2, id.vars=c("Nfert", "YrEst"))
spresp3$response<-as.character(lapply(strsplit(as.character(spresp3$variable), split="-"), "[", 1))
spresp3$species<-as.character(lapply(strsplit(as.character(spresp3$variable), split="-"), "[", 2))
spresp4<-cast(spresp3[,-3], Nfert + YrEst + species~response)

sprepsT5<-spresp4
sprepsT5$Trt<-"Grasses"
sprepsT8<-spresp4
sprepsT8$Trt<-"Grass/legumes"
sprepsT6<-spresp4
sprepsT6$Trt<-"Legumes"
sprepsALL<-rbind(sprepsT5,sprepsT8, sprepsT6)
splong<-unique(sprepsALL$species)
spshort<-c("A. gerardii", "E. canadensis", "P. virgatum","S. nutans", "A. canescens",
           "A. canadensis", "B. australis", "D. canadense", "P. purpurea")
sprepsALL$Nfert2<-factor(sprepsALL$Nfert, levels=c("0","1"),labels=c(expression(paste("0 kg N ", ha^-1)), expression(paste("67 kg N ", ha^-1))))
sprepsALL$species2<-factor(sprepsALL$species, levels=splong,labels=spshort)

ggplot(sprepsALL, aes(x=YrEst, y=Cover, color=species2, shape=species2, linetype=species2))+
  facet_grid(Trt~Nfert2, labeller = label_parsed)+
  geom_line()+
  geom_point()+
  geom_pointrange(aes(ymax = UCI, ymin=LCI, color=species2), linetype=1)+
  xlab(expression("Year since establishment " ))+
  ylab("Relative abundance (% cover)")+
  #scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(1,7))+
  scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999", "#009E73","gray30" , "#0072B2", "#D55E00", "#CC79A7","#F0E442"))+
  scale_shape_manual(values=c(3,15,0,2,1,17,16,4,18))+
  #scale_linetype_manual(values=c(1,2))+
  guides(color=guide_legend(ncol=2, title.theme=element_blank(), keywidth=.75, keyheight=.5))+
  theme(plot.title=element_text(size=12,face='bold', hjust=.02),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA),
        legend.position="top",
        legend.key=element_rect(fill='white', colour='white'),
        legend.key.size=unit(1, 'cm'),
        legend.text=element_text(face="italic", size=10),
        axis.line = element_line(color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank())

################
#################Special analysis of changes in sp. comp for trts 11 and 12
#################Have to alter the data set to get these to work. 
pdatx<-pdat[,c(1:6,53,13:43)]
pdatx$Liatris_aspera<-NA
outdat<-data.frame(Biomass=rep(NA, 14), Nfert=rep(c("0","1"),7), YrEst=rep(1:7, 2))
Trtmods<-vector("list",12)
Trtresponse<-vector("list",12)
iter<-10
newdat<-data.frame(Nfert=rep(c("0","1"),7), YrEst=rep(1:7, 2))
newdf<-data.frame(Nfert=factor(rep(c(0,1), each=7)), YrEst=1:7)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","gray30", "peachpuff")

for(j in c(11:12)){
  tt1<-droplevels(subset(pdatx, as.numeric(pdatx$Trt)==j))
  tt2<-Filter(function(x)!all(is.na(x)),tt1)
  numsp<-length(grep("_",colnames(tt2)))
  splist<-vector("list", numsp)
  datout<-vector("list", numsp)
  #need a way to exclude one species from Trt 12 (high-diversity),
  #Liatris (k = 18) always produces an error and shuts everything down.
  #Need a way to bypass this species without excluding it from beging so as to 
  #preserve the order of species i.e. when k = 31, that's Zizia
  for(k in 1:numsp){
    tt2$tab<-tt2[,7+k]
    mmod4<-lme(sqrt(tab)~Nfert*(YrEst+I(YrEst^2)+I(YrEst^3)), random=~1|Site/Plot,
               data=tt2, na.action=na.omit, control=iter2, method="ML")
    fmod<-stepAIC(mmod4, direction="backward")
    fmod<-update(fmod, .~., method="REML")
    splist[[k]]<-fmod
    rdat<-matrix(ncol=nrow(newdat), nrow=iter, data=NA)
    sst<-data.frame(matrix(nrow=3, ncol=14))
    if(length(fixef(fmod))>1){
      for(i in 1:iter) {
        subs<-sample(1:nrow(tt2), nrow(tt2), rep=T)
        newmod<-update(fmod, data=tt2[subs,])
        rdat[i,]<-predict(newmod, newdata=newdf, level=0)
        print(i/iter)
      }
      for(hh in 1:14){
        #rdat[rdat<0]<-0
        sst[1,hh]<-mean(rdat[,hh], na.rm=T)^2
        sst[2,hh]<-quantile(rdat[,hh], 0.975, na.rm=T)^2
        sst[3,hh]<-quantile(rdat[,hh], 0.025, na.rm=T)^2
      }
      datout[[k]]<-sst
    }else{
      datout[[k]]<-paste("No effect")
    }
  }
  Trtmods[[j]]<-splist
  Trtresponse[[j]]<-datout
}

Species number 5 has no effect
Trtresponse[[12]][[5]]<- data.frame(X1=rep(NA,3),X2=rep(NA,3),X3=rep(NA,3),X4=rep(NA,3),X5=rep(NA,3),X6=rep(NA,3),X7=rep(NA,3),X7=rep(NA,3),
                                    X9=rep(NA,3),X10=rep(NA,3),X11=rep(NA,3),X12=rep(NA,3),X13=rep(NA,3),X14=rep(NA,3))

spresp<-merge_all(Trtresponse[[12]])
spresp<-spresp[order(as.numeric(rownames(spresp))),]
tt1<-droplevels(subset(pdatx, as.numeric(pdatx$Trt)==12))
tt2<-Filter(function(x)!all(is.na(x)),tt1)
specs<-colnames(tt2)[grep("_",colnames(tt2))]
#spresp$species<-rep(colnames(tt2)[grep("_",colnames(tt2))], each=3)
#specs<-colnames(tt1)[grep("_",colnames(tt1))]
spresp2<-data.frame(t(spresp))
cnam<-rep(c("Cover","UCI","LCI"), length(specs))
cnam2<-rep(specs, each=3)
colnames(spresp2)<-paste(cnam,cnam2, sep='-')
spresp2$Nfert<-rep(c("0","1"), each=7)
spresp2$YrEst<-rep(c(1:7),2)
spresp3<-melt(spresp2, id.vars=c("Nfert", "YrEst"))
spresp3$response<-as.character(lapply(strsplit(as.character(spresp3$variable), split="-"), "[", 1))
spresp3$species<-as.character(lapply(strsplit(as.character(spresp3$variable), split="-"), "[", 2))
spresp4<-cast(spresp3[,-3], Nfert + YrEst + species~response)
#Subsetting spresp4 to only include species that are estimated to have relative
#abundance of 1% or greater during any Nfert/year combination.
spresp4$species<-as.factor(spresp4$species)
raresp<-levels(spresp4$species)[-c(1,2,3,5,6,7,10,11,12,13,16,17,18,21,22,25,27,26,28,29)]
spresp5<-droplevels(spresp4[spresp4$species %in% raresp,])
spresp5$Nfert<-as.factor(spresp5$Nfert)
splong<-unique(spresp5$species)
spshort<-c("A. canadensis", "D. canadense", "E. pallida", "G. lepidota", 
           "H. maximilianii", "M. fistulosa", "P. virgatum", "R. pinnata",
           "R. hirta", "Z. aurea")
spresp5$Nfert2<-factor(spresp5$Nfert, levels=c("0","1"), labels=c(expression(paste("0 kg N ", ha^-1)), expression(paste("67 kg N ", ha^-1))))
spresp5$species2<-factor(spresp5$species, levels=splong,labels=spshort)

sp24<-ggplot(spresp5, aes(x=YrEst, y=Cover, color=species2))+
  facet_grid(~Nfert2, labeller = label_parsed)+
  geom_line()+
  geom_point()+
  geom_pointrange(aes(ymax = UCI, ymin=LCI, color=species2))+
  xlab(expression("Year since establishment " ))+
  ylab("Relative abundance (% cover)")+
  scale_y_continuous(limits=c(0,40))+
  scale_x_continuous(limits=c(1,7))+
  scale_color_manual(values=cbPalette)+                   
  guides(color=guide_legend(title.theme=element_blank(), keywidth=.75, keyheight=.5))+
  ggtitle("(a)")+
  theme(plot.title=element_text(size=12,face='bold', hjust=.02),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        #legend.position="top",
        legend.key=element_rect(fill='white', colour='white'),
        legend.key.size=unit(.8, 'cm'),
        legend.text=element_text(face="italic", size=10),
        axis.line = element_line(color='black'),
        #axis.text.x=element_text(size=12, color='black'),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=12, color='black'),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank())

##################################
### same as above but fore 12 species mixture

spresp<-merge_all(Trtresponse[[11]])
spresp<-spresp[order(as.numeric(rownames(spresp))),]
tt1<-droplevels(subset(pdatx, as.numeric(pdatx$Trt)==11))
tt2<-Filter(function(x)!all(is.na(x)),tt1)
specs<-colnames(tt2)[grep("_",colnames(tt2))]
#spresp$species<-rep(colnames(tt2)[grep("_",colnames(tt2))], each=3)
#specs<-colnames(tt1)[grep("_",colnames(tt1))]
spresp2<-data.frame(t(spresp))
cnam<-rep(c("Cover","UCI","LCI"), length(specs))
cnam2<-rep(specs, each=3)
colnames(spresp2)<-paste(cnam,cnam2, sep='-')
spresp2$Nfert<-rep(c("0","1"), each=7)
spresp2$YrEst<-rep(c(1:7),2)
spresp3<-melt(spresp2, id.vars=c("Nfert", "YrEst"))
spresp3$response<-as.character(lapply(strsplit(as.character(spresp3$variable), split="-"), "[", 1))
spresp3$species<-as.character(lapply(strsplit(as.character(spresp3$variable), split="-"), "[", 2))
spresp4<-cast(spresp3[,-3], Nfert + YrEst + species~response)
#Subsetting spresp4 to only include species that are estimated to have relative
#abundance of 1% or greater during any Nfert/year combination.
spresp4$species<-as.factor(spresp4$species)
raresp<- levels(spresp4$species)[-c(1,3,11,14)]
spresp6<-droplevels(spresp4[spresp4$species %in% raresp,])

spresp6$Nfert<-as.factor(spresp6$Nfert)
splong2<-unique(spresp6$species)
spshort2<-c("A. gerardii", "A. canadensis", "B. australis", "D. canadense", "E. pallida",  
            "E. canadensis",  "H. maximilianii", "P. virgatum", "R. pinnata",
            "S. rigida")
spresp6$Nfert2<-factor(spresp6$Nfert, levels=c("0","1"), labels=c(expression(paste("0 kg N ", ha^-1)), expression(paste("67 kg N ", ha^-1))))
spresp6$species2<-factor(spresp6$species, levels=splong2,labels=spshort2)

sp12<-ggplot(spresp6, aes(x=YrEst, y=Cover, color=species2))+
  facet_grid(~Nfert2, labeller=label_parsed)+
  geom_line()+
  geom_point()+
  geom_pointrange(aes(ymax = UCI, ymin=LCI, color=species2))+
  xlab(expression("Year since establishment " ))+
  ylab("Relative abundance (% cover)")+
  #scale_y_continuous(limits=c(0,40))+
  scale_x_continuous(limits=c(1,7))+
  scale_color_manual(values=cbPalette)+   
  ggtitle("(b)")+
  guides(color=guide_legend(title.theme=element_blank(), keywidth=.75, keyheight=.5))+
  theme(plot.title=element_text(size=12,face='bold', hjust=.02),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        #legend.position="top",
        legend.text=element_text(face="italic", size=10),
        legend.key=element_rect(fill='white', colour='white'),
        legend.key.size=unit(.8, 'cm'),
        axis.line = element_line(color='black'),
        axis.text.x=element_text(size=12, color='black'),
        axis.text.y=element_text(size=12, color='black'),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank())
multiplot(sp24,sp12)

###
srmod<-lme(log10(TotalBiomass)~Nfert*YrEst*SR_P, random=~1+Nfert|Site/Year, method="REML",
            data=pdat2, na.action=na.omit, control=lmeControl(opt='optim'))

srmod2<-lme(log10(TotalBiomass)~Nfert*YrEst*SR_P, random=~1|Trt/Site/Year, method="REML",
           data=pdat2, na.action=na.omit, control=lmeControl(opt='optim'))

divmod<-lme(log10(TotalBiomass)~Nfert*YrEst*Div, random=~1+Nfert|Site/Year, method="REML",
            data=pdat2, na.action=na.omit, control=lmeControl(opt='optim'))

summary(lme(log10(TotalBiomass)~Nfert*SR_P, random=~1+Nfert|Site, method="REML",
           data=subset(pdat2, pdat2$YrEst==7),
           na.action=na.omit, control=lmeControl(opt='optim')))
pyldm1<-lme(logbiomass~Nfert*(YrEst+I(YrEst^2)+I(YrEst^3)), random=~1|Site/Plot,
            data=droplevels(subset(pdat2, as.numeric(pdat2$Trt)==j)),
            na.action=na.omit, control=iter2, method="ML")
pyldm2<-lme(logbiomass~Nfert*(YrEst+I(YrEst^2)), random=~1|Site/Plot,
            data=droplevels(subset(pdat2, as.numeric(pdat2$Trt)==j)),
            na.action=na.omit, control=iter2, method="ML")
pyldm3<-lme(logbiomass~Nfert*YrEst, random=~1|Site/Plot,
            data=droplevels(subset(pdat2, as.numeric(pdat2$Trt)==j)),
            na.action=na.omit, control=iter2, method="ML")
pyldm4<-lme(logbiomass~1, random=~1|Site/Plot,
            data=droplevels(subset(pdat2, as.numeric(pdat2$Trt)==j)),
            na.action=na.omit, control=iter2, method="ML")

AIC(pyldm1, pyldm2, pyldm3, pyldm4)
anova(pyldm3, pyldm4)
