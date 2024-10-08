#N-fixing project: Script version October 2024
#Author Linda Rosa Mueller


# Load the necessary libraries
library(dplyr);library(segmented);library(nlme);library(mgcv);library(gridExtra);library(betareg);
library(MASS);library(lme4);library(lmerTest);library(lsmeans);library(ggeffects);library(spdep);
library(ggplot2);library(effects);library(ncf);library(ape);library(sjPlot);library(gridExtra);
library(MuMIn);library(tidyverse);library(car); library(V8); library(rsq); library(DHARMa); library(grateful)

source("Master_modelling.R")  #change to functions?? 
options(na.action = "na.fail") #Change na. action


####DATASET####
#Load dataset from data cleaning script
gdat.ref <- readRDS("data/fulldata_for_analysis_2023.rds") %>%
  dplyr::select(c("entity_ID","entity_class","nfix","nfixno", "latitude","abs.lat","longitude", # Select columns
                  "landtype","status","presence","area","elev_range","precipitation","dist", "temperature","urbanland")) 

#### M1 BROAD PRESENCE MODEL####

#dataset for broad presence model including mainlands
gdat.ml.pres <- gdat.ref %>%
  dplyr::select(c("entity_ID","nfix","nfixno", "latitude","abs.lat","longitude", # Select columns
                  "landtype","status","presence","area","elev_range","precipitation","temperature"))%>%
  filter(!landtype=="other_island")%>%                                          #remove islands with undetermined type
  mutate(area = as.vector(log10((area)+.01)))%>%                                #log transform area
  drop_na()%>%
  mutate_at(c("abs.lat","area","elev_range","precipitation", "temperature"), scale) #scale all explanatory variables

test <- gdat.ml.pres%>%filter(landtype=="nonoceanic")%>%filter(status=="naturalized")%>%drop_na()

#### M1 Broad Presence (across islands and mainlands including native and naturalized floras)

broad.pres.model.full<- glm(presence~landtype*status+abs.lat+area+elev_range+precipitation+temperature, data =gdat.ml.pres, family= binomial(link ="logit"))
summary(broad.pres.model.full)

ref.land.stat <- lsmeans(broad.pres.model.full,pairwise~landtype*status, data= gdat.ml.pres, type="response")
ref.table.land.stat <- as.data.frame(ref.land.stat$lsmeans) 

#model including spatial correlation variable
rac <- Spat.cor.rep(broad.pres.model.full,gdat.ml.pres,2000)
broad.pres.model.rac.full <- glm(presence~landtype*status+abs.lat+area+elev_range+precipitation+temperature+rac, data = gdat.ml.pres, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(broad.pres.model.rac.full)

#M1 with only selected variables after stepwise regression
broad.pres.model<- glm(presence~landtype*status+abs.lat+area+elev_range, data =gdat.ml.pres, family= binomial(link ="logit"))
summary(broad.pres.model)

ref.land.stat <- lsmeans(broad.pres.model,pairwise~landtype*status, data= gdat.ml.pres, type="response")
ref.table.land.stat <- as.data.frame(ref.land.stat$lsmeans) 


##### Correlogram to test distance of spatial autocorrelation
mod <- broad.pres.model
dat<- gdat.ml.pres

#calculate residuals
resids.calc <- resid(mod)
sp.dat <- dat %>%
  mutate(resid = resids.calc)

#generate sp object
sp <- ncf::spline.correlog(x = as.numeric(sp.dat$latitude),
                           y = as.numeric(sp.dat$longitude),
                           z = as.numeric(sp.dat$resid),
                           xmax = 2000, resamp = 100, latlon=TRUE)
#save figure
png("figures/broadpresencecorrelogram.jpg", width = 10, height = 10, units = 'in', res = 300)
plot(sp)
dev.off()


#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(broad.pres.model,gdat.ml.pres,2000)
broad.pres.model.rac <- glm(presence~landtype*status+abs.lat+area+elev_range+rac, data = gdat.ml.pres, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(broad.pres.model.rac)

#check variance inflation factor (should be below 5 for all variables)
vif(broad.pres.model.rac)

#add rac variable to df to include in lsmeans.
gdat.ml.pres$rac <- rac

ref.land.stat.rac<-lsmeans(broad.pres.model.rac,pairwise~landtype*status, data = gdat.ml.pres, type="response")
ref.table.land.stat.rac<-as.data.frame(ref.land.stat.rac$lsmeans) 
ref.table.land.stat.rac.ntz <- ref.table.land.stat.rac%>%filter(!status=="native") #subset lsmeans to only include naturalized for plot 
ref.table.land.stat.rac.ntv <- ref.table.land.stat.rac%>%filter(!status=="naturalized") #subset lsmeans to only include native for plot

#contrasts
means <- emmeans(broad.pres.model.rac, ~landtype*status)
#look at means order to determine how to write contrasts:
means

#write contrasts:
contrasts <- list("mainland v all islands native" = c(-2,1,1,0,0,0),
                  "mainland v nonoceanic islands native" = c(-1,1,0,0,0,0),
                  "mainland v oceanic islands native" = c(-1,0,1,0,0,0),
                  "nonoceanic v oceanic islands native" = c(0,-1,1,0,0,0),
                  "mainland v all islands naturalized" = c(0,0,0,-2,1,1),
                  "mainland v oceanic islands naturalized" = c(0,0,0,-1,0,1),
                  "mainland v nonoceanic islands naturalized" = c(0,0,0,-1,1,0),
                  "nonoceanic v oceanic islands naturalized" = c(0,0,0,0,-1,1),
                  # is first contrast different native v non-native?
                  "mainland v all islands * status" = c(-2,1,1,2,-1,-1),
                  "mainland v oceanic islands * status" = c(-1,0,1,1,0,-1),
                  "mainland v nonoceanic islands * status" =  c(-1,1,0,1,-1,0),
                  "nonoceanic v oceanic island * status" = c(0,-1,1,0,1,-1),
                  "native v naturalized only on mainlands" =  c(-1,0,0,1,0,0),
                  "native v naturalized only on oceanic islands" =  c(0,0,-1,0,0,1),
                  "native v naturalized only on nonoceanic islands" =  c(0,-1,0,0,1,0),
                  "native v naturalized on both island types" =  c(0,-1,-1,0,1,1)
)

#extract results to df:
results <- lsmeans::contrast(means,contrasts)
results.df <- as.data.frame(results)
results.df


#subset data to only include naturalized for plot
gdat.ml.ntz<- gdat.ml.pres%>%filter(status=="naturalized")

#Create a custom color scale
colScale <- scale_colour_manual(values=c("coral"))
fillScale <- scale_fill_manual(values=c("coral"))

#plot the model results for presence of naturalized flora across islands and mainlands
broad.pres.ntz<- 
  ggplot(ref.table.land.stat.rac.ntz, aes(x = landtype, y = lsmean, fill = status, color = status))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) + 
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3)+
  geom_point(data = gdat.ml.ntz,aes(x=landtype,y=presence,color=status),position=position_jitterdodge(dodge.width = 1),size=2,alpha=0.3)+
  coord_cartesian(ylim=c(0,1))+
  colScale+
  fillScale+
  labs(y="Probability of N-fixing Plants")+
  xlab(" ")+
  theme_classic(base_size = 30) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  theme(legend.position = "top")+
  theme(legend.key = element_rect(colour =c("coral") ))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
  guides(color="none", fill = guide_legend(title="Floral status:",override.aes = list(shape = NA)))+
  theme(axis.text.y=element_text(size=50))+
  theme(axis.text.x=element_text(size=50))+
  ylim(0,1)

broad.pres.ntz

#Saving the plot as a png
png("figures/naturalized_broad_presence_mainland.jpg", width=10, height= 12, units='in', res=300)
broad.pres.ntz
dev.off()

#subset data to only include native for plot
gdat.ml.ntv<- gdat.ml.pres%>%filter(status=="native")

#Create a custom color scale
colScale <- scale_colour_manual(values=c("darkseagreen3"))
fillScale <- scale_fill_manual(values=c("darkseagreen3"))

#plot the model results for presence of native flora across islands and mainlands
broad.pres.ntv<- 
  ggplot(ref.table.land.stat.rac.ntv, aes(x = landtype, y = lsmean, fill = status, color = status))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) + 
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3)+
  geom_point(data = gdat.ml.ntv,aes(x=landtype,y=presence,color=status),position=position_jitterdodge(dodge.width = 1),size=2,alpha=0.3)+
  coord_cartesian(ylim=c(0,1))+
  colScale+
  fillScale+
  labs(y="Probability of N-fixing Plants")+
  xlab(" ")+
  theme_classic(base_size = 30) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  theme(legend.position = "top")+
  theme(legend.key = element_rect(colour =c("darkseagreen3") ))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
  guides(color="none", fill = guide_legend(title="Floral status:",override.aes = list(shape = NA)))+
  theme(axis.text.y=element_text(size=50))+
  theme(axis.text.x=element_text(size=50))+
  theme(axis.title.y=element_text(size=50))

broad.pres.ntv

#Saving the plot as a png
png("figures/native_broad_presence_mainland.jpg", width=10, height= 12, units='in', res=300)
broad.pres.ntv
dev.off()

####M2 BROAD PROPORTION MODEL####

#dataset for broad presence model including mainlands
gdat.ml.prop <- gdat.ref %>%
  dplyr::select(c("entity_ID","nfix","nfixno", "latitude","abs.lat","longitude", # Select columns
                  "landtype","status","presence","area","elev_range","precipitation", "temperature"))%>%
  filter(!landtype=="other_island")%>%                                          #remove islands with undetermined type
  mutate(area = as.vector(log10((area)+.01)))%>%                                #log transform area
  #find out rows that have no species counts
  mutate(species = nfix + nfixno) %>%
  filter(species > 0) %>%
  drop_na() %>%
  mutate_at(c("abs.lat","area","elev_range","precipitation","temperature"), scale) %>%  #scale all explanatory variables
  #remove outliers
  filter(precipitation < 3.5)

#### M2 Model Broad Proportion across islands and mainlands including native and naturalized flora 

broad.prop.model<- glmer(cbind(nfix,nfixno)~landtype*status+abs.lat+area+elev_range+precipitation + temperature+ (1|entity_ID), data =gdat.ml.prop, family= binomial(link ="logit"))
summary(broad.prop.model)

ref.land.stat <- lsmeans(broad.prop.model,pairwise~landtype*status, data= gdat.ml.prop, type="response")
ref.table.land.stat <- as.data.frame(ref.land.stat$lsmeans) 


##### Correlogram to test distance of spatial autocorrelation
mod <- broad.prop.model
dat<- gdat.ml.prop

#calc residuals
resids.calc <- resid(mod)
sp.dat <- dat %>%
  mutate(resid = resids.calc)

#generate sp object
sp <- ncf::spline.correlog(x = as.numeric(sp.dat$latitude),
                           y = as.numeric(sp.dat$longitude),
                           z = as.numeric(sp.dat$resid),
                           xmax = 4000, resamp = 100, latlon=TRUE)
#save figure
png("figures/broadpropcorrelogram.jpg", width = 10, height = 10, units = 'in', res = 300)
plot(sp)
dev.off()

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(broad.prop.model,gdat.ml.prop,2000)
broad.prop.model.rac <- glm(cbind(nfix,nfixno)~landtype*status+abs.lat+area+elev_range+precipitation+ temperature+rac, data = gdat.ml.prop, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(broad.prop.model.rac)

#test zero inflation: not a problem
testZeroInflation(broad.prop.model.rac)
#test dispersion: slightly overdispersed
testDispersion(broad.prop.model.rac)
#test vif: ok
vif(broad.prop.model.rac)

#check model assumptions
#check residuals for each variable separately
E2 <- resid(broad.prop.model.rac, type = "pearson")
#list of variables landtype, status, abs.lat, area, elev_range, precipitation
plot(E2 ~ precipitation, data = gdat.ml.prop)
plot(E2 ~ area, data = gdat.ml.prop)
plot(E2 ~ abs.lat, data = gdat.ml.prop)
plot(E2 ~ elev_range, data = gdat.ml.prop)

# try variables NOT in model:
plot(E2 ~ temperature, data = gdat.ml.prop)

#contrasts
means <- emmeans(broad.prop.model.rac, ~landtype*status)
#look at means order to determine how to write contrasts:
means

#write contrasts:
contrasts <- list("mainland v all islands native" = c(-2,1,1,0,0,0),
                  "mainland v nonoceanic islands native" = c(-1,1,0,0,0,0),
                  "mainland v oceanic islands native" = c(-1,0,1,0,0,0),
                  "nonoceanic v oceanic islands native" = c(0,-1,1,0,0,0),
                  "mainland v all islands naturalized" = c(0,0,0,-2,1,1),
                  "mainland v oceanic islands naturalized" = c(0,0,0,-1,0,1),
                  "mainland v nonoceanic islands naturalized" = c(0,0,0,-1,1,0),
                  "nonoceanic v oceanic islands naturalized" = c(0,0,0,0,-1,1),
                  # is first contrast different native v non-native?
                  "mainland v all islands * status" = c(-2,1,1,2,-1,-1),
                  "mainland v oceanic islands * status" = c(-1,0,1,1,0,-1),
                  "mainland v nonoceanic islands * status" =  c(-1,1,0,1,-1,0),
                  "nonoceanic v oceanic island * status" = c(0,-1,1,0,1,-1),
                  "native v naturalized only on mainlands" =  c(-1,0,0,1,0,0),
                  "native v naturalized only on oceanic islands" =  c(0,0,-1,0,0,1),
                  "native v naturalized only on nonoceanic islands" =  c(0,-1,0,0,1,0),
                  "native v naturalized on both island types" =  c(0,-1,-1,0,1,1)
)

#extract results to df:
results <- lsmeans::contrast(means,contrasts)
results.df <- as.data.frame(results)
results.df


####MODELS ONLY INCLUDING NATURALIZED FLORAS (M3,M4)#############################################

#Load datasets and filter out native
gdat.isl.ntz <- gdat.ref %>%
  dplyr::select(c("entity_ID","entity_class","nfix","nfixno", "latitude","abs.lat","longitude", # Select columns
                  "landtype","status","presence","area","dist","elev_range","precipitation", "temperature","urbanland"))%>%
  filter(status=="naturalized")%>%
  add_row(entity_ID = 0000, entity_class="Island", area=6)%>% #add breakpoint for testing of small island effect
  mutate(area = as.vector(log10((area)+.01)))%>%  #log10 transformation of area for models only; remove for figs 
  filter(landtype=="oceanic")%>%
  drop_na() %>%
  mutate_at(c("abs.lat","area","dist","elev_range","precipitation", "temperature","urbanland"), scale)#%>%
  
breakpoint_data <- gdat.isl.ntz%>%filter(entity_ID==0000) 
breakpoint_value <- breakpoint_data$area

####M3 Presence of naturalized N-fixing species on oceanic islands####

#285 objects
model.oc.pres.full<- glm(presence~abs.lat +area +dist +elev_range	+precipitation+temperature +urbanland	+area:dist +urbanland:dist, data=gdat.isl.ntz, family=binomial(link ="logit"))
summary(model.oc.pres.full)

##### Correlogram 
mod <- model.oc.pres.ntz
dat<- gdat.isl.ntz

#calc residuals
resids.calc <- resid(mod)
sp.dat <- dat %>%
  mutate(resid = resids.calc)

#generate sp object
sp <- ncf::spline.correlog(x = as.numeric(sp.dat$latitude),
                           y = as.numeric(sp.dat$longitude),
                           z = as.numeric(sp.dat$resid),
                           xmax = 1000, resamp = 100, latlon=TRUE)
#save figure
png("figures/ntzocprescorrelogram.jpg", width = 10, height = 10, units = 'in', res = 300)
plot(sp)
dev.off()

#spatial correlation
rac <- Spat.cor.rep(model.oc.pres.full,gdat.isl.ntz,2000)
model.oc.pres.rac.full <- glm(presence~abs.lat +area +dist +elev_range+precipitation+ temperature +urbanland+area:dist +urbanland:dist+rac, data = gdat.isl.ntz, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(model.oc.pres.rac.full)

#Model with selected variables after stepwise regression
model.oc.pres.ntz<- glm(presence~abs.lat +area+ temperature+area:dist, data=gdat.isl.ntz, family=binomial(link ="logit"))
summary(model.oc.pres.ntz)




#spatial correlation
rac <- Spat.cor.rep(model.oc.pres.ntz,gdat.isl.ntz,2000)
model.oc.pres.ntz.rac <- glm(presence~abs.lat +area+ temperature+area:dist+rac, data = gdat.isl.ntz, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(model.oc.pres.ntz.rac)

#calculate the VIF for each predictor variable in the model
vif(model.oc.pres.ntz.rac)

#calculating variable importance
rsq(model.oc.pres.ntz.rac)
prsq.oc.pres.ntz<-rsq.partial(model.oc.pres.ntz.rac, adj=FALSE)

#check model assumptions
simulationOutput <- simulateResiduals(fittedModel = model.oc.pres.ntz.rac, plot = F)
#access residuals:
residuals(simulationOutput)
#qqplot
plotQQunif(simulationOutput)
#residuals
plotResiduals(simulationOutput)
#overdispersion
testDispersion(simulationOutput)
#zero inflation
testZeroInflation(simulationOutput)

#check residuals for each variable separately
E2 <- resid(model.oc.pres.ntz.rac, type = "pearson")
#list of variables: abs.lat, area, dist, elev_range, precipitation, urbanland
plot(E2 ~ temperature, data = gdat.isl.ntz)

var.plot<- ggeffect(model.oc.pres.ntz.rac, terms=c("temperature"), type="re")
temp.plot<- plot(var.plot, colors="coral")

#Saving the plot as a png
png("figures/temperature_oc_pres_ntz_orange.jpg", width=10, height= 10, units='in', res=300)
temp.plot
dev.off()


#Figure area:distance interaction

#"#d7191c","#fdae61","#abd9e9","#2c7bb6
#"coral","deepskyblue3","cadetblue2","cornflowerblue" "sienna1
colScale <- scale_colour_manual(values =c ("coral","coral2","coral3","sandybrown"))
fillScale <- scale_fill_manual(values =c ("coral","coral2","coral3","coral"))

model.oc.pres.ntz.rac$area
sd(gdat.isl.ntz$area)

area.minmax <- data.frame(area = c(mean(gdat.isl.ntz$area)-sd(gdat.isl.ntz$area),mean(gdat.isl.ntz$area),mean(gdat.isl.ntz$area)+sd(gdat.isl.ntz$area)), size =c ("small", "medium","large")) 

dist.range <- with(gdat.isl.ntz, expand.grid(dist = seq(min(dist), max(dist), length = 1000))) 

ex.grid <- expand.grid(area = area.minmax$area, dist = dist.range$dist)

pred.dat <- ex.grid %>% mutate(rac = mean(rac), abs.lat= mean(gdat.isl.ntz$abs.lat), urbanland=mean(gdat.isl.ntz$urbanland))
pred <- predict.glm(model.oc.pres.ntz.rac,type="response",newdata = pred.dat,se=TRUE)

pdat <- cbind(ex.grid,pred) %>%
  left_join(area.minmax)

areadist_oc_pres_ntz_plot<- ggplot(pdat, aes(x = dist, y = fit, fill=size, color=size))+
  geom_line() +
  geom_ribbon(aes(ymin=fit-se.fit, ymax=fit+se.fit), alpha=0.7)+ 
  #geom_point(data=gdat.isl.ntz, mapping=aes(x=dist, y=presence),alpha=0.2, size=2)+ 
  theme_minimal()+
  labs(y="Probability of N-fixing plant species")+
  xlab("Distance [km]")+
  theme(axis.text.y=element_text(size=30))+
  theme(axis.text.x=element_text(size=30))+
  theme(axis.title.y=element_text(size=30))+
  theme(axis.title.x=element_text(size=30))+
  #coord_cartesian(ylim=c(0,1))+
  colScale+
  fillScale+
  ylim(0,1)+
  #theme(legend.title = element_text(size=30))+ #change legend title font size
  #theme(legend.text = element_text(size=30))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=34), #change legend title font size
        legend.text = element_text(size=30),#change legend text font size
        legend.position = c(0.2, 0.85)) #change legend position

areadist_oc_pres_ntz_plot

#Saving the plot as a png
png("figures/areadist4_oc_pres_ntz_orange.jpg", width=10, height= 10, units='in', res=300)
areadist_oc_pres_ntz_plot
dev.off()


####Naturalized Oceanic Proportion####
oceanic.prop<-gdat.isl.ntz%>%
  mutate(species = nfix + nfixno)%>%
  filter(species > 0)%>%
  filter(precipitation< 4)


names(oceanic.prop)
cont.var <- c("abs.lat", "elev_range","area","dist","precipitation", "temperature", "urbanland")
Mypairs(oceanic.prop[,cont.var]) # area and elevation range 0.72
MyMultipanel.ggp2(oceanic.prop, cont.var, vary = "presence", ylab = "presence") # pretty hard to see what's going on

#####MODEL naturalized oceanic proportion ####
#number of objects 116
model.oc.prop.full<- glm(cbind(nfix,nfixno)~abs.lat +area +dist +elev_range	+precipitation+temperature +urbanland	+area:dist +urbanland:dist, data=oceanic.prop, family=binomial(link ="logit"))
summary(model.oc.prop.full)

#spatial correlation
rac <- Spat.cor.rep(model.oc.prop.full,oceanic.prop,2000)
model.oc.prop.rac.full <- glm(cbind(nfix,nfixno)~abs.lat +area +dist +elev_range+precipitation+temperature +urbanland	+area:dist +urbanland:dist+rac, data = oceanic.prop, family = binomial(link ="logit"))
summary(model.oc.prop.rac.full)

#selected variables
model.oc.prop<- glm(cbind(nfix,nfixno)~abs.lat +area +dist +elev_range+urbanland, data=oceanic.prop, family=binomial(link ="logit"))
summary(model.oc.prop)


##### CORRELOGRAM###
mod <- model.oc.prop
dat<- oceanic.prop

#calc residuals
resids.calc <- resid(mod)
sp.dat <- dat %>%
  mutate(resid = resids.calc)

#generate sp object
sp <- ncf::spline.correlog(x = as.numeric(sp.dat$latitude),
                           y = as.numeric(sp.dat$longitude),
                           z = as.numeric(sp.dat$resid),
                           xmax = 1000, resamp = 100, latlon=TRUE)
#save figure
png("figures/ntzocpropcorrelogram.jpg", width = 10, height = 10, units = 'in', res = 300)
plot(sp)
dev.off()

#spatial correlation
rac <- Spat.cor.rep(model.oc.prop,oceanic.prop,2000)
model.oc.prop.rac <- glm(cbind(nfix,nfixno)~abs.lat +area +dist+elev_range +urbanland+rac, data = oceanic.prop, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(model.oc.prop.rac)

#check variance inflation factor
vif(model.oc.prop.rac)

#check variable importance 
rsq(model.oc.prop.rac)
prsq.oc.prop.ntz<-rsq.partial(model.oc.prop.rac, adj=FALSE)

var.plot<- ggeffect(model.oc.prop.rac, terms=c("dist","area"), type="re")
plot(var.plot)

var.plot<- ggeffect(model.oc.prop.rac, terms=c("urbanland"), type="re")
landuse.plot<- plot(var.plot, colors="coral")

#Saving the plot as a png
png("figures/humanlanduse_oc_pres_ntz_orange.jpg", width=10, height= 10, units='in', res=300)
landuse.plot
dev.off()

#check for overdispersion
disp_check(model.oc.prop.rac, oceanic.prop)

#check model assumptions
#access residuals:
residuals(model.oc.prop.rac)

var.plot<- ggeffect(model.oc.prop.rac, terms=c("urbanland"), type="re")
plot(var.plot)

#check residuals for each variable separately
E2 <- resid(model.oc.prop.rac, type = "pearson")
#list of variables: abs.lat, area, dist, elev_range, precipitation, urbanland
plot(E2 ~ dist, data = oceanic.prop)



####NATIVE VERSION########################################################
gdat.isl.ntv <- gdat.ref %>%
  dplyr::select(c("entity_ID","entity_class","nfix","nfixno", "latitude","abs.lat","longitude", # Select columns
                  "landtype","status","presence","area","dist","elev_range","precipitation", "temperature"))%>%
  filter(status=="native")%>%
  #filter(area>6)%>%
  add_row(entity_ID = 0000, entity_class="Island", area=6)%>% # breakpoint
  mutate(area = as.vector(log10((area)+.01)))%>%  #log10 transformation of area for models only; remove for figs 
  filter(landtype=="oceanic")%>%
  drop_na()%>%
  mutate_at(c("abs.lat","area","dist","elev_range","precipitation", "temperature"), scale)

####NATIVE OCEANIC PRESENCE####
oceanic.pres.ntv <- gdat.isl.ntv%>%
  filter(!entity_ID == "594")# one data point removed due to extreme residual outlier status (see below)

# look at vars
plot(oceanic.pres.ntv$presence~oceanic.pres.ntv$area) # increased area seems to predict greater presence
plot(oceanic.pres.ntv$presence~oceanic.pres.ntv$dist) # unsure
plot(oceanic.pres.ntv$presence~oceanic.pres.ntv$abs.lat) # not sure
plot(oceanic.pres.ntv$presence~oceanic.pres.ntv$precipitation) # increase?

#####MODEL native oceanic presence####
#number of objects 289
model.oc.pres.full <- glm(presence~abs.lat + area + dist +elev_range +precipitation+temperature +area:dist , data=oceanic.pres.ntv, family=binomial(link ="logit"))
summary(model.oc.pres.full)

#spatial correlation
rac <- Spat.cor.rep(model.oc.pres.full,oceanic.pres.ntv,2000)
model.oc.pres.rac.full <- glm(presence~abs.lat + area + dist +elev_range +precipitation+temperature +area:dist +rac, data = oceanic.pres.ntv, family = binomial(link ="logit")) 
summary(model.oc.pres.rac.full)

#SELECTED VARIABLES
model.oc.pres <- glm(presence~abs.lat + area + dist +precipitation, data=oceanic.pres.ntv, family=binomial(link ="logit"))
summary(model.oc.pres)

##### CORRELOGRAM###
mod <- model.oc.pres
dat<- oceanic.pres.ntv

#calc residuals
resids.calc <- resid(mod)
sp.dat <- dat %>%
  mutate(resid = resids.calc)

#generate sp object
sp <- ncf::spline.correlog(x = as.numeric(sp.dat$latitude),
                           y = as.numeric(sp.dat$longitude),
                           z = as.numeric(sp.dat$resid),
                           xmax = 1000, resamp = 100, latlon=TRUE)
#save figure
png("figures/ntvocprescorrelogram.jpg", width = 10, height = 10, units = 'in', res = 300)
plot(sp)
dev.off()


#spatial correlation
rac <- Spat.cor.rep(model.oc.pres,oceanic.pres.ntv,2000)
model.oc.pres.rac <- glm(presence~abs.lat + area + dist +precipitation +rac, data = oceanic.pres.ntv, family = binomial(link ="logit")) 
summary(model.oc.pres.rac)

#check variance inflation factor 
vif(model.oc.pres.rac)

# pinpoint the weird one: entity_ID = 594
outlier <- oceanic.pres.ntv %>% filter(abs.lat < -1) %>% filter(dist < 1.1 & dist > 0.9) %>% filter(area < 1.1 & area > 0.9)

#check variable importance
rsq(model.oc.pres.rac)
prsq.oc.pres.ntv<-rsq.partial(model.oc.pres.rac, adj=FALSE)

#check model assumptions
simulationOutput <- simulateResiduals(fittedModel = model.oc.pres.rac, plot = F)
#access residuals:
residuals(simulationOutput)
#qqplot
plotQQunif(simulationOutput)
#residuals
plotResiduals(simulationOutput)
#overdispersion
testDispersion(simulationOutput)
#zero inflation
testZeroInflation(simulationOutput)

#check residuals for each variable separately
E2 <- resid(model.oc.pres.rac, type = "pearson")
#list of variables: abs.lat, area, dist, elev_range, precipitation
plot(E2 ~ precipitation, data = oceanic.pres.ntv)
#outlier in precipitation

var.plot<- ggeffect(model.oc.pres.rac, terms=c("area"), type="re")
plot(var.plot)

var.plot<- ggeffect(model.oc.pres.rac, terms=c("precipitation"), type="re")
precip.ntv.plot<- plot(var.plot, colors="darkseagreen3")

#Saving the plot as a png
png("figures/precip_oc_pres_ntv_orange.jpg", width=10, height= 10, units='in', res=300)
precip.ntv.plot
dev.off()

#distance plot
new.dat.oc <- with(oceanic.pres.ntv, expand.grid(dist = seq(min(dist), max(dist), length = 1000))) %>%
  mutate(rac = mean(rac), abs.lat= mean(oceanic.pres.ntv$abs.lat), area=mean(oceanic.pres.ntv$area), elev_range=mean(oceanic.pres.ntv$elev_range), precipitation=mean(oceanic.pres.ntv$precipitation)) #CREATE range for area from min to max value

pred.ml <- predict(model.oc.pres.rac, newdata = new.dat.oc, type = "response", se = TRUE) %>% # create new dataframe where area changes and all other variables are kept on mean
  as.data.frame() %>%
  mutate(dist = new.dat.oc$dist, landtype = "oceanic")

oc_pres_native_dist_plot<-ggplot(pred.ml, aes(x = dist, y = fit))+
  geom_line(color="darkseagreen3") +
  geom_ribbon(aes(ymin=fit-se.fit, ymax=fit+se.fit), alpha=0.7,color="darkseagreen3", fill="darkseagreen3")+ 
  #geom_point(data=oceanic.pres, mapping=aes(x=dist, y=presence),alpha=0.2, size=2,color="coral", fill="coral")+ 
  theme_minimal()+
  labs(y="Probability of N-fixing plants")+
  xlab("Distance [km]")+  
  theme(axis.text.y=element_text(size=30))+
  theme(axis.text.x=element_text(size=30))+
  theme(axis.title.y=element_text(size=30))+
  theme(axis.title.x=element_text(size=30))

oc_pres_native_dist_plot

#Saving the plot as a png
png("figures/dist_oc_pres_native_points.jpg", width=10, height= 10, units='in', res=300)
oc_pres_native_dist_plot
dev.off()

####NATIVE OCEANIC PROPORTION ####
oceanic.prop.ntv <- gdat.isl.ntv%>%
  #remove outliers
  filter(!entity_ID == 675) %>%
  filter(!entity_ID == 921) %>%
  filter(!entity_ID == 11474)%>%
  mutate(species = nfix + nfixno)%>% #find out values that have no species
  filter(species > 0)%>%
  filter(precipitation<4)

dat <- oceanic.prop.ntv
names(oceanic.prop.ntv)
cont.var <- c("abs.lat", "elev_range","area","dist","precipitation","temperature")
Mypairs(dat[,cont.var]) # area and elevation range 0.72
MyMultipanel.ggp2(dat, cont.var, vary = "presence", ylab = "proportion") # pretty hard to see what's going on

#####MODEL native oceanic proportion####
#cbind(nfix,nfixno)~abs.lat +area +dist +elev_range +precipitation + area:dist, number of objects 287
model.oc.prop<- glm(cbind(nfix,nfixno)~abs.lat+area +dist +precipitation + area:dist, data=oceanic.prop.ntv, family=binomial(link ="logit"))
summary(model.oc.prop)

#spatial correlation
rac <- Spat.cor.rep(model.oc.prop,oceanic.prop.ntv,2000)
model.oc.prop.rac <- glm(cbind(nfix,nfixno)~abs.lat +area +dist  +precipitation + area:dist+rac, data = oceanic.prop.ntv, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(model.oc.prop.rac)

#rac <- Spat.cor.rep(model.oc.prop,oceanic.prop.ntv,2000)
#model.oc.prop.rac.testr2 <- glm(cbind(nfix,nfixno)~area +dist  +precipitation+temperature + area:dist+rac, data = oceanic.prop.ntv, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac

#variance inflation factor
vif(model.oc.prop.rac)

#variable importance
rsq(model.oc.prop.rac)
prsq.oc.prop.ntv<-rsq.partial(model.oc.prop.rac, adj=FALSE)

#check dispersion
disp_check(model.oc.prop.rac,oceanic.prop.ntv)

##### CORRELOGRAM###
mod <- model.oc.prop
dat<- oceanic.prop.ntv

#calc residuals
resids.calc <- resid(mod)
sp.dat <- dat %>%
  mutate(resid = resids.calc)

#generate sp object
sp <- ncf::spline.correlog(x = as.numeric(sp.dat$latitude),
                           y = as.numeric(sp.dat$longitude),
                           z = as.numeric(sp.dat$resid),
                           xmax = 2000, resamp = 100, latlon=TRUE)
#save figure
png("figures/ntvocpropcorrelogram.jpg", width = 10, height = 10, units = 'in', res = 300)
plot(sp)
dev.off()

#check model assumptions
simulationOutput <- simulateResiduals(fittedModel = model.oc.prop.rac, plot = F)
#access residuals:
residuals(simulationOutput)
#qqplot
plotQQunif(simulationOutput)
#residuals
plotResiduals(simulationOutput)
#overdispersion
testDispersion(simulationOutput)
#zero inflation
testZeroInflation(simulationOutput)

#check residuals for each variable separately
E2 <- resid(model.oc.prop.rac, type = "pearson")
#list of variables: abs.lat, area, dist, elev_range, precipitation
plot(E2 ~ abs.lat, data = oceanic.prop.ntv)
#outlier in precipitation


var.plot<- ggeffect(model.oc.prop.rac, terms=c("temperature"), type="re")
temp.ntv.plot<- plot(var.plot, colors="darkseagreen3")

#Saving the plot as a png
png("figures/temp_oc_prop_ntv_orange.jpg", width=10, height= 10, units='in', res=300)
temp.ntv.plot
dev.off()


var.plot<- ggeffect(model.oc.prop.rac, terms=c("precipitation"), type="re")
temp.ntv.plot<- plot(var.plot, colors="darkseagreen3")

#Saving the plot as a png
png("figures/temp_oc_prop_ntv_orange.jpg", width=10, height= 10, units='in', res=300)
temp.ntv.plot
dev.off()

#model validation
mod <- model.oc.prop.rac
dat <- oceanic.prop.ntv
resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$area)
indep_cat_plot(mod,dat,dat$dist)
indep_cat_plot(mod,dat,dat$precipitation)

# get outlier in residuals:
F2 <- fitted(mod)
dat$F2<- F2 > 0.13 # 675,921,11474

var.plot<- ggeffect(model.oc.prop.rac, terms=c("dist","area"), type="re")
areadist.ntv.plot<-plot(var.plot)

#Saving the plot as a png
png("figures/areadist_oc_prop_ntv_orange.jpg", width=10, height= 10, units='in', res=300)
areadist.ntv.plot
dev.off()



#Figure area:distance interaction

colScale <- scale_colour_manual(values =c ("darkseagreen1","darkseagreen3","darkseagreen4","darkseagreen"))
fillScale <- scale_fill_manual(values =c ("darkseagreen1","darkseagreen3","darkseagreen4","darkseagreen"))

model.oc.prop.rac$area
sd(oceanic.prop.ntv$area)

area.minmax <- data.frame(area = c(mean(oceanic.prop.ntv$area)-sd(oceanic.prop.ntv$area),mean(oceanic.prop.ntv$area),mean(oceanic.prop.ntv$area)+sd(oceanic.prop.ntv$area)), size =c ("small", "medium","large")) 

dist.range <- with(oceanic.prop.ntv, expand.grid(dist = seq(min(dist), max(dist), length = 1000))) 

ex.grid <- expand.grid(area = area.minmax$area, dist = dist.range$dist)

pred.dat <- ex.grid %>% mutate(rac = mean(rac), abs.lat= mean(oceanic.prop.ntv$abs.lat), precipitation=mean(oceanic.prop.ntv$precipitation), temperature=mean(oceanic.prop.ntv$temperature))
pred <- predict.glm(model.oc.prop.rac,type="response",newdata = pred.dat,se=TRUE)

pdat <- cbind(ex.grid,pred) %>%
  left_join(area.minmax)

areadist_oc_pres_ntv_plot<- ggplot(pdat, aes(x = dist, y = fit, fill=size, color=size))+
  geom_line() +
  geom_ribbon(aes(ymin=fit-se.fit, ymax=fit+se.fit), alpha=0.7)+ 
  #geom_point(data=oceanic.prop.ntv, mapping=aes(x=dist, y=presence),alpha=0.2, size=2)+ 
  theme_minimal()+
  labs(y="Probability of N-fixing plant species")+
  xlab("Distance [km]")+
  theme(axis.text.y=element_text(size=30))+
  theme(axis.text.x=element_text(size=30))+
  theme(axis.title.y=element_text(size=30))+
  theme(axis.title.x=element_text(size=30))+
  #coord_cartesian(ylim=c(0,1))+
  colScale+
  fillScale+
  ylim(0,0.1)+
  #theme(legend.title = element_text(size=30))+ #change legend title font size
  #theme(legend.text = element_text(size=30))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=34), #change legend title font size
        legend.text = element_text(size=30),#change legend text font size
        legend.position = c(0.2, 0.85)) #change legend position

areadist_oc_pres_ntv_plot

#Saving the plot as a png
png("figures/areadist_oc_pres_ntv.jpg", width=10, height= 10, units='in', res=300)
areadist_oc_pres_ntv_plot
dev.off()


####VARIABLE IMPORTANCE COMBINED FIGURES NATURALIZED AND NATIVE####

#Create a custom color scale
colScale <- scale_colour_manual(values=c("darkseagreen3","coral"))
fillScale <- scale_fill_manual(values=c("darkseagreen3","coral"))


#####OCEANIC PRESENCE####
pres.vrsq.ntz<- cbind(prsq.oc.pres.ntz$variable, prsq.oc.pres.ntz$partial.rsq)
pres.vrsq.ntz.zeros<- data.frame(V1 = c("dist", "precipitation","elev_range", "urbanland"),  # create vector for unsignificant variables with zero
                            V2 = c(0,0,0,0))
pres.vrsq.ntz.new<- rbind(pres.vrsq.ntz,pres.vrsq.ntz.zeros) #add unsignificant variables
status<- c("naturalized","naturalized","naturalized","naturalized","naturalized","naturalized","naturalized","naturalized","naturalized") #create column with status for each variable

pres.var.impo.ntz<- cbind(pres.vrsq.ntz.new,status) #add status column
pres.var.impo.ntz.norac <- data.frame(pres.var.impo.ntz%>%filter(!pres.var.impo.ntz$V1=="rac"))

#same for native
pres.vrsq.ntv<- cbind(prsq.oc.pres.ntv$variable, prsq.oc.pres.ntv$partial.rsq)
pres.vrsq.ntv.zeros<- data.frame(V1 = c("temperature","area:dist","elev_range", "urbanland"),
                            V2 = c(0,0,0,0))

pres.vrsq.ntv.new<- rbind(pres.vrsq.ntv,pres.vrsq.ntv.zeros)
status<- c("native","native","native","native","native","native","native","native","native")

pres.var.impo.ntv<- cbind(pres.vrsq.ntv.new, status)
pres.var.impo.ntv.norac <- data.frame(pres.var.impo.ntv%>%filter(!pres.var.impo.ntv$V1=="rac"))


pres.comb.var.impo<- rbind(pres.var.impo.ntv.norac, pres.var.impo.ntz.norac)#combine naturalized and native version

pres.comb.var.impo2<- as.data.frame(pres.comb.var.impo)%>%   #make sure all zeros are saved as numeric and round to 5 digits
  mutate(V2= round(as.numeric(V2), digits=5))

colnames(pres.comb.var.impo2) <- c("variable", "partialrsq","status")  #change headers

#change order of variables for plot and change variable names
pres.comb.var.impo3<- pres.comb.var.impo2 %>% 
  mutate(variable = fct_relevel(variable, 
                                "area","dist","area:dist","abs.lat","elev_range", "precipitation","temperature","urbanland"))

combined.plot.oc.pres<- ggplot(pres.comb.var.impo3, aes(x=variable, y=partialrsq, fill=status, color = status)) +
  geom_bar(stat="identity", position="dodge", alpha=0.6) +
  colScale+
  fillScale+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90 ,size=20), axis.text.y = element_text(angle = 90, size=20))+
  #coord_flip()+
  xlab("")+
  ylab("partial R squared")
  #facet_grid(rows=vars(variable))
combined.plot.oc.pres

#Saving the plot as a png
png("figures/var_imp_oc_pres_combined_n.jpg", width=10, height= 10, units='in', res=300)
combined.plot.oc.pres
dev.off()

####OCEANIC PROPORTION####
#prsq.oc.prop.ntz #prsq.oc.prop.ntv
vrsq.ntz<- cbind(prsq.oc.prop.ntz$variable, prsq.oc.prop.ntz$partial.rsq) #create dataframe for plot
vrsq.ntz.zeros<- data.frame(V1 = c("precipitation","temperature","area:dist"),
                            V2 = c(0,0,0))

vrsq.ntz.new<- rbind(vrsq.ntz,vrsq.ntz.zeros)
status<- c("naturalized","naturalized","naturalized","naturalized","naturalized","naturalized","naturalized","naturalized","naturalized") #create column with status for each variable

var.impo.ntz<- cbind(vrsq.ntz.new,status) #add status column
var.impo.ntz.norac <- data.frame(var.impo.ntz%>%filter(!var.impo.ntz$V1=="rac"))

#same for native
vrsq.ntv<- cbind(prsq.oc.prop.ntv$variable, prsq.oc.prop.ntv$partial.rsq)
vrsq.ntv.zeros<- data.frame(V1 = c("urbanland","elev_range","temperature"),
                            V2 = c(0,0,0))

vrsq.ntv.new<- rbind(vrsq.ntv,vrsq.ntv.zeros)
status<- c("native","native","native","native","native","native","native","native","native")

var.impo.ntv<- cbind(vrsq.ntv.new, status)
var.impo.ntv.norac <- data.frame(var.impo.ntv%>%filter(!var.impo.ntv$V1=="rac"))

comb.var.impo<- rbind(var.impo.ntv.norac, var.impo.ntz.norac)  #combine naturalized and native version
comb.var.impo2<- as.data.frame(comb.var.impo)%>%   #make sure all zeros are saved as numeric and round to 5 digits
  mutate(V2= round(as.numeric(V2), digits=5))

colnames(comb.var.impo2) <- c("variable", "partialrsq","status")  #change headers

#change order of variables for plot, arrange(variable) didn't change plot
comb.var.impo3<- comb.var.impo2 %>% 
  mutate(variable = fct_relevel(variable, 
                                "area","dist","area:dist","abs.lat","elev_range","precipitation","temperature","urbanland"))

combined.plot.oc.prop<- ggplot(comb.var.impo3, aes(x=variable, y=partialrsq, fill=status, color=status)) +
  geom_bar(stat="identity", position="dodge",alpha=0.6) +
  colScale+
  fillScale+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90,size=20), axis.text.y = element_text(angle = 0, size=20))+
  coord_flip()+
  xlab("")+
  ylab("partial R squared")
combined.plot.oc.prop

#Saving the plot as a png
png("figures/var_imp_oc_prop_combinedsc_october24.jpg", width=10, height= 10, units='in', res=300)
combined.plot.oc.prop
dev.off()


####NOTES####
library(partR2)
R2_BM <- partR2(model.oc.pres.rac, data = gdat.isl.ntz, R2_type = "marginal", nboot = 10)
R2_BM

##' What is the percentage of zeros?
100 * sum(oceanic.prop$nfix == 0, na.rm = TRUE) / nrow(gdat.isl.ntz)
100 * sum(oceanic.prop$nfixno == 0, na.rm = TRUE) / nrow(gdat.isl.ntz) 

####Package versions and citations####
# Get package versions to report in MS:
packageVersion(c("dplyr")) #1.1.4
packageVersion(c("segmented")) #2.0.2
packageVersion(c("nlme")) #3.1.163
packageVersion(c("mgcv")) #1.9.0
packageVersion(c("gridExtra"))#2.3
packageVersion(c("betareg")) #3.1.4
packageVersion(c("MASS")) # 7.3.60
packageVersion(c("lme4")) #1.1.35.1
packageVersion(c("lmerTest")) #3.1.3
packageVersion(c("lsmeans")) #2.30.0
packageVersion(c("ggeffects")) #1.4.0
packageVersion(c("ggplot2")) # 3.5.0
packageVersion(c("effects")) #4.2.2
packageVersion(c("ncf")) #1.3.2
packageVersion(c("ape")) #5.7.1
packageVersion(c("sjPlot")) #2.8.15
packageVersion(c("gridExtra")) #2.3
packageVersion(c("MuMIn")) #1.47.5
packageVersion(c("tidyverse")) # 2.0.0
packageVersion(c("car")) #3.1.2
packageVersion(c("V8")) #4.4.1
packageVersion(c("rsq")) # 2.6
packageVersion(c("Matrix")) #1.6.5


get_pkgs_info(pkgs = c("dplyr", "segmented","nlme","mgcv","gridExtra","betareg",
                       "MASS","lme4","lmerTest","lsmeans","ggeffects","spdep",
                       "ggplot2","effects","ncf","ape","sjPlot",
                       "MuMIn","tidyverse","car", "V8", "rsq", "DHARMa"), out.dir = getwd())

citation(c("dplyr"))
citation(c("segmented")) #2.0.2
citation(c("nlme")) #3.1.163
citation(c("mgcv")) #1.9.0
citation(c("gridExtra"))#2.3
citation(c("betareg")) #3.1.4
citation(c("MASS")) # 7.3.60
citation(c("lme4")) #1.1.35.1
citation(c("lmerTest")) #3.1.3
citation(c("lsmeans")) #2.30.0
citation(c("ggeffects")) #1.4.0
citation(c("ggplot2")) # 3.5.0
citation(c("effects")) #4.2.2
citation(c("ncf")) #1.3.2
citation(c("ape")) #5.7.1
citation(c("sjPlot")) #2.8.15
citation(c("MuMIn")) #1.47.5
citation(c("tidyverse")) # 2.0.0
citation(c("car")) #3.1.2
citation(c("V8")) #4.4.1
citation(c("rsq")) # 2.6
citation(c("Matrix"))

