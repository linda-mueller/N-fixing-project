#N-fixing project: Script version November 2024
#Author Linda Rosa Mueller

####Libraries and initial dataset####

# Load the necessary libraries
library(dplyr);library(segmented);library(nlme);library(mgcv);library(gridExtra);library(betareg);
library(MASS);library(lme4);library(lmerTest);library(lsmeans);library(ggeffects);library(spdep);
library(ggplot2);library(effects);library(ncf);library(ape);library(sjPlot);library(gridExtra);
library(MuMIn);library(tidyverse);library(car); library(V8); library(rsq); library(DHARMa); library(grateful);library(broom)

source("extra_functions.R")
options(na.action = "na.fail") #Change na. action

#Load dataset from data cleaning script
gdat.ref <- readRDS("data/fulldata_for_analysis_2024.rds") %>%
  dplyr::select(c("entity_ID","entity_class","nfix","nfixno", "latitude","abs.lat","longitude", # Select columns
                  "landtype","status","presence","area","elev_range","precipitation","dist", "temperature","urbanland")) 

#### M1 Broad Presence (across islands and mainlands including native and naturalized floras)####

#dataset for broad presence model including mainlands
gdat.ml.pres <- gdat.ref %>%
  dplyr::select(c("entity_ID","nfix","nfixno", "latitude","abs.lat","longitude", # Select columns
                  "landtype","status","presence","area","elev_range","precipitation","temperature"))%>%
  filter(!landtype=="other_island")%>%                                          #remove islands with undetermined type
  mutate(area = as.vector(log10((area)+.01)))%>%                                #log transform area
  drop_na()%>%
  mutate_at(c("abs.lat","area","elev_range","precipitation", "temperature"), scale) #scale all explanatory variables

#check correlations
dat <- gdat.ml.pres
names(gdat.ml.pres)
cont.var <- c("abs.lat", "elev_range","area","precipitation","temperature")
Mypairs(dat[,cont.var]) # area and elevation range 0.56, abs.lat and temperature -0.88


#### M1 Broad Presence (across islands and mainlands including native and naturalized floras)
broad.pres.model.full<- glm(presence~landtype*status+abs.lat+area+elev_range+precipitation+temperature, data =gdat.ml.pres, family= binomial(link ="logit"))
summary(broad.pres.model.full)

#ref.land.stat <- lsmeans(broad.pres.model.full,pairwise~landtype*status, data= gdat.ml.pres, type="response")
#ref.table.land.stat <- as.data.frame(ref.land.stat$lsmeans) 

#model including spatial correlation variable
rac <- Spat.cor.rep(broad.pres.model.full,gdat.ml.pres,2000)
broad.pres.model.rac.full <- glm(presence~landtype*status+abs.lat+area+elev_range+precipitation+temperature+rac, data = gdat.ml.pres, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(broad.pres.model.rac.full)

#M1 with only selected variables after stepwise regression
broad.pres.model<- glm(presence~landtype*status+area+elev_range+precipitation+temperature, data =gdat.ml.pres, family= binomial(link ="logit"))
summary(broad.pres.model)

ref.land.stat <- lsmeans(broad.pres.model,pairwise~landtype*status, data= gdat.ml.pres, type="response")
ref.table.land.stat <- as.data.frame(ref.land.stat$lsmeans) 

#Correlogram to test distance of spatial autocorrelation
correlogram(broad.pres.model, gdat.ml.pres, "figures/M1_broadpres_correlogram.jpg")

#model including spatial correlation variable
rac <- Spat.cor.rep(broad.pres.model,gdat.ml.pres,2000)
broad.pres.model.rac <- glm(presence~landtype*status+area+elev_range+precipitation+temperature+rac, data = gdat.ml.pres, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(broad.pres.model.rac)

#add rac variable to dataframe to include in lsmeans
gdat.ml.pres$rac <- rac

ref.land.stat.rac<-lsmeans(broad.pres.model.rac,pairwise~landtype*status, data = gdat.ml.pres, type="response")
ref.table.land.stat.rac<-as.data.frame(ref.land.stat.rac$lsmeans) 

#check variance inflation factor (should be below 5 for all variables)
vif(broad.pres.model.rac)

#check model assumptions
simulationOutput <- simulateResiduals(fittedModel = broad.pres.model.rac, plot = F)
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

#Naturalized Plot 
#subset data to only include naturalized for plot
gdat.ml.ntz<- gdat.ml.pres%>%filter(status=="naturalized")
ref.table.land.stat.rac.ntz <- ref.table.land.stat.rac%>%filter(!status=="native") #subset lsmeans to only include naturalized for plot 

#Create a custom color scale
colScale <- scale_colour_manual(values=c("coral"))
fillScale <- scale_fill_manual(values=c("coral"))

#plot the model results for M1 presence of naturalized flora across islands and mainlands
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
png("figures/M1_broad_presence_naturalized.jpg", width=10, height= 12, units='in', res=300)
broad.pres.ntz
dev.off()

#NATIVE PLOT  
#subset data to only include native species for plot
gdat.ml.ntv<- gdat.ml.pres%>%filter(status=="native")
ref.table.land.stat.rac.ntv <- ref.table.land.stat.rac%>%filter(!status=="naturalized") #subset lsmeans to only include native for plot

#Create a custom color scale
colScale <- scale_colour_manual(values=c("darkseagreen3"))
fillScale <- scale_fill_manual(values=c("darkseagreen3"))

#plot the model results for M1 presence of native flora across islands and mainlands
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
png("figures/M1_broad_presence_native.jpg", width=10, height= 12, units='in', res=300)
broad.pres.ntv
dev.off()


####M2 Broad Proportion (across islands and mainlands including native and naturalized floras)####

#dataset for M2 broad proportion model including mainlands
gdat.ml.prop <- gdat.ref %>%
  dplyr::select(c("entity_ID","nfix","nfixno", "latitude","abs.lat","longitude", # Select columns
                  "landtype","status","presence","area","elev_range","precipitation", "temperature"))%>%
  filter(!landtype=="other_island")%>%                                          #remove islands with undetermined type
  mutate(area = as.vector(log10((area)+.01)))%>%                                #log transform area
  mutate(species = nfix + nfixno) %>%   #find out rows that have no species counts
  filter(species > 0) %>%
  drop_na() %>%
  mutate_at(c("abs.lat","area","elev_range","precipitation","temperature"), scale) %>%  #scale all explanatory variables
  filter(precipitation < 3.5)   #remove outliers

ml<- gdat.ml.prop%>% filter(status=="naturalized")%>%filter(landtype== "mainland")
nocisl<- gdat.ml.prop%>%filter(status=="naturalized") %>%filter(landtype== "nonoceanic")
ocisl<- gdat.ml.prop%>%filter(status=="naturalized") %>%filter(landtype== "oceanic")


#check correlations
dat <- gdat.ml.prop
names(gdat.ml.prop)
cont.var <- c("abs.lat", "elev_range","area","precipitation","temperature")
Mypairs(dat[,cont.var]) # area and elevation range 0.54, abs.lat and temperature -0.88

#### M2 Model Broad Proportion across islands and mainlands including native and naturalized flora 
broad.prop.model<- glm(cbind(nfix,nfixno)~landtype*status+abs.lat+area+elev_range+precipitation+temperature, data =gdat.ml.prop, family= binomial(link ="logit"))
summary(broad.prop.model)

#ref.land.stat <- lsmeans(broad.prop.model,pairwise~landtype*status, data= gdat.ml.prop, type="response")
#ref.table.land.stat <- as.data.frame(ref.land.stat$lsmeans) 

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(broad.prop.model,gdat.ml.prop,2000)
broad.prop.model.rac <- glm(cbind(nfix,nfixno)~landtype*status+abs.lat+area+elev_range+precipitation+temperature+rac, data = gdat.ml.prop, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(broad.prop.model.rac)


#M2 with only selected variables after stepwise regression because vif value of abs.lat is too high
broad.prop.model<- glm(cbind(nfix,nfixno)~landtype*status+area+elev_range+precipitation+temperature, data =gdat.ml.prop, family= binomial(link ="logit"))
summary(broad.prop.model)

ref.land.stat <- lsmeans(broad.prop.model,pairwise~landtype*status, data= gdat.ml.prop, type="response")
ref.table.land.stat <- as.data.frame(ref.land.stat$lsmeans) 

#Correlogram to test distance of spatial autocorrelation
correlogram(broad.prop.model, gdat.ml.prop, "figures/M2_broadprop_correlogram.jpg")

#M2 including spatial correlation variable (rac)
rac <- Spat.cor.rep(broad.prop.model,gdat.ml.prop,2000)
broad.prop.model.rac <- glm(cbind(nfix,nfixno)~landtype*status+area+elev_range+precipitation+temperature+rac, data = gdat.ml.prop, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(broad.prop.model.rac)

#add rac to df to include in lsmeans.
gdat.ml.prop$rac <- rac

ref.land.stat.rac<-lsmeans(broad.prop.model.rac,pairwise~landtype*status, data = gdat.ml.prop, type="response")
ref.table.land.stat.rac<-as.data.frame(ref.land.stat.rac$lsmeans) 

#check variance inflation factor (should be below 5 for all variables)
vif(broad.prop.model.rac)

#check model assumptions
simulationOutput <- simulateResiduals(fittedModel = broad.prop.model.rac, plot = F)
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
E2 <- resid(broad.prop.model.rac, type = "pearson")
#list of variables landtype, status, area, elev_range, precipitation, temperature
plot(E2 ~ precipitation, data = gdat.ml.prop)
plot(E2 ~ area, data = gdat.ml.prop)
plot(E2 ~ temperature, data = gdat.ml.prop)
plot(E2 ~ elev_range, data = gdat.ml.prop)

# try variables NOT in model:
plot(E2 ~ abs.lat, data = gdat.ml.prop)

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


#Plot M2 Broad Proportion naturalized

#subset data to only include naturalized for plot
ref.table.land.stat.rac.ntz <- ref.table.land.stat.rac%>%filter(status=="naturalized") #subset lsmeans to only include naturalized for plot 

#Create a custom color scale
colScale <- scale_colour_manual(values=c("coral","coral","coral","coral"))
fillScale <- scale_fill_manual(values=c("coral","coral","coral","coral"))

gdat3cat<- gdat.ml.prop %>% mutate(propnfix=nfix/(nfix+nfixno))%>%filter(status=="naturalized") #only naturalized 

broadpropntz<- 
  ggplot(ref.table.land.stat.rac.ntz, aes(x = landtype, y = lsmean, fill = status, color = status))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) + 
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3)+
  geom_point(data = gdat3cat,aes(x=landtype,y=propnfix,color=landtype),position=position_jitterdodge(dodge.width = 1),size=2,alpha=0.3)+
  coord_cartesian(ylim=c(0,0.02))+
  colScale+
  fillScale+
  labs(y="Proportion N-fixing Plant Species")+
  xlab(" ")+
  #ggtitle("Model Proportion~ Landtype* Status only ntz")+
  theme_classic(base_size = 30) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  theme(legend.position = "top")+
  theme(legend.key = element_rect(colour ="coral" ))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
  guides(color="none", fill = guide_legend(title="Floral status:",override.aes = list(shape = NA,size = 10)))+
  theme(plot.title=element_text(size=20))

broadpropntz

#Saving the plot as a png
png("figures/M2_broad_proportion_naturalized.jpg", width=10, height= 10, units='in', res=300)
broadpropntz
dev.off()

#Native plot
ref.table.land.stat.rac.ntv <- ref.table.land.stat.rac%>%filter(status=="native") #subset lsmeans to only include naturalized for plot 

#Create a custom color scale
colScale <- scale_colour_manual(values=c("darkseagreen3","darkseagreen3"))
fillScale <- scale_fill_manual(values=c("darkseagreen3","darkseagreen3"))

gdat3cat<- gdat.ml.prop %>% mutate(propnfix=nfix/(nfix+nfixno))%>%filter(status=="native")

broadpropntv<- 
  ggplot(ref.table.land.stat.rac.ntv, aes(x = landtype, y = lsmean, fill = status, color = status))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) + 
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3)+
  geom_point(data = gdat3cat,aes(x=landtype,y=propnfix,color=status),position=position_jitterdodge(dodge.width = 1),size=2,alpha=0.3)+
  coord_cartesian(ylim=c(0,0.16))+
  colScale+
  fillScale+
  labs(y="Proportion N-fixing Plant Species")+
  xlab(" ")+
  #ggtitle("Model Proportion~ Landtype* Status")+
  theme_classic(base_size = 30) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  theme(legend.position = "top")+
  theme(legend.key = element_rect(colour =c("darkseagreen3","darkseagreen3") ))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
  guides(color="none", fill = guide_legend(title="Floral status:",override.aes = list(shape = NA)))+
  theme(axis.text.y=element_text(size=30))+
  theme(axis.text.x=element_text(size=30))                  # All font sizes

broadpropntv

#Saving the plot as a png
png("figures/M2_broad_proportion_native.jpg", width=10, height= 10, units='in', res=300)
broadpropntv
dev.off()


####Data only including NATURALIZED species on OCEANIC ISLANDS (M3,M4)####

#Load dataset and subset to only naturalized species and oceanic islands
gdat.isl.ntz <- gdat.ref %>%
  dplyr::select(c("entity_ID","entity_class","nfix","nfixno", "latitude","abs.lat","longitude", # Select columns
                  "landtype","status","presence","area","dist","elev_range","precipitation", "temperature","urbanland"))%>%
  filter(status=="naturalized")%>%
  mutate(area = as.vector(log10((area)+.01)))%>%  #log10 transformation of area for models only; remove for figures
  filter(landtype=="oceanic")%>%
  drop_na()%>%
  mutate_at(c("abs.lat","area","dist","elev_range","precipitation", "temperature","urbanland"), scale)#%>%  #scale all explanatory variables

#check correlations
dat <- gdat.isl.ntz
names(gdat.isl.ntz)
cont.var <- c("abs.lat", "area","elev_range","precipitation","temperature","dist", "urbanland")
Mypairs(dat[,cont.var]) # area and elevation range 0.72, abs.lat and temperature 0.88


####M3 Presence of naturalized N-fixing species on oceanic islands####

model.oc.pres.full<- glm(presence~abs.lat +elev_range	+precipitation+temperature	+area*dist +urbanland*dist, data=gdat.isl.ntz, family=binomial(link ="logit"))
summary(model.oc.pres.full)

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(model.oc.pres.full,gdat.isl.ntz,2000)
model.oc.pres.rac.full <- glm(presence~abs.lat +elev_range+precipitation+ temperature+area*dist +urbanland*dist+rac, data = gdat.isl.ntz, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(model.oc.pres.rac.full)

#model with selected variables after stepwise regression
model.oc.pres<- glm(presence~abs.lat+area+dist+area:dist, data=gdat.isl.ntz, family=binomial(link ="logit"))
summary(model.oc.pres)

#Correlogram to test distance of spatial autocorrelation
correlogram(model.oc.pres, gdat.isl.ntz, "figures/M3_ntzocpres_correlogram.jpg")

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(model.oc.pres,gdat.isl.ntz,2000)
model.oc.pres.rac<- glm(presence~ abs.lat+area+dist+area:dist+rac, data = gdat.isl.ntz, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(model.oc.pres.rac)

# convert the model output into a df
coef_dataM3 <- tidy(model.oc.pres.rac, conf.int = TRUE)

#check variance inflation factor 
vif(model.oc.pres.rac)

#calculating variable importance using Rsquared

library(glmtoolbox)#doesnt have a partial R2 function
adjR2(model.oc.pres.rac)



####using manually coded function for Partial R square and Rsquared
resp<-gdat.isl.ntz$presence
formula<- presence~ abs.lat+area+area:dist+rac
dat<-gdat.isl.ntz

#calculate R square of model
R2<- r2_glm(resp, formula, dat)

#for a single variable
#variable <- "area:dist"
#p2<- partial_r2_glm(resp, formula, dat, variable)

#trying partial R2 with loop for all variables
var_list<- c("abs.lat","area","area:dist","rac")
pr2.oc.pres.ntz<- partial_r2_glm_list(presence,formula,dat,var_list)

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
#list of variables: abs.lat, area, dist
plot(E2 ~ abs.lat, data = gdat.isl.ntz)

var.plot<- ggeffect(model.oc.pres.rac, terms=c("abs.lat"), type="re")
effect.plot<- plot(var.plot, colors="coral")
effect.plot


#Figure area:distance interaction
#create a color scale
colScale <- scale_colour_manual(values =c ("coral","coral2","coral3","sandybrown"))
fillScale <- scale_fill_manual(values =c ("coral","coral2","coral3","coral"))

model.oc.pres.rac$area
sd(gdat.isl.ntz$area)

area.minmax <- data.frame(area = c(mean(gdat.isl.ntz$area)-sd(gdat.isl.ntz$area),mean(gdat.isl.ntz$area),mean(gdat.isl.ntz$area)+sd(gdat.isl.ntz$area)), size =c ("small", "medium","large")) 

dist.range <- with(gdat.isl.ntz, expand.grid(dist = seq(min(dist), max(dist), length = 1000))) 

ex.grid <- expand.grid(area = area.minmax$area, dist = dist.range$dist)

pred.dat <- ex.grid %>% mutate(rac = mean(rac), abs.lat= mean(gdat.isl.ntz$abs.lat))
pred <- predict.glm(model.oc.pres.rac,type="response",newdata = pred.dat,se=TRUE)

pdat <- cbind(ex.grid,pred) %>%
  left_join(area.minmax)

#plot the area distance interaction
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
  colScale+
  fillScale+
  ylim(0,1)+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=34), #change legend title font size
        legend.text = element_text(size=30),#change legend text font size
        legend.position = c(0.2, 0.85)) #change legend position

areadist_oc_pres_ntz_plot

#Saving the plot as a png
png("figures/M3_areadist_interaction.jpg", width=10, height= 10, units='in', res=300)
areadist_oc_pres_ntz_plot
dev.off()


####M4 Proportion of naturalized N-fixing species on oceanic islands####

#load data and filter out all locations without species counts
oceanic.prop<-gdat.isl.ntz%>%
  mutate(species = nfix + nfixno)%>%
  filter(species >0)%>%
  filter(precipitation< 4) #take out outlier (only do this in scaled version for modelling remove when creating figures)

#check correlations
names(oceanic.prop)
cont.var <- c("abs.lat", "elev_range","area","dist","precipitation", "temperature", "urbanland")
Mypairs(oceanic.prop[,cont.var]) # area and elevation range 0.77, abs.lat temperature -0.87


#####M4  naturalized proportion on oceanic islands
model.oc.prop.full<- glm(cbind(nfix,nfixno)~abs.lat +area +dist +elev_range	+precipitation+temperature +urbanland	+area:dist +urbanland:dist, data=oceanic.prop, family=binomial(link ="logit"))
summary(model.oc.prop.full)

#model including spatial correlation variable
rac <- Spat.cor.rep(model.oc.prop.full,oceanic.prop,2000)
model.oc.prop.rac.full <- glm(cbind(nfix,nfixno)~abs.lat +area +dist +elev_range+precipitation+temperature +urbanland	+area:dist +urbanland:dist+rac, data = oceanic.prop, family = binomial(link ="logit"))
summary(model.oc.prop.rac.full)

#M4 with only selected variables after stepwise regression
model.oc.prop<- glm(cbind(nfix,nfixno)~abs.lat +area +dist, data=oceanic.prop, family=binomial(link ="logit"))
summary(model.oc.prop)

#Correlogram to test distance of spatial autocorrelation
correlogram(model.oc.prop, oceanic.prop, "figures/M4_ntzocprop_correlogram.jpg")

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(model.oc.prop,oceanic.prop,2000)
model.oc.prop.rac <- glm(cbind(nfix,nfixno)~abs.lat +area +dist+rac, data = oceanic.prop, family = binomial(link ="logit"))
summary(model.oc.prop.rac)


# convert the model output into a df for estimates plot
coef_dataM4 <- tidy(model.oc.prop.rac, conf.int = TRUE)

#check variance inflation factor (should be below 5 for all variables)
vif(model.oc.prop.rac)

#check variable importance 
library(glmtoolbox)#doesnt have a partial R2 function
adjR2(model.oc.prop.rac)

####using manually coded function for Partial R square and Rsquared
formula<- cbind(dat$nfix,dat$nfixno)~ abs.lat +area +dist+rac
dat<-oceanic.prop
resp<-cbind(dat$nfix,dat$nfixno)

#calculate R square of model
R2<- r2_glm(resp, formula, dat)

#for a single variable
variable <- "area"
p2<- partial_r2_glm(resp, formula, dat, variable)

#trying partial R2 with loop for all variables
var_list<- c("abs.lat","area","dist","rac")
pr2.oc.prop.ntz<- partial_r2_glm_list(resp,formula,dat,var_list)

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

#check for overdispersion
disp_check(model.oc.prop.rac, oceanic.prop)

#access residuals:
residuals(model.oc.prop.rac)

#check residuals for each variable separately
E2 <- resid(model.oc.prop.rac, type = "pearson")
#list of variables: abs.lat, area, dist, elev_range, precipitation, urbanland
plot(E2 ~ dist, data = oceanic.prop)

#test area distance effect
var.plot<- ggeffect(model.oc.prop.rac, terms=c("dist","area"), type="re")
plot(var.plot)

#plot effect of human land use
var.plot<- ggeffect(model.oc.prop.rac, terms=c("urbanland"), type="re")
landuse.plot<- plot(var.plot, colors="coral")
landuse.plot

#Saving the plot as a png
png("figures/M4_landuse_oc_pres_ntz.jpg", width=10, height= 10, units='in', res=300)
landuse.plot
dev.off()

####Data only including NATIVE species on OCEANIC ISLANDS (M5,M6)####

#Load dataset and subset to only native species and oceanic islands
gdat.isl.ntv <- gdat.ref %>%
  dplyr::select(c("entity_ID","entity_class","nfix","nfixno", "latitude","abs.lat","longitude", # Select columns
                  "landtype","status","presence","area","dist","elev_range","precipitation", "temperature"))%>%
  filter(status=="native")%>%
  #filter(area>6)%>%
  add_row(entity_ID = 0000, entity_class="Island", area=6)%>% # breakpoint
  mutate(area = as.vector(log10((area)+.01)))%>%  #log10 transformation of area for models only; remove for figs 
  filter(landtype=="oceanic")%>%
  drop_na()%>%
  mutate_at(c("abs.lat","area","dist","elev_range","precipitation", "temperature"), scale)  #scale all explanatory variables


####M5 presence of native N-fixing species on oceanic islands####

#create data subset for presence analysis
oceanic.pres.ntv <- gdat.isl.ntv%>%
  filter(!entity_ID == "594")# one data point removed due to extreme residual outlier status (see below)

#check correlations
dat <- oceanic.pres.ntv
names(oceanic.pres.ntv)
cont.var <- c("abs.lat", "area","elev_range","precipitation","temperature","dist")
Mypairs(dat[,cont.var]) # area and elevation range 0.73, abs.lat and temperature -0.89

#M5
model.oc.pres.full <- glm(presence~abs.lat + area + dist +elev_range +precipitation+temperature +area:dist , data=oceanic.pres.ntv, family=binomial(link ="logit"))
summary(model.oc.pres.full)

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(model.oc.pres.full,oceanic.pres.ntv,2000)
model.oc.pres.rac.full <- glm(presence~abs.lat + area + dist +elev_range +precipitation+temperature +area:dist +rac, data = oceanic.pres.ntv, family = binomial(link ="logit")) 
summary(model.oc.pres.rac.full)

#M5 with only selected variables after stepwise regression
model.oc.pres <- glm(presence~ area + dist +precipitation+temperature , data=oceanic.pres.ntv, family=binomial(link ="logit"))
summary(model.oc.pres)

#Correlogram to test distance of spatial autocorrelation
correlogram(model.oc.pres, oceanic.pres.ntv, "figures/M5_ntvocpres_correlogram.jpg")

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(model.oc.pres,oceanic.pres.ntv,2000)
model.oc.pres.rac <- glm(presence~ area + dist +precipitation+temperature +rac, data = oceanic.pres.ntv, family = binomial(link ="logit")) 
summary(model.oc.pres.rac)


# Tidy the model output
coef_dataM5 <- tidy(model.oc.pres.rac, conf.int = TRUE)

#check variance inflation factor (should be below 5 for all variables)
vif(model.oc.pres.rac)

#check variable importance
library(glmtoolbox)#doesnt have a partial R2 function
adjR2(model.oc.pres.rac)
####using manually coded function for Partial R square and Rsquared

formula<- presence~ area + dist +precipitation+temperature +rac
dat<-oceanic.pres.ntv
resp<-dat$presence


#calculate R square of model
R2<- r2_glm(resp, formula, dat)

#for a single variable
variable <- "area"
p2<- partial_r2_glm(resp, formula, dat, variable)

#trying partial R2 with loop for all variables
var_list<- c("area","dist","precipitation","temperature","rac")
pr2.oc.pres.ntv<- partial_r2_glm_list(resp,formula,dat,var_list)


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

# pinpoint the outlier: entity_ID = 594
outlier <- oceanic.pres.ntv %>% filter(abs.lat < -1) %>% filter(dist < 1.1 & dist > 0.9) %>% filter(area < 1.1 & area > 0.9)

var.plot<- ggeffect(model.oc.pres.rac, terms=c("area"), type="re")
plot(var.plot)

var.plot<- ggeffect(model.oc.pres.rac, terms=c("precipitation"), type="re")
precip.ntv.plot<- plot(var.plot, colors="darkseagreen3")
precip.ntv.plot

#Saving the plot as a png
png("figures/M5_precip_oc_pres_ntv.jpg", width=10, height= 10, units='in', res=300)
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
png("figures/M5_dist_oc_pres_ntv.jpg", width=10, height= 10, units='in', res=300)
oc_pres_native_dist_plot
dev.off()

####M6 proportion of native N-fixing species on oceanic islands ####

#subset data for proportion analysis
oceanic.prop.ntv <- gdat.isl.ntv%>%
  #remove outliers
  filter(!entity_ID == 675) %>%
  filter(!entity_ID == 921) %>%
  filter(!entity_ID == 11474)%>%
  mutate(species = nfix + nfixno)%>% #find out values that have no species counts
  filter(species > 0)%>%
  filter(precipitation<4)

#check correlations
dat <- oceanic.prop.ntv
names(oceanic.prop.ntv)
cont.var <- c("abs.lat", "elev_range","area","dist","precipitation","temperature")
Mypairs(dat[,cont.var]) # area and elevation range 0.73, abs.lat and temperature -0.89


#####M6 Model proportion of native N-fixing species on oceanic islands
model.oc.prop.full<- glm(cbind(nfix,nfixno)~abs.lat+area +dist +elev_range+precipitation +temperature+ area:dist, data=oceanic.prop.ntv, family=binomial(link ="logit"))
summary(model.oc.prop.full)

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(model.oc.prop.full,oceanic.prop.ntv,2000)
model.oc.prop.rac.full <- glm(cbind(nfix,nfixno)~abs.lat +area +dist+elev_range +precipitation+temperature + area:dist+rac, data = oceanic.prop.ntv, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(model.oc.prop.rac.full)

#selected variables
model.oc.prop<- glm(cbind(nfix,nfixno)~abs.lat+area +dist+precipitation +temperature+ area:dist, data=oceanic.prop.ntv, family=binomial(link ="logit"))
summary(model.oc.prop)

#Correlogram to test distance of spatial autocorrelation
correlogram(model.oc.prop, oceanic.prop.ntv, "figures/M6_ntvocprop_correlogram.jpg")

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(model.oc.prop,oceanic.prop.ntv,2000)
model.oc.prop.rac <- glm(cbind(nfix,nfixno)~abs.lat +area +dist+precipitation+temperature + area:dist+rac, data = oceanic.prop.ntv, family = binomial(link ="logit")) #has to be exactly the same as the model but with +rac
summary(model.oc.prop.rac)


# convert the model output into a df for estimates plot
coef_dataM6 <- tidy(model.oc.prop.rac, conf.int = TRUE)


#check variance inflation factor (should be below 5 for all variables)
vif(model.oc.prop.rac)

#variable importance using Rsquared and partial Rsquared

#check variable importance
library(glmtoolbox)#doesnt have a partial R2 function
adjR2(model.oc.prop.rac)

rsq(model.oc.prop.rac)
rsq.partial(model.oc.prop.rac, adj=FALSE)

####using manually coded function for Partial R square and Rsquared

formula<- (cbind(nfix,nfixno))~abs.lat+area+dist+precipitation+temperature+area:dist+rac
dat<-oceanic.prop.ntv
resp<- cbind(dat$nfix,dat$nfixno)

#calculate R square of model
R2<- r2_glm(resp, formula, dat)

#for a single variable
variable <- "area"
p2<- partial_r2_glm(resp, formula, dat, variable)

#trying partial R2 with loop for all variables
var_list<- c("abs.lat","area","dist","precipitation","temperature","area:dist","rac")
pr2.oc.prop.ntv<- partial_r2_glm_list(resp,formula,dat,var_list)


#check model assumptions
#check dispersion
disp_check(model.oc.prop.rac,oceanic.prop.ntv)

simulationOutput <- simulateResiduals(fittedModel = model.oc.prop.rac.full, plot = F)
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

#var.plot<- ggeffect(model.oc.prop.rac, terms=c("temperature"), type="re")
#temp.ntv.plot<- plot(var.plot, colors="darkseagreen3")

#Saving the plot as a png
#png("figures/temp_oc_prop_ntv_green.jpg", width=10, height= 10, units='in', res=300)
#temp.ntv.plot
#dev.off()


var.plot<- ggeffect(model.oc.prop.rac, terms=c("precipitation"), type="re")
precip.ntv.plot<- plot(var.plot, colors="darkseagreen3")

#Saving the plot as a png
png("figures/M6_precip_oc_prop_ntv.jpg", width=10, height= 10, units='in', res=300)
precip.ntv.plot
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
png("figures/M6_areadist_oc_prop_ntv.jpg", width=10, height= 10, units='in', res=300)
areadist.ntv.plot
dev.off()


#Figure native area:distance interaction

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
png("figures/M6_areadist_interaction.jpg", width=10, height= 10, units='in', res=300)
areadist_oc_pres_ntv_plot
dev.off()


#####FIGURES: Plot of Estimates #####

###Presence Models M3 and M5

#add distinction between models
coef_dataM3$model <- "Model 3"
coef_dataM5$model <- "Model 5"

#combine estimates of both models into one df
combined_data <- rbind(coef_dataM3, coef_dataM5)
# Specify the order of variables
combined_data$term <- factor(combined_data$term, levels = c("rac","temperature", "precipitation","abs.lat","area:dist","dist","area", "(Intercept)")) 
combined_data<- combined_data%>% filter(!term=="(Intercept)")#filter out intercept because variables are scaled 

# Create the combined plot
oc_pres_est_plot<-
  ggplot(combined_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, size=1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  scale_color_manual(values = c("Model 3" = "coral", "Model 5" = "darkseagreen3")) +  
  labs(title = "Coefficient Estimates for Presence on oceanic islands (M3 & M5)") +
  theme_minimal() +
  coord_flip()+
  theme(axis.text.x = element_text(size = 14))

oc_pres_est_plot

#Saving the plot as a png
png("figures/M3M5_estimates_oc_pres.jpg", width=10, height= 9, units='in', res=300)
oc_pres_est_plot
dev.off()


###Proportion Models M4 and M6

#add distinction between models
coef_dataM4$model <- "Model 4"
coef_dataM6$model <- "Model 6"

#combine estimates of both models into one df
combined_data <- rbind(coef_dataM4, coef_dataM6)
# Specify the order of variables
combined_data$term <- factor(combined_data$term, levels = c("rac","temperature", "precipitation","abs.lat","area:dist","dist","area", "(Intercept)")) 
combined_data<- combined_data%>% filter(!term=="(Intercept)")
  
# Create the combined plot
oc_prop_est_plot<-
  ggplot(combined_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, size=1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  scale_color_manual(values = c("Model 4" = "coral", "Model 6" = "darkseagreen3")) +  
  labs(title = "Coefficient Estimates for Proportion on oceanic islands (M4 & M6)") +
  theme_minimal() +
  coord_flip()+
  theme(axis.text.x = element_text(size = 14))

oc_prop_est_plot

#Saving the plot as a png
png("figures/M4M6_estimates_oc_prop.jpg", width=10, height= 9, units='in', res=300)
oc_prop_est_plot
dev.off()


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

