#################################################################################################################################
############################################   -----  DATA PREPARATION   -----   ############################################

## load data
data2018 <- read.csv("data/data2018_update.csv", h=T,sep=",")
head(data2018)
summary(data2018)
data2018 <- plyr::rename(data2018,c("district.x"="district"))

## add a column with the proportion of mutiled individuals
data2018$prop_mutil <- data2018$nbindiv_mutil/data2018$nbindiv_checked

## add a column with individual indentification ('ID')
data2018$ID <- c(1:nrow(data2018))

## save new table
write.csv(data2018,"data/data2018_ready.csv", row.names = F)


#################################################################################################################################
############################################   -----  PRELIMINARY ANALYSES  -----   ############################################

data.tab <- read.csv("data/data2018_ready.csv", h=T)
drop_var <- c("extleft_fing","medleft_fing","intleft_fing","backleft_fing","extright_fing","medright_fing","intright_fing","backright_fing","left_tot","right_tot")
data.tab <- data.tab[,!colnames(data.tab) %in% drop_var]
head(data.tab)

## check collinearity
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.5/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor*2)
  # text(0.5, 0.5, txt, cex = cex.cor * r)
}

## scatterplots - collinearity among variables
pdf(file="figs/collinearity_explvar.pdf",width=20,height=12)
pairs(~ humpop_dens + human_flow + hrs_openmarkets + bakeries_dens + hairdres_dens + greenspace_dens + noise_pollut + airqual_pca1 + waste_dens, pch=16,lower.panel=panel.smooth,upper.panel=panel.cor,cex=1,bg="light blue", cex.labels = 1.8,cex.axis=2,font.labels=2, data=data.tab)
dev.off()

## VIF analysis
library(usdm)
list_var <- c("humpop_dens","human_flow","hrs_openmarkets","bakeries_dens","hairdres_dens","greenspace_dens", "noise_pollut","airqual_pca1","waste_dens")
subdata <- data.tab[,colnames(data.tab) %in% list_var]
vif.out <- vif(subdata)
vif.out
write.csv(vif.out,"outputs/vif_results.csv")
## >> no collinearity issues

#at the scale of triris sector (when analysing prop of individuals mutilated / triris)
subdata.triris <- data.tab[,c("district","sect_triris","humpop_dens","human_flow", "openmarkets","hrs_openmarkets","bakeries_dens","hairdres_dens","greenspace_dens", "noise_pollut","airqual_pca1","waste_dens")]
subdata.triris <- subdata.triris[!duplicated(subdata.triris$sect_triris),]
list_var <- c("humpop_dens","human_flow", "openmarkets","hrs_openmarkets","bakeries_dens","hairdres_dens","greenspace_dens", "noise_pollut","airqual_pca1","waste_dens")
subdata <- subdata.triris[,colnames(subdata.triris) %in% list_var]
subdata$openmarkets <- as.numeric(subdata$openmarkets)
vif.out <- vifcor(subdata)
vif.out


#################################################################################################################################
############################################   -----  STATISTICAL ANALYSES   -----   ############################################

library(glmmTMB)
library(rsq)
library(modEvA)
library(sjPlot)
source("Rcodes/utilities.R")

#load data
data.df <- read.csv("data/data2018_ready.csv", h=T)
data.df <- data.df[complete.cases(data.df),]
head(data.df)

#set dummy variables as factors
data.df$morph <- as.factor(data.df$morph)
data.df$openmarkets <- as.factor(data.df$openmarkets)
data.df$foodmarkets <- as.factor(data.df$foodmarkets)

#create new DF in which continuous variables are scaled (to be used with multivariate models)
std_var <- c("humpop_dens","human_flow","hairdres_dens","bakeries_dens","greenspace_dens","noise_pollut", "waste_dens","hrs_openmarkets","hrsarea_openmarkets")
data_std.df <- data.df
data_std.df[,colnames(data_std.df) %in% std_var] <- scale(data_std.df[,colnames(data_std.df) %in% std_var])  

#### i. Test effect of individual variability (phenotype) on the presence of mutilation ####

# ## using GLM
# glm_pheno1.mutil <- glm(mutilation ~ age , data = data.df, family = "binomial")
# summary(glm_pheno1.mutil)
# glm_pheno2.mutil <- glm(mutilation ~ colour_simplif, data = data.df, family = "binomial")
# summary(glm_pheno2.mutil)
# glm_pheno3.mutil <- glm(mutilation ~ morph, data = data.df, family = "binomial")
# summary(glm_pheno3.mutil)
# glm_pheno.mutil <- glm(mutilation ~ age-1 + colour_simplif-1 + morph-1, data = data.df, family = "binomial")
# summary(glm_pheno.mutil)
# Anov.pheno_mutil <- Anova(glm_pheno.mutil) #variance partitioning
# Anov.pheno_mutil
# rsq.n(glm_pheno.mutil) #Nagelekerke r-squared
# 
# varPart(A = rsq.n(glm_pheno1.mutil), B = rsq.n(glm_pheno2.mutil), C = rsq.n(glm_pheno3.mutil), AB = rsq.n(glm(mutilation ~ age + colour_simplif , data = data.df, family = "binomial")), BC = rsq.n(glm(mutilation ~ colour_simplif + morph, data = data.df, family = "binomial")), AC = rsq.n(glm(mutilation ~ age + morph , data = data.df, family = "binomial")), ABC = rsq.n(glm_pheno.mutil), A.name = "Age", B.name = "Colour (simplif)", C.name = "Morph", main = "Phenotype Variability", plot.unexpl = TRUE, plot.digits = 5)

## using GLMM (sampling block effect treated as random)

#test random intercept (block effect)
mdl0.mutil <- glm(mutilation ~ 1 , data = data.df, family = binomial)
mdl01.mutil <- glmmTMB::glmmTMB(mutilation ~ 1 +(1|ID), data = data.df, family = binomial)
mdl02.mutil <- glmmTMB::glmmTMB(mutilation ~ 1 +(1|district/sect_triris), data = data.df, family = binomial)
mdl03.mutil <- glmmTMB::glmmTMB(mutilation ~ 1 +(1|district/sect_triris) +(1|ID), data = data.df, family = binomial)

# #test explanatory variables separately and then together
# glmm_pheno1.mutil <- glmmTMB(mutilation ~ age + (1|district/sect_triris), data = data.df, family = binomial)
# summary(glmm_pheno1.mutil) 
# glmm_pheno2.mutil <- glmmTMB(mutilation ~ colour_simplif + (1|district/sect_triris), data = data.df, family = binomial)
# summary(glmm_pheno2.mutil) 
# glmm_pheno3.mutil <- glmmTMB(mutilation ~ morph + (1|district/sect_triris), data = data.df, family = binomial)
# summary(glmm_pheno3.mutil) 
# glmm_pheno4.mutil <- glmmTMB(mutilation ~ age + colour_simplif + (1|district/sect_triris), data = data.df, family = binomial)
# summary(glmm_pheno4.mutil) 
# glmm_pheno.mutil <- glmmTMB(mutilation ~ age + colour_simplif + morph + (1|district/sect_triris), data = data.df, family = binomial)
# summary(glmm_pheno.mutil)
# #likelihood ratio test (using anova)
# anova(glmm_pheno4.mutil, glmm_pheno.mutil) #morph should not be included in further analyses
# anova(glmm_pheno1.mutil, glmm_pheno4.mutil); anova(glmm_pheno2.mutil, glmm_pheno4.mutil) #age and colour should be further taken into account


### Based on previous analyses on a separated dataset, Fred found an age effect on the mutilation occurrence (signif. diff. between juveniles and adults; juveniles still have  all their toe), but no effect of color morph or plumage menalism.
### for further analyses, we consider only observations on adults.
adult.df <- subset(data.df, age == "A")
adult_std.df <- subset(data_std.df, age == "A")

#### ii. Test effect of the environment on the presence of mutilation / the amount of mutilations ####
## we first run univariate models, to test each environmental predictor separately and potential non linear effect.
## we then, run multivariate models. AIC is used to confirm the integration of non-linear terms in the model.


## Response variable: PRESENCE of mutilations ('mutilation') ##
out_mutilation.df <- matrix(NA,nrow=0,ncol=5)
alloutlinear_mutilation.df <- matrix(NA,nrow=0,ncol=5)
alloutquadrat_mutilation.df <- matrix(NA,nrow=0,ncol=5)

##effect of human density
glmm_env1.mutil <- glmmTMB(mutilation ~ humpop_dens + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env1.mutil)
sjstats::r2(glmm_env1.mutil) ##   Marginal R2: 0.001; Conditional R2: 0.036
glmm_env11.mutil <- glmmTMB(mutilation ~ humpop_dens + I(humpop_dens^2) + (1|district/sect_triris), data = adult_std.df, family = binomial) ## need to use scaled variable otherwise it doesn't converge...
summary(glmm_env11.mutil) #test done with polynomial terms as well, but no effect either
#save significant results (or by default, linear effect model)
out_mutilation.df <- rbind(out_mutilation.df,c(rownames(summary(glmm_env1.mutil)$coefficients$cond)[2],summary(glmm_env1.mutil)$coefficients$cond[2,]))
# out_mutilation.df <- rbind(out_mutilation.df,summary(glmm_env1.mutil)$coefficients$cond)
#save all outputs
alloutlinear_mutilation.df <- rbind(alloutlinear_mutilation.df,c(rownames(summary(glmm_env1.mutil)$coefficients$cond)[2],summary(glmm_env1.mutil)$coefficients$cond[2,]))
alloutquadrat_mutilation.df <- rbind(alloutquadrat_mutilation.df,cbind(rownames(summary(glmm_env11.mutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env11.mutil)$coefficients$cond[2:3,])))

##effect of flow of passers
glmm_env2.mutil <- glmmTMB(mutilation ~ human_flow + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env2.mutil) 
glmm_env21.mutil <- glmmTMB(mutilation ~ human_flow + I(human_flow^2) + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env21.mutil) #test done with polynomial terms as well, but no effect either
#save significant results (or by default, linear effect model)
out_mutilation.df <- rbind(out_mutilation.df,c(rownames(summary(glmm_env2.mutil)$coefficients$cond)[2],summary(glmm_env2.mutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_mutilation.df <- rbind(alloutlinear_mutilation.df,c(rownames(summary(glmm_env2.mutil)$coefficients$cond)[2],summary(glmm_env2.mutil)$coefficients$cond[2,]))
alloutquadrat_mutilation.df <- rbind(alloutquadrat_mutilation.df,cbind(rownames(summary(glmm_env21.mutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env21.mutil)$coefficients$cond[2:3,])))

##effect of open markets / food markets (open+closed+organic markets)
glmm_env3.mutil <- glmmTMB(mutilation ~ openmarkets + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env3.mutil)
#save significant results (or by default, linear effect model)
out_mutilation.df <- rbind(out_mutilation.df,c(rownames(summary(glmm_env3.mutil)$coefficients$cond)[2],summary(glmm_env3.mutil)$coefficients$cond[2,]))
# glmm_env3.mutil <- glmmTMB(mutilation ~ foodmarkets + (1|district/sect_triris), data = adult_std.df, family = binomial)
# summary(glmm_env3.mutil)
# out_mutilation.df <- rbind(out_mutilation.df,c(rownames(summary(glmm_env3.mutil)$coefficients$cond)[2],summary(glmm_env3.mutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_mutilation.df <- rbind(alloutlinear_mutilation.df,c(rownames(summary(glmm_env3.mutil)$coefficients$cond)[2],summary(glmm_env3.mutil)$coefficients$cond[2,]))

#effect of hours of opening food markets (proportional to the block area)
glmm_env31.mutil <- glmmTMB(mutilation ~ hrs_openmarkets + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env31.mutil)
glmm_env32.mutil <- glmmTMB(mutilation ~ hrs_openmarkets + I(hrs_openmarkets^2) + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env32.mutil)  #test done with polynomial terms as well, but no effect either
#save significant results (or by default, linear effect model)
out_mutilation.df <- rbind(out_mutilation.df,c(rownames(summary(glmm_env31.mutil)$coefficients$cond)[2],summary(glmm_env31.mutil)$coefficients$cond[2,]))
# glmm_env3.mutil <- glmmTMB(mutilation ~ hrsarea_totopenhrs + (1|district/sect_triris), data = adult_std.df, family = binomial)
# summary(glmm_env3.mutil)
# out_mutilation.df <- rbind(out_mutilation.df,c(rownames(summary(glmm_env3.mutil)$coefficients$cond)[2],summary(glmm_env3.mutil)$coefficients$cond[2,]))
# glmm_env3.mutil <- glmmTMB(mutilation ~  poly(hrsarea_totopenhrs,2) + (1|district/sect_triris), data = adult_std.df, family = binomial)
# summary(glmm_env3.mutil)  #test done with polynomial terms as well, but no effect either
#save all outputs
alloutlinear_mutilation.df <- rbind(alloutlinear_mutilation.df,c(rownames(summary(glmm_env31.mutil)$coefficients$cond)[2],summary(glmm_env31.mutil)$coefficients$cond[2,]))
alloutquadrat_mutilation.df <- rbind(alloutquadrat_mutilation.df,cbind(rownames(summary(glmm_env32.mutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env32.mutil)$coefficients$cond[2:3,])))

#effect of bakeries density
glmm_env4.mutil <- glmmTMB(mutilation ~  bakeries_dens + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env4.mutil)  
glmm_env41.mutil <- glmmTMB(mutilation ~  bakeries_dens + I(bakeries_dens^2) + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env41.mutil)  #test done with polynomial terms as well, but no effect either
#save significant results (or by default, linear effect model)
out_mutilation.df <- rbind(out_mutilation.df,c(rownames(summary(glmm_env4.mutil)$coefficients$cond)[2],summary(glmm_env4.mutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_mutilation.df <- rbind(alloutlinear_mutilation.df,c(rownames(summary(glmm_env4.mutil)$coefficients$cond)[2],summary(glmm_env4.mutil)$coefficients$cond[2,]))
alloutquadrat_mutilation.df <- rbind(alloutquadrat_mutilation.df,cbind(rownames(summary(glmm_env41.mutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env41.mutil)$coefficients$cond[2:3,])))

#effect of greenspace density
glmm_env5.mutil <- glmmTMB(mutilation ~ greenspace_dens + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env5.mutil)
glmm_env51.mutil <- glmmTMB(mutilation ~ greenspace_dens+ I(greenspace_dens^2) + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env51.mutil) #test done with polynomial terms as well, but no effect either
#save significant results (or by default, linear effect model)
out_mutilation.df <- rbind(out_mutilation.df,c(rownames(summary(glmm_env5.mutil)$coefficients$cond)[2],summary(glmm_env5.mutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_mutilation.df <- rbind(alloutlinear_mutilation.df,c(rownames(summary(glmm_env5.mutil)$coefficients$cond)[2],summary(glmm_env5.mutil)$coefficients$cond[2,]))
alloutquadrat_mutilation.df <- rbind(alloutquadrat_mutilation.df,cbind(rownames(summary(glmm_env51.mutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env51.mutil)$coefficients$cond[2:3,])))

#effect of noise pollution
glmm_env6.mutil <- glmmTMB(mutilation ~ noise_pollut + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env6.mutil) #SIGNIFICANT! (linear effect only) [noise_pollut  0.21384    0.08117   2.635  0.00842 **]
glmm_env61.mutil <- glmmTMB(mutilation ~ noise_pollut+ I(noise_pollut^2) + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env61.mutil) #test done with polynomial terms as well, but linear effect only
#save significant results (or by default, linear effect model)
out_mutilation.df <- rbind(out_mutilation.df,c(rownames(summary(glmm_env6.mutil)$coefficients$cond)[2],summary(glmm_env6.mutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_mutilation.df <- rbind(alloutlinear_mutilation.df,c(rownames(summary(glmm_env6.mutil)$coefficients$cond)[2],summary(glmm_env6.mutil)$coefficients$cond[2,]))
alloutquadrat_mutilation.df <- rbind(alloutquadrat_mutilation.df,cbind(rownames(summary(glmm_env61.mutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env61.mutil)$coefficients$cond[2:3,])))
#plot effects
noisepollut.mutil <- glmmTMB(mutilation ~ noise_pollut + (1|district/sect_triris), data = adult.df, family = binomial)
pdf(file="figs/effect_noisepollution_probmutilation.pdf",width=8,height=8)
plot_model(noisepollut.mutil,type="pred", pred.type = "fe",terms = "noise_pollut",allow.new.levels=TRUE, axis.lim=list(c(min(adult.df$noise_pollut),max(adult.df$noise_pollut)),c(0,1)))
dev.off()


#effect of gas pollution (air quality index)
glmm_env7.mutil <- glmmTMB(mutilation ~ airqual_pca1 + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env7.mutil) #SIGNIFICANT! [airqual_pca1  0.12006    0.05034   2.385   0.0171 * ]
glmm_env71.mutil <- glmmTMB(mutilation ~ airqual_pca1+ I(airqual_pca1^2) + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env71.mutil) #test done with polynomial terms as well, but linear effect only
#save significant results (or by default, linear effect model)
out_mutilation.df <- rbind(out_mutilation.df,c(rownames(summary(glmm_env7.mutil)$coefficients$cond)[2],summary(glmm_env7.mutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_mutilation.df <- rbind(alloutlinear_mutilation.df,c(rownames(summary(glmm_env7.mutil)$coefficients$cond)[2],summary(glmm_env7.mutil)$coefficients$cond[2,]))
alloutquadrat_mutilation.df <- rbind(alloutquadrat_mutilation.df,cbind(rownames(summary(glmm_env71.mutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env71.mutil)$coefficients$cond[2:3,])))
#plot effects
airqual.mutil <- glmmTMB(mutilation ~ airqual_pca1 + (1|district/sect_triris), data = adult.df, family = binomial)
pdf(file="figs/effect_airquality_probmutilation.pdf",width=8,height=8)
plot_model(airqual.mutil,type="pred", pred.type = "fe",terms = "airqual_pca1",allow.new.levels=TRUE, axis.lim=list(c(min(adult.df$airqual_pca1),max(adult.df$airqual_pca1)),c(0,1)))
dev.off()

#effect of waste amount and waste density
glmm_env8.mutil <- glmmTMB(mutilation ~ waste_dens + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env8.mutil)
glmm_env81.mutil <- glmmTMB(mutilation ~ waste_dens+ I(waste_dens^2) + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env81.mutil)
#save significant results (or by default, linear effect model)
out_mutilation.df <- rbind(out_mutilation.df,c(rownames(summary(glmm_env8.mutil)$coefficients$cond)[2],summary(glmm_env8.mutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_mutilation.df <- rbind(alloutlinear_mutilation.df,c(rownames(summary(glmm_env8.mutil)$coefficients$cond)[2],summary(glmm_env8.mutil)$coefficients$cond[2,]))
alloutquadrat_mutilation.df <- rbind(alloutquadrat_mutilation.df,cbind(rownames(summary(glmm_env81.mutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env81.mutil)$coefficients$cond[2:3,])))

# ##effect of fabric manufucture density
# glmm_env9.mutil <- glmmTMB(mutilation ~  fabricmanuf_dens + (1|district/sect_triris), data = adult.df, family = binomial)
# summary(glmm_env9.mutil)  
# out_mutilation.df <- rbind(out_mutilation.df,c(rownames(summary(glmm_env9.mutil)$coefficients$cond)[2],summary(glmm_env9.mutil)$coefficients$cond[2,]))
# glmm_env9.mutil <- glmmTMB(mutilation ~  poly(fabricmanuf_dens,2) + (1|district/sect_triris), data = adult.df, family = binomial)
# summary(glmm_env9.mutil)  #polynomial term (negatively??) significant (diff AIC<2)
# out_mutilation.df <- rbind(out_mutilation.df,cbind(rownames(summary(glmm_env9.mutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env9.mutil)$coefficients$cond[2:3,])))

##effect of hairdressers density
glmm_env10.mutil <- glmmTMB(mutilation ~  hairdres_dens + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env10.mutil)  #close from significant
glmm_env101.mutil <- glmmTMB(mutilation ~  hairdres_dens + I(hairdres_dens^2) + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env101.mutil)  #close from significant - quadratic effect 
#save significant results (or by default, linear effect model)
out_mutilation.df <- rbind(out_mutilation.df,c(rownames(summary(glmm_env10.mutil)$coefficients$cond)[2],summary(glmm_env10.mutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_mutilation.df <- rbind(alloutlinear_mutilation.df,c(rownames(summary(glmm_env10.mutil)$coefficients$cond)[2],summary(glmm_env10.mutil)$coefficients$cond[2,]))
alloutquadrat_mutilation.df <- rbind(alloutquadrat_mutilation.df,cbind(rownames(summary(glmm_env101.mutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env101.mutil)$coefficients$cond[2:3,])))
#plot effects
hairdres.mutil <- glmmTMB(mutilation ~ hairdres_dens + (1|district/sect_triris), data = adult.df, family = binomial)
pdf(file="figs/effect_hairdressers_lineareffect_probmutilation.pdf",width=8,height=8)
plot_model(hairdres.mutil,type="pred", pred.type = "fe",terms = "hairdres_dens",allow.new.levels=TRUE, axis.lim=list(c(min(adult.df$hairdres_dens),max(adult.df$hairdres_dens)),c(0,1)))
dev.off()
hairdres.mutil <- glmmTMB(mutilation ~ hairdres_dens + I(hairdres_dens^2)+ (1|district/sect_triris), data = adult.df, family = binomial)
pdf(file="figs/effect_hairdressers_quadeffect_probmutilation.pdf",width=8,height=8)
plot_model(hairdres.mutil,type="pred", pred.type = "fe",terms = "hairdres_dens",allow.new.levels=TRUE, axis.lim=list(c(min(adult.df$hairdres_dens),max(adult.df$hairdres_dens)),c(0,1)))
dev.off()


#save outputs of univariate models
colnames(out_mutilation.df) <- c("Variable","Estimate","SE","z-value","P-value")
write.csv(out_mutilation.df,"outputs/outputs_univariate_mdls_mutilation.csv")
colnames(alloutlinear_mutilation.df) <- c("Variable","Estimate","SE","z-value","P-value")
write.csv(alloutlinear_mutilation.df,"outputs/outputs_univariate_alllinearmdls_mutilation.csv")
colnames(alloutquadrat_mutilation.df) <- c("Variable","Estimate","SE","z-value","P-value")
write.csv(alloutquadrat_mutilation.df,"outputs/outputs_univariate_allquadratmdls_mutilation.csv")


##all effects in a single model (and accounting for effect of indiv variability) - use DF with scaled variables
glmm_env.mutil <- glmmTMB(mutilation ~ humpop_dens + human_flow + openmarkets + hrs_openmarkets + bakeries_dens + hairdres_dens + noise_pollut + greenspace_dens + waste_dens + airqual_pca1 + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env.mutil) #noise pollution and airquality still significant
sjstats::r2(glmm_env.mutil)
#look at collinearity
vif(glmm_env.mutil$frame[,c(2:3,5:11)])
#run new model after dropping waste_dens
glmm_env.mutil <- glmmTMB(mutilation ~ humpop_dens + human_flow + openmarkets + hrs_openmarkets + bakeries_dens + hairdres_dens + noise_pollut + greenspace_dens + airqual_pca1 + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env.mutil) #noise pollution and airquality still significant
sjstats::r2(glmm_env.mutil)
#look at collinearity
vif(glmm_env.mutil$frame[,c(2:3,5:10)])
#save outputs
out_multivar_mutil.df <- as.data.frame(summary(glmm_env.mutil)$coefficients$cond)
write.csv(out_multivar_mutil.df,"outputs/outputs_full_multivariate_mdl_mutilation.csv")


#reduced model with only significant variables or almost significant
glmm_env_red.mutil <- glmmTMB(mutilation ~ noise_pollut + airqual_pca1 + hairdres_dens + I(hairdres_dens^2) + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env_red.mutil) #noise pollution and air quality  still significant. hairdessers density close from significant
# aov(glmm_env.mutil)
out_multivar_mutil.df <- as.data.frame(summary(glmm_env_red.mutil)$coefficients$cond)
write.csv(out_multivar_mutil.df,"outputs/outputs_reduced_multivariate_mdl_mutilation.csv")
sjstats::r2(glmm_env_red.mutil)  ## Marginal R2: 0.031; Conditional R2: 0.037

#without hairdressers density
glmm_env_red2.mutil <- glmmTMB(mutilation ~ noise_pollut + airqual_pca1 + (1|district/sect_triris), data = adult_std.df, family = binomial)
summary(glmm_env_red2.mutil) #noise pollution and air quality  still significant. hairdessers density close from significant
# aov(glmm_env.mutil)
out_multivar_mutil.df <- as.data.frame(summary(glmm_env_red2.mutil)$coefficients$cond)
write.csv(out_multivar_mutil.df,"outputs/outputs_reduced2_multivariate_mdl_mutilation.csv")
sjstats::r2(glmm_env_red2.mutil)




########################################################
## Response variable: NUMBER of mutilations ('total') ##

#subset data - keep only records of mutilations
data_nbmutil.df <- subset(adult.df, mutilation == 1)
summary(data_nbmutil.df)
data_nbmutil_std.df <- subset(adult_std.df, mutilation == 1)
summary(data_nbmutil_std.df)

out_nbmutilation.df <- matrix(NA,nrow=0,ncol=5)
alloutlinear_nbmutilation.df <- matrix(NA,nrow=0,ncol=5)
alloutquadrat_nbmutilation.df <- matrix(NA,nrow=0,ncol=5)

#effect of human density
glmm_env1.nbmutil <- glmmTMB(total ~ humpop_dens + (1|district/sect_triris),  data = data_nbmutil_std.df, family = poisson)
summary(glmm_env1.nbmutil)
glmm_env11.nbmutil <- glmmTMB(total ~ humpop_dens + I(humpop_dens^2) + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env11.nbmutil) #SIGNIFICANT - polynomial term (diff AIC >2)
#save significant results (or by default, linear effect model)
out_nbmutilation.df <-  rbind(out_nbmutilation.df,cbind(rownames(summary(glmm_env11.nbmutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env11.nbmutil)$coefficients$cond[2:3,])))
#save all outputs
alloutlinear_nbmutilation.df <- rbind(alloutlinear_nbmutilation.df,c(rownames(summary(glmm_env1.nbmutil)$coefficients$cond)[2],summary(glmm_env1.nbmutil)$coefficients$cond[2,]))
alloutquadrat_nbmutilation.df <- rbind(alloutquadrat_nbmutilation.df,cbind(rownames(summary(glmm_env11.nbmutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env11.nbmutil)$coefficients$cond[2:3,])))
sjstats::r2(glmm_env11.nbmutil)
#plot effects
humpop.nbmutil <- glmmTMB(total ~ poly(humpop_dens,2)+ (1|district/sect_triris), na.action = na.omit, data = data_nbmutil.df, family = poisson)
pdf(file="figs/effect_humanpop_quadeffect_probmutilation.pdf",width=8,height=8)
plot_model(humpop.nbmutil, type="pred", pred.type = "fe", terms = "humpop_dens",allow.new.levels=TRUE)
dev.off()


#effect of flow of passers
glmm_env2.nbmutil <- glmmTMB(total ~ human_flow + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env2.nbmutil)
glmm_env21.nbmutil <- glmmTMB(total ~ human_flow + I(human_flow^2) + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env21.nbmutil)
#save significant results (or by default, linear effect model)
out_nbmutilation.df <- rbind(out_nbmutilation.df,c(rownames(summary(glmm_env2.nbmutil)$coefficients$cond)[2],summary(glmm_env2.nbmutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_nbmutilation.df <- rbind(alloutlinear_nbmutilation.df,c(rownames(summary(glmm_env2.nbmutil)$coefficients$cond)[2],summary(glmm_env2.nbmutil)$coefficients$cond[2,]))
alloutquadrat_nbmutilation.df <- rbind(alloutquadrat_nbmutilation.df,cbind(rownames(summary(glmm_env21.nbmutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env21.nbmutil)$coefficients$cond[2:3,])))

#effect of presence of open/food markets
glmm_env3.nbmutil <- glmmTMB(total ~ openmarkets  + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env3.nbmutil)
#save significant results (or by default, linear effect model)
out_nbmutilation.df <- rbind(out_nbmutilation.df,c(rownames(summary(glmm_env3.nbmutil)$coefficients$cond)[2],summary(glmm_env3.nbmutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_nbmutilation.df <- rbind(alloutlinear_nbmutilation.df,c(rownames(summary(glmm_env3.nbmutil)$coefficients$cond)[2],summary(glmm_env3.nbmutil)$coefficients$cond[2,]))

# glmm_env3.nbmutil <- glmmTMB(total ~ foodmarkets  + (1|district/sect_triris), data = data_nbmutil.df, family = poisson)
# summary(glmm_env3.nbmutil)
#effect of nb of opening hours of open/food markets
glmm_env31.nbmutil <- glmmTMB(total ~ hrs_openmarkets  + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env31.nbmutil)
glmm_env32.nbmutil <- glmmTMB(total ~ hrs_openmarkets+I(hrs_openmarkets^2)  + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env32.nbmutil)
#save significant results (or by default, linear effect model)
out_nbmutilation.df <- rbind(out_nbmutilation.df,c(rownames(summary(glmm_env31.nbmutil)$coefficients$cond)[2],summary(glmm_env31.nbmutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_nbmutilation.df <- rbind(alloutlinear_nbmutilation.df,c(rownames(summary(glmm_env31.nbmutil)$coefficients$cond)[2],summary(glmm_env31.nbmutil)$coefficients$cond[2,]))
alloutquadrat_nbmutilation.df <- rbind(alloutquadrat_nbmutilation.df,cbind(rownames(summary(glmm_env32.nbmutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env32.nbmutil)$coefficients$cond[2:3,])))


#effect of bakeries density
glmm_env4.nbmutil <- glmmTMB(total ~ bakeries_dens + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env4.nbmutil)
glmm_env41.nbmutil <- glmmTMB(total ~ bakeries_dens+I(bakeries_dens^2) + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env41.nbmutil)
#save significant results (or by default, linear effect model)
out_nbmutilation.df <- rbind(out_nbmutilation.df,c(rownames(summary(glmm_env4.nbmutil)$coefficients$cond)[2],summary(glmm_env4.nbmutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_nbmutilation.df <- rbind(alloutlinear_nbmutilation.df,c(rownames(summary(glmm_env4.nbmutil)$coefficients$cond)[2],summary(glmm_env4.nbmutil)$coefficients$cond[2,]))
alloutquadrat_nbmutilation.df <- rbind(alloutquadrat_nbmutilation.df,cbind(rownames(summary(glmm_env41.nbmutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env41.nbmutil)$coefficients$cond[2:3,])))

#effect of greenspace density
glmm_env5.nbmutil <- glmmTMB(total ~  greenspace_dens + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env5.nbmutil) #close from significant
glmm_env51.nbmutil <- glmmTMB(total ~ greenspace_dens+ I(greenspace_dens^2) + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env51.nbmutil)
#save significant results (or by default, linear effect model)
out_nbmutilation.df <- rbind(out_nbmutilation.df,c(rownames(summary(glmm_env5.nbmutil)$coefficients$cond)[2],summary(glmm_env5.nbmutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_nbmutilation.df <- rbind(alloutlinear_nbmutilation.df,c(rownames(summary(glmm_env5.nbmutil)$coefficients$cond)[2],summary(glmm_env5.nbmutil)$coefficients$cond[2,]))
alloutquadrat_nbmutilation.df <- rbind(alloutquadrat_nbmutilation.df,cbind(rownames(summary(glmm_env51.nbmutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env51.nbmutil)$coefficients$cond[2:3,])))

#effect of noise pollution
glmm_env6.nbmutil <- glmmTMB(total ~  noise_pollut + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env6.nbmutil)
glmm_env61.nbmutil <- glmmTMB(total ~  noise_pollut+I(noise_pollut^2) + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env61.nbmutil)
#save significant results (or by default, linear effect model)
out_nbmutilation.df <- rbind(out_nbmutilation.df,c(rownames(summary(glmm_env6.nbmutil)$coefficients$cond)[2],summary(glmm_env6.nbmutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_nbmutilation.df <- rbind(alloutlinear_nbmutilation.df,c(rownames(summary(glmm_env6.nbmutil)$coefficients$cond)[2],summary(glmm_env6.nbmutil)$coefficients$cond[2,]))
alloutquadrat_nbmutilation.df <- rbind(alloutquadrat_nbmutilation.df,cbind(rownames(summary(glmm_env61.nbmutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env61.nbmutil)$coefficients$cond[2:3,])))

#effect of air pollution (air quality index)
glmm_env7.nbmutil <- glmmTMB(total ~ airqual_pca1 + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env7.nbmutil) #close from significant
glmm_env71.nbmutil <- glmmTMB(total ~  airqual_pca1+I(airqual_pca1^2) + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env71.nbmutil)
#save significant results (or by default, linear effect model)
out_nbmutilation.df <- rbind(out_nbmutilation.df,c(rownames(summary(glmm_env7.nbmutil)$coefficients$cond)[2],summary(glmm_env7.nbmutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_nbmutilation.df <- rbind(alloutlinear_nbmutilation.df,c(rownames(summary(glmm_env7.nbmutil)$coefficients$cond)[2],summary(glmm_env7.nbmutil)$coefficients$cond[2,]))
alloutquadrat_nbmutilation.df <- rbind(alloutquadrat_nbmutilation.df,cbind(rownames(summary(glmm_env71.nbmutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env71.nbmutil)$coefficients$cond[2:3,])))

#effect of waste density
glmm_env8.nbmutil <- glmmTMB(total ~ waste_dens + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env8.nbmutil)
glmm_env81.nbmutil <- glmmTMB(total ~ waste_dens + I(waste_dens^2) + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env81.nbmutil)
#save significant results (or by default, linear effect model)
out_nbmutilation.df <- rbind(out_nbmutilation.df,c(rownames(summary(glmm_env8.nbmutil)$coefficients$cond)[2],summary(glmm_env8.nbmutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_nbmutilation.df <- rbind(alloutlinear_nbmutilation.df,c(rownames(summary(glmm_env8.nbmutil)$coefficients$cond)[2],summary(glmm_env8.nbmutil)$coefficients$cond[2,]))
alloutquadrat_nbmutilation.df <- rbind(alloutquadrat_nbmutilation.df,cbind(rownames(summary(glmm_env81.nbmutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env81.nbmutil)$coefficients$cond[2:3,])))


##effect of hairdressers density
glmm_env10.nbmutil <- glmmTMB(total ~  hairdres_dens + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env10.nbmutil)  
glmm_env101.nbmutil <- glmmTMB(total ~ hairdres_dens+ I(hairdres_dens^2) + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env101.nbmutil) 
#save significant results (or by default, linear effect model)
out_nbmutilation.df <- rbind(out_nbmutilation.df,c(rownames(summary(glmm_env10.nbmutil)$coefficients$cond)[2],summary(glmm_env10.nbmutil)$coefficients$cond[2,]))
#save all outputs
alloutlinear_nbmutilation.df <- rbind(alloutlinear_nbmutilation.df,c(rownames(summary(glmm_env10.nbmutil)$coefficients$cond)[2],summary(glmm_env10.nbmutil)$coefficients$cond[2,]))
alloutquadrat_nbmutilation.df <- rbind(alloutquadrat_nbmutilation.df,cbind(rownames(summary(glmm_env101.nbmutil)$coefficients$cond)[2:3],as.matrix(summary(glmm_env101.nbmutil)$coefficients$cond[2:3,])))


#save outputs of univariate models
colnames(out_nbmutilation.df) <- c("Variable","Estimate","SE","z-value","P-value")
write.csv(out_nbmutilation.df,"outputs/outputs_univariate_mdls_nbmutilation.csv")
colnames(alloutlinear_nbmutilation.df) <- c("Variable","Estimate","SE","z-value","P-value")
write.csv(alloutlinear_nbmutilation.df,"outputs/outputs_univariate_alllinearmdls_nbmutilation.csv")
colnames(alloutquadrat_nbmutilation.df) <- c("Variable","Estimate","SE","z-value","P-value")
write.csv(alloutquadrat_nbmutilation.df,"outputs/outputs_univariate_allquadratmdls_nbmutilation.csv")


#all effects in a single model (no age effect - not relevant)
#colour_simplif was then removed because of the disproportionate number of individual of each morph
glmm_env.nbmutil <- glmmTMB(total ~ humpop_dens + I(humpop_dens^2) + human_flow + openmarkets  + hrs_openmarkets + bakeries_dens + greenspace_dens + hairdres_dens + noise_pollut + waste_dens + airqual_pca1 + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env.nbmutil)
sjstats::r2(glmm_env.nbmutil) #Conditional R2: 0.061
my_rsq(glmm_env.nbmutil)
#look at collinearity
vif(glmm_env.nbmutil$frame[,c(2:4,6:12)])
#run model after removing waste density (correlated to humanpop_dens)
glmm_env1.nbmutil <- glmmTMB(total ~ humpop_dens + I(humpop_dens^2) + human_flow + openmarkets  + hrs_openmarkets + bakeries_dens + greenspace_dens + hairdres_dens + noise_pollut + airqual_pca1 + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env1.nbmutil)
sjstats::r2(glmm_env1.nbmutil) ##   Marginal R2: 0.061; Conditional R2: 0.061
my_rsq(glmm_env1.nbmutil) 
#look at collinearity
vifstep(glmm_env1.nbmutil$frame[,c(2:4,6:11)],th=3)
#run model after removing hairdressers density (correlated to bakeries density???)
glmm_env2.nbmutil <- glmmTMB(total ~ humpop_dens + I(humpop_dens^2) + human_flow + openmarkets + hrs_openmarkets + bakeries_dens + greenspace_dens + noise_pollut + airqual_pca1 + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env2.nbmutil)
sjstats::r2(glmm_env2.nbmutil) ##   Marginal R2: 0.061; Conditional R2: 0.061
my_rsq(glmm_env2.nbmutil) 
#look at collinearity
vifstep(glmm_env2.nbmutil$frame[,c(2:4,6:10)])
#save outputs
out_multivar_nbmutil.df <- as.data.frame(summary(glmm_env2.nbmutil)$coefficients$cond)
colnames(out_multivar_nbmutil.df) <- c("Estimate","SE","z-value","P-value")
write.csv(out_multivar_nbmutil.df,"outputs/outputs_full_multivariate_mdl_nbmutilation.csv")


## reduced model with only significant variables
glmm_env3.nbmutil <- glmmTMB::glmmTMB(total ~  humpop_dens + I(humpop_dens^2) + greenspace_dens +  airqual_pca1 + (1|district/sect_triris), data = data_nbmutil_std.df, family = poisson)
summary(glmm_env3.nbmutil)
sjstats::r2(glmm_env3.nbmutil) ##   Marginal R2: 0.059; Conditional R2: 0.059
out_multivar_nbmutil.df <- as.data.frame(summary(glmm_env3.nbmutil)$coefficients$cond)
colnames(out_multivar_nbmutil.df) <- c("Estimate","SE","z-value","P-value")
write.csv(out_multivar_nbmutil.df,"outputs/outputs_reduced_multivariate_mdl_nbmutilation.csv")




# ## Response variable: proportion of individuals with mutilations (per triris sector) - among adults only ##
# 
# #subset data
# #juveniles should be removed (they usually have no mutilation...)
# newtab <- adult.df
# newtab$indiv <- 1
# mutil.df <- aggregate.data.frame(newtab[,c(17,41)],list(newtab$sect_triris),sum)
# colnames(mutil.df) <- c("sect_triris", "nb_mutil","nb_checked")
# mutil.df$prop_mutil <- mutil.df$nb_mutil/mutil.df$nb_checked
# rownames(mutil.df) <- mutil.df$sect_triris
# data_propmutil.df <- adult.df[,c("district","sect_triris","humpop_dens","human_flow","openmarkets","bakeries_dens","greenspace_dens","noise_pollut","GHG","PM10","waste_dens")]
# data_propmutil.df <- data_propmutil.df[!duplicated(data_propmutil.df$sect_triris),]
# data_propmutil.df <- data.frame(data_propmutil.df, prop_mutil = mutil.df[as.character(data_propmutil.df$sect_triris),4])
# 
# #using a beta regression
# library(betareg)
# 
# betar_env1.propmutil <- betareg(prop_mutil ~  humpop_dens, data = data_propmutil.df)
# summary(betar_env1.propmutil)
# AIC(betar_env1.propmutil)
# 
# betar_env2.propmutil <- betareg(prop_mutil ~ poly(human_flow,3), data = data_propmutil.df)
# summary(betar_env2.propmutil)
# AIC(betar_env2.propmutil)
# 
# betar_env3.propmutil <- betareg(prop_mutil ~ openmarkets, data = data_propmutil.df)
# summary(betar_env3.propmutil)
# AIC(betar_env3.propmutil)
# 
# betar_env4.propmutil <- betareg(prop_mutil ~ poly(bakeries_dens,2), data = data_propmutil.df)
# summary(betar_env4.propmutil)
# AIC(betar_env4.propmutil)
# 
# betar_env5.propmutil <- betareg(prop_mutil ~ greenspace_dens, data = data_propmutil.df)
# summary(betar_env5.propmutil)
# AIC(betar_env5.propmutil)
# 
# betar_env6.propmutil <- betareg(prop_mutil ~ noise_pollut, data = data_propmutil.df)
# summary(betar_env6.propmutil)
# AIC(betar_env6.propmutil)
# 
# betar_env7.propmutil <- betareg(prop_mutil ~ GHG, data = data_propmutil.df)
# summary(betar_env7.propmutil)
# AIC(betar_env7.propmutil)
# 
# betar_env8.propmutil <- betareg(prop_mutil ~ waste_dens, data = data_propmutil.df)
# summary(betar_env8.propmutil)
# AIC(betar_env8.propmutil)
# 
# betar_env.propmutil <- betareg(prop_mutil ~ humpop_dens + poly(human_flow,2) + openmarkets + poly(bakeries_dens,2) + greenspace_dens + noise_pollut + GHG + waste_dens, data = data_propmutil.df)
# summary(betar_env.propmutil)
# AIC(betar_env.propmutil)
# 
# 
# #using a logistic regression with district as a random effect
# glmm_env.propmutil <- glmmTMB(prop_mutil ~ humpop_dens + poly(human_flow,2) + openmarkets + poly(bakeries_dens,2) + greenspace_dens + noise_pollut + GHG + waste_dens + (1|district), data = data_propmutil.df, family=list(family="beta",link="logit"))
# summary(glmm_env.propmutil)
