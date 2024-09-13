#Beeby et al. "Climate and nutrition drive gut microbiome variation in a fruit-specialist primate."
#
#Nina Beeby (ninabeeby@gmail.com) and Timothy H Webster (timothy.h.webster@utah.edu), corresponding authors
#
#
#START OF SCRIPT
#Outline code for analyses
#
#Packages required
library(readr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(vegan)
library(lme4)
library(MASS)
library(fitdistrplus)
library(stats4)
library(zoo)
library(car)
library(nlme)
library(lmerTest)
library(padr)
library(purrr)
library(plotrix)
library(MetBrewer)
library(forcats)
library(MuMIn)
library(MCMCglmm)
library(ggforce)
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#FIGURE 1
#
#Select required columns
Diversity_Indices <- GM_Climate_Nutrition_Diversity %>%
  dplyr::select(monthnum, shannon_entropy, faith_pd)
#Convert to long format
Diversity_Indices <- Diversity_Indices %>%
  pivot_longer(!monthnum, names_to = "DiversityIndex", values_to = "count")
#Plot
MonthlyIndices <- Diversity_Indices %>%
  ggplot(., aes(x=factor(monthnum, level=c('1','2','3','4','5','6','7','8','9','10','11','12')), y=count)) +
  geom_point(aes(color=DiversityIndex)) +
  geom_smooth(aes(group=DiversityIndex, color=DiversityIndex), method='loess', formula=y~x) +
  xlab("Month") + ylab("Diversity Value") + labs(color="Diversity Index") + theme_classic() +
  scale_color_manual(values = c("shannon_entropy" = "#17154f", "faith_pd" = "#b0799a"))
MonthlyIndices
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#TABLE 1 prep
#Test for linearity (for each of shannon and faith, repeat for each variable)
shannon_diet_ndf <-ggplot(GM_Climate_Nutrition, aes(x=shannon_entropy, y=intake_ndf_30_day)) +
  geom_point() + geom_smooth() + xlab("Shannon") + ylab("NDF Intake") + stat_cor(method = "pearson")
#Test nutrient intakes for multicollinearity (removed g_total_food and intake_ap_30 day for VIFs < 5)
nutrients <- GM_Climate_Nutrition %>%
  dplyr::select(shannon_entropy, temp_30_day, rainfall_30_day, intake_ndf_30_day, intake_fat_30_day, intake_tnc_30_day, prop_fruit_30_day, diet_div_30_day)
model_nutrients <- lm(shannon_entropy ~ ., data=nutrients)  
summary(model_nutrients)
vif(model_nutrients) 
#Test random effects for multicollinearity (no issues)
nutrients_random <- GM_Climate_Nutrition %>%
  dplyr::select(shannon_entropy, sex, library, batch)
model_nutrients_random <- lm(shannon_entropy ~ ., data=nutrients_random)
summary(model_nutrients_random)
#
#TABLE 1
#Models predicting Diversity Indices (plus check residuals for normality)
#Shannon
Shannon_model = lmer(formula = shannon_entropy ~ temp_30_day + rainfall_30_day + intake_ndf_30_day + intake_tnc_30_day + intake_fat_30_day + prop_fruit_30_day + diet_div_30_day + (1|batch),
                     data=GM_Climate_Nutrition_Diversity, na.action=na.fail)
summary(Shannon_model)
res <- resid(Shannon_model)
plot(fitted(Shannon_model), res)
abline(0,0)
qqnorm(res)
qqline(res)
plot(density(res))
Shannon_Models <- dredge(Shannon_model, rank = "AICc")
View(Shannon_Models)
summary(dredge(Shannon_model)) #model selection 
summary(model.avg(get.models(Shannon_Models, subset = TRUE))) #model averaging
#Faith
Faith_model = lmer(formula = faith_pd ~ temp_30_day + rainfall_30_day + intake_ndf_30_day + intake_tnc_30_day + intake_fat_30_day + prop_fruit_30_day + diet_div_30_day + (1|batch),
                   data=GM_Climate_Nutrition_Diversity, na.action=na.fail)
summary(Faith_model)
res1 <- resid(Faith_model)
plot(fitted(Faith_model), res1)
abline(0,0)
qqnorm(res1)
qqline(res1)
plot(density(res1))
Faith_Models <- dredge(Faith_model, rank = "AICc") #model selection
View(Faith_Models)
summary(model.avg(get.models(Faith_Models, subset = TRUE))) #model averaging
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#FIGURE 2
#Characterising community composition at phylum level
Phylum_plot <- GM_Climate_Nutrition %>%
  dplyr::filter(!Phylum=="NA") %>% dplyr::filter(!Family=="Mitochondria") %>% dplyr::filter(!Family=="Chloroplast") %>% dplyr::filter(!Family=="uncultured") %>%
  ggplot(., aes(x = reorder(SampleID, Number), y = count, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill", colour="transparent") + ylab("Abundance") + xlab("Sample") +
  scale_fill_manual(values=met.brewer("Renoir", 24)) + theme_bw(base_size = 14) + guides(fill=guide_legend(ncol=2)) + theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), strip.background = element_rect(fill="white"))
Phylum_plot  
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#TABLE 2 prep
#remove all predictor variable rows with NAs for MCMCglmm
GLMM_GM_Dataset <- GM_Climate_Nutrition %>%
  filter(!intake_ndf_30_day=="NA") %>%
  filter(!diet_div_30_day=="NA") %>%
  filter(!prop_fruit_30_day=="NA")
View(GLMM_GM_Dataset)
#Scale predictor variables, if necessary
GLMM_GM_Dataset_Scaled <- GLMM_GM_Dataset %>%
  mutate(temp_scaled=scale(temp_30_day)) %>% mutate(rain_scaled=scale(rainfall_30_day)) %>% mutate(ndf_scaled=scale(intake_ndf_30_day)) %>% mutate(tnc_scaled=scale(intake_tnc_30_day)) %>%
  mutate(fat_scaled=scale(intake_fat_30_day)) %>% mutate(fruit_scaled=scale(prop_fruit_30_day)) %>% mutate(dietdiv_scaled=scale(diet_div_30_day))
View(GLMM_GM_Dataset_Scaled)
#
#TABLE 2
#Models predicting microbial relative abundances using MCMCglmm
###library size not included, as SampleID captures this###
mf <- 10 #multiplier for model iterations and chains 
Priors <- list(R = list(V = diag(1), nu = 2), #Specify priors for each random effect
                G = list(G1 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100),
                         G2 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100),
                         G3 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100),
		        G4 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100)))
GLMM <- MCMCglmm(count ~ temp_scaled + rain_scaled + ndf_scaled + tnc_scaled + fat_scaled + fruit_scaled + dietdiv_scaled, 
                  random= ~ batch + ASV_ID + individualID + SampleID, 
                  prior=Priors,
                  data=GLMM_GM_Dataset_Scaled,
                  verbose=T, pr=T, pl=T, # will save posteriors for random effect levels & model predictions for each iteration 
                  family = "poisson",
                  nitt = 13000*mf,
                  thin = 10*mf,burnin=3000*mf)
summary(GLMM)
summary(GLMM$Sol)
plot(GLMM$VCV)
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#FIGURE 3, FIGURE 4, SUPPLEMENTARY FIGURE 2 prep
#Discretize predictor variables
GLMM_GM_Dataset_Discrete <- GLMM_GM_Dataset %>%
  mutate(temp_discrete = ifelse(temp_30_day <  15, "Low", "High")) %>% mutate(rain_discrete = ifelse(rainfall_30_day <  400, "Low", "High")) %>% mutate(ndf_discrete = ifelse(intake_ndf_30_day <  50, "Low", "High")) %>% mutate(tnc_discrete = ifelse(intake_tnc_30_day <  50, "Low", "High")) %>%
  mutate(fat_discrete = ifelse(intake_fat_30_day <  10, "Low", "High")) %>% mutate(fruit_discrete = ifelse(prop_fruit_30_day <  50, "Low", "High")) %>% mutate(dietdiv_discrete = ifelse(diet_div_30_day <  1, "Low", "High"))
View(GLMM_GM_Dataset_Discrete)
#
#FIGURE 3, FIGURE 4, SUPPLEMENTARY FIGURE 2
#Discretized variable models predicting microbial relative abundances using MCMCglmm
#include interaction between ASV_ID and each fixed effect variable (to later extract posterior distributions)
mf <- 10 #multiplier for model iterations and chains 
Priors2 <- list(R = list(V = diag(1), nu = 2), #Specify priors - 1x "G" for each random effect (including interaction terms)
                G = list(G1 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100),
                         G2 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100),
                         G3 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100),
                         G4 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100),
                         G5 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100),
                         G6 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100),
                         G7 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100),
                         G8 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100),
                         G9 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100),
                         G10 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100),
		        G11 = list(V = diag(1), nu = 2, alpha.mu = 0, alpha.V = diag(1)*100)))
GLMM2 <- MCMCglmm(count ~ temp_discrete + rain_discrete + ndf_discrete + tnc_discrete + fat_discrete + fruit_discrete + dietdiv_discrete, 
                  random= ~ batch + ASV_ID + individualID + SampleID + ASV_ID:temp_discrete + ASV_ID:rain_discrete + ASV_ID:ndf_discrete + ASV_ID:tnc_discrete + ASV_ID:fat_discrete + ASV_ID:fruit_discrete + ASV_ID:dietdiv_discrete, 
                  prior=Priors2,
                  data=GLMM_GM_Dataset_Discrete,
                  verbose=T, pr=T, pl=T, # will save posteriors for random effect levels & model predictions for each iteration 
                  family = "poisson",
                  nitt = 13000*mf,
                  thin = 10*mf,burnin=3000*mf)
summary(GLMM2)
summary(GLMM2$Sol)
plot(GLMM2$VCV)
#Pull out posterior distributions for groups of interest (done at phylum, family, and genus levels)
Solutions<-GLMM2$Sol[, 1:dim(GLMM2$Z)[2]]
SolutionsDF<- Solutions %>% as.data.frame()
#Run function below to apply to posterior distributions for each ASV in each group of predictor, calcs diffs between groups
RanefPostComp <- function(GLMM_GM_Dataset_Discrete, comps, taxa, var){
  require(tidyverse)
  require(MCMCglmm)
  not_all_na <- function(x) any(!is.na(x))
  group <- GLMM_GM_Dataset_Discrete[, var]
  taxa <- GLMM_GM_Dataset_Discrete[, taxa]
  GLMM_GM_Dataset_Discrete<- GLMM_GM_Dataset_Discrete %>% filter(term!="(Intercept)") %>% 
    mutate(taxa=as.factor(taxa))
  LevelComps <-list()
  Levels=levels(GLMM_GM_Dataset_Discrete[, "taxa"])
  for(x in 1:length(Levels)){
    post1 <- GLMM_GM_Dataset_Discrete %>% subset(taxa==Levels[x] & group==comps[1]) %>% pull(posterior)
    post2 <- GLMM_GM_Dataset_Discrete %>% subset(taxa==Levels[x] & group==comps[2]) %>% pull(posterior)
    if(length(post2>0 ) & length(post1>0)){
      postComp <- post2-post1 %>% as.matrix() 
    }
    else{
      postComp <- rep(NA, 1000)
    }
    LevelComps[[x]] <-postComp
  }
  LevelCompsFin <- bind_cols(LevelComps) %>% as.matrix()
  colnames(LevelCompsFin) <- paste(Levels, comps[[2]], sep=".")
  LevelCompsFin %>% as.data.frame() %>% select_if(not_all_na) -> GLMM_GM_Dataset_Discrete 
  
  df1 <- apply(GLMM_GM_Dataset_Discrete,2, function(a) mean(a)[1])
  df1l <- apply(GLMM_GM_Dataset_Discrete,2,function(a)  HPDinterval(as.mcmc(a))[1])
  df1u <- apply(GLMM_GM_Dataset_Discrete,2,function(a)  HPDinterval(as.mcmc(a))[2])
  df1v <- apply(GLMM_GM_Dataset_Discrete,2,function(a)  var(a)[1])
  
  df2 <- cbind(df1, df1l, df1u, df1v)
  colnames(df2) <- c("mean", "lower", "upper", "var")
  df2 <- df2 %>% as.data.frame() %>% 
    rownames_to_column(var="term") %>% 
    mutate(taxa=sapply(str_split(term, "[[.]]"), '[', 1),
           effect=var, 
           nonZeroCI= ifelse(lower<0&upper<0, "Y", 
                             ifelse(lower>0&upper>0, "Y", "N"))) 
  colnames(df2)[which(names(df2) == "taxa")] <- paste(taxa)
  return(df2)
  
}
###Apply function to each predictor
###e.g., NDF (for all following code)
SolKeepNDF <- SolutionsDF %>% dplyr::select(contains("ASV_ID:ndf_discrete")) %>% colnames()
SolutionsNDF<- SolutionsDF %>% dplyr::select(SolKeepNDF) 
SolutionsHighNDF <- SolutionsNDF %>% 
  gather(key="level", value="posterior") %>% as.data.frame() %>% 
  separate(level, into=c("term", "ASV_ID", "ndf_discrete"), sep="[.]")
ASVByNDF <- RanefPostComp(SolutionsHighNDF, comps=c("Low", "High"), taxa="ASV_ID", var="ndf_discrete")
ASVByNDF <- ASVByNDF %>%
  dplyr::rename(ASV_ID = "0001e08442abc4b10626e0fdd0f9219d") #add back in dropped column name for ASV_ID
View(ASVByNDF)
#Prep data frame for visualization
Taxo <- GLMM_GM_Dataset_Discrete %>% 
  filter(!duplicated(ASV_ID)) %>%
  dplyr::select(ASV_ID, Phylum, Family, Genus) %>%
  as.data.frame()
Prev <- GLMM_GM_Dataset_Discrete %>% group_by(ASV_ID) %>% 
  summarise(prevASV_ID=length(count[count>0])/ 153) %>%  # n samples 
  arrange(-prevASV_ID)
TaxoDf <- left_join(Prev, Taxo) 
TaxNDF <- merge(ASVByNDF, TaxoDf, by="ASV_ID", all.x=TRUE)
View(TaxNDF)
#Identify taxa with strongest effects
keepPhy <- TaxNDF %>% group_by(Phylum) %>% na.omit() %>% summarise(nASV_ID=length(ASV_ID)) %>% filter(nASV_ID>1)
keepPhy <- keepPhy %>% pull(Phylum) %>% as.character() # retain Phylum with > 1 ASV_ID for visualisation 
TopNDFASV_ID <- TaxNDF %>% filter(nonZeroCI=="Y") # pull out taxa effects with CI not spanning zero 
Top50 <- TaxNDF %>% arrange(-mean) %>% filter(Phylum %in% keepPhy) %>% dplyr::select(ASV_ID, Phylum, lower, upper, mean) %>% filter(lower>0) %>%  #make sure errors don't span zero 
  top_n(50) %>% mutate(group="real")
Bottom50 <- TaxNDF %>% arrange(-mean) %>% filter(Phylum %in% keepPhy) %>% dplyr::select(ASV_ID, Phylum, lower, upper, mean) %>%  filter(upper<0.1) %>% 
  top_n(-50) %>% mutate(group="real")
#Create filler categories for axis breaks for forest plot
DummiesNDF<-paste0("dummy", rep(1:10))
GroupNDF<-rep("dummy", 10)
DummyNDF<- cbind(DummiesNDF, GroupNDF) %>% as.data.frame()
colnames(DummyNDF)<- c("ASV_ID", "group")
select50<- bind_rows(Top50, DummyNDF, Bottom50) 
select_x_NDF<- select50 %>% pull(ASV_ID)
#Make Forest Plot
NDFTopBottom <- select50 %>% filter(!is.na(Phylum)) %>% 
  ggplot(aes(x=reorder(ASV_ID, mean), y=mean, colour=Phylum)) +
  geom_point(position=position_dodge(w=0.4), size=1.5) + geom_errorbar(position=position_dodge(w=0.5), aes(ymin=lower,ymax=upper), size=0.9,width=0.4) + geom_segment(aes(y=0, yend=0, x=0, xend=37), colour="black", size=0.5) +
  theme_classic(base_size = 14) + labs(x=NULL) +
  scale_colour_manual(values = c("Actinobacteriota"= "#17154f", "Bacteroidota"= "#2f357c", "Cyanobacteria"= "#b0799a", "Firmicutes"= "#f6b3b0", "Proteobacteria"= "#e48171", "Spirochaetota"= "#e69b00", "Synergistota"= "#f5bb50", "Verrucomicrobiota"= "#ada43b")) + #*always need to select colors that match other plot*
  scale_x_discrete(limits=rev(select_x_NDF)) + 
  coord_flip() +
  theme(strip.text.x = element_text(size = 0.8), strip.background =element_rect(fill="white"), axis.text.y = element_blank(), axis.ticks.y=element_blank())             
NDFTopBottom
#
#Load percentages data frame to plot strongest positiver and strongest negative shifts for each predictor variable
#E.g., NDF
NDFplotpositive <- Diff_Phylum %>%
  filter(Predictor=="NDF") %>% filter(Rank=="Positive Shifts") %>%
  ggplot(., aes(y = Percent, x = reorder(Phylum, -Percent))) + geom_bar(aes(fill = Phylum), stat = "identity", position = position_dodge()) +
  scale_fill_manual(values=c("Actinobacteriota"= "#17154f", "Bacteroidota"= "#2f357c", "Cyanobacteria"= "#b0799a", "Firmicutes"= "#f6b3b0", "Proteobacteria"= "#e48171", "Spirochaetota"= "#e69b00", "Synergistota"= "#f5bb50", "Verrucomicrobiota"= "#ada43b")) + labs(y = "% Top Taxa", x = NULL) + ylim(0,100) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none", axis.title.x = element_blank())
NDFplotpositive
NDFplotnegative <- Diff_Phylum %>%
  filter(Predictor=="NDF") %>% filter(Rank=="Negative Shifts") %>%
  ggplot(., aes(y = Percent, x = reorder(Phylum, -Percent))) + geom_bar(aes(fill = Phylum), stat = "identity", position = position_dodge()) +
  scale_fill_manual(values=c("Actinobacteriota"= "#17154f", "Bacteroidota"= "#2f357c", "Cyanobacteria"= "#b0799a", "Firmicutes"= "#f6b3b0", "Proteobacteria"= "#e48171", "Spirochaetota"= "#e69b00", "Synergistota"= "#f5bb50", "Verrucomicrobiota"= "#ada43b")) + labs(y = "% Top Taxa", x = NULL) + ylim(0,100) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none", axis.title.x = element_blank())
NDFplotnegative
#*Repeat for other taxonomic levels*
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#SUPPLEMENTARY FIGURE 1
#Plots showing data spread for each variable, with threshold for discretization (repeat for each plot)
#e.g., NDF: 50 kcal
ndf_plot <- GM_Climate_Nutrition_reordered %>%
  ggplot(., aes(x = Date, y = intake_ndf_30_day)) + 
  geom_point() + labs(y="30 Day Mean NDF Intake") + theme_classic() + theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank())
ndf_plot
#Combine plots
correlations <- ggarrange(temp_plot, rain_plot, fruit_plot, ndf_plot, tnc_plot, fat_plot, nrow=3, ncol=2)
correlations
#
#
#END OF SCRIPT
