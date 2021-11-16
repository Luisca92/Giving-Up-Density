################################## SET WORKING DIRECTORY ################################# 
setwd("~/Desktop/R") # Set the working directory
Treatments <- read.csv("Data/Treatments.csv")

###################################### LOAD PACKAGES ##################################### 
library(Rmisc) # Necessary for summarySE
library(dplyr) # Necessary for group_by
library(reshape) # Necessary for melt & cast
library(ggplot2) # Necessary for ggplots
library(MASS) # Necessary for glm.nb
library(ggpubr) # Necessary for ggarrange
library(conover.test) # Necessary for conover.test
library(Hmisc) # 
library(lme4) # Necessary for lmer
library(car) # Necessary for Levene's Test
asinTransform <- function(p) { asin(sqrt(p)) } # Function for arcsin square root transformation.
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-(sd(x)/sqrt(60))
  ymax <- m+(sd(x)/sqrt(60))
  return(c(y=m,ymin=ymin,ymax=ymax))
}
##################################### DATAFRAME SETUP #################################### 
# GUD Dataframe
GUD <- read.csv("Data/GUD.csv")
GUD <- melt(GUD, id=c("Plot","Tray","Treatment"))
colnames(GUD) <- c("Plot", "Tray", "Treatment", "Date", "GUD")

Temp <- aggregate(FernPlots[, 4], list(FernPlots$TSubplot), mean)
Temp2 <- aggregate(FernPlots[, 34], list(FernPlots$TSubplot), mean)
colnames(Temp) <- c("Subplot", "Fronds") # Rename columns
Temp$Canopy <- Temp2[,2]
Temp$Plot <- round(Temp$Subplot)
GUD$Fronds <- with(Temp, Fronds[match(GUD$Tray, Subplot)]) # Match Fronds by Subplot 
GUD$Treatment <- as.factor(GUD$Treatment)
GUD$Treatment <-revalue(GUD$Treatment, c("Control"="Nat. Succ."))
GUD$Date <- as.factor(GUD$Date)
GUD$Tray <- as.factor(GUD$Tray)
GUD$Plot <- as.factor(GUD$Plot)

# Fern Dataframe
Fern <- read.csv("Data/Fern.csv") 
GrassLianaPlots <- subset(Fern, Fern$Grass.... > 49 | Fern$Liana.... > 49) # Create dataframe with grass or liana dominated quadrats.
Fern$TSubplot <- as.numeric(Fern$Plot)+(as.numeric(Fern$Subplot)/10) # Create unique codename for each subplot.
Fern$TQuadrat <- Fern$Plot+(Fern$Quadrat*0.01) # Create unique codename for each quadrat
Fern$Canopy <- with(Fern, (Fern$Canopy.Cover+Fern$Canopy.Cover..2)/2, na.rm = TRUE) # Create variable based on canopy cover measures.
FernPlots <- Fern[!(Fern$TQuadrat %in% GrassLianaPlots$TQuadrat),] # Exclude the grass/liana dominated quadrats.

# Camera Dataframe
Cam <- read.csv("Data/CamAct.csv")
Cam <- as.data.frame(t(Cam))
Cam$Time <- row.names(Cam)
colnames(Cam) <- c("All", "Ate", "Forest", "Animal", "Wind", "Nat. Succ.", "Time") # Rename columns
Cam <- Cam[-1,]
row.names(Cam) <- 1:49
Cam <- Cam[-49,]

Cam$Time <- levels(Cam$Time) <- c("0:00","0:30","1:00","1:30","2:00","2:30","3:00","3:30","4:00","4:30","5:00","5:30","6:00","6:30","7:00","7:30","8:00","8:30","9:00","9:30","10:00","10:30","11:00","11:30","12:00","12:30","13:00","13:30","14:00","14:30","15:00","15:30","16:00","16:30","17:00","17:30","18:00","18:30","19:00","19:30","20:00","20:30","21:00","21:30","22:00","22:30","23:00","23:30")

Cam2 <- Cam[,-1] # With Total Ate
Cam3 <- Cam2[,-1] # Only treatments
Cam4 <- melt(Cam3, id=c("Time")) # For GG
colnames(Cam4) <- c("Time", "Habitat", "Frequency") # Rename columns
Cam4$Frequency <- as.numeric(Cam4$Frequency)
Cam4$Time <- as.chron(Cam4$Time)

Cam4$hms <- format(Cam4$Time, format = "%H:%M")
Cam4$hms2 <- as.POSIXlt(Cam4$hms, format = "%H:%M")
Cam4$hms2 <- format(Cam4$hms2, format = "%H:%M")

###################################### Setup Fern x GUD ###################################### 
Temp3 <- aggregate(Temp[, 2:3], list(Temp$Plot), mean, na.rm=TRUE) # This is summary of fern and canopy per plot. 
colnames(Temp3) <- c("Plot", "Subplot", "Fronds") # Rename columns
# Lets do another way where we only retain the subplots used in the GUD study. 
# We'll need to get summaries per subplot. 
Temp4 <- Temp

Temp4$GUD <- with(GUD, GUD[match(Temp4$Subplot, Tray)]) # Match GUDs by Subplot 
Temp4$Canopy <- with(Temp, Canopy[match(Temp4$Subplot, Subplot)])

GUDX <- Temp4[complete.cases(Temp4), ] # Excellent. So this has it per subplot, but frond density is low in 1.1. 
# Dont forget zero-inflated data. 

GUDX$Treatment <- with(Treatments, Treatment[match(GUDX$Plot, Plot)]) # Match Treatment by Plot
Temp3$Treatment <- with(Treatments, Treatment[match(Temp3$Plot, Plot)]) # Match Treatment by Plot

G2 <- GUDX # GUDX has the subplot values.
G3 <- Temp3 # has plot averages for the subplot

##################################### GUD X TREATMENT #################################### 
# New Chris Strat
# You can rank order the data (GUDs) 
# from least to greatest within each treatment, and then run an ANOVA on the ranked data. 

# we arcsin square-root-transformed proportional GUDs (final food density/initial food density). 
# However, following this transformation, Levene’s test of heterogeneity of error variance was highly significant,
# indicating that error variances were unequal. To increase confidence in the analysis, 
# we then repeated the analysis on rank-transformed GUDs. GUD$Treatment <- as.factor(GUD$Treatment)

# Releveling
GUD$Treatment <- relevel(GUD$Treatment, ref = "Wind") # Relevel 
GUD$Treatment <- relevel(GUD$Treatment, ref = "Animal") # Relevel 
GUD$Treatment <- relevel(GUD$Treatment, ref = "Forest") # Relevel 

GUD$Prop <- GUD$GUD/15 # (final food density/initial food density). 
GUD$Prop2 <- pmin(GUD$Prop, 1) # This makes it so that anything above the proportion is turned to 1. 

leveneTest(GUD$GUD, GUD$Treatment) # F = 9.3423, df = 3, p = 7.354e-06
leveneTest(GUD$Prop, GUD$Treatment) # F = 9.3423, df = 3, p = 7.354e-06

GUD$Prop2Trans <- asinTransform(GUD$Prop2) # ArcSinSquareRoot Transformation
GUD$Prop2Trans2 <- asinTransform(GUD$Prop) # ArcSinSquareRoot Transformation With NAs

leveneTest(GUD$Prop2Trans, GUD$Treatment) # F = 11.119, df = 3, p = 7.473e-07
leveneTest(GUD$Prop2Trans2, GUD$Treatment) # F = 8.7472, df = 3, p = 1.998e-05

# Time to Rank Them
GUD$Prop2TransRank <- rank(GUD$Prop2Trans) # Ranked! 
GUD$Prop2TransRank2 <- rank(GUD$Prop2Trans2) # Ranked! with NAs

leveneTest(GUD$Prop2TransRank, GUD$Treatment) # F = 3.9429, df = 3, p =  0.009027
leveneTest(GUD$Prop2TransRank2, GUD$Treatment) # F = 4.0107, df = 3, p =  0.008252 

# Still significant Levene's
GUD$Treatment <- relevel(GUD$Treatment, ref = "Forest") # Relevel 
lmRank <- lm(Prop2TransRank ~ Treatment, data = GUD)
summary(lmRank)

# F-statistic: 15.53 on 3 and 236 DF,  p-value: 2.964e-09
# Animal x Control: -1.548   0.1229 
# Animal x Forest: 4.815 2.63e-06 ***
# Animal x Wind: 2.394   0.0174 * 
# Wind x Control: -3.943 0.000106 ***
# Wind x Forest: 2.420 0.016257 *
# Forest x Control: -6.363 1.02e-09 ***

GUD$Treatment <- relevel(GUD$Treatment, ref = "Forest") # Relevel 
lmRank2 <- lm(Prop2TransRank2 ~ Treatment, data = GUD) # This makes the NA ranks go higher.
summary(lmRank2)
# F-statistic: 15.02 on 3 and 236 DF,  p-value: 5.607e-09
# Animal x Control: -1.594   0.1122 
# Animal x Forest: 4.650 5.52e-06 ***
# Animal x Wind: 2.419   0.0163 *  
# Wind x Control: -4.013 8.05e-05 ***
# Wind x Forest: 2.231   0.0266 * 
# Forest x Control: -6.245 1.96e-09 ***

# Time for a mixed-effects option
lmer1 <- lmer(GUD ~ Treatment + (1|Tray), data = GUD)  # Normal distribution. 

summary(lmer1) # Tray     (Intercept) 4.671
anova(lmer1) # F-statistic: 1.2605 on 3 and 20 DF,  p-value: 0.3147

GUD$Treatment <- relevel(GUD$Treatment, ref = "Forest") # Relevel 

lmer2 <- lmer(Prop2TransRank ~ Treatment + (1|Tray), data = GUD)  
anova(lmer2) # F-statistic: 3.4872 on 3 and 20 DF,  p-value: 0.0349*
summary(lmer2) # Tray     (Intercept) 1504     38.78 residual = 2709

1504+2709 # 4213
1504/4213*100 # 35.70 % of the variance is explained by tray/plot. 
# Animal x Control: -0.734   0.4717
# Animal x Forest: -2.281  0.03363 * 
# Animal x Wind:   1.134   0.2700  
# Wind x Control: -1.868   0.0765 . 
# Wind x Forest: -1.147  0.26500 
# Forest x Control: -3.015  0.00684 ** 

lmer3 <- lmer(Prop2TransRank ~ Treatment + (1|Tray) + (1|Date), data = GUD)  
anova(lmer3) # F statistic  3.4872 on 3 and 20 DF, p-value = 0.0349 *
summary(lmer3)
1565.8 + 620.8 + 2088.5 # 4275.1
1565.8 / 4275.1 * 100 #  36.626 % Plot
620.8 / 4275.1 * 100 # 14.521 %

# Animal x Control: -0.734   0.4717
# Animal x Forest: -2.281  0.03363 * 
# Animal x Wind:   1.134   0.2700  
# Wind x Control: -1.868   0.0765 . 
# Wind x Forest: -1.147  0.26500 
# Forest x Control: -3.015  0.00684 ** 

# Seeds often registered values above 15 g when not consumed. Seeds imbibed moisture from the air. 
# Here we cap all average GUDs with 15 g cap. 
GUD$GUD2 <- pmin(GUD$GUD, 15)
# Releveling
GUD$Treatment <- relevel(GUD$Treatment, ref = "Wind") # Relevel 
GUD$Treatment <- relevel(GUD$Treatment, ref = "Animal") # Relevel 
GUD$Treatment <- relevel(GUD$Treatment, ref = "Forest") # Relevel 

gd <- GUD %>% 
  group_by(Treatment) %>% 
  summarise(GUD2 = mean(GUD2))
gd

Fig2A <- ggplot(GUD, aes(x = Treatment, y = GUD2, color = Treatment, fill = Treatment)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  ggrepel::geom_text_repel(aes(label = Plot), color = "black", size = 2, segment.color = "grey") +
  geom_point() +
  guides(color = "none", fill = "none") +
  scale_y_continuous(breaks=seq(0,15,5), limits = c(0,18)) +
  scale_fill_manual(values=c("peachpuff4","#9999CC", "#66CC99", "#CC6666")) +
  scale_color_manual(values=c("peachpuff4","#9999CC", "#66CC99", "#CC6666")) +
  theme_bw(base_size = 15) +
  labs(title = "", x = "Habitat", y = "Giving-Up Density (g)")

Fig2A


###################################### GUD X Fern ###################################### 
# Match the other GUDs
GUDX$Prop2TransRank <- with(GUD, Prop2TransRank[match(GUDX$Subplot, Tray)]) # Match Treatment by plot. 

GUDX$GUD2 <- pmin(GUDX$GUD, 15)
shapiro.test(GUDX$Fronds) # W = 0.88676144, p-value = 0.03397521 not Normal
shapiro.test(GUDX$GUD2) # W = 0.5252275, p-value = 1.337515e-06

hist(sqrt(GUDX$Fronds))
shapiro.test(sqrt(GUDX$Fronds)) # W = 0.97702474, p-value = 0.9142472
# Square Rooting the Fronds helps with normality.


GUDX2 <- droplevels(subset(GUDX, GUDX$Plot != "1")) # Exclude Plot 1
# Exclude Plot 1 as that subplot is mostly grass?

Fig2B <- ggplot(GUDX, aes(x=Fronds, y=GUD2)) + 
  geom_point() +
  geom_smooth(method=lm, se = FALSE, linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = Plot), color = "black", size = 2, segment.color = "grey") +
  xlab(expression(Frond~Density~"m"^{2})) + 
  scale_y_continuous(breaks=seq(0,15,5), limits = c(0,16)) +
  ylab("Giving-Up Density (g)") + 
  theme_bw(base_size = 15) +theme(legend.position="bottom")
Fig2B

Fig3 <- ggplot(GUDX2, aes(x=Fronds, y=GUD2)) + 
  geom_point() +
  geom_smooth(method=lm, se = FALSE, linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = Plot), color = "black", size = 2, segment.color = "grey") +
  xlab(expression(Frond~Density~"m"^{2})) + 
  scale_y_continuous(breaks=seq(0,15,5), limits = c(0,16)) +
  ylab("Giving-Up Density (g)") + 
  theme_bw(base_size = 15) +theme(legend.position="bottom")
Fig3

inter<-interaction(GUDX2$Treatment, GUDX2$Fronds)
kruskal.test(GUD2 ~ inter, data = GUDX2)

# Kruskal-Wallis chi-squared = 16, df = 16, p-value = 0.4529608
lm1 <- lm(GUD2 ~ Fronds * Treatment, data = GUDX2)
anova(lm1) # Fronds:Treatment  F =  2.99940, p = 0.091277 .
summary(lm1) # F-statistic: 1.477594 on 5 and 11 DF,  p-value: 0.2733157
# Multiple R-squared:  0.4017827,	Adjusted R-squared:  0.1298658 

lm2 <- lm(Prop2TransRank ~ Fronds * Treatment, data = GUDX2)
anova(lm2) # F =  2.70672 0.069738 . Marginally sig ANCOVA
summary(lm2) # F-statistic: 3.186984 on 5 and 164 DF,  p-value: 0.008985387. SIG model.

# THIS
WindMod <- lm(Prop2TransRank ~ Fronds , data=droplevels(subset(GUDX2, GUDX2$Treatment == "Wind")))  # build linear regression model on full data
AnimMod <- lm(Prop2TransRank ~ Fronds , data=droplevels(subset(GUDX2, GUDX2$Treatment == "Animal")))  # build linear regression model on full data
ContMod <- lm(Prop2TransRank ~ Fronds , data=droplevels(subset(GUDX2, GUDX2$Treatment == "Nat. Succ.")))  # build linear regression model on full data

summary(WindMod) # F-statistic: 0.5759923 on 1 and 58 DF,  p-value: 0.4509609
summary(AnimMod) # F-statistic: 0.8321731 on 1 and 58 DF,  p-value: 0.3654212
summary(ContMod) # F-statistic: 7.887673 on 1 and 48 DF,  p-value: 0.007176919

lm3 <- lm(Prop2TransRank ~ Fronds * Treatment, data = GUDX2)

# Square Rooting the Fronds 

lm3 <- lm(Prop2TransRank ~ sqrt(Fronds) * Treatment, data = GUDX)
summary(lm3) # F-statistic: 2.903552 on 5 and 164 DF,  p-value: 0.01537837,  R-squared:  0.08132389,	Adjusted R-squared:  0.05331547 
anova(lm3) # F 1.65682, p = 0.193920  Here's the money. 

lm4 <- lm(Prop2TransRank ~ sqrt(Fronds), data = GUDX2)
summary(lm4) #F = 0.00566, p = 0.94101
anova(lm4)

# Reported test. Treatments pooled. 
lm5 <- lm(Prop2TransRank ~ sqrt(Fronds), data = GUDX)
lm6 <- lm(Prop2TransRank ~ Fronds, data = GUDX)
summary(lm5) # F-statistic: 0.0717 on 1 and 16 DF,  p-value: 0.792, Adjusted R-squared:  -0.0578 
summary(lm6) # F-statistic: 0.0717 on 1 and 16 DF,  p-value: 0.792, Adjusted R-squared:  -0.0538

cor(GUDX$Prop2TransRank ~ GUDX$Fronds)

Fig2A
# FIG2B
Fig2C

ggarrange(Fig2A, Fig2B, # Main Fig.
          align='h',labels=c('(A)','(B)'),
          common.legend = F, nrow = 1, ncol = 2)

###################################### GUD X Canopy ###################################### 

GUD3 <- aggregate(GUD2$GUD, by=list(Plot=GUD2$Plot, Subplot=GUD2$Tray, Canopy=GUD2$Canopy, Fronds=GUD2$Fronds, Treatment=GUD2$Treatment), FUN= mean, na.rm = "TRUE")
GUD3$Subplot <- as.factor(GUD3$Subplot)
colnames(GUD3) <- c("Plot","Subplot", "Canopy","Fronds", "Treatment", "GUD") # Rename columns

GUD4 <- droplevels(subset(GUD3, GUD3$Treatment != "Forest")) # Exclude Forest
GUD4$Treatment <- relevel(GUD4$Treatment, ref = "Wind") 
GUD4$Treatment <- relevel(GUD4$Treatment, ref = "Animal") 

ggplot(GUD4, aes(x=Fronds, y=GUD, color=Treatment, fill=Treatment)) + 
  geom_point() +
  geom_smooth(method=lm, se = TRUE) +
  xlab(expression(Average~Frond~Density~"1m"^{2})) + 
  ylab("Giving-Up Density") + 
  scale_color_manual(values=c("#CC6666", "#9999CC", "#66CC99")) +
  scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99")) +
  theme_minimal(base_size = 25) +theme(legend.position="bottom")

ggplot(GUD4, aes(x=Fronds, y=GUD)) + 
  geom_point() +
  geom_smooth(method=lm, se = TRUE) +
  xlab(expression(Average~Frond~Density~"1m"^{2})) + 
  ylab("Giving-Up Density") + 
  scale_color_manual(values=c("#CC6666", "#9999CC", "#66CC99")) +
  scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99")) +
  theme_minimal(base_size = 25) +theme(legend.position="bottom")

# We need plot summaries... otherwise this is pseudo replicated. 
summary(GUD4$Plot)

GUD5 <- cast(GUD4, Plot~Treatment, mean)

ggplot(GUD4, aes(x=Canopy, y=GUD, color=Treatment, fill=Treatment)) + 
  geom_point() +
  geom_smooth(method=lm, se = TRUE) +
  xlab("Canopy Cover (%)") + 
  ylab("Giving-Up Density") + 
  scale_color_manual(values=c("#CC6666", "#9999CC", "#66CC99")) +
  scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99")) +
  theme_minimal(base_size = 25) +theme(legend.position="bottom")

ggplot(GUD4, aes(x=Canopy, y=GUD)) + 
  geom_point() +
  geom_smooth(method=lm, se = TRUE) +
  xlab("Canopy Cover (%)") + 
  ylab("Giving-Up Density") + 
  scale_color_manual(values=c("#CC6666", "#9999CC", "#66CC99")) +
  scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99")) +
  theme_minimal(base_size = 25) +theme(legend.position="bottom")

M1 <- lm(GUD ~ Fronds * Treatment, data = GUD4)
anova(M1); summary(M1) # No interaction

M2 <- lm(GUD ~ Canopy * Treatment, data = GUD4)
anova(M2); summary(M2) # No interaction

AnimMod <- lm(GUD ~ Fronds, data=droplevels(subset(GUD4, GUD4$Treatment == "Animal")))  # build linear regression model on full data
anova(AnimMod);summary(AnimMod) # F-statistic: 1.547 on 1 and 4 DF,  p-value: 0.2815, Multiple R-squared:  0.2789,	Adjusted R-squared:  0.09864 
AnimMod2 <- lm(GUD ~ Canopy, data=droplevels(subset(GUD4, GUD4$Treatment == "Animal")))  # build linear regression model on full data
anova(AnimMod2);summary(AnimMod2) # F-statistic: 0.5455 on 1 and 4 DF,  p-value: 0.5012, Multiple R-squared:   0.12,	Adjusted R-squared:   -0.1 

##########################################################################################
####################################### PHOTOGRAPHS ###################################### 
##########################################################################################

ggplot(Cam4, aes(x=hms2, y=Frequency)) + 
  geom_area(stat = "bin") +
  xlab("Habitat") +
  ylab("Frequency") +   
  theme_bw(base_size = 25) + theme(legend.position="none") 

ggplot(Cam4, aes(x=hms, y=Frequency, fill=Habitat)) + 
  geom_area() 

ggplot(Cam4, aes(x= hms2, y= Frequency, group = Habitat, fill = Habitat)) +
  xlab("Time of Day") + ylab("Frequency") +
  geom_area(position = "stack") + 
  theme_bw(base_size = 20) + theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("peachpuff4","#66CC99", "#9999CC", "#CC6666")) 

ggplot(Cam4, aes(hms2, Frequency, group = Habitat)) +
  geom_area(aes(fill = Habitat), alpha = .4) +
  xlab("Time of Day") + ylab("Photographed Rodents") +
  theme_bw(base_size = 15) + theme(legend.position="bottom") +
  scale_y_continuous(limits = c(0,6), expand = c(0, 0)) + 
  geom_line(aes(group = Habitat), position = "stack") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("peachpuff4","#66CC99", "#9999CC", "#CC6666")) 

ggplot(Cam4, aes(hms2, Frequency, group = Habitat)) +
  geom_area(aes(fill = Habitat), position = "stack") +
  xlab("Time of Day") + ylab("Frequency") +
  theme_bw(base_size = 15) + theme(legend.position="bottom", legend.direction = 'horizontal') +
  scale_y_continuous(limits = c(0,6), expand = c(0, 0)) + 
  geom_line(aes(group = Habitat), position = "stack") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("peachpuff4","#66CC99", "#9999CC", "#CC6666")) + 
  guides(group = guide_legend(title.hjust = 0.5)) 

guides(color = guide_legend(title.position = "top", 
                            # hjust = 0.5 centres the title horizontally
                            title.hjust = 0.5,
                            label.position = "bottom")) 

ggplot(Cam4, aes(hms2, Frequency, group = Habitat)) +
  geom_area(aes(fill = Habitat), position = "stack") +
  xlab("Time of Day") + ylab("Frequency") +
  theme_bw(base_size = 15) + theme(legend.position="bottom") +
  scale_y_continuous(limits = c(0,6), expand = c(0, 0)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("#66CC99","#66CC99", "#66CC99", "#66CC99")) 


# GUD x Canopy

Fern <- read.csv("Data/Fern.csv") # Check Excel File for deleted row in Plot 19.
Fern$TSubplot <- as.numeric(Fern$Plot)+(as.numeric(Fern$Subplot)/10) # Create unique codename for each subplot.
Fern$TQuadrat <- Fern$Plot+(Fern$Quadrat*0.01) # Create unique codename for each quadrat
Fern$Treatment <- with(Slope, Treatment[match(Fern$Plot, Plot)]) # Match Treatment by plot. 
Fern$Mean.Frond <- as.numeric(as.character(Fern$Mean.Frond)) # Make numeric.
Fern$Canopy <- with(Fern, (Fern$Canopy.Cover+Fern$Canopy.Cover..2)/2) # Create variable based on canopy cover measures.
Fern$Treatment <- relevel(Fern$Treatment, ref = "Control") # Relevel 
Fern$Plot <- as.factor(Fern$Plot) # Need to turn into factor variables. 
Fern$Subplot <- as.factor(Fern$Subplot) # Need to turn into factor variables. 
Fern$Quadrat <- as.factor(Fern$Quadrat) # Need to turn into factor variables. 


ggplot(Cam4, aes(hms2, Frequency)) +
  geom_area(aes(fill = Habitat))

ggplot(Cam4, aes(x = hms2, fill = Habitat)) +
  geom_bar(position = "fill")

ggplot(Cam4, aes(factor(hms2), fill = factor(Habitat))) + geom_bar()

####################################### Eliminating Capped ###################################### 

GUD5 <- GUD2[which(GUD2$GUD < 15),]

GUD5$Treatment <- relevel(GUD5$Treatment, ref = "Wind") 
GUD5$Treatment <- relevel(GUD5$Treatment, ref = "Animal") 
GUD5$Treatment <- relevel(GUD5$Treatment, ref = "Forest") 
GUD5$Date <- as.factor(GUD5$Date)

gd <- GUD5 %>% 
  group_by(Treatment) %>% 
  summarise(GUD = mean(GUD))
gd

ggplot(GUD5, aes(x = Treatment, y = GUD, color = Treatment, fill = Treatment)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  ggrepel::geom_text_repel(aes(label = Tray), color = "black", size = 2.5, segment.color = "grey") +
  geom_point() +
  scale_fill_manual(values=c("peachpuff4","#66CC99", "#9999CC", "#CC6666", "khaki3")) +
  scale_color_manual(values=c("peachpuff4","#66CC99", "#9999CC", "#CC6666", "khaki3")) +
  guides(color = "none", fill = "none") +
  theme_bw(base_size = 15) +
  labs(title = "", x = "Habitat", y = "Giving-Up Density (GUD)")

ggplot(GUD5, aes(x = Treatment, y = GUD, color = Treatment, fill = Treatment)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  ggrepel::geom_text_repel(aes(label = Plot), color = "black", size = 2.5, segment.color = "grey") +
  geom_point() +
  scale_fill_manual(values=c("peachpuff4","#66CC99", "#9999CC", "#CC6666", "khaki3")) +
  scale_color_manual(values=c("peachpuff4","#66CC99", "#9999CC", "#CC6666", "khaki3")) +
  guides(color = "none", fill = "none") +
  theme_bw(base_size = 15) +
  labs(title = "", x = "Habitat", y = "Giving-Up Density (GUD)")

ggplot(GUD5, aes(x = Treatment, y = GUD, color = Treatment, fill = Treatment)) +# By date
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  ggrepel::geom_text_repel(aes(label = Date), color = "black", size = 2.5, segment.color = "grey") +
  geom_point() +
  scale_fill_manual(values=c("peachpuff4","#66CC99", "#9999CC", "#CC6666", "khaki3")) +
  scale_color_manual(values=c("peachpuff4","#66CC99", "#9999CC", "#CC6666", "khaki3")) +
  guides(color = "none", fill = "none") +
  theme_bw(base_size = 15) +
  labs(title = "", x = "Habitat", y = "Giving-Up Density (GUD)")

GUD5$Date <- relevel(GUD5$Date, ref = "X9.Jun") # Relevel 
GUD5$Date <- relevel(GUD5$Date, ref = "X8.Jun") # Relevel 
GUD5$Date <- relevel(GUD5$Date, ref = "X7.Jun") # Relevel 
GUD5$Date <- relevel(GUD5$Date, ref = "X6.Jun") # Relevel 
GUD5$Date <- relevel(GUD5$Date, ref = "X30.May") # Relevel 
GUD5$Date <- relevel(GUD5$Date, ref = "X29.May") # Relevel 
GUD5$Date <- relevel(GUD5$Date, ref = "X28.May") # Relevel 
GUD5$Date <- relevel(GUD5$Date, ref = "X27.May") # Relevel 

tgc <- summarySE(GUD5, measurevar="GUD", groupvars=c("Treatment", "Date"))

tgc <- summarySE(GUD5, measurevar="GUD", groupvars=c("Plot"))


ggplot(tgc, aes(x=Date, y=GUD, fill=Date)) + 
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin=GUD-se, ymax=GUD+se), width=.1) +
  guides(fill = FALSE) +
  labs(title = "", x = "Habitat", y = "Giving-Up Density (GUD)") +
  theme_bw(base_size = 15)

gd <- GUD5 %>% 
  group_by(Date) %>% 
  summarise(GUD = mean(GUD))
gd

ggplot(GUD5, aes(x = Date, y = GUD, color = Date, fill = Date)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  ggrepel::geom_text_repel(aes(label = Plot), color = "black", size = 2.5, segment.color = "grey") +
  geom_point() +
  guides(color = "none", fill = "none") +
  theme_bw(base_size = 15) +
  labs(title = "", x = "Date", y = "Giving-Up Density (GUD)")

ggplot(GUD5, aes(x = Date, y = GUD, fill = Date)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  ggrepel::geom_text_repel(aes(label = Plot), color = "black", size = 2.5, segment.color = "grey") +
  geom_point(aes(color = Treatment)) +
  guides(color = "none", fill = "none") +
  theme_bw(base_size = 15) +
  labs(title = "", x = "Date", y = "Giving-Up Density (GUD)")

ggplot(GUD5, aes(x=Treatment, y=GUD, fill = Treatment)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  xlab("Habitat") +
  ylab("Giving-Up Density (g)") + 
  theme(axis.title.x = element_text(size = 18, color = "white")) + 
  theme(axis.text.x = element_text(size = 7.5, color = "black", angle = 45, hjust = 1)) + 
  theme(axis.text.y = element_text(size = 7.5, color = "black")) +
  scale_fill_manual(values=c("peachpuff4","#66CC99", "#9999CC", "#CC6666", "khaki3")) +
  theme_minimal(base_size = 25) +theme(legend.position="none")

