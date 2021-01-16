# Radiological-anatomical distribution patterns of Jakob-Creutzfeldt disease in relation to functional cortical systems
Diploma thesis

The R scripts and Matlab code used for this thesis are found below.

**Note:** This is my first attempt at using Matlab and R code for research purposes. As a medical student, there is no teaching at all in programming. Therefore, many more experienced programmers will probably find several inconsistencies or even mistakes in my code. Thus, I would be more than happy to receive feedback and tips to improve my coding skills.

Please note, that the following codes cannot run as one script on RStudio. A couple of times .csv files were exported in the middle of a script and reformatted before they could be reread by the script!

## Extraktion of XYZ coordinates from spm12 volumes
After the whole data preprocessing, registration and segmentation steps were done via spm12 and MRIcroGL, the XYZ coordinates of the segmentations were extracted via MATLAB. 
In order for this to work, spm12 must be installed on MATLAB before.

This is the code used to extract the XYZ coordinates in MATLAB:
```
%% load data
 [DataFileImg, DataDirectoryImg] = uigetfile('*.nii','DWI series','MultiSelect', 'off');
% nii_img = load_untouch_nii(fullfile(DataDirectoryImg, DataFileImg));
% images = single(nii_img.img(:,:,:,1:end));
% [lines columns slices ImageNumber] = size(images);
% header = nii_img.hdr;

%%
volumeInfo=spm_vol( [DataDirectoryImg, DataFileImg]);
[intensityValues, xyzCoordinates]= spm_read_vols(volumeInfo);
linInt=find(intensityValues~=0);
xyzCoordinates=transpose(xyzCoordinates);
XYZmask=xyzCoordinates(linInt,:);
writematrix(XYZmask,'Matlab.csv');

%%
spm 
```
After executing the code, a matrix of coordinates should appear in Matlab. You now may export them as an .csv file.

## Label4MRI
In this step, the coordinates are labelled automatically to the AAL atlas or Brodmann areas defined in MRIcroGL.
This R package was created by yunshiuan and can be found here: https://github.com/yunshiuan/label4MRI

I used the following R script:
```
library("devtools")
library(label4MRI)
library(dplyr)

#Get data
m <- read.csv(".csv") #insert scan or patient ID

#Assign areas
Results <- t(mapply(
  FUN = mni_to_region_name, x = m$x, y = m$y, z = m$z, template = c("aal", "ba")))
write.csv(Results,".csv")
```
Export the labelled areas to another .csv file for further processing.

## Assignments to phylogenetic areas
I have assigned labelled areas to phylogenetic areas in my thesis using the following script:
As one will notice, the statistical methods and plots on which my results are based can be found later in the code.
```
library(dplyr)
library(ggplot2)
library(ggthemes)
library(extrafont)
library(showtext)
library(devtools)
showtext_auto() 

#Transform datasheet
m <- read.csv(".csv") #insert dataframe name
data <- filter(m, aal.distance == 0)
x <- table(data$aal.label)
count <- data %>% count(aal.label)
write.csv(count, "freq_of_areas.csv")

#Assign to phylogenetic areas
m <- read.csv("freq_of_areas.csv")
Iso <- filter(m, aal.label %in% c(
  "Precentral_L",
  "Precentral_R",
  "Frontal_Sup_L",
  "Frontal_Sup_R",
  "Frontal_Sup_Orb_L",
  "Frontal_Sup_Orb_R",
  "Frontal_Mid_L",
  "Frontal_Mid_R",
  "Frontal_Mid_Orb_L",
  "Frontal_Mid_Orb_R",
  "Frontal_Inf_Oper_L",
  "Frontal_Inf_Oper_R",
  "Frontal_Inf_Tri_L",
  "Frontal_Inf_Tri_R",
  "Frontal_Inf_Orb_L",
  "Frontal_Inf_Orb_R",
  "Rolandic_Oper_L",
  "Rolandic_Oper_R",
  "Supp_Motor_Area_L",
  "Supp_Motor_Area_R",
  "Frontal_Sup_Medial_L",
  "Frontal_Sup_Medial_R",
  "Frontal_Med_Orb_L",
  "Frontal_Med_Orb_R",
  "Rectus_L",
  "Rectus_R",
  "Insula_L",
  "Insula_R",
  "Amygdala_L",
  "Amygdala_R",
  "Calcarine_L",
  "Calcarine_R",
  "Cuneus_L",
  "Cuneus_R",
  "Lingual_L",
  "Lingual_R",
  "Occipital_Sup_L",
  "Occipital_Sup_R",
  "Occipital_Mid_L",
  "Occipital_Mid_R",
  "Occipital_Inf_L",
  "Occipital_Inf_R",
  "Fusiform_L",
  "Fusiform_R",
  "Postcentral_L",
  "Postcentral_R",
  "Parietal_Sup_L",
  "Parietal_Sup_R",
  "Parietal_Inf_L",
  "Parietal_Inf_R",
  "SupraMarginal_L",
  "SupraMarginal_R",
  "Angular_L",
  "Angular_R",
  "Precuneus_L",
  "Precuneus_R",
  "Paracentral_Lobule_L",
  "Paracentral_Lobule_R",
  "Caudate_L",
  "Caudate_R",
  "Putamen_L",
  "Putamen_R",
  "Pallidum_L",
  "Pallidum_R",
  "Thalamus_L",
  "Thalamus_R",
  "Heschl_L",
  "Heschl_R",
  "Temporal_Sup_L",
  "Temporal_Sup_R",
  "Temporal_Pole_Sup_L",
  "Temporal_Pole_Sup_R",
  "Temporal_Mid_L",
  "Temporal_Mid_R",
  "Temporal_Pole_Mid_L",
  "Temporal_Pole_Mid_R",
  "Temporal_Inf_L",
  "Temporal_Inf_R",
  "Cerebelum_Crus1_L",
  "Cerebelum_Crus1_R",
  "Cerebelum_Crus2_L",
  "Cerebelum_Crus2_R",
  "Cerebelum_3_L",
  "Cerebelum_3_R",
  "Cerebelum_4_5_L",
  "Cerebelum_4_5_R",
  "Cerebelum_6_L",
  "Cerebelum_6_R",
  "Cerebelum_7b_L",
  "Cerebelum_7b_R",
  "Cerebelum_8_L",
  "Cerebelum_8_R",
  "Cerebelum_9_L",
  "Cerebelum_9_R",
  "Cerebelum_10_L",
  "Cerebelum_10_R"
)
    )
Iso <- sum(Iso$n)

#Assign to phylogenetic areas
Meso <- filter(m, aal.label %in% c(
  "Cingulum_Ant_L",
  "Cingulum_Ant_R",
  "Cingulum_Mid_L",
  "Cingulum_Mid_R",
  "Cingulum_Post_L",
  "Cingulum_Post_R",
  "ParaHippocampal_L",
  "ParaHippocampal_R"
)
      )
Meso <- sum(Meso$n)

Allo <- filter(m, aal.label %in% c(
  "Olfactory_L",
  "Olfactory_R",
  "Hippocampus_L",
  "Hippocampus_R"
)
    )
Allo <- sum(Allo$n)

IsoR <- filter(m, aal.label %in% c(
  "Precentral_R",
  "Frontal_Sup_R",
  "Frontal_Sup_Orb_R",
  "Frontal_Mid_R",
  "Frontal_Mid_Orb_R",
  "Frontal_Inf_Oper_R",
  "Frontal_Inf_Tri_R",
  "Frontal_Inf_Orb_R",
  "Rolandic_Oper_R",
  "Supp_Motor_Area_R",
  "Frontal_Sup_Medial_R",
  "Frontal_Med_Orb_R",
  "Rectus_R",
  "Insula_R",
  "Amygdala_R",
  "Calcarine_R",
  "Cuneus_R",
  "Lingual_R",
  "Occipital_Sup_R",
  "Occipital_Mid_R",
  "Occipital_Inf_R",
  "Fusiform_R",
  "Postcentral_R",
  "Parietal_Sup_R",
  "Parietal_Inf_R",
  "SupraMarginal_R",
  "Angular_R",
  "Precuneus_R",
  "Paracentral_Lobule_R",
  "Caudate_R",
  "Putamen_R",
  "Pallidum_R",
  "Thalamus_R",
  "Heschl_R",
  "Temporal_Sup_R",
  "Temporal_Pole_Sup_R",
  "Temporal_Mid_R",
  "Temporal_Pole_Mid_R",
  "Temporal_Inf_R",
  "Cerebelum_Crus1_R",
  "Cerebelum_Crus2_R",
  "Cerebelum_3_R",
  "Cerebelum_4_5_R",
  "Cerebelum_6_R",
  "Cerebelum_7b_R",
  "Cerebelum_8_R",
  "Cerebelum_9_R",
  "Cerebelum_10_R"
))
IsoR <- sum(IsoR$n)

MesoR <- filter(m, aal.label %in% c(
  "Cingulum_Ant_R",
  "Cingulum_Mid_R",
  "Cingulum_Post_R",
  "ParaHippocampal_R"
)
    )
MesoR <- sum(MesoR$n)


AlloR <- filter(m, aal.label %in% c(
  "Olfactory_R",
  "Hippocampus_R"
)
    )
AlloR <- sum(AlloR$n)

IsoL <- filter(m, aal.label %in% c(
"Precentral_L",
"Frontal_Sup_L",
"Frontal_Sup_Orb_L",
"Frontal_Mid_L",
"Frontal_Mid_Orb_L",
"Frontal_Inf_Oper_L",
"Frontal_Inf_Tri_L",
"Frontal_Inf_Orb_L",
"Rolandic_Oper_L",
"Supp_Motor_Area_L",
"Frontal_Sup_Medial_L",
"Frontal_Med_Orb_L",
"Rectus_L",
"Insula_L",
"Amygdala_L",
"Calcarine_L",
"Cuneus_L",
"Lingual_L",
"Occipital_Sup_L",
"Occipital_Mid_L",
"Occipital_Inf_L",
"Fusiform_L",
"Postcentral_L",
"Parietal_Sup_L",
"Parietal_Inf_L",
"SupraMarginal_L",
"Angular_L",
"Precuneus_L",
"Paracentral_Lobule_L",
"Caudate_L",
"Putamen_L",
"Pallidum_L",
"Thalamus_L",
"Heschl_L",
"Temporal_Sup_L",
"Temporal_Pole_Sup_L",
"Temporal_Mid_L",
"Temporal_Pole_Mid_L",
"Temporal_Inf_L",
"Cerebelum_Crus1_L",
"Cerebelum_Crus2_L",
"Cerebelum_3_L",
"Cerebelum_4_5_L",
"Cerebelum_6_L",
"Cerebelum_7b_L",
"Cerebelum_8_L",
"Cerebelum_9_L",
"Cerebelum_10_L"
)
      )
IsoL <- sum(IsoL$n)

MesoL <- filter(m, aal.label %in% c(
  "Cingulum_Ant_L",
  "Cingulum_Mid_L",
  "Cingulum_Post_L",
  "ParaHippocampal_L"
)
      )
MesoL <- sum(MesoL$n)

AlloL <- filter(m, aal.label %in% c(
  "Olfactory_L",
  "Hippocampus_L"
)
      )
AlloL <- sum(AlloL$n)

r <- data.frame(Phylogenetic_areas = c(
  "Isocortex", "Isocortex", "Mesocortex", "Mesocortex", "Allocortex", "Allocortex"), 
  Volume = c(IsoL, IsoR, MesoL, MesoR, AlloL, AlloR), 
  Side = c("Left", "Right", "Left", "Right", "Left", "Right")) 
write.csv(r, "Auswertung.csv")  

#Plots after some reformatting in Auswertung.csv
r <- read.csv("Auswertung.csv") #insert dataframe name
r <- r[1:6, ]
p <- ggplot(data=r, aes(x = Phylogenetic_areas, y = Volume_total, fill = Side)) + 
    geom_bar(stat="identity", position=position_dodge())
p2 <- p + scale_x_discrete(limits=c("Isocortex", "Mesocortex", "Allocortex")) +
  theme_tufte(base_size=13,  base_family="Myriad Pro") +
  ylim(0, 50000)
p2


#Student's t-test and mean (SD)
adx <- 1
ady <- 2
x <- c(r$Vol_1[adx], r$Vol_2[adx], r$Vol_3[adx], r$Vol_4[adx], r$Vol_5[adx],r$Vol_6[adx],r$Vol_7[adx],r$Vol_8[adx],r$Vol_9[adx],r$Vol_10[adx])
y <- c(r$Vol_1[ady], r$Vol_2[ady], r$Vol_3[ady], r$Vol_4[ady], r$Vol_5[ady],r$Vol_6[ady],r$Vol_7[ady],r$Vol_8[ady],r$Vol_9[ady],r$Vol_10[ady])
t.test(x, y)
b <- x+y
sd(b)
mean(b)

#Consider percentage of Allocortex/Mesocortex
m2 <- read.csv("Auswertung.csv")
Iso <- m2$Volume_total[1:2]
Iso <- sum(Iso)
Meso <- (m2$Volume_total[3:4])
Meso <- sum(Meso)*17
Allo <- (m2$Volume_total[5:6])
Allo <- sum(Allo)*8.5
Meso <- Meso/Iso
Allo <- Allo/Iso
Iso <- Iso/Iso

rat <- data.frame(
  Phylogenetic = c("Isocortex", "Mesocortex", "Allocortex"), 
  Volume = c(Iso, Meso, Allo)
    )

iso2 <- m2$Volume_total2[1:20]
meso2 <- (m2$Volume_total2[21:40])*17
meso2 <- meso2/iso2
allo2 <- (m2$Volume_total2[41:60])*8.5
allo2 <- allo2/iso2
mean(allo2)
sd(allo2)


#Plot 
pp <- ggplot(data=rat, aes(x = Phylogenetic, y = Volume)) + 
  geom_bar(stat="identity", position=position_dodge(), fill = "powderblue")+
  geom_errorbar(aes(ymin=Volume-1.33281, ymax=Volume+1.33281), width=.2,
                position=position_dodge(.9))
pp2 <- pp + scale_x_discrete(limits=c("Isocortex", "Mesocortex", "Allocortex")) +
  theme_tufte(base_size=13,  base_family="Myriad Pro") +
  ylim(0, 3.5)
pp2
```
## Definition of long-range fibre tracts and assignment to Brodmann areas
I have allocated brain lesions on MRI automatically to long-range fibre tracts by defining them using the following script:
As one will notice, the statistical methods and plots on which my results are based can be found later in the code.

```
library(dplyr)
library(ggplot2)
library(ggthemes)
library(extrafont)
library(showtext)
library(devtools)
showtext_auto() 

#Assign Brodmann to tracts
m <- read.csv("Tracts.csv") #insert dataframe with Brodmann labels
data <- filter(m, ba.distance == 0)
x <- table(data$ba.label)
co <- data %>% count(ba.label)

SLF <- filter(co, ba.label %in% c(
  "Left-BA31",
  "Left-BA39",
  "Left-BA40",
  "Left-BA44",
  "Left-BA46",
  "Left-BA6",
  "Left-BA7",
  "Left-BA8",
  "Left-BA9",
  "Left-SensoryAssoc (5)",
  "Right-BA31",
  "Right-BA39",
  "Right-BA40",
  "Right-BA44",
  "Right-BA46",
  "Right-BA6",
  "Right-BA7",
  "Right-BA8",
  "Right-BA9",
  "Right-SensoryAssoc (5)"
)
)

ILF <- filter(co, ba.label %in% c(
  "Left-BA19",
  "Left-BA7",
  "Left-Parahip (36)",
  "Left-PrimVisual (17)",
  "Left-VisualAssoc (18)",
  "Right-BA19",
  "Right-BA7",
  "Right-Parahip (36)",
  "Right-PrimVisual (17)",
  "Right-VisualAssoc (18)"
)
)

UF <- filter(co, ba.label %in% c(
  "Left-BA10",
  "Left-BA24",
  "Left-BA25",
  "Left-BA32",
  "Left-BA38",
  "Left-BA45",
  "Left-BA46",
  "Left-BA47",
  "Left-Insula (13)",
  "Left-Parahip (36)",
  "Right-BA10",
  "Right-BA24",
  "Right-BA25",
  "Right-BA32",
  "Right-BA38",
  "Right-BA45",
  "Right-BA46",
  "Right-BA47",
  "Right-Insula (13)",
  "Right-Parahip (36)"
)
)
AF <- filter(co, ba.label %in% c(
  "Left-BA22",
  "Left-BA45",
  "Left-BA47",
  "Left-BA8",
  "Left-BA9",
  "Left-Insula (13)",
  "Left-Parahip (36)",
  "Left-PrimAuditory (41)",
  "Right-BA22",
  "Right-BA45",
  "Right-BA47",
  "Right-BA8",
  "Right-BA9",
  "Right-Insula (13)",
  "Right-Parahip (36)",
  "Right-PrimAuditory (41)"
)
)

MdLF <- filter(co, ba.label %in% c(
  "Left-BA19",
  "Left-BA7",
  "Left-PrimVisual (17)",
  "Left-VisualAssoc (18)",
  "Right-BA19",
  "Right-BA38",
  "Right-BA7",
  "Right-PrimVisual (17)",
  "Right-VisualAssoc (18)"
)
)

CB <- filter(co, ba.label %in% c(
  "Left-BA19",
  "Left-BA25",
  "Left-BA30",
  "Left-BA31",
  "Left-BA32",
  "Left-BA34",
  "Left-BA39",
  "Left-BA46",
  "Left-BA8",
  "Left-BA9",
  "Left-Insula (13)",
  "Left-Parahip (36)",
  "Right-BA11",
  "Right-BA25",
  "Right-BA30",
  "Right-BA31",
  "Right-BA32",
  "Right-BA34",
  "Right-BA39",
  "Right-BA46",
  "Right-BA8",
  "Right-BA9",
  "Right-Insula (13)",
  "Right-Parahip (36)"
)
)

ND <- filter(co, ba.label %in% c(
  "Left-Amygdala (53)",
  "Left-BA20",
  "Left-BA21",
  "Left-BA23",
  "Left-Caudate (48)",
  "Left-Fusiform (37)",
  "Left-GlobPal (51)",
  "Left-Hippocampus (54)",
  "Left-Hypothalamus (55)",
  "Left-NucAccumb (52)",
  "Left-PrimMotor (4)",
  "Left-PrimSensory (1)",
  "Left-Putamen (49)",
  "Left-Thalamus (50)",
  "Right-Amygdala (53)",
  "Right-BA20",
  "Right-BA23",
  "Right-Caudate (48)",
  "Right-Fusiform (37)",
  "Right-GlobPal (51)",
  "Right-Hippocampus (54)",
  "Right-Hypothalamus (55)",
  "Right-NucAccumb (52)",
  "Right-PrimMotor (4)",
  "Right-PrimSensory (1)",
  "Right-Putamen (49)",
  "Right-Thalamus (50)"
)
)

#Get accumulated volumes of tracts

SLF <- sum(SLF$n)
ILF <- sum(ILF$n)
UF <- sum(UF$n)
AF <- sum(AF$n)
MdLF <- sum(MdLF$n)
CB <- sum(CB$n)
ND <- sum(ND$n)

m <- data.frame(Tracts=c(
  "SLF","ILF","UF","AF","MdLF","CB","ND"), 
  Volume=c(SLF,ILF,UF,AF,MdLF,CB,ND)
)

# Plot volumes
pp <- ggplot(data=m, aes(x = Tracts, y = Volume)) + geom_bar(stat="identity", position=position_dodge())
pp2 <- pp + scale_x_discrete(limits=c("SLF","ILF","UF","AF","MdLF","CB","ND"))
pp2

#Choose patient from table Tracts.csv
Pat <- m$Vol_10
TV <- m$Voxel_total[10]
SLF <- Pat[1]
ILF <- Pat[2]
UF <- Pat[3]
AF <- Pat[4]
MdLF <- Pat[5]
CB <- Pat[6]
ND <- Pat[7]

#Get Ratios Tracts/Total_volume
SLF <- SLF/TV
ILF <- ILF/TV
UF <- UF/TV
AF <- AF/TV
MdLF <- MdLF/TV
CB <- CB/TV
ND <- ND/TV

m3 <- data.frame(Tracts=c(
  "SLF","ILF","UF","AF","MdLF","CB","ND"), 
  Percentage=c(SLF,ILF,UF,AF,MdLF,CB,ND)
)

sd(m$Percentages_MdLF)

#Barplots
ppp <- ggplot(data=m, aes(x = Tracts, y = Percentage_mean)) + 
    geom_bar(stat="identity", fill = "powderblue", position=position_dodge()) + 
    theme_tufte(base_size=13,  base_family="Myriad Pro") +
  geom_rangeframe ()
ppp2 <- ppp + scale_x_discrete(limits=c("SLF","ILF","UF","AF","MdLF","CB","ND")) + 
  ylim(0,1)

ppp2


#Means, SD of Percentage of sJCD lesions
m <- read.csv("Widni.csv")
data <- filter(m, ba.distance == 0)
x <- table(data$ba.label)
co <- data %>% count(ba.label)

no_as <- filter(co, ba.label %in% c(
  "Left-Amygdala (53)",
  "Left-BA20",
  "Left-BA21",
  "Left-BA23",
  "Left-Caudate (48)",
  "Left-Fusiform (37)",
  "Left-GlobPal (51)",
  "Left-Hippocampus (54)",
  "Left-Hypothalamus (55)",
  "Left-NucAccumb (52)",
  "Left-PrimMotor (4)",
  "Left-PrimSensory (1)",
  "Left-Putamen (49)",
  "Left-Thalamus (50)",
  "Right-Amygdala (53)",
  "Right-BA20",
  "Right-BA23",
  "Right-Caudate (48)",
  "Right-Fusiform (37)",
  "Right-GlobPal (51)",
  "Right-Hippocampus (54)",
  "Right-Hypothalamus (55)",
  "Right-NucAccumb (52)",
  "Right-PrimMotor (4)",
  "Right-PrimSensory (1)",
  "Right-Putamen (49)",
  "Right-Thalamus (50)"
)
    )

as <- filter(co, ba.label %in% c(
  "Left-BA10",
  "Left-BA11",
  "Left-BA19",
  "Left-BA21",
  "Left-BA22",
  "Left-BA24",
  "Left-BA25",
  "Left-BA30",
  "Left-BA31",
  "Left-BA32",
  "Left-BA34",
  "Left-BA38",
  "Left-BA39",
  "Left-BA40",
  "Left-BA44",
  "Left-BA45",
  "Left-BA46",
  "Left-BA47",
  "Left-BA6",
  "Left-BA7",
  "Left-BA8",
  "Left-BA9",
  "Left-Insula (13)",
  "Left-Parahip (36)",
  "Left-PrimAuditory (41)",
  "Left-PrimVisual (17)",
  "Left-SensoryAssoc (5)",
  "Left-VisualAssoc (18)",
  "Right-BA10",
  "Right-BA11",
  "Right-BA19",
  "Right-BA21",
  "Right-BA22",
  "Right-BA24",
  "Right-BA25",
  "Right-BA30",
  "Right-BA31",
  "Right-BA32",
  "Right-BA34",
  "Right-BA38",
  "Right-BA39",
  "Right-BA40",
  "Right-BA44",
  "Right-BA45",
  "Right-BA46",
  "Right-BA47",
  "Right-BA6",
  "Right-BA7",
  "Right-BA8",
  "Right-BA9",
  "Right-Insula (13)",
  "Right-Parahip (36)",
  "Right-PrimAuditory (41)",
  "Right-PrimVisual (17)",
  "Right-SensoryAssoc (5)",
  "Right-VisualAssoc (18)"
)
    )

as <- sum(as$n)
no_as <- sum(no_as$n)
k <- data.frame(As = as, No_as = no_as)
k

m <- read.csv("Tracts.csv")
Assoc <- m$Assoc
Nassoc <- m$No_assoc
ass <- mean(Assoc)
nass <- mean(Nassoc)

Ratioass <- ass/(ass+nass)
Rationass <- 1 - Ratioass

#Bar chart Percentage
val <- data.frame(Tract = c("Tracts", "ND"), Percentage = c(Ratioass, Rationass))
Lesions <- c("No association cortices", "Assigned to association cortices")
bp<- ggplot(val, aes(x="Percentage", y=Percentage, fill = Lesions))+
  geom_bar(width = 1, stat = "identity") +
  theme_tufte(base_size=13,  base_family="Myriad Pro") +
  geom_rangeframe(y = c(0,1))
bp + scale_fill_manual(values=c("powderblue", "tomato"))
bp

#sd
Ass <- m$Assoc
Nass <- m$No_assoc
Ass <- Ass/(Ass+Nass)
Nass <- Nass/(Ass+Nass)
sd(Ass)
sd(Nass)

```
## Definition of cortical functional regions 
In the following R script I have assigned JCD lesions to cortical systems which are based on the dual origin concept.
The statistical methods and plots on which my results are based can be found later in the code.

```
library(dplyr)
library(ggplot2)
library(ggthemes)
library(extrafont)
library(showtext)
library(devtools)
library(label4MRI)
library(psych)
showtext_auto() 

#Transform datasheet
m <- read.csv("Podr.csv")
data <- filter(m, ba.distance == 0)
x <- table(data$ba.label)
count <- data %>% count(ba.label)

# sum(count$n)    Sum absolute volumes
# z <- data.frame(Name = c(1,10), Volume_total = c("9237","1500","12127","9180","13135","7422","12331","21497","12796","8816"))
# write.csv(z, "Volumes absolute.csv")    create dataframe


#Assign Auditory system
Au_core_L <- filter(count, ba.label %in%
  "Left-PrimAuditory (41)")
Au_core_L <- sum(Au_core_L$n)
Au_core_R <- filter(count, ba.label %in% 
  "Right-PrimAuditory (41)")
Au_core_R <- sum(Au_core_R$n)
Au_belt_L <- filter(count, ba.label %in% c(
  "Left-BA10",
  "Left-BA22"
  ))
Au_belt_L <- sum(Au_belt_L$n)
Au_belt_R <- filter(count, ba.label %in% c(
  "Right-BA10",
  "Right-BA22"
  ))
Au_belt_R <- sum(Au_belt_R$n)
Au_root_L <- filter(count, ba.label %in% c(
  "Left-BA38",
  "Left-BA22"
  ))
Au_root_L <- sum(Au_root_L$n)
Au_root_R <- filter(count, ba.label %in% c(
  "Right-BA38",
  "Left-BA22"
  ))
Au_root_R <- sum(Au_root_R$n)
Au_MM_L <- filter(count, ba.label %in% c(
  "Left-PrimAuditory (41)",
  "Left-Insula (13)",
  "Left-BA23",
  "Left-BA30",
  "Left-Parahip (36)"
  ))
Au_MM_L <- sum(Au_MM_L$n)
Au_MM_R <- filter(count, ba.label %in% c(
  "Right-PrimAuditory (41)",
  "Right-Insula (13)",
  "Right-BA23",
  "Right-BA30",
  "Right-Parahip (36)"
  ))
Au_MM_R <- sum(Au_MM_R$n)
Au_PL_L <- filter(count, ba.label %in% c(
  "Left-BA23",
  "Left-BA30",
  "Left-Parahip (36)",
  "Right-Insula (13)"
  ))
Au_PL_L <- sum(Au_PL_L$n)
Au_PL_R <- filter(count, ba.label %in% c(
  "Right-Insula (13)",
  "Right-BA23",
  "Right-BA30",
  "Right-Parahip (36)"
  ))
Au_PL_R <- sum(Au_PL_R$n)
Au_PF_L <- filter(count, ba.label %in% c(
  "Left-Insula (13)",
  "Left-BA11",
  "Left-BA10",
  "Left-BA25",
  "Left-BA47",
  "Left-BA8",
  "Left-BA9",
  "Left-BA46",
  "Left-BA45"
  ))
Au_PF_L <- sum(Au_PF_L$n)
Au_PF_R <- filter(count, ba.label %in% c(
  "Right-Insula (13)",
  "Right-BA11",
  "Right-BA10",
  "Right-BA25",
  "Right-BA47",
  "Right-BA8",
  "Right-BA9",
  "Right-BA46",
  "Right-BA45"
  ))
Au_PF_R <- sum(Au_PF_R$n)
#Assign Visual system
Vi_Core_L <- filter(count, ba.label %in% c(
  "Left-PrimVisual (17)",
  "Left-VisualAssoc (18)",
  "Left-BA19"
  ))
Vi_Core_L <- sum(Vi_Core_L$n)
Vi_Core_R <- filter(count, ba.label %in% c(
  "Right-PrimVisual (17)",
  "Right-VisualAssoc (18)",
  "Right-BA19"
  ))
Vi_Core_R <- sum(Vi_Core_R$n)
Vi_belt_L <- filter(count, ba.label %in% c(
  "Left-BA21",
  "Left-BA22",
  "Left-VisualAssoc (18)",
  "Left-BA19",
  "Left-BA40"
  ))
Vi_belt_L <- sum(Vi_belt_L$n)
Vi_belt_R <- filter(count, ba.label %in% c(
  "Right-BA21",
  "Right-BA22",
  "Right-VisualAssoc (18)",
  "Right-BA19",
  "Right-BA40"
  ))
Vi_belt_R <- sum(Vi_belt_R$n)
Vi_root_L <- filter(count, ba.label %in% c(
  "Left-VisualAssoc (18)",
  "Left-BA19"
  ))
Vi_root_L <- sum(Vi_root_L$n)
Vi_root_R <- filter(count, ba.label %in% c(
  "Right-VisualAssoc (18)",
  "Right-BA19"
  ))
Vi_root_R <- sum(Vi_root_R$n)
Vi_MM_L <- filter(count, ba.label %in% c(
  "Left-BA7",
  "Left-Insula (13)",
  "Left-BA23",
  "Left-BA30",
  "Left-Parahip (36)",
  "Left-PrimAuditory (41)"
  ))
Vi_MM_L <- sum(Vi_MM_L$n)
Vi_MM_R <- filter(count, ba.label %in% c(
  "Right-BA7",
  "Right-Insula (13)",
  "Right-BA23",
  "Right-BA30",
  "Right-Parahip (36)",
  "Right-PrimAuditory (41)"
  ))
Vi_MM_R <- sum(Vi_MM_R$n)
Vi_PL_L <- filter(count, ba.label %in% c(
  "Left-BA23",
  "Left-BA30",
  "Left-Parahip (36)",
  "Right-Insula (13)"
  ))
Vi_PL_L <- sum(Vi_PL_L$n)
Vi_PL_R <- filter(count, ba.label %in% c(
  "Right-Insula (13)",
  "Right-BA23",
  "Right-BA30",
  "Right-Parahip (36)"
  ))
Vi_PL_R <- sum(Vi_PL_R$n)
Vi_PF_L <- filter(count, ba.label %in% c(
  "Left-BA8",
  "Left-BA45",
  "Left-BA47",
  "Left-BA11",
  "Left-BA6",
  "Left-BA9",
  "Left-BA46"
  ))
Vi_PF_L <- sum(Vi_PF_L$n)
Vi_PF_R <- filter(count, ba.label %in% c(
  "Right-BA8",
  "Right-BA45",
  "Right-BA47",
  "Right-BA11",
  "Right-BA6",
  "Right-BA9",
  "Right-BA46"
  ))
Vi_PF_R <- sum(Vi_PF_R$n)
#Assign Somatosensory system
SS_Core_L <- filter(count, ba.label %in%
  "Left-PrimSensory (1)")
SS_Core_L <- sum(SS_Core_L$n)
SS_Core_R <- filter(count, ba.label %in% 
  "Right-PrimSensory (1)")
SS_Core_R <-sum(SS_Core_R$n)
SS_root_L <- filter(count, ba.label %in%
  "Left-SensoryAssoc (5)")
SS_root_L <- sum(SS_root_L$n)
SS_root_R <- filter(count, ba.label %in%
  "Right-SensoryAssoc (5)")
SS_root_R <- sum(SS_root_R$n)
SS_belt_L <- filter(count, ba.label %in% c(
  "Left-PrimSensory (1)",
  "Left-SensoryAssoc (5)",
  "Left-BA40"
  ))
SS_belt_L <- sum(SS_belt_L$n)
SS_belt_R <- filter(count, ba.label %in% c(
  "Right-PrimSensory (1)",
  "Right-SensoryAssoc (5)",
  "Right-BA40"
  ))
SS_belt_R <- sum(SS_belt_R$n)
SS_MM_L <- filter(count, ba.label %in% c(
  "Left-BA39",
  "Left-PrimAuditory (41)",
  "Left-Insula (13)",
  "Left-BA23",
  "Left-BA30",
  "Left-Parahip (36)"
  ))
SS_MM_L <- sum(SS_MM_L$n)
SS_MM_R <- filter(count, ba.label %in% c(
  "Right-BA39",
  "Right-PrimAuditory (41)",
  "Right-Insula (13)",
  "Right-BA23",
  "Right-BA30",
  "Right-Parahip (36)"
  ))
SS_MM_R <- sum(SS_MM_R$n)
SS_PL_L <- filter(count, ba.label %in% c(
  "Left-BA23",
  "Left-BA30",
  "Left-Parahip (36)",
  "Left-BA11",
  "Left-Insula (13)"
  ))
SS_PL_L <- sum(SS_PL_L$n)
SS_PL_R <- filter(count, ba.label %in% c(
  "Right-BA23",
  "Right-BA30",
  "Right-Parahip (36)",
  "Right-BA11",
  "Right-Insula (13)"
  ))
SS_PL_R <- sum(SS_PL_R$n)
SS_PF_L <- filter(count, ba.label %in% c(
  "Left-BA47",
  "Left-BA9",
  "Left-BA46"
  ))
SS_PF_L <- sum(SS_PF_L$n)
SS_PF_R <- filter(count, ba.label %in% c(
  "Right-BA47",
  "Right-BA9",
  "Right-BA46"
  ))
SS_PF_R <- sum(SS_PF_R$n)
#Assign Motor System
Mo_Core_L <- filter(count, ba.label %in%
  "Left-PrimMotor (4)")
Mo_Core_L <-sum(Mo_Core_L$n)
Mo_Core_R <- filter(count, ba.label %in%
  "Right-PrimMotor (4)")
Mo_Core_R <- sum(Mo_Core_R$n)
Mo_root_L <- filter(count, ba.label %in% c(
  "Left-BA6",
  "Left-PrimMotor (4)"
  ))
Mo_root_L <- sum(Mo_root_L$n)
Mo_root_R <- filter(count, ba.label %in% c(
  "Right-BA6",
  "Right-PrimMotor (4)"
  ))
Mo_root_R <- sum(Mo_root_R$n)
Mo_belt_L <- filter(count, ba.label %in% c(
  "Left-BA6",
  "Left-SensoryAssoc (5)",
  "Left-BA7",
  "Left-BA39",
  "Left-BA40"
))
Mo_belt_L <- sum(Mo_belt_L$n)
Mo_belt_R <- filter(count, ba.label %in% c(
  "Right-BA6",
  "Right-SensoryAssoc (5)",
  "Right-BA7",
  "Right-BA39",
  "Right-BA40"
  ))
Mo_belt_R <-sum(Mo_belt_R$n)
Mo_MM_L <- filter(count, ba.label %in% c(
  "Left-Insula (13)",
  "Left-PrimAuditory (41)",
  "Left-BA31",
  "Left-BA39",
  "Left-BA40"
  ))
Mo_MM_L <- sum(Mo_MM_L$n)
Mo_MM_R <- filter(count, ba.label %in% c(
  "Right-Insula (13)",
  "Right-PrimAuditory (41)",
  "Right-BA31",
  "Right-BA39",
  "Right-BA40"
  ))
Mo_MM_R <- sum(Mo_MM_R$n)
Mo_PL_L <- filter(count, ba.label %in% c(
  "Left-Insula (13)",
  "Left-BA23",
  "Left-BA30",
  "Left-BA31"
  ))
Mo_PL_L <- sum(Mo_PL_L$n)
Mo_PL_R <- filter(count,ba.label %in% c(
  "Right-Insula (13)",
  "Right-BA23",
  "Right-BA30",
  "Right-BA31"
  ))
Mo_PL_R <- sum(Mo_PL_R$n)
Mo_PF_L <- filter(count, ba.label %in% c(
  "Left-BA46",
  "Left-BA9",
  "Left-BA8",
  "Left-BA45",
  "Left-BA47",
  "Left-BA11"
  ))
Mo_PF_L <- sum(Mo_PF_L$n)
Mo_PF_R <- filter(count, ba.label %in% c(
  "Right-BA46",
  "Right-BA9",
  "Right-BA8",
  "Right-BA45",
  "Right-BA47",
  "Right-BA11"
  ))
Mo_PF_R <- sum(Mo_PF_R$n)

#create data.frame
m_left1 <- data.frame(    #change number 1 to 10 after m_left
  System = c(
  "A_core","V_core", "S_core", "M_core",
  "A_belt","V_belt", "S_belt", "M_belt",
  "A_root", "V_root", "S_root", "M_root",
  "A_MM", "V_MM", "S_MM", "M_MM",
  "A_PF", "V_PF", "S_PF", "M_PF",
  "A_PL", "V_PL", "S_PL", "M_PL"
  ), Volumes = c(
    Au_core_L, Vi_Core_L, SS_Core_L, Mo_Core_L, 
    Au_belt_L, Vi_belt_L, SS_belt_L, Mo_belt_L,
    Au_root_L, Vi_root_L, SS_root_L, Mo_root_L,
    Au_MM_L, Vi_MM_L, SS_MM_L, Mo_MM_L,
    Au_PL_L, Vi_PL_L, SS_PL_L, Mo_PL_L,
    Au_PF_L, Vi_PF_L, SS_PF_L, Mo_PF_L
    ))
left <- data.frame(
  m_left1, 
  m_left2, 
  m_left3,
  m_left4,
  m_left5,
  m_left6,
  m_left7,
  m_left8,
  m_left9,
  m_left10)
write.csv(left, "Left Systems.csv")

m_right1 <- data.frame(  #change number 1 to 10 after m_right
  System = c(
  "A_core","V_core", "S_core", "M_core",
  "A_belt","V_belt", "S_belt", "M_belt",
  "A_root", "V_root", "S_root", "M_root",
  "A_MM", "V_MM", "S_MM", "M_MM",
  "A_PF", "V_PF", "S_PF", "M_PF",
  "A_PL", "V_PL", "S_PL", "M_PL"
  ), Volumes = c(
    Au_core_R, Vi_Core_R, SS_Core_R, Mo_Core_R,
    Au_belt_R, Vi_belt_R, SS_belt_R, Mo_belt_R,
    Au_root_R, Vi_root_R, SS_root_R, Mo_root_R,
    Au_MM_R, Vi_MM_R, SS_MM_R, Mo_MM_R,
    Au_PL_R, Vi_PL_R, SS_PL_R, Mo_PL_R,
    Au_PF_R, Vi_PF_R, SS_PF_R, Mo_PF_R
    ))
right <- data.frame(
  m_right1,
  m_right2,
  m_right3, 
  m_right4, 
  m_right5, 
  m_right6,
  m_right7,
  m_right8,
  m_right9,
  m_right10)
write.csv(right, "Right Systems.csv")


#Bar plot absolute
n <- read.csv("Right Systems.csv") # or "Left Systems.csv"
#reorder
n$Subsystem <- factor(n$Subsystem, levels=c("core", "belt", "root", "MM", "PL", "PF"))

pp <- ggplot(data=n, aes(fill = Subsystem, x = System, y = Volumes_avg)) + 
  geom_bar(stat="identity", position= "dodge") +
  geom_errorbar(aes(ymin=Volumes_avg, ymax=Volumes_avg+SD, color=factor(Subsystem)), 
                size=.1,
                width=0,                    # Width of the error bars
                position=position_dodge(.9)) + 
  scale_color_manual("Subsystems", values=c(
    "#3d9093","#00bfc4","#66d9dc","#938e46","#d87c54","#853d93"))+
  scale_fill_manual("Subsystems", values= c(
    "core" = "#3d9093", "belt" = "#00bfc4", "root" = "#66d9dc", "MM"="#938e46", "PL"="#d87c54", "PF"="#853d93"))+
  theme_tufte(base_size=13,  base_family="Myriad Pro")+
  ylim(0, 4500)
pp

# Bar plot percentages
n <- read.csv("Right Systems.csv") # or "Left Systems.csv"
#reorder
n$Subsystem <- factor(n$Subsystem, levels=c("core", "belt", "root", "MM", "PL", "PF"))

pp2 <- ggplot(data=n, aes(fill = Subsystem, x = System, y = Percentage_avg)) + 
  geom_bar(stat="identity", position= "dodge") +
  geom_errorbar(aes(ymin=Percentage_avg, ymax=Percentage_avg+SD_perc, color=factor(Subsystem)), 
                size=.1,
                width=0,                    # Width of the error bars
                position=position_dodge(.9)) + 
  scale_color_manual("Subsystems", values=c("#3d9093","#00bfc4","#66d9dc","#938e46","#d87c54","#853d93"))+
  scale_fill_manual("Subsystems", values= c("core" = "#3d9093", "belt" = "#00bfc4", "root" = "#66d9dc", "MM"="#938e46", "PL"="#d87c54", "PF"="#853d93"))+
  theme_tufte(base_size=13,  base_family="Myriad Pro")+
  ylim(0, 0.4)
pp2

# Statistics
b <- read.csv("Both_sides.csv")

au <- data.frame(Side = c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), 
                 Volumes = c(4062,6953,0,7923,3156,2053,3634,777,8061,8658,750,5,3159,3456,5310,3913,3288,4378,6334,5000))

Mo <- data.frame(Side = c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
                 Volumes = c(2798,5756,0,5793,2397,1087,4579,1034,7136,11804,474,0,3597,2172,6359,5245,5320,2263,6503,3593))

SS <- data.frame(Side = c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
                 Volumes = c(2012,2987,0,3621,1513,663,2443,535,4127,5927,210,0,1945,1630,3778,3192,2973,1693,3258,2636))

Vi <- data.frame(Side = c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
                 Volumes = c(5419,4022,0,5644,2263,1200,6904,806,5338,10460,245,0,2444,2109,4229,7099,5736,1820,6821,8523))

describeBy(b$Volumes_avg, b$Side)                 
t.test(Vi$Volumes~Vi$Side, var.equal = FALSE, alternative ="two.sided")


au <- c(4062,6953,0,7923,3156,2053,3634,777,8061,8658,750,5,3159,3456,5310,3913,3288,4378,6334,5000)
Mo <- c(2798,5756,0,5793,2397,1087,4579,1034,7136,11804,474,0,3597,2172,6359,5245,5320,2263,6503,3593)
SS <- c(2012,2987,0,3621,1513,663,2443,535,4127,5927,210,0,1945,1630,3778,3192,2973,1693,3258,2636)
Vi <- c(5419,4022,0,5644,2263,1200,6904,806,5338,10460,245,0,2444,2109,4229,7099,5736,1820,6821,8523)

au <- data.frame(System = 1, Volume = au)
Mo <- data.frame(System = 2, Volume = Mo)
SS <- data.frame(System = 3, Volume = SS)
Vi <- data.frame(System = 4, Volume = Vi)

z <- rbind(au, Mo)
y <- rbind(SS, Vi)
x <- rbind(z, y)

describeBy(x$Volume, x$System)

stat <- aov(x$Volume~x$System)
hist(rstandard(stat))
kruskal.test(x$Volume~x$System)

#The following dataframes were copied from the former .csv file
core <- c(51,0,0,16,0,0,39,1,16,51,891,9,0,267,207,0,1492,0,116,902,52,8,0,51,75,69,37,15,83,283,15,7,0,40,88,51,2,41,109,157,
          0,0,19,0,239,189,82,19,56,135,0,0,92,80,23,1305,832,0,700,2038,0,0,51,18,320,114,100,44,79,24,0,0,121,20,188,37,42,30,48,35)
belt <- c(896,1634,0,1673,810,523,1117,23,1884,1673,1839,115,0,456,352,224,2382,41,985,2439,195,121,0,192,133,133,449,67,475,1089,478,861,0,803,522,254,1675,276,1463,3631,
          270,0,706,732,1051,989,468,1103,1211,1068,0,0,565,229,1034,2494,1259,187,1490,2484,0,0,468,70,1122,774,456,228,535,118,55,0,908,294,1660,1841,1858,294,1773,652)
root <- c(355,0,0,27,93,100,301,0,363,492,853,3,0,260,207,0,1448,0,116,873,0,1,0,0,0,0,0,11,15,82,109,157,0,458,340,182,25,178,635,1284,
          0	,0	,27	,153	,135	,286	,0	,288	,464	,417,
          0	,0	,91,	80,	23	,1264,	832,	0	,680	,1687,
          0	,0	,0	,0	,17,	0,	52,	5,	34,	0,
          0,	0	,354,	86,	725,	199,	418,	140,	696,	93)
MM <- c(533,	480,	0,	905,	289,	62,	258,	174,	862,	1015,
        562	,549,	0	,915	,289	,62	,552	,244,	931,	1667,
        745	,1009,	0	,1139	,501	,121	,1204,	191,	1338,	2061,
        917	,1088	,0	,929	,551	,185	,1703	,231	,1741	,2831,
        11,	0	,340	,335	,689	,577	,801	,382	,635	,748,
        53	,0	,352	,390	,773	,826	,1175,	382,	914,	767,
        24	,0	,586	,456,	926,	1347,	1553	,382	,1025	,1229,
        120	,0	,1043,	508,	1819,	1939	,1563,	602,	1502,	1223)
PL <- c(1800,	4362	,0	,4437,	1623,	938,	1624,	330	,4285	,4531,
        847,	2869	,0	,2881,	867,	484	,735	,272	,2539	,3683,
        538,	1368,	0	,1350,	515,	278	,534,	78,	1349,	1424,
        753,	2719,	0	,2463	,615	,353	,712,	135,	2013,	2556,
        458,	5,	1746,	1901	,2746,	1484,	1218,	2223	,3389,	2019,
        181,	0	,1023	,995	,1926,	822,	919	,888	,2458,	934,
        175,	0,	519,	751,	943,	566,	90	,649,	979,	652,
        181,	0	,790,	929,	1389,	660,	543,	778,	1810	,876)
PF <- c(427,	477,	0	,865	,341	,430,	295	,249	,651,	896,
        427,	477,	0	,865,	341,	430,	295	,249	,651,	896,
        482,	480,	0,	889,	289,	62,	219,	173,	867,	988,
        526	,924	,0	,1100,	281,	62	,462,	173,	1175	,1345,
        11,	0	,321,	335,	450,	388,	719	,363,	579,	613,
        11	,0	,321	,335,	450,	388,	719,	363,	579,	613,
        11	,0	,321,	335,	450,	391,	722,	385,	606,	613,
        118,	0	,381	,335,	578,	569,	896,	419,	674	,714)
core <- data.frame(Sub = 1, Volume = core)
belt <- data.frame(Sub = 2, Volume = belt)
root <- data.frame(Sub = 3, Volume = root)
MM <- data.frame(Sub = 4, Volume = MM)
PL <- data.frame(Sub = 5, Volume = PL)
PF <- data.frame(Sub = 6, Volume = PF)

a <- rbind(core, belt)
b <- rbind(root, MM)
c <- rbind(PL, PF)
x <- rbind(a, b)
y <- rbind(x, c)

#statistics
stat <- aov(y$Volume~y$Sub)

hist(rstandard(stat))
summary(stat)
plot(stat, 2)
kruskal.test(y$Volume~y$Sub)
pairwise.wilcox.test(y$Volume, y$Sub, paired = FALSE, p.adjust.method = "bonferroni")

d <- c(37, 45, 415, 33, 35)
median(d)
```
