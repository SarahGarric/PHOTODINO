
# Libraires ---------------------------------------------------------------

library(tidyverse)
library(readxl)
library(lubridate)
library(phytotools)
library(gridExtra)
library(r2symbols)
library(forcats)
library(cowplot)
library(phytotools)
library(ggprism)
library(scales)
library(rstatix)
library(ggpubr)

# Graphical parameters ----------------------------------------------------

Taille_Point = 6
Taille_Point_Shift = 4
Taille_Titre = 17
Taille_Axe = 15
Taille_Strip = 17
Taille_Titre = 22
Taille_Axe_Shift = 18
Couleur ="#98a1b3"
Taille_Ligne=1
Taille_SD=0.4
Taille_Vline = 0.8
Taille_Hline = 0.9
Taille_Rectangle = 1.5
Taille_Ticks = .2
Epaisseur_Ticks = 1
Taille_Ticks_shift = .3
Epaisseur_Ticks_shift = 1.2
Couleur_Bordure = "gray47"
Couleur_Titre = "gray37"

Gris_Clair = "gray80"
jaune_clair <- "#ffa14d"
or <- "#f0d063"


# Courbes PE --------------------------------------------------------------



files_names <- list.files("Data/Courbes_PE")
files_names_full <- list.files("Data/Courbes_PE", full.names = T)
files <- map_dfr(files_names_full, .f=read_delim, 
                 delim = ";", escape_double = FALSE, locale = locale(), 
                 trim_ws = TRUE, .id = "source" )
nb_files <- length(files_names)
data_names <- vector("list",length=nb_files)

for (i in 1 : nb_files) {
  data_names[i] <- strsplit(files_names[i], split=".CSV")
}



Courbes_PE <- tibble()


for (i in 1:nb_files) {
  Fichier <- subset(files, source==i)
  Fichier$Type <-  as.character(data_names[[i]])
  Courbes_PE <- bind_rows(Courbes_PE, Fichier)
}

Courbes_PE <- Courbes_PE[,c(1:9, 12,17,27,37)]
Courbes_PE <- na.omit(Courbes_PE)

Courbes_PE <- separate(Courbes_PE, Type, c(NA, NA,"Espece", "Food", "Jour", "Replicat" ))

Courbes_PE$Yield <- (Courbes_PE$`Fm'3`-Courbes_PE$F3)/Courbes_PE$`Fm'3`
Courbes_PE$rETR <- (Courbes_PE$PAR*Courbes_PE$Yield)

files_names <- list.files("Data/Sigma")
files_names_full <- list.files("Data/Sigma", full.names = T)
files <- map_dfr(files_names_full, .f=read_delim, 
                 delim = ";", escape_double = FALSE, locale = locale(), 
                 trim_ws = TRUE, .id = "source" )
nb_files <- length(files_names)
data_names <- vector("list",length=nb_files)


for (i in 1 : nb_files) {
  data_names[i] <- strsplit(files_names[i], split=".csv")
}



Sigma <- tibble()

for (i in 1:nb_files) {
  Fichier <- subset(files, source==i)
  Fichier$Type <-  as.character(data_names[[i]])
  Sigma <- bind_rows(Sigma, Fichier)
}

Sigma <- Sigma[,c(7,9, 13,16)]
Sigma <- na.omit(Sigma)

Sigma <- separate(Sigma[,3:4], Type,  c(NA,"Espece", "Food", "Jour", "Replicat" ))


Fv_Fm <- read_excel("Data/Fv_Fm.xlsx")

Courbes_PE <- merge(Courbes_PE, Sigma, c( "Espece", "Food", "Jour", "Replicat"))
Courbes_PE <- merge(Courbes_PE, Fv_Fm, c( "Espece", "Food", "Jour", "Replicat" ))


Courbes_PE$Sigma <- as.numeric(Courbes_PE$Sigma)


Courbes_PE$PAR_II <- Courbes_PE$PAR*Courbes_PE$Sigma*0.6022

Courbes_PE$ETR_II <- (Courbes_PE$PAR_II*Courbes_PE$Yield)/Courbes_PE$Fv_Fm

## Model calcul --------------

PE_model <- tibble()
PE_curves <- tibble()
PE_Parameters <- tibble()

for (i in c(6)) {
  for (k in c(0)) {
    for (j in c("A")) {
      PE <- subset(Courbes_PE, Espece=="P" & Jour==i & Food==k & Replicat==j)
      data <- fitEP(PE$PAR, PE$ETR_II, normalize = F, 
                    lowerlim = c(4.8, 0, 900), upperlim = c(10, 800, 1100))
      PE$model <- PE$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE$PAR+(1/data$alpha[1]))
      ETRmax=data$ps[1]
      Ek=ETRmax/data$alpha[1]
      alpha=data$alpha[1]
      eopt=data$eopt[1]
      Fv_Fm_mean=mean(PE$Fv_Fm)
      PE_Model <- tibble(Espece="P", Jour=i , Food=k ,  Replicat=j, PAR=c(1:934))
      PE_Model$model = PE_Model$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE_Model$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE_Model$PAR+(1/data$alpha[1]))
      data <- tibble(Espece="P", Jour=i , Food=k ,  Replicat=j,
                     alpha = alpha,
                     ETRmax=ETRmax,
                     Ek=Ek,
                     eopt=eopt,
                     Fv_Fm_mean=Fv_Fm_mean)
      PE_model <- bind_rows(PE_model , PE_Model)
      PE_Parameters <- bind_rows(PE_Parameters, data)
      PE_curves <- bind_rows(PE_curves, PE)
    }
  }
}



for (i in c(6)) {
  for (k in c(0)) {
    for (j in c( "B")) {
      PE <- subset(Courbes_PE, Espece=="P" & Jour==i & Food==k & Replicat==j)
      data <- fitEP(PE$PAR, PE$ETR_II, normalize = F, lowerlim = c(4.7, 0, 1000), upperlim = c(10, 900, 1200))
      PE$model <- PE$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE$PAR+(1/data$alpha[1]))
      ETRmax=data$ps[1]
      Ek=ETRmax/data$alpha[1]
      alpha=data$alpha[1]
      eopt=data$eopt[1]
      Fv_Fm_mean=mean(PE$Fv_Fm)
      PE_Model <- tibble(Espece="P", Jour=i , Food=k ,  Replicat=j, PAR=c(1:934))
      PE_Model$model = PE_Model$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE_Model$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE_Model$PAR+(1/data$alpha[1]))
      data <- tibble(Espece="P", Jour=i , Food=k ,  Replicat=j,
                     alpha = alpha,
                     ETRmax=ETRmax,
                     Ek=Ek,
                     eopt=eopt,
                     Fv_Fm_mean=Fv_Fm_mean)
      PE_model <- bind_rows(PE_model , PE_Model)
      PE_Parameters <- bind_rows(PE_Parameters, data)
      PE_curves <- bind_rows(PE_curves, PE)
    }
  }
}


for (i in c(6)) {
  for (k in c(0)) {
    for (j in c("C")) {
      PE <- subset(Courbes_PE, Espece=="P" & Jour==i & Food==k & Replicat==j)
      data <- fitEP(PE$PAR, PE$ETR_II, normalize = F, lowerlim = c(4.6, 0, 1000), upperlim = c(10, 900, 1300))
      PE$model <- PE$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE$PAR+(1/data$alpha[1]))
      ETRmax=data$ps[1]
      Ek=ETRmax/data$alpha[1]
      alpha=data$alpha[1]
      eopt=data$eopt[1]
      Fv_Fm_mean=mean(PE$Fv_Fm)
      PE_Model <- tibble(Espece="P", Jour=i , Food=k ,  Replicat=j, PAR=c(1:934))
      PE_Model$model = PE_Model$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE_Model$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE_Model$PAR+(1/data$alpha[1]))
      data <- tibble(Espece="P", Jour=i , Food=k ,  Replicat=j,
                     alpha = alpha,
                     ETRmax=ETRmax,
                     Ek=Ek,
                     eopt=eopt,
                     Fv_Fm_mean=Fv_Fm_mean)
      PE_model <- bind_rows(PE_model , PE_Model)
      PE_Parameters <- bind_rows(PE_Parameters, data)
      PE_curves <- bind_rows(PE_curves, PE)
    }
  }
}



for (i in c(6)) {
  for (k in c(0)) {
    for (j in c("D")) {
      PE <- subset(Courbes_PE, Espece=="P" & Jour==i & Food==k & Replicat==j)
      data <- fitEP(PE$PAR, PE$ETR_II, normalize = F, lowerlim = c(4.9, 000, 1000), upperlim = c(10, 950, 1500))
      PE$model <- PE$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE$PAR+(1/data$alpha[1]))
      ETRmax=data$ps[1]
      Ek=ETRmax/data$alpha[1]
      alpha=data$alpha[1]
      eopt=data$eopt[1]
      Fv_Fm_mean=mean(PE$Fv_Fm)
      PE_Model <- tibble(Espece="P", Jour=i , Food=k ,  Replicat=j, PAR=c(1:934))
      PE_Model$model = PE_Model$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE_Model$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE_Model$PAR+(1/data$alpha[1]))
      data <- tibble(Espece="P", Jour=i , Food=k ,  Replicat=j,
                     alpha = alpha,
                     ETRmax=ETRmax,
                     Ek=Ek,
                     eopt=eopt,
                     Fv_Fm_mean=Fv_Fm_mean)
      PE_model <- bind_rows(PE_model , PE_Model)
      PE_Parameters <- bind_rows(PE_Parameters, data)
      PE_curves <- bind_rows(PE_curves, PE)
    }
  }
}

ggplot()+
  geom_point(data=subset(PE_curves, Jour==6 & Replicat =="D"), aes(x=PAR, y=ETR_II, group=interaction(Replicat, Jour), colour=Replicat))+
  geom_line(data=subset(PE_model, Jour==6 & Replicat =="D"), aes(x=PAR, y=model, group=interaction(Replicat, Jour), colour=Replicat))



PE_Parameters

ggplot()+
  geom_point(data=subset(PE_curves, Jour==6), aes(x=PAR, y=ETR_II, group=interaction(Replicat, Jour), colour=Replicat))+
  geom_line(data=subset(PE_model, Jour==6), aes(x=PAR, y=model, group=interaction(Replicat, Jour), colour=Replicat))

ggplot()+
  geom_point(data=subset(PE_curves, Jour==6), aes(x=PAR, y=ETR_II, group=interaction(Replicat, Jour), colour=Replicat))+
  geom_line(data=subset(PE_model, Jour==6), aes(x=PAR, y=model, group=interaction(Replicat, Jour), colour=Replicat))+
  xlim(0,100)+
  ylim(0,300)

for (i in c(0)) {
  for (k in c(0)) {
    for (j in c("A", "B","C","D")) {
      PE <- subset(Courbes_PE, Espece=="P" & Jour==i & Food==k & Replicat==j)
      data <- fitEP(PE$PAR, PE$ETR_II, normalize = F, lowerlim = c(3, 0, 0), upperlim = c(5, 1600, 1100), 
                    fitmethod=c("Nelder-Mead"))
      PE$model <- PE$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE$PAR+(1/data$alpha[1]))
      ETRmax=data$ps[1]
      Ek=ETRmax/data$alpha[1]
      alpha=data$alpha[1]
      eopt=data$eopt[1]
      Fv_Fm_mean=mean(PE$Fv_Fm)
      PE_Model <- tibble(Espece="P", Jour=i , Food=k ,  Replicat=j, PAR=c(1:934))
      PE_Model$model = PE_Model$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE_Model$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE_Model$PAR+(1/data$alpha[1]))
      data <- tibble(Espece="P", Jour=i , Food=k ,  Replicat=j,
                     alpha = alpha,
                     ETRmax=ETRmax,
                     Ek=Ek,
                     eopt=eopt,
                     Fv_Fm_mean=Fv_Fm_mean)
      PE_model <- bind_rows(PE_model , PE_Model)
      PE_Parameters <- bind_rows(PE_Parameters, data)
      PE_curves <- bind_rows(PE_curves, PE)
    }
  }
}





for (i in c(0)) {
  for (k in c(0,40,200)) {
    for (j in c("A", "B","C")) {
      PE <- subset(Courbes_PE, Espece=="M" & Jour==i & Food==k & Replicat==j)
      data <- fitEP(PE$PAR, PE$ETR_II, normalize = FALSE)
      PE$model <- PE$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE$PAR+(1/data$alpha[1]))
      ETRmax=data$ps[1]
      Ek=ETRmax/data$alpha[1]
      alpha=data$alpha[1]
      eopt=data$eopt[1]
      Fv_Fm_mean=mean(PE$Fv_Fm)
      PE$Food = "0"
      PE_Model <- tibble(Espece="M", Jour=i , Food=0 ,  Replicat=j, PAR=c(1:934))
      PE_Model$model = PE_Model$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE_Model$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE_Model$PAR+(1/data$alpha[1]))
      data <- tibble(Espece="M", Jour=i , Food=0 ,  Replicat=j,
                     alpha = alpha,
                     ETRmax=ETRmax,
                     eopt=eopt,
                     Ek=Ek,
                     Fv_Fm_mean=Fv_Fm_mean)
      PE_model <- bind_rows(PE_model , PE_Model)
      PE_Parameters <- bind_rows(PE_Parameters, data)
      PE_curves <- bind_rows(PE_curves, PE)
    }
  }
}

PE_Parameters




for (i in c(7)) {
  for (k in c(0,40,200)) {
    for (j in c("A", "B","C")) {
      PE <- subset(Courbes_PE, Espece=="M" & Jour==i & Food==k & Replicat==j)
      data <- fitEP(PE$PAR, PE$ETR_II, normalize = FALSE)
      PE$model <- PE$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE$PAR+(1/data$alpha[1]))
      ETRmax=data$ps[1]
      Ek=ETRmax/data$alpha[1]
      alpha=data$alpha[1]
      eopt=data$eopt[1]
      Fv_Fm_mean=mean(PE$Fv_Fm)
      PE_Model <- tibble(Espece="M", Jour=i , Food=k ,  Replicat=j, PAR=c(1:934))
      PE_Model$model = PE_Model$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE_Model$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE_Model$PAR+(1/data$alpha[1]))
      data <- tibble(Espece="M", Jour=i , Food=k ,  Replicat=j,
                     alpha = alpha,
                     ETRmax=ETRmax,
                     eopt=eopt,
                     Ek=Ek,
                     Fv_Fm_mean=Fv_Fm_mean)
      PE_model <- bind_rows(PE_model , PE_Model)
      PE_Parameters <- bind_rows(PE_Parameters, data)
      PE_curves <- bind_rows(PE_curves, PE)
    }
  }
}

ggplot()+
  geom_point(data=subset(PE_curves, Espece=="M" & Jour==0), aes(x=PAR, y=ETR_II, group=interaction(Replicat, Jour, Food), colour=Replicat))+
  geom_line(data=subset(PE_model,  Espece=="M" & Jour==0), aes(x=PAR, y=model, group=interaction(Replicat, Jour, Food), colour=Replicat))

ggplot()+
  geom_point(data=subset(PE_curves, Espece=="M" & Jour==7 & Food==200), aes(x=PAR, y=ETR_II, group=interaction(Replicat, Jour), colour=Replicat))+
  geom_line(data=subset(PE_model, Espece=="M" & Jour==7 & Food==200), aes(x=PAR, y=model, group=interaction(Replicat, Jour), colour=Replicat))+
  xlim(0,100)+
  ylim(0,300)


for (i in c(0)) {
  for (k in c(0,40,200)) {
    for (j in c("A", "B","C")) {
      PE <- subset(Courbes_PE, Espece=="D" & Jour==i & Food==k & Replicat==j)
      data <- fitEP(PE$PAR, PE$ETR_II, normalize = FALSE)
      PE$model <- PE$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE$PAR+(1/data$alpha[1]))
      ETRmax=data$ps[1]
      Ek=ETRmax/data$alpha[1]
      alpha=data$alpha[1]
      eopt=data$eopt[1]
      Fv_Fm_mean=mean(PE$Fv_Fm)
      PE$Food = "0"
      PE_Model <- tibble(Espece="D", Jour=i , Food=0 ,  Replicat=j, PAR=c(1:1650))
      PE_Model$model = PE_Model$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE_Model$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE_Model$PAR+(1/data$alpha[1]))
      data <- tibble(Espece="D", Jour=i , Food=0 ,  Replicat=j,
                     alpha = alpha,
                     eopt=eopt,
                     ETRmax=ETRmax,
                     Ek=Ek,
                     Fv_Fm_mean=Fv_Fm_mean)
      PE_model <- bind_rows(PE_model , PE_Model)
      PE_Parameters <- bind_rows(PE_Parameters, data)
      PE_curves <- bind_rows(PE_curves, PE)
    }
  }
}


for (i in c(8)) {
  for (k in c(0,40,200)) {
    for (j in c("A", "B","C")) {
      PE <- subset(Courbes_PE, Espece=="D" & Jour==i & Food==k & Replicat==j)
      data <- fitEP(PE$PAR, PE$ETR_II, normalize = FALSE)
      PE$model <- PE$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE$PAR+(1/data$alpha[1]))
      ETRmax=data$ps[1]
      Ek=ETRmax/data$alpha[1]
      alpha=data$alpha[1]
      eopt=data$eopt[1]
      Fv_Fm_mean=mean(PE$Fv_Fm)
      PE_Model <- tibble(Espece="D", Jour=i , Food=k ,  Replicat=j, PAR=c(1:1650))
      PE_Model$model = PE_Model$PAR/((1/(data$alpha[1]*data$eopt[1]^2))*PE_Model$PAR^2+(1/data$ps[1]-2/(data$alpha[1]*data$eopt[1]))*PE_Model$PAR+(1/data$alpha[1]))
      data <- tibble(Espece="D", Jour=i , Food=k ,  Replicat=j,
                     alpha = alpha,
                     eopt=eopt,
                     ETRmax=ETRmax,
                     Ek=Ek,
                     Fv_Fm_mean=Fv_Fm_mean)
      PE_model <- bind_rows(PE_model , PE_Model)
      PE_Parameters <- bind_rows(PE_Parameters, data)
      PE_curves <- bind_rows(PE_curves, PE)
    }
  }
}

ggplot()+
  geom_point(data=subset(PE_curves, Espece=="D" & Jour==8), aes(x=PAR, y=ETR_II, group=interaction( Jour, Food), colour=Food))

ggplot()+
  geom_point(data=subset(PE_curves, Espece=="D" & Jour==8 & Food==200), aes(x=PAR, y=ETR_II, group=interaction(Replicat, Jour), colour=Replicat))+
  geom_line(data=subset(PE_model, Espece=="D" & Jour==8 & Food==200), aes(x=PAR, y=model, group=interaction(Replicat, Jour), colour=Replicat))+
  xlim(0,100)+
  ylim(0,300)

## Model mean -----------------

PE_curves_mean <- summarise(group_by(PE_curves, Espece, Food, Jour, PAR),
                       ETR_II_mean=mean(ETR_II),
                       ETR_II_sd=sd(ETR_II))

PE_curves_mean$state <-PE_curves_mean$Jour
PE_curves_mean$state <- str_replace(PE_curves_mean$state, "0", "T0")
PE_curves_mean$state <- str_replace(PE_curves_mean$state, "6", "Tf")
PE_curves_mean$state <- str_replace(PE_curves_mean$state, "7", "Tf")
PE_curves_mean$state <- str_replace(PE_curves_mean$state, "8", "Tf")
PE_curves_mean$Food <- as.numeric(PE_curves_mean$Food)
PE_curves_mean$Food <- as.factor(PE_curves_mean$Food)

PE_model <- summarise(group_by(PE_model, Espece, Food, Jour, PAR),
                      model_mean=mean(model),
                      model_sd=sd(model))

PE_model$state <-PE_model$Jour
PE_model$state <- str_replace(PE_model$state, "0", "T0")
PE_model$state <- str_replace(PE_model$state, "6", "Tf")
PE_model$state <- str_replace(PE_model$state, "7", "Tf")
PE_model$state <- str_replace(PE_model$state, "8", "Tf")

PE_model$Food <- as.factor(PE_model$Food)

PE_Parameters$Food2 <- PE_Parameters$Food
PE_Parameters$Food2[PE_Parameters$Jour==0] <- 1

Parameters_mean1 <- summarise(group_by(PE_Parameters, Espece, Food, Jour),
                             mean_alpha = mean(alpha),
                             mean_eopt = mean(eopt),
                             mean_ETRmax=mean(ETRmax),
                             mean_Ek=mean(Ek),
                             mean_Fv_Fm=mean(Fv_Fm_mean),
                             sd_alpha = sd(alpha),
                             sd_eopt = sd(eopt),
                             sd_ETRmax=sd(ETRmax),
                             sd_Ek=sd(Ek),
                             sd_Fv_Fm=sd(Fv_Fm_mean))

Parameters_mean <- summarise(group_by(PE_Parameters, Espece, Food2, Jour),
                              mean_alpha = mean(alpha),
                              mean_eopt = mean(eopt),
                              mean_ETRmax=mean(ETRmax),
                              sd_eopt = sd(eopt),
                              mean_Ek=mean(Ek),
                              mean_Fv_Fm=mean(Fv_Fm_mean),
                              sd_alpha = sd(alpha),
                              sd_ETRmax=sd(ETRmax),
                              sd_Ek=sd(Ek),
                              sd_Fv_Fm=sd(Fv_Fm_mean))

Parameters_mean$state <-Parameters_mean$Jour
Parameters_mean$state <- str_replace(Parameters_mean$state, "0", "T0")
Parameters_mean$state <- str_replace(Parameters_mean$state, "6", "Tf")
Parameters_mean$state <- str_replace(Parameters_mean$state, "7", "Tf")
Parameters_mean$state <- str_replace(Parameters_mean$state, "8", "Tf")

#PE_model$state <- factor(PE_model$state, levels=c("Initial", "Final"))
#PE_curves_mean$state <- factor(PE_curves_mean$state, levels=c("Initial", "Final"))
Parameters_mean$state <- factor(Parameters_mean$state, levels=c("T0", "Tf"))
# PE_model$Espece <- factor(PE_model$Espece, levels=c("P", "M", "D"))
# PE_curves_mean$Espece <- factor(PE_curves_mean$Espece, levels=c("P", "M", "D"))

PE_model$state <- factor(PE_model$state, levels=c("T0", "Tf")
                         )
PE_curves_mean$state <- factor(PE_curves_mean$state, levels=c("T0", "Tf"))

PE_model$Espece[PE_model$Espece == "P"] <- "T. amphioxeia"
PE_model$Espece[PE_model$Espece == "M"] <- "M. rubrum"
PE_model$Espece[PE_model$Espece == "D"] <- "D. acuminata"

PE_model$Espece <- factor(PE_model$Espece, levels=c("T. amphioxeia",
                                                                "M. rubrum", 
                                                                "D. acuminata"))

PE_curves_mean$Espece[PE_curves_mean$Espece == "P"] <- "T. amphioxeia"
PE_curves_mean$Espece[PE_curves_mean$Espece == "M"] <- "M. rubrum"
PE_curves_mean$Espece[PE_curves_mean$Espece == "D"] <- "D. acuminata"

PE_curves_mean$Espece <- factor(PE_curves_mean$Espece, levels=c("T. amphioxeia",
                                                                  "M. rubrum", 
                                                                  "D. acuminata"))

Parameters_mean$Espece[Parameters_mean$Espece == "P"] <- "T. amphioxeia"
Parameters_mean$Espece[Parameters_mean$Espece == "M"] <- "M. rubrum"
Parameters_mean$Espece[Parameters_mean$Espece == "D"] <- "D. acuminata"

Parameters_mean$Espece <- factor(Parameters_mean$Espece, levels=c("T. amphioxeia",
                                                                  "M. rubrum", 
                                                                  "D. acuminata"))
Parameters_mean$Food2 <- factor(Parameters_mean$Food2, 
                                levels=c(1, 0, 200,40))

PE_plot <- merge(PE_curves_mean, PE_model, by=c("Espece", "Jour", "Food", "PAR", "state"), all=T)

PE_plot_P <- subset(PE_plot, Espece=="T. amphioxeia")

PE_plot_P$Food <- "-1"

PE_plot <- bind_rows(PE_plot_P, subset(PE_plot, Espece!="P. prolonga"))


Courbes_PE_Plagio_I <- ggplot(data = subset(PE_plot, Espece=="T. amphioxeia" &  state=="T0"), aes(x=PAR, group=Food, shape=Food))+
  geom_errorbar(aes(  y=ETR_II_mean, ymin=ETR_II_mean-ETR_II_sd, ymax=ETR_II_mean+ETR_II_sd), colour="gray50")+
  scale_x_continuous(limits = c(0,1020),
                     breaks=seq(0,1000, 200),
                     minor_breaks = seq(0,1000,50))+
  scale_y_continuous(limits = c(0,1100),
                     breaks=seq(0,1100, 500),
                     minor_breaks = seq(0,1100,100))+
  guides(x="prism_minor",y="prism_minor")+
  geom_line( aes( y=model_mean), size=Taille_Ligne, colour="gray50")+
  geom_point(aes( y=ETR_II_mean), size=Taille_Point-2, colour="gray50")+
  facet_grid(Espece~state, space = "free_x")+
  scale_shape_manual(values=c(19,15,18,17), 
                     labels=c(expression(bold(paste("Unfed"))),
                              expression(bold(paste("Unfed"))), 
                              expression(bold(paste('HL prey'))),
                              expression(bold(paste('LL prey')))))+
  #xlim(0,1050)+ 
  theme_linedraw()+
  theme(axis.ticks.length = unit(Taille_Ticks_shift, 'cm'),         
        axis.ticks = element_line(size=Epaisseur_Ticks_shift, colour=Couleur_Bordure),
        panel.border  = element_rect( size = Taille_Rectangle, fill=NA, colour = Couleur_Bordure),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid = element_blank(),
        axis.title = element_text(colour=Couleur_Bordure,face="bold",size=Taille_Titre), 
        axis.text = element_text(colour=Couleur_Bordure,size=Taille_Axe),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.key.width = unit(3.5, "line"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent'),
        legend.text = element_text(size=Taille_Axe-3,face="bold", colour=Couleur_Bordure,
                                   margin = unit(c(1, 0, 0, 0), "mm")),
        strip.background = element_blank(),
        strip.text = element_text(colour = Couleur_Bordure, face="bold", size=Taille_Axe),
        strip.text.y = element_blank(), #element_text(face = "bold.italic"),
        strip.placement = "top")+
  annotate(geom = 'segment', y = Inf, yend = Inf, color = Couleur_Bordure, x = -Inf, xend = Inf, size = Taille_Rectangle+.5) +
  annotate(geom = 'segment', y = Inf, yend = -Inf, color = Couleur_Bordure, x = -Inf, xend = -Inf, size = Taille_Rectangle+.5) + 
  #Labels
  labs(#x=expression(paste("irradiance (?mol de photons ",  m^{-2}, s^{-1}, ')')),
    x=expression(bold(paste("Irradiance (µmol photons ",  "m"^{"-2"}, "s"^{"-1"}, ')'))),
    y=expression(bold(paste("ETRII (e- ",  "PSII"^{"-1"}, ')'))))
Courbes_PE_Plagio_I

Courbes_PE_Plagio_F <- ggplot(data = subset(PE_plot, Espece=="T. amphioxeia" &  state=="Tf"), aes(x=PAR, group=Food, shape=Food))+
  geom_errorbar(aes(  y=ETR_II_mean, ymin=ETR_II_mean-ETR_II_sd, ymax=ETR_II_mean+ETR_II_sd), colour="gray50")+
  scale_x_continuous(limits = c(0,1020),
                     breaks=seq(0,1000, 200),
                     minor_breaks = seq(0,1000,50))+
  scale_y_continuous(limits = c(0,1100),
                     breaks=seq(0,1100, 500),
                     minor_breaks = seq(0,1100,100))+
  guides(x="prism_minor",y="prism_minor")+
  geom_line( aes( y=model_mean), size=Taille_Ligne, colour="gray50")+
  geom_point(aes( y=ETR_II_mean), size=Taille_Point-2, colour="gray50")+
  facet_grid(Espece~state, space = "free_x")+
  scale_shape_manual(values=c(19,15,17,18), 
                     labels=c(expression(bold(paste("Unfed"))),
                              expression(bold(paste("Unfed"))), 
                              expression(bold(paste('HL prey'))),
                              expression(bold(paste('LL prey')))))+
  #xlim(0,1050)+ 
  theme_linedraw()+
  theme(axis.ticks.length = unit(Taille_Ticks_shift, 'cm'),         
        axis.ticks = element_line(size=Epaisseur_Ticks_shift, colour=Couleur_Bordure),
        panel.border  = element_rect( size = Taille_Rectangle, fill=NA, colour = Couleur_Bordure),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid = element_blank(),
        axis.title = element_text(colour=Couleur_Bordure,face="bold",size=Taille_Titre), 
        axis.text = element_text(colour=Couleur_Bordure,size=Taille_Axe),
        axis.text.y=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none", #c(0.38,0.42),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.key.width = unit(3.5, "line"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent'),
        legend.text = element_text(size=Taille_Axe,face="bold", colour=Couleur_Bordure,
                                   margin = unit(c(1, 0, 0, 0), "mm")),
        strip.background = element_blank(),
        strip.text = element_text(colour = Couleur_Bordure, face="bold", size=Taille_Titre),
        strip.text.y = element_text(face = "bold.italic"),
        strip.placement = "top")+
  annotate(geom = 'segment', y = Inf, yend = Inf, color = Couleur_Bordure, x = -Inf, xend = Inf, size = Taille_Rectangle+.5) +
  annotate(geom = 'segment', y = Inf, yend = -Inf, color = Couleur_Bordure, x = -Inf, xend = -Inf, size = Taille_Rectangle+.5) + 
  #Labels
  labs(#x=expression(paste("irradiance (?mol de photons ",  m^{-2}, s^{-1}, ')')),
    x=expression(bold(paste("Irradiance (µmol photons ",  "m"^{"-2"}, "s"^{"-1"}, ')'))),
    y=expression(bold(paste("ETRII (e- ",  "PSII"^{"-1"}, ')'))))
Courbes_PE_Plagio_F

Courbes_PE_Meso_I <- ggplot(data = subset(PE_plot, Espece=="M. rubrum" &  state=="T0"), aes(x=PAR, group=Food, shape=Food))+
  geom_errorbar(aes(  y=ETR_II_mean, ymin=ETR_II_mean-ETR_II_sd, ymax=ETR_II_mean+ETR_II_sd), colour="gray50")+
  scale_x_continuous(limits = c(0,1020),
                     breaks=seq(0,1000, 200),
                     minor_breaks = seq(0,1000,50))+
  scale_y_continuous(limits = c(0,1100),
                     breaks=seq(0,1100, 500),
                     minor_breaks = seq(0,1100,100))+
  guides(x="prism_minor",y="prism_minor")+
  geom_line( aes( y=model_mean), size=Taille_Ligne, colour="gray50")+
  geom_point(aes( y=ETR_II_mean), size=Taille_Point-2, colour="gray50")+
  facet_grid(Espece~state, space = "free_x")+
  scale_shape_manual(values=c(19))+
  #xlim(0,1050)+ 
  theme_linedraw()+
  theme(axis.ticks.length = unit(Taille_Ticks_shift, 'cm'),         
        axis.ticks = element_line(size=Epaisseur_Ticks_shift, colour=Couleur_Bordure),
        panel.border  = element_rect( size = Taille_Rectangle, fill=NA, colour = Couleur_Bordure),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid = element_blank(),
        axis.title = element_text(colour=Couleur_Bordure,face="bold",size=Taille_Titre), 
        axis.text = element_text(colour=Couleur_Bordure,size=Taille_Axe),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.key.width = unit(3.5, "line"),
        legend.spacing.x = unit(-0.2, 'cm'),
        legend.key = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent'),
        legend.text = element_text(size=Taille_Axe,face="bold", colour=Couleur_Bordure,
                                   margin = unit(c(0, 0, 0, 0), "mm")),
        strip.background = element_blank(),
        strip.text = element_text(colour = Couleur_Bordure, face="bold", size=Taille_Axe),
        strip.text.y = element_blank(), #element_text(face = "bold.italic"),
        strip.text.x = element_blank(),
        strip.placement = "top")+
  annotate(geom = 'segment', y = Inf, yend = Inf, color = Couleur_Bordure, x = -Inf, xend = Inf, size = Taille_Rectangle+.5) +
  annotate(geom = 'segment', y = Inf, yend = -Inf, color = Couleur_Bordure, x = -Inf, xend = -Inf, size = Taille_Rectangle+.5) + 
  #Labels
  labs(#x=expression(paste("irradiance (?mol de photons ",  m^{-2}, s^{-1}, ')')),
    x=expression(bold(paste("Irradiance (µmol photons ",  "m"^{"-2"}, "s"^{"-1"}, ')'))),
    y=expression(bold(paste("ETRII (e- ",  "PSII"^{"-1"}, ')'))))
Courbes_PE_Meso_I

Courbes_PE_Meso_F <- ggplot(data = subset(PE_plot, Espece=="M. rubrum" &  state=="Tf"), aes(x=PAR, group=Food, colour=Food, shape=Food))+
  geom_errorbar(aes(  y=ETR_II_mean, ymin=ETR_II_mean-ETR_II_sd, ymax=ETR_II_mean+ETR_II_sd), colour="gray50")+
  scale_x_continuous(limits = c(0,1020),
                     breaks=seq(0,1000, 200),
                     minor_breaks = seq(0,1000,50))+
  scale_y_continuous(limits = c(0,1100),
                     breaks=seq(0,1100, 500),
                     minor_breaks = seq(0,1100,100))+
  scale_colour_manual(values = c(Gris_Clair, or, jaune_clair), 
                      labels=c(expression(bold(paste("Unfed"))), 
                               expression(bold(paste('HL prey'))),
                               expression(bold(paste('LL prey')))))+
  guides(x="prism_minor",y="prism_minor")+
  geom_line( aes( y=model_mean), size=Taille_Ligne)+
  geom_point(aes( y=ETR_II_mean), size=Taille_Point-2)+
  facet_grid(Espece~state, space = "free_x")+
  scale_shape_manual(values=c(15,17,18), 
                     labels=c(expression(bold(paste("Unfed"))), 
                              expression(bold(paste('HL prey'))),
                              expression(bold(paste('LL prey')))))+
  #xlim(0,1050)+ 
  theme_linedraw()+
  theme(axis.ticks.length = unit(Taille_Ticks_shift, 'cm'),         
        axis.ticks = element_line(size=Epaisseur_Ticks_shift, colour=Couleur_Bordure),
        panel.border  = element_rect( size = Taille_Rectangle, fill=NA, colour = Couleur_Bordure),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        axis.title = element_text(colour=Couleur_Bordure,face="bold",size=Taille_Titre), 
        axis.text = element_text(colour=Couleur_Bordure,size=Taille_Axe),
        axis.text.y=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.position = c(0.78,0.21),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.key.width = unit(3.5, "line"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent'),
        legend.text = element_text(size=Taille_Axe-1,face="bold", colour=Couleur_Bordure,
                                   margin = unit(c(1, 0, 0, 0), "mm")),
        strip.background = element_blank(),
        strip.text = element_text(colour = Couleur_Bordure, face="bold", size=Taille_Titre),
        strip.text.y = element_text(face = "bold.italic"),
        strip.text.x = element_blank(),
        strip.placement = "top")+
  annotate(geom = 'segment', y = Inf, yend = Inf, color = Couleur_Bordure, x = -Inf, xend = Inf, size = Taille_Rectangle+.5) +
  annotate(geom = 'segment', y = Inf, yend = -Inf, color = Couleur_Bordure, x = -Inf, xend = -Inf, size = Taille_Rectangle+.5) + 
  #Labels
  labs(#x=expression(paste("irradiance (?mol de photons ",  m^{-2}, s^{-1}, ')')),
    x=expression(bold(paste("Irradiance (µmol photons ",  "m"^{"-2"}, "s"^{"-1"}, ')'))),
    y=expression(bold(paste("ETRII (e- ",  "PSII"^{"-1"}, ')'))))
Courbes_PE_Meso_F

Courbes_PE_Dino_I <- ggplot(data = subset(PE_plot, Espece=="D. acuminata" &  state=="T0"), aes(x=PAR, group=Food, shape=Food))+
  geom_errorbar(aes(  y=ETR_II_mean, ymin=ETR_II_mean-ETR_II_sd, ymax=ETR_II_mean+ETR_II_sd), colour="gray50")+
  scale_x_continuous(limits = c(0,1020),
                     breaks=seq(0,1000, 200),
                     minor_breaks = seq(0,1000,50))+
  scale_y_continuous(limits = c(0,1100),
                     breaks=seq(0,1100, 500),
                     minor_breaks = seq(0,1100,100))+
  guides(x="prism_minor",y="prism_minor")+
  geom_line( aes( y=model_mean), size=Taille_Ligne, colour="gray50")+
  geom_point(aes( y=ETR_II_mean), size=Taille_Point-2, colour="gray50")+
  facet_grid(Espece~state, space = "free_x")+
  scale_shape_manual(values=c(19),
                     labels=c(expression(bold(paste("Unfed"))),
                              expression(bold(paste("Unfed"))),
                              expression(bold(paste('HL fed'))),
                              expression(bold(paste('LL fed')))))+
  theme_linedraw()+
  theme(axis.ticks.length = unit(Taille_Ticks_shift, 'cm'),         
        axis.ticks = element_line(size=Epaisseur_Ticks_shift, colour=Couleur_Bordure),
        panel.border  = element_rect( size = Taille_Rectangle, fill=NA, colour = Couleur_Bordure),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid = element_blank(),
        axis.title = element_text(colour=Couleur_Bordure,face="bold",size=Taille_Titre), 
        axis.text = element_text(colour=Couleur_Bordure,size=Taille_Axe),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.key.width = unit(3.5, "line"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent'),
        legend.text = element_text(size=Taille_Axe-3,face="bold", colour=Couleur_Bordure,
                                   margin = unit(c(1, 0, 0, 0), "mm")),
        strip.background = element_blank(),
        strip.text = element_text(colour = Couleur_Bordure, face="bold", size=Taille_Axe),
        strip.text.y = element_blank(), #element_text(face = "bold.italic"),
        strip.text.x = element_blank(),
        strip.placement = "top")+
  annotate(geom = 'segment', y = Inf, yend = Inf, color = Couleur_Bordure, x = -Inf, xend = Inf, size = Taille_Rectangle+.5) +
  annotate(geom = 'segment', y = Inf, yend = -Inf, color = Couleur_Bordure, x = -Inf, xend = -Inf, size = Taille_Rectangle+.5) + 
  #Labels
  labs(#x=expression(paste("irradiance (?mol de photons ",  m^{-2}, s^{-1}, ')')),
    x=expression(bold(paste("Irradiance (µmol photons ",  "m"^{"-2"}, "s"^{"-1"}, ')'))),
    y=expression(bold(paste("ETRII (e- ",  "PSII"^{"-1"}, ')'))))
Courbes_PE_Dino_I

Courbes_PE_Dino_F <- ggplot(data = subset(PE_plot, Espece=="D. acuminata" &  state=="Tf"), aes(x=PAR, group=Food, colour=Food, shape=Food))+
  geom_errorbar(aes(  y=ETR_II_mean, ymin=ETR_II_mean-ETR_II_sd, ymax=ETR_II_mean+ETR_II_sd), colour="gray50")+
  scale_x_continuous(limits = c(0,1020),
                     breaks=seq(0,1000, 200),
                     minor_breaks = seq(0,1000,50))+
  scale_y_continuous(limits = c(0,1100),
                     breaks=seq(0,1100, 500),
                     minor_breaks = seq(0,1100,100))+
  scale_colour_manual(values = c(Gris_Clair,  jaune_clair,or))+
  guides(x="prism_minor",y="prism_minor")+
  geom_line( aes( y=model_mean), size=Taille_Ligne)+
  geom_point(aes( y=ETR_II_mean), size=Taille_Point-2)+
  facet_grid(Espece~state, space = "free_x")+
  scale_shape_manual(values=c(15,18,17), 
                     labels=c(expression(bold(paste("Unfed"))),
                              expression(bold(paste("Unfed"))), 
                              expression(bold(paste('HL fed'))),
                              expression(bold(paste('LL fed')))))+
  #xlim(0,1050)+ 
  theme_linedraw()+
  theme(axis.ticks.length = unit(Taille_Ticks_shift, 'cm'),         
        axis.ticks = element_line(size=Epaisseur_Ticks_shift, colour=Couleur_Bordure),
        panel.border  = element_rect( size = Taille_Rectangle, fill=NA, colour = Couleur_Bordure),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid = element_blank(),
        axis.title = element_text(colour=Couleur_Bordure,face="bold",size=Taille_Titre), 
        axis.text = element_text(colour=Couleur_Bordure,size=Taille_Axe),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none", #c(0.38,0.42),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.key.width = unit(3.5, "line"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent'),
        legend.text = element_text(size=Taille_Axe-3,face="bold", colour=Couleur_Bordure,
                                   margin = unit(c(1, 0, 0, 0), "mm")),
        strip.background = element_blank(),
        strip.text = element_text(colour = Couleur_Bordure, face="bold", size=Taille_Titre),
        strip.text.y = element_text(face = "bold.italic"),
        strip.text.x = element_blank(),
        strip.placement = "top")+
  annotate(geom = 'segment', y = Inf, yend = Inf, color = Couleur_Bordure, x = -Inf, xend = Inf, size = Taille_Rectangle+.5) +
  annotate(geom = 'segment', y = Inf, yend = -Inf, color = Couleur_Bordure, x = -Inf, xend = -Inf, size = Taille_Rectangle+.5) + 
  #Labels
  labs(#x=expression(paste("irradiance (?mol de photons ",  m^{-2}, s^{-1}, ')')),
    x=expression(bold(paste("Irradiance (µmol photons ",  "m"^{"-2"}, "s"^{"-1"}, ')'))),
    y=expression(bold(paste("ETRII (e- ",  "PSII"^{"-1"}, ')'))))
Courbes_PE_Dino_F

## Plot PE curves -----

Combine_Courbes_PE <- plot_grid(Courbes_PE_Plagio_I, NULL, Courbes_PE_Plagio_F,
                                NULL,NULL,NULL,
                                Courbes_PE_Meso_I, NULL, Courbes_PE_Meso_F,
                                NULL,NULL,NULL,
                                Courbes_PE_Dino_I, NULL, Courbes_PE_Dino_F,
                                ncol = 3,
                                nrow = 5,
                                rel_widths = c(3,-0.4,3),
                                rel_heights = c(3,-0.5,3,-0.5,3),
                                align = "hv",
                                axis = "tblr",
                                labels = c("","","",
                                           "","","",
                                           "","","",
                                           "","","",
                                           "","",""
                                           ),
                                label_colour = Couleur_Bordure,
                                label_size = Taille_Titre+4,
                                label_x=0.19,
                                label_y = 0.90
                                )

ggsave(Combine_Courbes_PE, filename="Figures/Supp_3.tiff", dpi=300, height = 11, width = 12)


  