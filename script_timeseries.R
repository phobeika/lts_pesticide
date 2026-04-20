# Projet Time-Series
# 2A Paul Hobeika et Sara Chikhi

# Packages
library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(ggthemes)
library(tseries)
library(forecast)
library(astsa)
library(portes)
library(zoo)

# Partie 1 : Données ----
# Import des données
ind <- read_csv2("donnees/valeurs_mensuelles.csv", col_names = TRUE) %>% 
  slice(4:437) %>% # on supprime les premières lignes qui ne correspondent pas à des données
  select(-Codes) %>% # on supprime aussi la colonne "Codes" qui est toujours égale à 'A'
  rename(Date = 'Libellé', 
         Valeur = "Indice CVS-CJO de la production industrielle (base 100 en 2021) - Fabrication de pesticides et d'autres produits agrochimiques (NAF rév. 2, niveau groupe, poste 20.2)") %>% 
  mutate(Date = ym(Date), # on met la variable de période au format 'Date'
         Valeur = as.numeric(Valeur)) # et l'indice en nombre plutôt qu'en caractères

# Graphe de la série
graphe_brut <- ind %>% 
  arrange(Date) %>% 
  ggplot(aes(x = Date, y = Valeur)) +
  geom_line() +
  theme_minimal()
ggsave(filename = "donneesbrutes.png", device = "png")

# On fait une transformation de Box-Cox pour diminuer la variation de la
# variance en fonction du temps
lambda <- BoxCox.lambda(ind$Valeur)
ind_bc <- BoxCox(ind$Valeur, lambda)
plot(ind$Date, ind_bc, type = "l")

# Puis on prend la "first difference" pour supprimer les variations "lentes"
ind_diff1 <- diff(ind$Valeur)
plot(ind$Date[-1], ind_diff1, type = "l") # Ça semble relativement stationnaire

# on vérifie la stationnarité à l'aide de la fonction d'autocorrélation
acf(ind_diff1, main = "Autocorrélation après transformation de Box-Cox")

# on observe que l'autocorrélation chute rapidement à 0 
# c'est signe de stationnarité

# nous effectuons enfin le test de racine unitaire de Phillips-Perron
# c'est un test non paramétrique (pas besoin de spécifier un niveau d'autocorrélation)
pp.test(ind_diff1)
# l'hypothèse nulle est rejetée : il n'y a pas de racine unité
# la série est stationnaire

# nous pouvons appliquer un modèle ARMA(p,q)

# Partie II : Modèles ARMA ----

# le plot d'ACF suggère qu'il faudra utiliser les deux dernières valeurs pour 
# pour prédire la suivante ; soit p=2, éventuellement p=1

# on s'intéresse à l'autocorrélation partielle pour définir q
pacf(ind_diff1, main = "Autocorrélation partielle après transformation de Box-Cox")

# le graph suggère un q entre 2 et 5
# on s'intéresse aux statistiques des différentes spécifications

for (i in c(1,2)) {
  for (j in c(2:5)) {
   (sarima(ind_diff1, i, 1,j, details=FALSE))$ttable
  }
}

# on exclut ARMA(2,3), ARMA(2,4) et ARMA(2,5) qui contiennent des valeurs manquantes
# ce qui signifie qu'il y a redondance
# ARMA (1,3), ARMA(1,5) et ARMA(2,2) comportent des coeff non significatifs
# un ARMA(1,2) ou ARMA (1,4) pourrait convenir 
# (modèles pour lesquels tous les coefficients sont significatifs)

arma12 <- sarima(ind_diff1,1,1,2)
arma14 <- sarima(ind_diff1,1,1,4)

# les résidus présentent des autocorrélations très similaires

# test de Ljung-Box (ou test du porte-manteau)
LjungBox(arma12$fit)
LjungBox(arma14$fit)

# Le test de Ljung Box rejette le modèle ARMA(1,4) 
# lorsqu'il prend en compte 5 autocorrélations et 25 autocorrélations
# tandis que ARMA(1,2) n'est jamais rejeté par le test

# Comparons les AIC et BIC de chacun des modèles
AIC(arma12$fit, arma14$fit)
BIC(arma12$fit, arma14$fit)

# Les tests AIC et BIC sélectionnent tous deux le modèle ARMA(1,2)

# Partie III : Prévision ----
# prévision avec région de confiance à 95 %
ind_diff1_zoo <- zoo(ind_diff1, ind$Date)
ind_diff1_date <- ts(coredata(ind_diff1_zoo), frequency = 12) 

model <- Arima(ind_diff1_date, order = c(1, 0, 2))
prediction <- forecast(model, h = 2, level = 95)

autoplot(prediction) +
  # scale_x_continuous(labels = function(x) {
  #   format(seq(tail(ind$Date,-50)[1], by = "6 months", length.out=length(x)), "%Y-%m")
  # })+
  ggtitle("Prévision de la série à t+2, région de confiance à 95%") +
  xlab("Temps") + ylab("Xt") +
  theme_minimal() +
  coord_cartesian(xlim = c(tail(time(ind_diff1_date),50)[1], tail(time(prediction$mean), 2)[2]))
#REPRENDRE


