# Projet Time-Series
# 2A Paul Hobeika et Sara Chikhi

# Packages
library(dplyr)
library(readr)
library(lubridate)
library(fUnitRoots)
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
ind_zoo <- zoo(ind_bc, order.by = ind$Date)
ind_diff <- diff(ind_zoo, lag=1)
# plot(ind$Date[-1], ind_diff1, type = "l") # Ça semble relativement stationnaire

plot(cbind(ind_zoo, ind_diff))

# on prépare la série sans NA
ind_diff1 <- coredata(ind_diff)
ind_diff1 <- x[!is.na(ind_diff1)]

# Test de racine unitaire ----
# on fait une régression linéaire pour savoir comment spécifier notre test ADF
summary(lm(ind_diff1 ~ ind$Date[-1]))
# la constante est faible et graphiquement on n'observe pas de tendance
adf <- adfTest(ind_diff1, lags = 0, type = "nc")

# ------------------------------------------------------------------------------
# Test de Ljung-Box
# ------------------------------------------------------------------------------
# Cette fonction calcule les p-values du test de Ljung-Box
# pour les retards l = 1, ..., k.
#
# H_0 : rho_1 = rho_2 = ... = rho_l = 0,
# où l désigne le retard testé.
# Autrement dit, il n'y a pas d'autocorrélation résiduelle
# jusqu'au retard l considéré.
#
# Arguments :
# - series : série sur laquelle le test est appliqué (des résidus, en général)
# - k      : nombre maximal de retards testés (e.g. 24 pour les data mensuelles)
# - fitdf  : nombre de degrés de liberté à retirer
#            (le nombre de paramètres estimés par le modèle qui a donné ces résidus)
#
# La fonction renvoie une matrice à deux colonnes :
# - lag  : le retard considéré
# - pval : la p-value associée au test

Qtests <- function(series, k, fitdf = 0) {
  pvals <- matrix(NA, nrow = k, ncol = 2)
  colnames(pvals) <- c("lag", "pval")
  
  for (l in 1:k) {
    pvals[l, "lag"] <- l
    
    if (l <= fitdf) {
      pvals[l, "pval"] <- NA
    } else {
      pvals[l, "pval"] <- Box.test(
        series,
        lag = l,
        type = "Ljung-Box",
        fitdf = fitdf
      )$p.value
    }
  }
  
  return(pvals)
}

# Vérification de la validité du test ADF ----
# On teste alors les résidus pour vérifier l'absence
# d'autocorrélation résiduelle
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))

# les résidus apparaissent autocorrélés
# on ajoute des retards jusqu'à obtenir des résidus non autocorrélés
# On ajoute des retards jusqu'à obtenir des résidus décorrélés
series <- ind_diff1
kmax <- 24
adftype <- "nc"

adfTest_valid <- function(series, kmax, adftype) {
  k <- 0
  is_white_noise <- FALSE
  
  # Tant que les résidus restent autocorrélés
  # et que le nombre maximal de retards n'est pas dépassé
  while (is_white_noise == FALSE && k <= kmax) {
    # Message de suivi
    cat(paste0("ADF with ", k, " lags: residuals OK? "))
    
    # Calcul du test ADF avec k retards
    adf <- adfTest(series, lags = k, type = adftype)
    
    # Récupération des p-values des tests de Ljung-Box
    # sur les résidus de la régression ADF
    pvals <- Qtests(
      adf@test$lm$residuals,
      24,
      fitdf = length(adf@test$lm$coefficients)
    )[, 2]
    
    # Si aucune p-value n'est inférieure à 5 %,
    # on considère que les résidus sont décorrélés
    if (sum(pvals < 0.05, na.rm = TRUE) == 0) {
      is_white_noise <- 1
      cat("OK\n")
    } else {
      cat("nope\n")
      k <- k + 1
    }
  }
  
  # Si aucun nombre de retards <= kmax ne convient,
  # le dernier test obtenu ne doit pas être interprété
  if (is_white_noise == FALSE) {
    warning("Aucun lag <= kmax ne permet d'obtenir des résidus décorrélés.")
  }
  
  return(adf)
}
adf <- adfTest_valid(ind_diff1, 24, adftype = "nc")

Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
# aucun test de Ljung Box ne rejette l'hypothèse nulle
# on peut interprêter le test
adf
# L'hypothèse de racine unitaire est rejetée

# nous pouvons appliquer un modèle ARMA(p,q)

# Partie II : Modèles ARMA ----
# Méthodologie Box-Jenkins

# on s'intéresse aux autocorrélations pour déterminer q
acf(ind_diff1, main = "Autocorrélation partielle après transformation de Box-Cox")
# le plot des autocorrélations montre une autocorrélation positive pour
# les deux dernières valeurs; cela indique un q_max = 2

# on s'intéresse à l'autocorrélation partielle pour définir p
pacf(ind_diff1, main = "Autocorrélation partielle après transformation de Box-Cox")
# le graph suggère un p_max = 4
# on s'intéresse aux statistiques des différentes spécifications

# Test sur les paramètres ----
for (i in c(0:4)) {
  for (j in c(0:2)) {
   (sarima(ind_diff1, i, 0,j, details=FALSE))$ttable
  }
}

# les modèles pour lesquels tous les coefficients sont significatifs sont
# MA(1), MA(2), AR(1), AR(2), AR(3), AR(4) ARMA(1,1), ARMA(2,2)

# Test sur les résidus ----
# ces tests vont permettre d'écarter les modèles pas assez riches

for (i in c(0:4)) {
  for (j in c(0:2)) {
    tmp <- capture.output(
      model <- sarima(ind_diff1, i, 0, j, details = FALSE)$fit
    )
    # model <- invisible(sarima(ind_diff1, i, 0,j, details=FALSE))$fit
    print(paste0("ARMA(", i,",", j,")"))
    print(LjungBox(model))
  }
}

# on rejette MA(1), AR(1), AR(2), AR(3), AR(4)
# ce sont tous les modèles pour lesquel, sur au moins 1 niveau de lag
# la p-valeur est faible et on rejette H0 (H0 = pas d'autocorrélation entre 
# les résidus jusqu'au lag considéré)
# ce sont des modèles pour lesquels on a des autocorrélations entre les résidus

# Finalement, on garde MA(2), ARMA(1,1), ARMA(2,2)

ma2 <- sarima(ind_diff1, 0, 0, 2)
arma11 <- sarima(ind_diff1, 1, 0, 1)
arma22 <- sarima(ind_diff1, 2, 0, 2)

# Les deux graphiques montrent des comportements similaires des résidus
# les autocorrélations sont significativement nulles

# le graphique "Normal Q-Q Plot of Std Residuals" montre qu'ils
# sont distribués de façon gaussienne

# Pour les deux modèles, les résidus ont bien une allure de white noise

# Critères d'information ----
# Les critères d'information permettent de comparer différents modèles
# en sanctionnant les modèles trop complexes
# ils vont peut-être nous permettre de trancher entre MA(2) et ARMA(1,1)
AIC(ma2$fit, arma11$fit, arma22$fit)
BIC(ma2$fit, arma11$fit, arma22$fit)

# Les tests AIC et BIC sont les plus faibles pour le modèle ARMA(2,2)
# C'est donc le modèle qui semble le plus adapté

# Partie III : Prévision ----

model <- Arima(ind_diff1, order = c(2, 0, 2))
prediction <- forecast(model, h = 2, level = 95)

# graphique 
autoplot(prediction) +
  scale_x_continuous(labels = function(x) {
     format(seq(tail(ind$Date,-50)[1], by = "6 months", length.out=length(x)), "%Y-%m")
   })+
  ggtitle("Prévision de la série à t+2, région de confiance à 95%") +
  xlab("Temps") + ylab("Xt") +
  theme_minimal() +
  coord_cartesian(xlim = c(tail(time(ind_diff1),50)[1], tail(time(prediction$mean), 2)[2]))