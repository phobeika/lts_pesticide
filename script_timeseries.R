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
ggsave(filename = "donneesbrutes.png", plot = ind, device = "png")



# Partie I
#=========

# 1. Que représente la série choisie ?


# On fait une transformation de Box-Cox pour diminuer la variations de la
# variance en fonction du temps
lambda <- BoxCox.lambda(ind$Valeur)
ind_bc <- BoxCox(ind$Valeur, lambda)
plot(ind$Date, ind$Valeur, type = "l")

# Puis on prend la "first difference" pour supprimer les variations "lentes"
ind_diff1 <- diff(ind$Valeur)
plot(ind$Date[-1], ind_diff1, type = "l") # Ça semble relativement stationnaire

