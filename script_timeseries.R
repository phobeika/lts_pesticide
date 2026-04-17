# Projet Time-Series
# 2A Paul Hobeika et Sara Chikhi

# pacakges
library(dplyr)
library(readr)
library(lubridate)
# import des données
ind <- read_csv2("donnees/valeurs_mensuelles.csv", col_names = TRUE) %>% 
  slice(4:437) %>% # on supprime les premières lignes qui ne correspondent pas à des données
  select(-Codes) %>% # on supprime aussi la colonne "Codes" qui est toujours égale à 'A'
  rename(Période = 'Libellé', Indice = "Indice CVS-CJO de la production industrielle (base 100 en 2021) - Fabrication de pesticides et d'autres produits agrochimiques (NAF rév. 2, niveau groupe, poste 20.2)") %>% 
  mutate(Période = ym(Période))

# Graphe de la série



# Partie I
#=========

# 1. Que représente la série choisie ?



