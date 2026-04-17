# Projet Time-Series
# 2A Paul Hobeika et Sara Chikhi

# Packages
library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(ggthemes)

# Import des données
ind <- read_csv2("donnees/valeurs_mensuelles.csv", col_names = TRUE) %>% 
  slice(4:437) %>% # on supprime les premières lignes qui ne correspondent pas à des données
  select(-Codes) %>% # on supprime aussi la colonne "Codes" qui est toujours égale à 'A'
  rename(Date = 'Libellé', 
         Valeur = "Indice CVS-CJO de la production industrielle (base 100 en 2021) - Fabrication de pesticides et d'autres produits agrochimiques (NAF rév. 2, niveau groupe, poste 20.2)") %>% 
  mutate(Date = ym(Date), # on met la variable de période au format 'Date'
         Valeur = as.numeric(Valeur)) # et l'indice en nombre plutôt qu'en caractères


# Graphe de la série
ind %>% 
  arrange(Date) %>% 
  ggplot(aes(x = Date, y = Valeur)) +
  geom_line() +
  theme_minimal() %>% 
  ggsave("donnes_brutes.png")


# Partie I
#=========

# 1. Que représente la série choisie ?



