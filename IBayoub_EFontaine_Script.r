install.packages("xts")
install.packages("tseries")
install.packages("urca")
install.packages("forecast")

library(forecast)
library(urca)
library(tseries)
library(xts)
library(zoo)

#Option pour voir les graphiques dans VS Code (si supporté)
options(vsc.plot = TRUE)


# 1) Import et formatage des données


raw_data <- read.csv("/home/onyxia/work/Time_Series_Project/valeurs_mensuelles.csv", sep = ";", skip = 4, header = FALSE, stringsAsFactors = FALSE)

#Ne conserver que les lignes où la première colonne ressemble à une date AAAA-MM
raw_data <- raw_data[grepl("^\\d{4}-\\d{2}", raw_data$V1), ]

#Renommer les colonnes
colnames(raw_data) <- c("Période", "Valeur", "Code")

#Transformer la colonne Période en date au format standard (1er du mois)
raw_data$Date <- as.Date(paste0(raw_data$Période, "-01"))

#Convertir la colonne Valeur en numérique (remplace la virgule par un point si besoin)
raw_data$Valeur <- as.numeric(gsub(",", ".", raw_data$Valeur))

#Réorganiser les colonnes
data_clean <- raw_data[, c("Date", "Valeur", "Code")]

#Afficher les premières lignes
head(data_clean)


# 2) Première approche de la série


#Affiche la tendance : Serie0
png("/home/onyxia/work/Time_Series_Project/Graphes_ACF/Serie0_brute.png", width = 800, height = 600)

plot(data_clean$Date, data_clean$Valeur, type = "l", col = "black",
     main = "Série brute – Indice brut de production",
     xlab = "Time", ylab = "X_t")

dev.off()


# Observons le comportement de la série: elle semble 
# suivre une tendance quadratique avec une saisonnalité (annuelle?).
# Elle semble également hétéroscedastique. 


#Appliquons alors un transformation de Box-Cox (à cause de la tendance 
#quadratique, hétéroscédasticité, positivité et sauts) pour corriger ce caractère 
#hétéroscédastique et la linéariser.

# Conversion de ta série en ts 
ts_data <- ts(data_clean$Valeur, start = c(1989, 1), frequency = 12) 

# Calcul du lambda optimal
lambda <- BoxCox.lambda(ts_data)

# Application de la transformation Box-Cox
ts_boxcox <- BoxCox(ts_data, lambda)

# Visualisation: Serie1
png("/home/onyxia/work/Time_Series_Project/Graphes_ACF/Serie1_BoxCox.png", width = 800, height = 600)

plot(ts_boxcox, main = paste("Série après transformation Box-Cox (lambda =", round(lambda, 2), ")"), ylab = "G_t")

dev.off()


#La transformation a réduit l’hétéroscédasticité 
#(la variance est plus constante qu’avant).
#Toutefois, fa forme en U subsiste : on voit encore une tendance 
#non linéaire (pas stationnaire en moyenne).
#Aucun motif saisonnier clair n’apparaît à l’œil nu, 
#mais ça reste à tester formellement.

#Avant toute différenciation supplémentaire, regardon l'ACF de la série
#pour détecter : une tendance résiduelle et une saisonnalité.


#ACF1

png("/home/onyxia/work/Time_Series_Project/Graphes_ACF/ACF1_BoxCox.png", width = 800, height = 600)

acf(ts_boxcox, main = "ACF – Série après transformation BoxCox")

dev.off()

test_kpss <- kpss.test(ts_boxcox)
test_kpss$statistic
# statistique = 1.0537 > 0.216, le seuil critique de 1 % 
# La p-value est donc inférieure à 1 %, ce qui conduit à rejeter l’hypothèse nulle de stationnarité.


#Avec l'ACF, nous observons que la première barre à lag 1 
#est très élevée (~0.9), les valeurs décroissent lentement, 
#de manière régulière. Enfin, il n'y a pas de pic net à un 
#lag spécifique, mais un effet persistant.

#Nous pouvons en conclure que c'est le signe d’une tendance non stationnaire.
#L’ACF ne coupe pas brutalement : elle décroît lentement, ce qui 
#indique que la série transformée n’est pas encore stationnaire.
#Il n’y a pas de motif saisonnier net (comme des pics tous les 12 lags), 
#donc pas de saisonnalité apparente à ce stade.

#Poursuivons donc la transformation de notre série pour obtenir 
#une série stationnaire: différention la.

#Différenciation à l'ordre 1 de la série: Serie2

ts_boxcox_diff <- diff(ts_boxcox)


png("/home/onyxia/work/Time_Series_Project/Graphes_ACF/Serie2_Diff.png", width = 800, height = 600)

plot(ts_boxcox_diff, main = "Série après transformation BoxCox et différenciation", ylab = "Z_t")


dev.off()


#ACF de la série différenciée: ACF2

png("/home/onyxia/work/Time_Series_Project/Graphes_ACF/ACF2_Diff.png", width = 800, height = 600)

acf(ts_boxcox_diff, main = "ACF de la série différenciée")

dev.off()



boucle_adf_ljung <- function(serie, max_k = 10) {
  cat("Test ADF + Ljung-Box sur les résidus\n")
  serie <- na.omit(serie)
  n <- length(serie)
  
  for (k in 0:max_k) {
    # Test ADF (ur.df) avec k retards
    test_adf <- ur.df(serie, type = "drift", lags = k)
    stat_adf <- test_adf@teststat[1]
    seuil_adf <- test_adf@cval[1, "5pct"]

    # Construction des variables pour régression linéaire
    y <- diff(serie)[(k+1):(n-1)]
    z_lag1 <- serie[(k+1):(n-1)]
    
    if (k >= 1) {
      z_diff_lags <- embed(diff(serie), k + 1)[, -1]
      fit <- lm(y ~ z_lag1 + z_diff_lags)
    } else {
      fit <- lm(y ~ z_lag1)
    }

    res <- residuals(fit)
    pval_ljung <- Box.test(res, lag = 10, type = "Ljung-Box")$p.value

    cat(sprintf("k = %d | ADF stat = %.3f | seuil 5%% = %.3f | Ljung-Box p = %.4f",
                k, stat_adf, seuil_adf, pval_ljung), "\n")
    
    if (pval_ljung > 0.05 && stat_adf < seuil_adf) {
      cat(sprintf(" Meilleur k trouvé : %d (résidus non autocorrélés, ADF significatif)\n", k))
      break
    }
  }
}


boucle_adf_ljung(ts_boxcox_diff)  # ou ts_dG ou ts_ddG selon la série testée

#Graphiquement, la série semble osciller autour d’une moyenne 
#constante, avec une variance plus stable.
#L’amplitude est homogène (entre -10 et +10).

#Toutefois, l'ACF met en valeur des pics d’autocorrélation 
#notables à plusieurs lags, bien en dehors des bandes de confiance.
#Ces pics reviennent de façon régulière (environ tous les 12 lags).

#Il y aurait donc une structure saisonnière résiduelle, 
#probablement annuelle (tous les 12 mois).

#Appliquons alors une différenciation saisonniaire.

ts_boxcox_diff_seas <- diff(ts_boxcox_diff, lag = 12)


png("/home/onyxia/work/Time_Series_Project/Graphes_ACF/Serie3_Saisons.png", width = 800, height = 600)

plot(ts_boxcox_diff_seas, main = "Série après la différenciation saisonnière",ylab="Y_t")

dev.off()

png("/home/onyxia/work/Time_Series_Project/Graphes_ACF/ACF3_Saisons.png", width = 800, height = 600)

acf(ts_boxcox_diff_seas, main = "ACF de la série saisonnièrement différenciée")


dev.off()

png("/home/onyxia/work/Time_Series_Project/Graphes_ACF/PACF3_Saisons.png", width = 800, height = 600)
pacf(ts_boxcox_diff_seas, main = "PACF de la série saisonnièrement différenciée")
dev.off()
#La plupart des autocorrélations sont dans les bandes de confiance.
#Un seul pic légèrement au-dessus des bandes à lag 13.
#Les autres lags (y compris multiples de 12) sont faibles.
#Ce pic n’invalide pas la stationnarité, mais indique qu’il faut 
#inclure un terme saisonnier AR ou MA dans notre modèle.

#Maintenant, faisons un test ADF pour valider l'hypothèse de 
#stationnarité de notre série différenciée.


# Test ADF sans tendance (car la série est différenciée)
boucle_adf_ljung(ts_ddG)

# Le résultat nous confirme que la série est bien stationnaire.

#Nous pouvons dès lors procéder à l’identification du modèle ARIMA.


png("/home/onyxia/work/Time_Series_Project/Graphes_ACF/ACF_identification_model.png", width = 800, height = 600)
acf(ts_boxcox_diff_seas, main = "")
dev.off()

png("/home/onyxia/work/Time_Series_Project/Graphes_ACF/PACF_identification_model.png", width = 800, height = 600)
pacf(ts_boxcox_diff_seas, lag.max = 48, main = "")
dev.off()




# D'après ce que nous observons les combinaisons de paramètres candidates à tester se situent dans l’espace suivant : 
# p∈{0,1,2,3,4,5}, q∈{0,1}, P∈{0,1}, Q∈{0,1} avec d=1, D=1, s=12.

library(forecast)

# Série stationnaire
fit_auto <- auto.arima(ts_boxcox_diff_seas, stationary = TRUE)

summary(fit_auto)

# Ajuster un modèle ARIMA(0,0,1)(0,0,1)[12] sur une série déjà différenciée (stationnaire)
fit <- Arima(ts_boxcox_diff_seas, 
             order = c(0, 0, 1), 
             seasonal = list(order = c(0, 0, 0), period = 12), 
             include.mean = FALSE)

# Résumé du modèle
print("Résumé du modèle :")
print(summary(fit))

# Extraire les coefficients
coefs <- coef(fit)
print("Coefficients estimés :")
print(coefs)

# Erreurs standards des coefficients
erreurs <- sqrt(diag(vcov(fit)))
print("Erreurs standards :")
print(erreurs)

# Statistiques t
t_stat <- coefs / erreurs
print("Statistiques t :")
print(t_stat)

# Valeurs p
p_values <- 2 * (1 - pnorm(abs(t_stat)))
print("Valeurs p :")
print(p_values)

# Test de Ljung-Box sur les résidus
ljung_pval <- Box.test(residuals(fit), lag = 6, type = "Ljung-Box")$p.value
print("p-value du test de Ljung-Box (lag = 6) :")
print(ljung_pval)


# Afficher le résumé du modèle

#ARIMA(2,0,7)(0,0,1)[12] with zero mean 

#Coefficients:
#          ar1      ar2     ma1     ma2      ma3      ma4      ma5     ma6      ma7     sma1
#       -1.1462   -0.9795   0.5644  0.2832  -0.6809  -0.0228  -0.0347  0.0735  -0.0151  -0.7317
#s.e.   0.0169   0.0144  0.0537  0.0583   0.0592   0.0665   0.0655  0.0687   0.0535   0.0392

#sigma^2 = 1.188:  log likelihood = -617.34
#AIC=1256.67   AICc=1257.34   BIC=1300.85

#Training set error measures:
#                     ME    RMSE       MAE       MPE     MAPE      MASE         ACF1
#Training set 0.05395005 1.07679 0.7402778 -242.6238 667.4359 0.4307027 -0.004582152

#png("/home/onyxia/work/Time_Series_Project/Graphes_ACF/residuals.png", width = 800, height = 600)

checkresiduals(fit)
dev.off()


# Ljung-Box test

#data:  Residuals from ARIMA(2,0,7)(0,0,1)[12] with zero mean
#Q* = 22.193, df = 14, p-value = 0.07473

#Model df: 10.   Total lags used: 24

Box.test(residuals(fit), lag = 20, type = "Ljung-Box")

#Box-Ljung test
#data:  residuals(fit)
#X-squared = 19.414, df = 20, p-value = 0.4951

pred <- forecast(fit, h = 2, level = 95)

library(ggplot2)

# Affichage de la prévision avec région de confiance

png("/home/onyxia/work/Time_Series_Project/Graphes_ACF/prevision.png", width = 800, height = 600)

autoplot(pred) +
  ggtitle("Prévision de la série Xt avec région de confiance à 95%") +
  xlab("Temps") + ylab("Xt") +
  coord_cartesian(xlim = c(tail(time(ts_boxcox_diff_seas), 50)[1], tail(time(pred$mean), 2)[2]))

dev.off()


# Définir les modèles ARIMA à tester
modeles <- list(
  list(order = c(0,0,1), seasonal = c(0,0,0)),
  list(order = c(1,0,1), seasonal = c(0,0,0)),
  list(order = c(2,0,1), seasonal = c(0,0,0)),
  list(order = c(3,0,1), seasonal = c(0,0,0)),
  list(order = c(4,0,1), seasonal = c(0,0,0)),
  
  list(order = c(0,0,2), seasonal = c(0,0,0)),
  list(order = c(1,0,2), seasonal = c(0,0,0)),
  list(order = c(3,0,2), seasonal = c(0,0,0)),
  list(order = c(4,0,2), seasonal = c(0,0,0)),
  
  list(order = c(2,0,1), seasonal = c(0,0,1)),
  list(order = c(3,0,1), seasonal = c(0,0,1)),
  list(order = c(4,0,1), seasonal = c(0,0,1)),
  
  list(order = c(2,0,1), seasonal = c(1,0,1)),
  list(order = c(3,0,1), seasonal = c(1,0,1)),
  list(order = c(4,0,1), seasonal = c(1,0,1))
)

# Stocker les résultats
resultats <- data.frame()

# Boucle sur chaque modèle
for (i in seq_along(modeles)) {
  mod <- modeles[[i]]
  
  # Ajuster le modèle
  fit <- tryCatch({
    Arima(ts_boxcox_diff_seas,
          order = mod$order,
          seasonal = list(order = mod$seasonal, period = 12),
          include.mean = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(fit)) {
    coefs <- coef(fit)
    erreurs <- sqrt(diag(vcov(fit)))
    t_stat <- coefs / erreurs
    p_values <- 2 * (1 - pnorm(abs(t_stat)))
    ljung_test <- Box.test(residuals(fit), lag = 6, type = "Ljung-Box")
    
    for (j in seq_along(coefs)) {
      resultats <- rbind(resultats, data.frame(
        Modele = paste0("ARIMA(", paste(mod$order, collapse = ","), ")(", 
                        paste(mod$seasonal, collapse = ","), ")[12]"),
        Coef = names(coefs)[j],
        Estimate = coefs[j],
        Std_Error = erreurs[j],
        t_stat = t_stat[j],
        p_value = p_values[j],
        Ljung_Box_Stat = ljung_test$statistic,
        Ljung_Box_pval = ljung_test$p.value,
        AIC = AIC(fit),
        BIC = BIC(fit)
      ))
    }
  }
}

# Afficher le tableau final
print(resultats)

# Optionnel : sauvegarde en CSV
# write.csv(resultats, "resultats_arima.csv", row.names = FALSE)



