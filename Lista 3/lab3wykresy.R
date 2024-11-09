n <- 5000
p <- 950
X <- matrix(rnorm(n*p, 0, 1/sqrt(5000)), ncol=p, nrow=n)
beta <- c(rep(3.5, 20), rep(0, p-20))
epsilon <- rnorm(n)
Y <- as.vector(X%*%beta+epsilon)
epsilon_gwiazdka <- rnorm(n)

zad_1 <- function(k){
  X_k <- X[, 1:k]
  beta_est <- as.vector(solve(t(X_k)%*%X_k)%*%t(X_k)%*%Y)
  Y_daszek <- as.vector(X_k%*%beta_est)
  RSS <- sum((Y_daszek-Y)^2)
  PE <- sum((X_k%*%beta[1:k]+epsilon_gwiazdka-X_k%*%beta_est)^2)
  PE_1 <- RSS+2*1*k
  PE_2 <- RSS+2*k*(RSS/(n-k))
  diag_M <- diag(X_k%*%solve(t(X_k)%*%X_k)%*%t(X_k)) #X(X'X)^(-1)X')
  CV <- sum(((Y-Y_daszek)/(1-diag_M))^2)
  return(c(abs(PE-PE_1), abs(PE-PE_2), abs(PE-CV)))
}
zad_1_vec <- Vectorize(zad_1)
k_vec <- c(10, 20, 30, 100, 500, 950)
zad_1_odp <- data.frame(t(zad_1_vec(k_vec)))
kolumny_1 <- c("PE_1", "PE_2", "CV")
colnames(zad_1_odp) = kolumny_1

PE_1_boxplot <- data.frame(matrix(0, ncol=6, nrow=100))
colnames(PE_1_boxplot) <- as.character(k_vec)
PE_2_boxplot <- data.frame(matrix(0, ncol=6, nrow=100))
colnames(PE_2_boxplot) <- as.character(k_vec)
CV_boxplot <- data.frame(matrix(0, ncol=6, nrow=100))
colnames(CV_boxplot) <- as.character(k_vec)

# Populate boxplot data frames
for(i in 1:100){
  epsilon_gwiazdka <- rnorm(n)
  values <- zad_1_vec(k_vec)
  PE_1_boxplot[i, ] <- values[1, ]
  PE_2_boxplot[i, ] <- values[2, ]
  CV_boxplot[i, ] <- values[3, ]
}

# Reshape data for ggplot
PE_1_data <- melt(PE_1_boxplot)
PE_2_data <- melt(PE_2_boxplot)
CV_data <- melt(CV_boxplot)

# Create boxplots using ggplot2
b1 <- ggplot(PE_1_data, aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(title = "|PE-PE1|", x = "Ilość zmiennych", y = "|PE-PE1|") +
  theme_minimal()

b2 <- ggplot(PE_2_data, aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(title = "|PE-PE2|", x = "Ilość zmiennych", y = "|PE-PE2|") +
  theme_minimal()

b3<- ggplot(CV_data, aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(title = "|PE-CV|", x = "Ilość zmiennych", y = "|PE-CV|") +
  theme_minimal()

CV_data_subset <- filter(CV_data, CV_data$variable %in% c(10,20,30,100,500))
b4 <- ggplot(CV_data, aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(title = "|PE-CV| bez p = 950", x = "Ilość zmiennych", y = "|PE-CV|") +
  theme_minimal()
grid.arrange(b1,b2,b3, ncol = 3)




## Zadanie 2 powtorzenia 100 do csv (dlugo sie robi)

#library(reshape)
#library(gridExtra)
set.seed(1)

aic_fun <- function(dane){
  wyn_aic=fast_forward(dane, crit='aic')
  ind_aic=sort(as.numeric(wyn_aic$model))
  TD_aic=sum(ind_aic<=20)
  FD_aic=sum(ind_aic>20)
  Power_aic=TD_aic/20
  FDP_aic=FD_aic/max(length(ind_aic), 1)
  X_p <- X_k[, ind_aic]
  beta_est_p <- as.vector(solve(t(X_p)%*%X_p)%*%t(X_p)%*%Y)
  Y_daszek_p <- as.vector(X_p%*%beta_est_p)
  SE_aic <- sum((Y_daszek_p-X%*%beta)^2)
  aic <- c(TD_aic, FD_aic, Power_aic, FDP_aic, SE_aic)
  return(aic)
}

bic_fun <- function(dane){
  wyn_bic=fast_forward(dane, crit='bic')
  ind_bic=sort(as.numeric(wyn_bic$model))
  TD_bic=sum(ind_bic<=20)
  FD_bic=sum(ind_bic>20)
  Power_bic=TD_bic/20
  FDP_bic=FD_bic/max(length(ind_bic), 1)
  X_p <- X_k[, ind_bic]
  beta_est_p <- as.vector(solve(t(X_p)%*%X_p)%*%t(X_p)%*%Y)
  Y_daszek_p <- as.vector(X_p%*%beta_est_p)
  SE_bic <- sum((Y_daszek_p-X%*%beta)^2)
  bic <- c(TD_bic, FD_bic, Power_bic, FDP_bic, SE_bic)
  return(bic)
}

ric_fun <- function(dane){
  wyn_ric=fast_forward(dane, crit='mbic', const=sqrt(n))
  ind_ric=sort(as.numeric(wyn_ric$model))
  TD_ric=sum(ind_ric<=20)
  FD_ric=sum(ind_ric>20)
  Power_ric=TD_ric/20
  FDP_ric=FD_ric/max(length(ind_ric), 1)
  X_p <- X_k[, ind_ric]
  beta_est_p <- as.vector(solve(t(X_p)%*%X_p)%*%t(X_p)%*%Y)
  Y_daszek_p <- as.vector(X_p%*%beta_est_p)
  SE_ric <- sum((Y_daszek_p-X%*%beta)^2)
  ric <- c(TD_ric, FD_ric, Power_ric, FDP_ric, SE_ric)
  return(ric)
}

mbic_fun <- function(dane){
  wyn_mbic=fast_forward(dane, crit='mbic')
  ind_mbic=sort(as.numeric(wyn_mbic$model))
  TD_mbic=sum(ind_mbic<=20)
  FD_mbic=sum(ind_mbic>20)
  Power_mbic=TD_mbic/20
  FDP_mbic=FD_mbic/max(length(ind_mbic), 1)
  X_p <- X_k[, ind_mbic]
  beta_est_p <- as.vector(solve(t(X_p)%*%X_p)%*%t(X_p)%*%Y)
  Y_daszek_p <- as.vector(X_p%*%beta_est_p)
  SE_mbic <- sum((Y_daszek_p-X%*%beta)^2)
  mbic <- c(TD_mbic, FD_mbic, Power_mbic, FDP_mbic, SE_mbic)
  return(mbic)
}

mbic2_fun <- function(dane){
  wyn_mbic2=fast_forward(dane, crit='mbic2')
  ind_mbic2=sort(as.numeric(wyn_mbic2$model))
  TD_mbic2=sum(ind_mbic2<=20)
  FD_mbic2=sum(ind_mbic2>20)
  Power_mbic2=TD_mbic2/20
  FDP_mbic2=FD_mbic2/max(length(ind_mbic2), 1)
  X_p <- X_k[, ind_mbic2]
  beta_est_p <- as.vector(solve(t(X_p)%*%X_p)%*%t(X_p)%*%Y)
  Y_daszek_p <- as.vector(X_p%*%beta_est_p)
  SE_mbic2 <- sum((Y_daszek_p-X%*%beta)^2)
  mbic2 <- c(TD_mbic2, FD_mbic2, Power_mbic2, FDP_mbic2, SE_mbic2)
  return(mbic2)
}


n <- 5000
p <- 950
X <- matrix(rnorm(n*p, 0, 1/sqrt(5000)), ncol=p, nrow=n)
beta <- c(rep(3.5, 20), rep(0, p-20))
epsilon <- rnorm(n)
Y <- as.vector(X%*%beta+epsilon)
epsilon_gwiazdka <- rnorm(n)

zad_2 <- function(k){
  X_k = X[,1:k]
  d <- prepare_data(Y, X_k)
  AIC <- round(aic_fun(d),3)
  BIC <- round(bic_fun(d),3)
  RIC <- round(ric_fun(d),3)
  mBIC <- round(mbic_fun(d),3)
  mBIC2 <- round(mbic2_fun(d),3)
  wyniki <- t(cbind(AIC, BIC, RIC, mBIC, mBIC2))
  return(wyniki)
}


rep = 100
zad_3 <- function(k){
  wynik_moc <- matrix(NA, nrow = 5, ncol = rep)
  wynik_FDR<- matrix(NA, nrow = 5, ncol = rep)
  wynik_SE<- matrix(NA, nrow = 5, ncol = rep)
  for(i in 1:rep){
    epsilon <- rnorm(n)
    Y <- X %*% beta + epsilon
    X_k <- X[,1:k]
    wynik_moc[,i] <- zad_2(k)[,3]
    wynik_FDR[,i] <- zad_2(k)[,4]
    wynik_SE[,i] <- zad_2(k)[,5]
  }
  result <- cbind(rowMeans(wynik_moc),rowMeans(wynik_FDR), rowMeans(wynik_SE))
  return(result)
  
}

#X_k = X[,1:50]
#zad_3_50 <- zad_3(50)
X_k = X[,1:100]
zad_3_100 <- zad_3(100)
X_k = X[,1:200]
zad_3_200 <- zad_3(200)
X_k = X[,1:500]
zad_3_500 <- zad_3(500)
X_k = X[,1:950]
zad_3_950 <- zad_3(950)

zad_3_wynik <- cbind(zad_3_100,zad_3_200,zad_3_500,zad_3_950)
colnames(zad_3_wynik) <- c( "Power", "FDR", "SE", "Power", "FDR", "SE","Power", "FDR", "SE","Power", "FDR", "SE")


rownames(zad_3_wynik) <- c("AIC", "BIC", "RIC", "mBIC", "mBIC2")



write.csv(zad_3_wynik, "wynikii_v4_100_5000.csv")


5.292/(sqrt(24.01)*sqrt(29.16))
2.254/(sqrt(24.01)*sqrt(5.29))
2.484/(sqrt(29.16)*sqrt(5.29))
