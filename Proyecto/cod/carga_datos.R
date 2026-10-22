library("quantmod")
library("ggplot2")

# Cargar series de acciones
# SPY: S&P 500
# QQQ: NASDAQ-100
# DIA: Dow Jones
# IWM: Russel 2000
series <- c("SPY", "QQQ", "DIA", "IWM")
getSymbols(series, src = "yahoo", from = "2017-01-01", to = Sys.Date())

# Calcular y graficar los rendimientos
# SPY:
r_SPY <- na.omit(diff(log(Ad(SPY))))
df_SPY <- data.frame(fecha = index(r_SPY),
                 rend  = as.numeric(r_SPY))
rends_SPY <- ggplot(df_SPY, aes(x = rend)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, fill = "steelblue",
                 color = "white", alpha = 0.8) +
  geom_density(color = "red", linewidth = 1) +
  labs(title = "Distribución de los rendimientos logarítmicos diarios de SPY",
       x = "Rendimiento logarítmico", y = "Densidad") +
  theme_minimal()
ggsave(filename = "figs/rends_SPY.png",
       plot     = rends_SPY,
       width    = 8, height = 5,
       dpi      = 300)
# QQQ:
r_QQQ <- na.omit(diff(log(Ad(QQQ))))
df_QQQ <- data.frame(fecha = index(r_QQQ),
                     rend  = as.numeric(r_QQQ))
rends_QQQ <- ggplot(df_QQQ, aes(x = rend)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, fill = "steelblue",
                 color = "white", alpha = 0.8) +
  geom_density(color = "red", linewidth = 1) +
  labs(title = "Distribución de los rendimientos logarítmicos diarios de QQQ",
       x = "Rendimiento logarítmico", y = "Densidad") +
  theme_minimal()
ggsave(filename = "figs/rends_QQQ.png",
       plot     = rends_QQQ,
       width    = 8, height = 5,
       dpi      = 300)
# DIA:
r_DIA <- na.omit(diff(log(Ad(DIA))))
df_DIA <- data.frame(fecha = index(r_DIA),
                     rend  = as.numeric(r_DIA))
rends_DIA <- ggplot(df_DIA, aes(x = rend)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, fill = "steelblue",
                 color = "white", alpha = 0.8) +
  geom_density(color = "red", linewidth = 1) +
  labs(title = "Distribución de los rendimientos logarítmicos diarios de DIA",
       x = "Rendimiento logarítmico", y = "Densidad") +
  theme_minimal()
ggsave(filename = "figs/rends_DIA.png",
       plot     = rends_DIA,
       width    = 8, height = 5,
       dpi      = 300)
# IWM:
r_IWM <- na.omit(diff(log(Ad(IWM))))
df_IWM <- data.frame(fecha = index(r_IWM),
                     rend  = as.numeric(r_IWM))
rends_IWM <- ggplot(df_IWM, aes(x = rend)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, fill = "steelblue",
                 color = "white", alpha = 0.8) +
  geom_density(color = "red", linewidth = 1) +
  labs(title = "Distribución de los rendimientos logarítmicos diarios de IWM",
       x = "Rendimiento logarítmico", y = "Densidad") +
  theme_minimal()
ggsave(filename = "figs/rends_IWM.png",
       plot     = rends_IWM,
       width    = 8, height = 5,
       dpi      = 300)