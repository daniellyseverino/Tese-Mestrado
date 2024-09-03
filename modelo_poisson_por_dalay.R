# Autora: Danielly Santos Sevrino

# Estimacao por delay - Familia Poisson para as contagens de casos

rm(list = ls())
options(OutDec = ",") 

# Pacotes ======================================================================

library(rstan)
library(tidyverse)

# Dados para modelagem em T = 25 ===============================================

dados_dengue = readr::read_rds("dados/dengueData.RDS")
# filtrando apenas as primeiras 25 semanas (35 semanas = primeira onda)
dados_dengue = dados_dengue[1:25,]

T = 25
D = dim(dados_dengue)[2]

# criando os NAs
dados_dengue[outer(1:T, 0:(D - 1), FUN = "+") > T] <- NA
dados_NA = matrix(nrow = 35-T, ncol = D) %>% as.data.frame()
colnames(dados_NA) = paste0("d",0:(D-1))
dados_dengue = rbind(dados_dengue, dados_NA)


dados_dengue_longo <- dados_dengue %>% 
  rownames_to_column(var = "t") %>% 
  pivot_longer(cols = 2:(D+1), names_to = "d", values_to = "n_td") %>% 
  mutate(delay = d,
         delay = factor(delay, levels = c("d0","d1","d2","d3","d4","d5","d6","d7","d8","d9","d10")),
         d = str_remove_all(d,"d"),
         d = as.numeric(d)) %>% 
  group_by(d) %>% 
  mutate(t = rep(1:(T + (35-T)))) %>% 
  ungroup()

dados_dengue$N_t = rowSums(dados_dengue, na.rm = T)

dados_dengue_longo_completo = dados_dengue_longo %>% 
  mutate(d = d + 1) %>% 
  filter(!is.na(n_td)) %>% 
  arrange(d)


# Estimativas independentes por delay ==========================================

# Estimativas por coluna para d = 0, ..., D

# modelo_stan =  stan_model("Stan/modelLogisticbyDelay.stan")

stanLogistic = function(stanModel, y){
  
  y = y[!is.na(y)]
  n = length(y)
  
  # -> Stan configuration
  warmup = 1000
  chains = 1
  thin = 1
  sample_size = 10000
  number_interations = warmup + thin*sample_size
  
  # -> Preparing model
  
  params = c("mu", "a", "b", "c", "f")
  
  output = sampling(stanModel,
                    data = list(n = n, y = y),
                    iter = number_interations,
                    warmup = warmup,
                    chains = chains,
                    pars = params,
                    verbose = FALSE)
  
  # -> Extracting samples
  rstan::extract(output)
  
}

# para cada tempo do dalay i, vamos ter 35(T) - i estimativas de mu com mil repeticoes
# temos 10 mil estimativas de a,b,c,f

# estimativas_por_delay = apply(dados_dengue, 2, function(y) stanLogistic(modelo_stan, y = y))
# saveRDS(estimativas_por_delay, "Estimativas/Familia Poisson/estimativas_por_delay.rds",version = 2)
estimativas_por_delay = readRDS("Estimativas/Familia Poisson/estimativas_por_delay.rds")

# abcf_por_dalay = lapply(estimativas_por_delay, function(x) c(mean(x$a), mean(x$b), mean(x$c), mean(x$f)) )
# saveRDS(abcf_por_dalay, "Estimativas/Familia Poisson//abcf_por_dalay.rds", version = 2)
abcf_por_dalay = readRDS("Estimativas/Familia Poisson/abcf_por_dalay.rds")

# estimativas_theta_isolado = stanLogistic(modelo_stan, y = rowSums(dados_dengue, na.rm = TRUE))
# saveRDS(estimativas_theta_isolado, "Estimativas/Familia Poisson/estimativas_theta_isolado.rds",version = 2)
estimativas_theta_isolado = readRDS("Estimativas/Familia Poisson/estimativas_theta_isolado.rds")

ts.plot(estimativas_theta_isolado$a)
ts.plot(estimativas_theta_isolado$b)
ts.plot(estimativas_theta_isolado$c)
ts.plot(estimativas_theta_isolado$f)
