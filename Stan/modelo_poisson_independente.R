# Autora: Danielly Santos Sevrino

# Estimacao com estrutura independente - Familia Poisson para as contagens de casos

rm(list = ls())
options(OutDec = ",") 

# Pacotes ======================================================================

library(rstan)
library(tidyverse)

# Funcoes ======================================================================

genLog = function(t, a, b, c, f, logScale = TRUE){
  logV = log(f)+log(a)+log(c)-(c*t)-(f+1)*log( b+exp(-c*t) )
  if (logScale){
    return(logV);
  } else {
    return(exp(logV));
  }
}

# Carregando estimativas por delay =============================================

estimativas_por_delay = readRDS("Estimativas/Familia Poisson/estimativas_por_delay.rds")
abcf_por_dalay = readRDS("Estimativas/Familia Poisson/abcf_por_dalay.rds")
estimativas_theta_isolado = readRDS("Estimativas/Familia Poisson/estimativas_theta_isolado.rds")

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


# Modelo =======================================================================

# separando o delay 0 (novo 1)
td_dalay_1 = (dados_dengue_longo_completo %>% filter(d == 1)) %>% select(t,d)
n_dalay_1 = (dados_dengue_longo_completo %>% filter(d == 1)) %>% 
  select(n_td) %>% 
  mutate(n_td = as.numeric(n_td))

td_dalay_k = dados_dengue_longo_completo %>% filter(d != 1) %>% select(t,d)
n_dalay_k = dados_dengue_longo_completo %>% filter(d != 1) %>%
  select(n_td) %>% 
  mutate(n_td = as.numeric(n_td))

dados_stan = list(
  nk = n_dalay_1$n_td, n_k = n_dalay_k$n_td,
  T = 35, D = D,
  Tk = td_dalay_1$t, Dk = td_dalay_1$d,
  T_k =  td_dalay_k$t,  D_k =  td_dalay_k$d,
  qk = nrow(n_dalay_1), q_k = nrow(n_dalay_k)
)

modelo_stan_indep =  stan_model("Stan/modelLogistic.stan")

warmup = 5000
chains = 1
thin = 1
sample_size = 10000
number_interations = warmup + thin*sample_size

params = c("lambda",#"n_k"
           "alpha", "a_alpha", "b_alpha", "c_alpha", "f_alpha",
           "theta", "a_theta", "b_theta", "c_theta", "f_theta",
           "psi")

a_theta_inicial = mean(estimativas_theta_isolado$a) #191,7288
b_theta_inicial = mean(estimativas_theta_isolado$b) #0,002143543
c_theta_inicial = mean(estimativas_theta_isolado$c) #0,3907127
f_theta_inicial = mean(estimativas_theta_isolado$f) #1,000318

theta_t_inicial = genLog(t = 1:35, a = a_theta_inicial, 
                         b = b_theta_inicial, 
                         c = c_theta_inicial, 
                         f = f_theta_inicial, logScale = FALSE)

a_alpha_inicial = list()
b_alpha_inicial = list()
c_alpha_inicial = list()
f_alpha_inicial =list()
alfha_t_d_inicial = list()
for (d in 1:(D-1)) {
  
  a_alpha_inicial[[d]] = abcf_por_dalay[[d]][1]
  b_alpha_inicial[[d]] = abcf_por_dalay[[d]][2]
  c_alpha_inicial[[d]] = abcf_por_dalay[[d]][3]
  f_alpha_inicial[[d]] = abcf_por_dalay[[d]][4]
  
  
  alfha_t_d_inicial[[d]] = genLog(t = 1:35, a = abcf_por_dalay[[d]][1],
                                  b = abcf_por_dalay[[d]][2], 
                                  c = abcf_por_dalay[[d]][3], 
                                  f = abcf_por_dalay[[d]][4], logScale = TRUE)
}

alfha_t_d_inicial = matrix(unlist(alfha_t_d_inicial), ncol = 10)
a_alpha_inicial = unlist(a_alpha_inicial)
b_alpha_inicial = unlist(b_alpha_inicial)
c_alpha_inicial = unlist(c_alpha_inicial)
f_alpha_inicial = unlist(f_alpha_inicial)


parametros_iniciais = list(
  list(
    
    # alpha = alfha_t_d_inicial,
# 
#     a_alpha = a_alpha_inicial,
#     b_alpha = b_alpha_inicial,
#     c_alpha = c_alpha_inicial,
#     f_alpha = f_alpha_inicial,
    
    theta = theta_t_inicial,
    
    a_theta = a_theta_inicial,
    b_theta = b_theta_inicial,
    c_theta = c_theta_inicial,
    f_theta = f_theta_inicial,
    
    psi = n_dalay_1$n_td
  )
)

output_modelo_indep = rstan::sampling(modelo_stan_indep,
                                      data = dados_stan,
                                      iter = number_interations,
                                      warmup = warmup,
                                      chains = chains,
                                      pars = params,
                                      init = parametros_iniciais,
                                      verbose = FALSE)

estimativas_indep = rstan::extract(output_modelo_indep)
saveRDS(estimativas_indep, "Estimativas/Familia Poisson/estimativas_estrutura_independente.rds", version = 2)
estimativas_indep = readRDS("Estimativas/Familia Poisson/estimativas_estrutura_independente.rds")

ts.plot(estimativas_indep$a_theta)
ts.plot(estimativas_indep$b_theta)
ts.plot(estimativas_indep$c_theta)
ts.plot(estimativas_indep$f_theta)
acf(estimativas_indep$a_theta)
acf(estimativas_indep$b_theta)
acf(estimativas_indep$c_theta)
acf(estimativas_indep$f_theta)