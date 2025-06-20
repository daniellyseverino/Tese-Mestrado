#http://tabnet.datasus.gov.br/cgi/deftohtm.exe?sinannet/cnv/denguebbr.def

rm(list = ls())
library(tidyverse)
library(highcharter)

# Dados Dengue =================================================================

dengue23 = readxl::read_xlsx("Dados/Dengue/DataDengue2023.xlsx")
dengue24 = readxl::read_xlsx("Dados/Dengue/DataDengue2024.xlsx")

data_dengue = dengue23 %>% bind_rows(dengue24)

# Manipulacao dos Dados ========================================================

data_dengue = data_dengue %>% 
  pivot_longer(cols = 4:55, names_to = "t2", values_to = "n") %>% 
  dplyr::rename(t = `Semana epidem. 1º Sintomas(s)`, N = Total) %>% 
  mutate(
    t = as.numeric(gsub("[^0-9]", "", t)),
    t2 = as.numeric(gsub("[^0-9]", "", t2)),
    d = t2 - t,
    N = ifelse(N == "-"|is.na(N), 0, as.numeric(N)),
    n = ifelse(n == "-"|is.na(n), 0, as.numeric(n)),
    t = ifelse(Ano == 2023, t, t + 52)
  ) %>% 
  filter(d >= 0, UF != "Sudeste") %>% 
  select(-t2) %>% 
  filter(t >= 36 & t <= 96) %>% 
  mutate(t = t - 35)

data_dengue_br = data_dengue %>% 
  filter(UF == "Brasil") %>% 
  distinct(t, d, Ano, n, N)

data_dengue = data_dengue %>% 
  filter(UF != "Brasil") %>% 
  mutate(s = as.numeric(as.factor(UF)))

saveRDS(data_dengue_br, 'Dados/data_dengue_br.rds', version = 2)
saveRDS(data_dengue, 'Dados/data_dengue.rds', version = 2)

