

# D. TEMP-DRIVEN CHANGE IN ABUNDANCE SPECTRA =======================================
# for-loop to take random samples for distributions representing activation energies
# Then compare that to a projection with a constant temperature
sim <- 1:100

t <- c()
tt <- c()
groj <- c()
growth <- c()
data_list_s <- list()

for (i in sim) {
  
  t <- params@species_params
  
  t$ea_met <- ea$met[i]
  t$ea_int <- ea$int[i]
  t$ea_mor <- ea$mor[i]
  
  tt <- MizerParams(t, 
                    ea_gro = ea$gro[i],
                    ea_car = ea$car[i],
                    kappa_ben = kappa_ben,
                    kappa = kappa,
                    w_bb_cutoff = w_bb_cutoff,
                    w_pp_cutoff = w_pp_cutoff,
                    r_pp = r_pp,
                    r_bb = r_bb)
  
  proj <- project(tt, 
                  dt = dt,
                  effort = projectEffort_m,
                  temperature = projectTemp$temperature,
                  diet_steps = 10,
                  t_max = t_max,
                  t_ref = 10)   
  
  spec <- getSpectra(proj)
  
  spec$ea_met <- proj@params@species_params$ea_met[1]
  spec$ea_mor <- proj@params@species_params$ea_mor[1]
  spec$ea_int <- proj@params@species_params$ea_int[1]
  
  spec$ea_gro <- proj@params@ea_gro
  spec$ea_car <- proj@params@ea_car
  
  spec$sim <- i
  
  data_list_s[[i]] <- spec
  
}

str(data_list_s)

big_spectra_data <- dplyr::bind_rows(data_list_s)

str(big_spectra_data)

# Now 


#**** Plot size spectra ============================================================
# Test plotting
ggplot(big_spectra_data, aes(w, n, group = factor(sim))) +
  geom_line(size = 0.5, alpha = 0.1, color = "grey50") +
  labs(y = "Size (g)") +
  facet_wrap(~species, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_linetype_manual(values = c("solid", "twodash")) + 
  guides(color = FALSE) +
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  scale_y_log10() +
  scale_x_log10() +
  NULL



# Plot separately
plotSpec %>% filter(n > 0) %>% 
  ggplot(., aes(w, n, linetype = Scenario, color = Scenario)) + 
  geom_line(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = col[c(2,4)]) +
  labs(y = "Warm/Constant") +
  facet_wrap(~species, scales = "free") +
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  scale_y_log10() +
  scale_x_log10()  +
  NULL





# Plot ratio instead
spec_df_wide <- data.frame(w = conSpec$w,
                           n_diff = warmSpec$n / conSpec$n,
                           Species = conSpec$species)

spec_df_wide %>% filter(w < max(vbgedat$Weight_g)) %>% # Remove the large size-classes
  ggplot(., aes(w, n_diff, color = Species)) + 
  geom_hline(yintercept = 1, col = "red") +
  geom_line(size = 1.5, alpha = 0.8) +
  scale_color_viridis(discrete = TRUE) +
  labs(y = "Warm/Constant", x = "W (g)") +
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  scale_y_log10() +
  scale_x_log10()  +
  NULL









# E. VARY EFFORT GIVEN WARMING/NO WARMING ==========================================
# Here I should add a layer showing effort during calibration period instead - with and without warming
# I will do 5 sets of effort: 0.75*FMSY, 1.25FMSY, 1.5FMSY, 1.75*FSMY
# The control is no warming + FMSY

# FSMY effort
tail(projectEffort_ct, 1)[1]

# In this scenario, replace FMSY with the effort in the calibration

# 0.75FMSY
projectEffort_ct_075FMSY <- data.frame(projectEffort_ct)
projectEffort_ct_075FMSY[c(80:117), 1] <- tail(projectEffort_ct, 1)[1] * 0.75
projectEffort_ct_075FMSY[c(80:117), 2] <- tail(projectEffort_ct, 1)[2] * 0.75
projectEffort_ct_075FMSY[c(80:117), 3] <- tail(projectEffort_ct, 1)[3] * 0.75

# FMSY
projectEffort_ct_FMSY <- data.frame(projectEffort_ct)

# 1.25FMSY
projectEffort_ct_125FMSY <- data.frame(projectEffort_ct)
projectEffort_ct_125FMSY[c(80:117), 1] <- tail(projectEffort_ct, 1)[1] * 1.25
projectEffort_ct_125FMSY[c(80:117), 2] <- tail(projectEffort_ct, 1)[2] * 1.25
projectEffort_ct_125FMSY[c(80:117), 3] <- tail(projectEffort_ct, 1)[3] * 1.25

# 1.5FMSY
projectEffort_ct_15FMSY <- data.frame(projectEffort_ct)
projectEffort_ct_15FMSY[c(80:117), 1] <- tail(projectEffort_ct, 1)[1] * 1.5
projectEffort_ct_15FMSY[c(80:117), 2] <- tail(projectEffort_ct, 1)[2] * 1.5
projectEffort_ct_15FMSY[c(80:117), 3] <- tail(projectEffort_ct, 1)[3] * 1.5

# 1.75FMSY
projectEffort_ct_175FMSY <- data.frame(projectEffort_ct)
projectEffort_ct_175FMSY[c(80:117), 1] <- tail(projectEffort_ct, 1)[1] * 1.75
projectEffort_ct_175FMSY[c(80:117), 2] <- tail(projectEffort_ct, 1)[2] * 1.75
projectEffort_ct_175FMSY[c(80:117), 3] <- tail(projectEffort_ct, 1)[3] * 1.75

# Constant temperature models
m_cons_075 <- project(params3_upd, 
                      dt = 0.1,
                      effort = as.matrix(projectEffort_ct_075FMSY),
                      temperature = rep(10, nrow(projectEffort_ct_075FMSY)),
                      diet_steps = 10,
                      t_max = t_max,
                      t_ref = 10) 

m_cons_FMSY <- project(params3_upd, 
                       dt = 0.1,
                       effort = as.matrix(projectEffort_ct_FMSY),
                       temperature = rep(10, nrow(projectEffort_ct_FMSY)),
                       diet_steps = 10,
                       t_max = t_max,
                       t_ref = 10) 

m_cons_125 <- project(params3_upd, 
                      dt = 0.1,
                      effort = as.matrix(projectEffort_ct_125FMSY),
                      temperature = rep(10, nrow(projectEffort_ct_125FMSY)),
                      diet_steps = 10,
                      t_max = t_max,
                      t_ref = 10) 

m_cons_15 <- project(params3_upd, 
                     dt = 0.1,
                     effort = as.matrix(projectEffort_ct_15FMSY),
                     temperature = rep(10, nrow(projectEffort_ct_15FMSY)),
                     diet_steps = 10,
                     t_max = t_max,
                     t_ref = 10) 

m_cons_175 <- project(params3_upd, 
                      dt = 0.1,
                      effort = as.matrix(projectEffort_ct_175FMSY),
                      temperature = rep(10, nrow(projectEffort_ct_175FMSY)),
                      diet_steps = 10,
                      t_max = t_max,
                      t_ref = 10) 

# vbgecal temperature models
m_warm_075 <- project(params3_upd, 
                      dt = 0.1,
                      effort = as.matrix(projectEffort_ct_075FMSY),
                      temperature = projectTemp$temperature,
                      diet_steps = 10,
                      t_max = t_max,
                      t_ref = 10) 

m_warm_FMSY <- project(params3_upd, 
                       dt = 0.1,
                       effort = as.matrix(projectEffort_ct_FMSY),
                       temperature = projectTemp$temperature,
                       diet_steps = 10,
                       t_max = t_max,
                       t_ref = 10) 

m_warm_125 <- project(params3_upd, 
                      dt = 0.1,
                      effort = as.matrix(projectEffort_ct_125FMSY),
                      temperature = projectTemp$temperature,
                      diet_steps = 10,
                      t_max = t_max,
                      t_ref = 10) 

m_warm_15 <- project(params3_upd, 
                     dt = 0.1,
                     effort = as.matrix(projectEffort_ct_15FMSY),
                     temperature = projectTemp$temperature,
                     diet_steps = 10,
                     t_max = t_max,
                     t_ref = 10) 

m_warm_175 <- project(params3_upd, 
                      dt = 0.1,
                      effort = as.matrix(projectEffort_ct_175FMSY),
                      temperature = projectTemp$temperature,
                      diet_steps = 10,
                      t_max = t_max,
                      t_ref = 10) 


# Extract size-spectra from all models
# Starting with reference:
ref_spec <- getSpectra(m_cons_FMSY)

# Constant temperature
spec_m_cons_075 <- getSpectra(m_cons_075)
spec_m_cons_075$Fm <- "0.75*FMSY"
spec_m_cons_075$Temperature <- "Constant"
spec_m_cons_075$n_rel <- spec_m_cons_075$n / ref_spec$n

spec_m_cons_125 <- getSpectra(m_cons_125)
spec_m_cons_125$Fm <- "1.25*FMSY"
spec_m_cons_125$Temperature <- "Constant"
spec_m_cons_125$n_rel <- spec_m_cons_125$n / ref_spec$n

spec_m_cons_15 <- getSpectra(m_cons_15)
spec_m_cons_15$Fm <- "1.50*FMSY"
spec_m_cons_15$Temperature <- "Constant"
spec_m_cons_15$n_rel <- spec_m_cons_15$n / ref_spec$n

spec_m_cons_175 <- getSpectra(m_cons_175)
spec_m_cons_175$Fm <- "1.75*FMSY"
spec_m_cons_175$Temperature <- "Constant"
spec_m_cons_175$n_rel <- spec_m_cons_175$n / ref_spec$n

# Varying temperature
spec_m_warm_075 <- getSpectra(m_warm_075)
spec_m_warm_075$Fm <- "0.75*FMSY"
spec_m_warm_075$Temperature <- "Warming"
spec_m_warm_075$n_rel <- spec_m_warm_075$n / ref_spec$n

spec_m_warm_FMSY <- getSpectra(m_warm_FMSY)
spec_m_warm_FMSY$Fm <- "1.00*FMSY"
spec_m_warm_FMSY$Temperature <- "Warming"
spec_m_warm_FMSY$n_rel <- spec_m_warm_FMSY$n / ref_spec$n

spec_m_warm_125 <- getSpectra(m_warm_125)
spec_m_warm_125$Fm <- "1.25*FMSY"
spec_m_warm_125$Temperature <- "Warming"
spec_m_warm_125$n_rel <- spec_m_warm_125$n / ref_spec$n

spec_m_warm_15 <- getSpectra(m_warm_15)
spec_m_warm_15$Fm <- "1.50*FMSY"
spec_m_warm_15$Temperature <- "Warming"
spec_m_warm_15$n_rel <- spec_m_warm_15$n / ref_spec$n

spec_m_warm_175 <- getSpectra(m_warm_175)
spec_m_warm_175$Fm <- "1.75*FMSY"
spec_m_warm_175$Temperature <- "Warming"
spec_m_warm_175$n_rel <- spec_m_warm_175$n / ref_spec$n

# Rbind all data
warmFishSpectra <- rbind(spec_m_warm_075,
                         spec_m_warm_FMSY,
                         spec_m_warm_125,
                         spec_m_warm_15,
                         spec_m_warm_175,
                         spec_m_cons_075,
                         spec_m_cons_125,
                         spec_m_cons_15,
                         spec_m_cons_175)

warmFishSpectra %>% filter(w < max(vbgedat$Weight_g) & n_rel > 0) %>% # Remove the large size-classes
  ggplot(., aes(w, n_rel, color = Fm)) + 
  facet_grid(Temperature~species) +
  geom_hline(yintercept = 1, col = "black", linetype = "dashed") +
  geom_line(size = 1.7, alpha = 0.8) +
  scale_color_viridis(discrete = TRUE) +
  labs(y = "Warm/Constant", x = "W (g)") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  ylim(0.75, 1.35) +
  scale_x_log10() +
  NULL

# Different layout
warmFishSpectra %>% 
  filter(w < max(vbgedat$Weight_g) & n_rel > 0 & Fm %in% c("0.75*FMSY",
                                                           #"1.25*FMSY",
                                                           "1.50*FMSY")) %>% # Remove the large size-classes
  ggplot(., aes(w, n_rel, color = Fm, linetype = Temperature)) + 
  facet_wrap(~species, scales = "free_x") +
  geom_hline(yintercept = 1, col = "red", size = 0.5) +
  geom_line(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = col[c(1,4)]) +
  scale_linetype_manual(values = c("solid", "twodash")) +
  labs(y = "# / FMSY & Constant Temp.", x = "W (g)") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  ylim(0.75, 1.35) +
  scale_x_log10() +
  NULL


# R color brewer colors:
warmFishSpectra %>% 
  filter(w < max(vbgedat$Weight_g) & n_rel > 0 & Fm %in% c("0.75*FMSY",
                                                           #"1.25*FMSY",
                                                           "1.50*FMSY")) %>% # Remove the large size-classes
  ggplot(., aes(w, n_rel, color = Temperature, linetype = Fm)) + 
  facet_wrap(~species, scales = "free_x") +
  geom_hline(yintercept = 1, color = "black", size = 0.7) +
  geom_line(size = 1.4, alpha = 0.8) +
  scale_color_manual(values = c(brewer.pal(n = 8, name = "RdBu")[8], brewer.pal(n = 8, name = "RdBu")[1])) +
  scale_linetype_manual(values = c("solid", "twodash")) +
  labs(y = "# / FMSY & Constant Temp.", x = "W (g)") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  ylim(0.75, 1.35) +
  scale_x_log10() +
  NULL

# Interesting that fishing causes the accumulation of cod around 100g. Plot diet and compare 0.75*FMSY and 1.5*FMSY

#**** Plot Yield at 2050 ===========================================================
# 3. Yield. How about a simple value? Three columns (species), Three rows(F), two lines(warming or constant?)
# Create data frames of all models

# Warming
Yield_m_warm_075 <- data.frame(getYield(m_warm_075))[117, ]
Yield_m_warm_075$Temperature <- "Warming"
Yield_m_warm_075$Fm <- "0.75*FMSY"

Yield_m_warm_FMSY <- data.frame(getYield(m_warm_FMSY))[117, ]
Yield_m_warm_FMSY$Temperature <- "Warming"
Yield_m_warm_FMSY$Fm <- "1.00*FMSY"

Yield_m_warm_125 <- data.frame(getYield(m_warm_125))[117, ]
Yield_m_warm_125$Temperature <- "Warming"
Yield_m_warm_125$Fm <- "1.25*FMSY"

Yield_m_warm_15 <- data.frame(getYield(m_warm_15))[117, ]
Yield_m_warm_15$Temperature <- "Warming"
Yield_m_warm_15$Fm <- "1.50*FMSY"

Yield_m_warm_175 <- data.frame(getYield(m_warm_175))[117, ]
Yield_m_warm_175$Temperature <- "Warming"
Yield_m_warm_175$Fm <- "1.75*FMSY"

# Constant temperature
Yield_m_cons_075 <- data.frame(getYield(m_cons_075))[117, ]
Yield_m_cons_075$Temperature <- "Constant"
Yield_m_cons_075$Fm <- "0.75*FMSY"

Yield_m_cons_FMSY <- data.frame(getYield(m_cons_FMSY))[117, ]
Yield_m_cons_FMSY$Temperature <- "Constant"
Yield_m_cons_FMSY$Fm <- "1.00*FMSY"

Yield_m_cons_125 <- data.frame(getYield(m_cons_125))[117, ]
Yield_m_cons_125$Temperature <- "Constant"
Yield_m_cons_125$Fm <- "1.25*FMSY"

Yield_m_cons_15 <- data.frame(getYield(m_cons_15))[117, ]
Yield_m_cons_15$Temperature <- "Constant"
Yield_m_cons_15$Fm <- "1.50*FMSY"

Yield_m_cons_175 <- data.frame(getYield(m_cons_175))[117, ]
Yield_m_cons_175$Temperature <- "Constant"
Yield_m_cons_175$Fm <- "1.75*FMSY"

# Combine to single data frame
allYield <- rbind(Yield_m_cons_075,
                  Yield_m_cons_FMSY,
                  Yield_m_cons_125,
                  Yield_m_cons_15,
                  Yield_m_cons_175,
                  Yield_m_warm_075,
                  Yield_m_warm_FMSY,
                  Yield_m_warm_125,
                  Yield_m_warm_15,
                  Yield_m_warm_175)

# Make data frame long
allYield_long <- allYield %>% gather(Species, Yield, 1:3)

# Create relative yield (to Yield_m_cons_FMSY)
Yield_m_cons_FMSY

allYield_long$relYield <- 0

allYield_long$relYield <- ifelse(allYield_long$Species == "Cod",
                                 allYield_long$Yield / Yield_m_cons_FMSY$Cod,
                                 allYield_long$relYield)

allYield_long$relYield <- ifelse(allYield_long$Species == "Sprat",
                                 allYield_long$Yield / Yield_m_cons_FMSY$Sprat,
                                 allYield_long$relYield)

allYield_long$relYield <- ifelse(allYield_long$Species == "Herring",
                                 allYield_long$Yield / Yield_m_cons_FMSY$Herring,
                                 allYield_long$relYield)


ggplot(allYield_long, aes(Species, relYield, color = Temperature, shape = Fm)) +
  geom_jitter(size = 4, alpha = 0.6, width = 0.1, height = 0, stroke = 1.7) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~Species, scales = "free") +
  scale_color_manual(values = c(brewer.pal(n = 8, name = "RdBu")[8], brewer.pal(n = 8, name = "RdBu")[1])) +
  labs(y = "Relative Yield", x = "Species") +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 3/4) +
  NULL




