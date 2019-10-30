## IN PREP

# Borrow code from calibration.
# Do variations of the following analysis:

# F. Vary effort given warming scenarios
# G. Test time (temperature)-varying kappa

# Some preliminary code for that is below


#**** Plot growth rates ============================================================
# Here I need to extract the data that goes into the the plotGrowthCurves function.
# Copy and edit the function to only extract data. Do that for multiple models and plot

# 1. Changes in growth: calibration period relative to warming. 

# Define new function
getGrowth <- function(object, species,
                      max_age = 20, percentage = FALSE, print_it = TRUE) {
  if (is(object, "MizerSim")) {
    sim <- object
    if (missing(species)) {
      species <- dimnames(sim@n)$sp
    }
    # reorder list of species to coincide with order in sim
    idx <- which(dimnames(sim@n)$sp %in% species)
    species <- dimnames(sim@n)$sp[idx]
    age <- seq(0, max_age, length.out = 50)
    ws <- array(dim = c(length(species), length(age)),
                dimnames = list("Species" = species, "Age" = age))
    g <- getEGrowth(sim@params, sim@n[dim(sim@n)[1], , ], 
                    sim@n_pp[dim(sim@n)[1], ], sim@n_bb[dim(sim@n)[1], ], sim@n_aa[dim(sim@n)[1], ], sim@intTempScalar[,,1], sim@metTempScalar[,,1]) #AA
    for (j in 1:length(species)) {
      i <- idx[j]
      g_fn <- stats::approxfun(sim@params@w, g[i, ])
      myodefun <- function(t, state, parameters){
        return(list(g_fn(state)))
      }
      ws[j, ] <- deSolve::ode(y = sim@params@species_params$w_min[i],
                              times = age, func = myodefun)[, 2]
      if (percentage) {
        ws[j, ] <- ws[j, ] / sim@params@species_params$w_inf[i] * 100
      }
    }
    plot_dat <- reshape2::melt(ws)
    return(plot_dat) # Added this to extract growth data in data.frame
  }    
}

# Test if correct:
# growth_warm <- getGrowth(m6)
# plotGrowthCurves(m6, species = "Cod") + 
#   geom_line(data = subset(growth_warm, Species == "Cod"), aes(Age, value), color = "red", linetype = "dashed")

# Get growth data for all scenarios to be compared
# With projected temperature (until 2050)
growth_warm <- getGrowth(m6)
growth_warm$Scenario <- "Warm"

# With constant temperature (until 2050)
growth_con <- getGrowth(m5)
growth_con$Scenario <- "Constant"

# Only calibration scenario
# growth_cal <- getGrowth(m4)
# growth_cal$Scenario <- "Calibration"

# Combine
growth_df <- rbind(growth_con, growth_warm)

# Plot together
growth_df %>% filter(Age < 15) %>% 
  ggplot(., aes(Age, value, color = Scenario, linetype = Scenario)) + 
  geom_line(size = 1.5, alpha = 0.8) +
  labs(y = "Size (g)") +
  facet_wrap(~Species, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = c(brewer.pal(n = 8, name = "RdBu")[8], brewer.pal(n = 8, name = "RdBu")[1])) +
  scale_linetype_manual(values = c("solid", "twodash")) + 
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  NULL

# Test plotting % instead:
growth_warm <- getGrowth(m6)
growth_warm$Scenario <- "Warm"

# With constant temperature (until 2050)
growth_con <- getGrowth(m5)
growth_con$Scenario <- "Constant"

# Now with percentage instead
growth_df_wide <- data.frame(Age = growth_warm$Age,
                             value = growth_warm$value / growth_con$value,
                             Species = growth_warm$Species)

growth_df_wide %>% filter(Age < 15 & Age > 0) %>% 
  ggplot(., aes(Age, value)) + 
  geom_hline(yintercept = 1, size = 0.5, alpha = 0.8, color = "red") +
  geom_line(size = 1.5, alpha = 0.8, color = col[2]) +
  labs(y = "Relative Size (Warm/Constant)") +
  facet_wrap(~Species, scales = "free_y") +
  scale_y_continuous(limits = c(0.85, 1.05), expand = c(0, 0)) + 
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 3/4) +
  NULL


#**** Plot size spectra ============================================================
# Here I need to extract the data that goes into the the plotSpectra function.
# Copy and edit the function to only extract data. Do that for multiple models and plot

# 2. Changes in size-spectra relative to calibration period. Fig. 4 from Blanchard et al 2012, but lines don’t represent countries but species. Extract numbers at size. Plot difference. Check the @n for spectra, don't know which one for yield

# Plot difference in spectra between warming and constant. 

# Create getSpectra function..
getSpectra <- function(object){
  time_range  <- max(as.numeric(dimnames(object@n)$time))
  time_elements <- get_time_elements(object, time_range)
  n <- apply(object@n[time_elements, , ,drop = FALSE], c(2, 3), mean)
  species <- balticParams$species
  n <- n[as.character(dimnames(n)[[1]]) %in% species, , drop = FALSE]
  power <- 1
  n <- sweep(n, 2, params@w^power, "*")
  specDat <- data.frame(w = rep(as.numeric(dimnames(n)$w), length(species)),
                        n = c(as.numeric(n[1, ]), 
                              as.numeric(n[2, ]),
                              as.numeric(n[3, ])),
                        species = rep(species, each = 100))
  return(specDat)
}

# Testing it's correct
# plotSpectra(m6, plankton = F, benthos = F, algae = F) + 
#   geom_line(data = subset(getSpectra(m6), n > 0), aes(w, n, group = species), linetype = 2, color = "red") +
#   NULL


# Compare spectra in m5 and m6
warmSpec <- getSpectra(m6)
warmSpec$Scenario <- "Warming"

conSpec <- getSpectra(m5)
conSpec$Scenario <- "Constant"

plotSpec <- rbind(conSpec, warmSpec)

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


# F. VARY EFFORT GIVEN WARMING/NO WARMING ==========================================
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



# G. TEST EXTREME KAPPA AND LAMBA VALUES ===========================================
# Project with temperatue and effort varying through time
# This is equivalent to the m6-model with vbgecal temperature (centered to 10) and
# historical effort

params4_upd <- params3_upd

params4_upd@lambda <- 2.13
params4_upd@r_pp <- 200

# 1. How can we change r_pp when it's not in @params? Only hardwired changes?

# 2. Is r * (aw^b) the same as r*a*w^b ?

w <- seq(1, 10, 1)
a <- 0.1
b <- 0.7
r <- 4

r*a*w^b
r*(a*w^b)


params_test <- MizerParams(params4_upd@species_params,
                           kappa_ben = kappa_ben,
                           kappa = kappa,
                           w_bb_cutoff = w_bb_cutoff,
                           w_pp_cutoff = w_pp_cutoff,
                           r_pp = 4,
                           r_bb = r_bb)

params_test2 <- MizerParams(params4_upd@species_params,
                            kappa_ben = kappa_ben,
                            kappa = kappa,
                            w_bb_cutoff = w_bb_cutoff,
                            w_pp_cutoff = w_pp_cutoff,
                            r_pp = 999999,
                            r_bb = r_bb)

m7a <- project(params_test2, 
               dt = 0.1,
               effort = projectEffort_ct,
               temperature = projectTemp$temperature,
               diet_steps = 10,
               t_max = t_max,
               t_ref = 10) 

m7b <- project(params_test, 
               dt = 0.1,
               effort = projectEffort_ct,
               temperature = projectTemp$temperature,
               diet_steps = 10,
               t_max = t_max,
               t_ref = 10) 


plotBiomass(m7a)
plotBiomass(m7b)

tail(getYield(m7a))
tail(getYield(m7b))

str(m7@params)
str(m7@params@kappa)

m7@params@species_params
m7@params@kappa
m7@params@lambda
m7@params@rr_pp
m7@params@r_pp

# Here kappa is 6 and lambda is 2.13
# How much are kappa & lambda predicted to change according to Barnes?

# First delta-temp in time series
max(projectTemp$temperature) - min(projectTemp$temperature)

# How much does kappa and lambda change from default when delta T is 2.5?
# Lambda

# Lambda = -1.175 -0.002*T

df <- data.frame(temp  = c(10, 12),
                 m     = seq(1e-2, to = 1, length.out = 40),
                 log_a = c(9.8, 9.7),# 9.8-(2*0.045)
                 b     = c(-1.1, -1.104)) # -1.1-(0.002*2)


df$B <- df$log_a + log10(df$m)*df$b

ggplot(df, aes(m, B, color = factor(temp))) +
  geom_line() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  NULL


