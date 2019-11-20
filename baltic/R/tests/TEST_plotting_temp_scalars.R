# C. PLOT TEMP SCALARS =============================================================



script <- getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/R/functions/tempFun.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

# Read in temperature data
temp_datRCP8.5 <- read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/mizer-rewiring/rewire-temp/baltic/data/Climate/Test_RCP8.5_from_graph.csv"), sep = ";", stringsAsFactors = FALSE)

# Calculate average temperature by year
tempDat <- data.frame(temp_datRCP8.5 %>% 
                        dplyr::filter(year < 2050.5) %>% 
                        dplyr::mutate(Year.r = factor(round(year, digits = 0))) %>% 
                        dplyr::group_by(Year.r) %>% 
                        dplyr::summarize(mean_temp = mean(rel_temp_70.99Ave)) %>% 
                        dplyr::mutate(Year.num = as.numeric(as.character(Year.r))))

t_ref <- 9.57

tempDat$mean_temp_scaled <- tempDat$mean_temp + t_ref

#temp <- seq(10, 12, 0.1)
temp <- seq(from = min(tempDat$mean_temp_scaled),
            to = max(tempDat$mean_temp_scaled),
            by = 0.5)

dat <- data.frame(expand.grid(temp = temp,
                              ea_met = ea$met))

dat$ea_mor <- ea$mor
dat$ea_int <- ea$int
dat$ea_gro <- ea$gro
dat$ea_car <- ea$car

data_list <- list()

for(i in 1:5) {
  
  scal <- data.frame(scal = as.numeric(tempFun(temperature = dat$temp, 
                                               t_ref = 10, 
                                               Ea = dat[, i + 1], 
                                               c_a = 0, 
                                               w = 10)),
                     temp = dat$temp,
                     rate = factor(rep(colnames(dat)[i + 1], length(dat$temp))),
                     id = i)
  
  data_list[[i]] <- scal
  
}

str(data_list)
data_list[[1]]

big_dat <- data_list %>% map_df(~as.data.frame(., stringsAsFactors = FALSE))

unique(big_dat$rate)

head(big_dat, 20)

# Plot
col <- rev(RColorBrewer::brewer.pal("Dark2", n = 5))

big_dat$rate <- factor(big_dat$rate, levels = c("Maximum\nconsumption rate", "Metabolic rate", "Background\nmortality rate", 
                                                "Resource growth rate", "Resource\ncarrying capacity"))

ggplot(big_dat, aes(x = temp, scal, group = temp)) + 
  geom_line(color = "black", linetype = "dashed") +
  facet_wrap(~rate, scales = "free") +
  theme_classic() +
  guides(fill = FALSE) +
  labs(x = expression(paste("Temperature [", degree*C, "]")),
       y = "Rate scalars") +
  NULL

#ggsave("baltic/figures/supp/random_rate_scalar.pdf", plot = last_plot(), width = 19, height = 19, units = "cm")