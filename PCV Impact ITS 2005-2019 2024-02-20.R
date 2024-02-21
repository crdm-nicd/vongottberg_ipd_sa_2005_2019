#-----------------------------------#
#### PCV Impact on IPD 2005-2019 ####
#-----------------------------------#

#Created by: Jackie Kleynhans
#First created: 21 November 2023
#Purpose: Interrupted time series analysis to assess reduction in IPD post PCV introduction

#Set working directory
setwd("M:/Anne/PCV Impact 2005-2019/GitHub")

#Import libraries ----
library(readr)
library(dplyr)
library(tidyr)
library(janitor)
library(MASS)
library(lubridate)
library(ggplot2)
library(epiR)
library(writexl)
library(tibble)
library(cowplot)

#Clear workspace
rm(list=ls())

#Import data ----
an_imp_ipd <- read_csv("Annual imputed IPD counts 2024-02-20.csv")
qr_imp_ipd <- read_csv("Quarterly imputed IPD counts 2024-02-20.csv")

#Call ITS function
source('functions its 2024-02-20.R')

#Set seed
set.seed(1111)

#Model options (uncomment one)
model_to_run <- "Annual Poisson"
#model_to_run <- "Annual Negative Binomial"
#model_to_run <- "Quarterly Poisson"
#model_to_run <- "Quarterly Poisson Cos Season"
#model_to_run <- "Quarterly Negative Binomial Cos Season"
#model_to_run <- "Quarterly Negative Binomial"

if(grepl("Annual", model_to_run)){
  input_data <- an_imp_ipd
}
if(grepl("Quarterly", model_to_run)){
  input_data <- qr_imp_ipd
}
if(grepl("Poisson", model_to_run)){
  distribution <- 'pois'
}
if(grepl("Binomial", model_to_run)){
  distribution <- 'negbin'
}
if(grepl("Cos", model_to_run)){
  other_co_variates <- c('time_index', 'cos3')
}else{other_co_variates <- c('time_index')}

model_desc <- gsub(" ", "_", model_to_run)

#List of all models to run
dslist <- unique(input_data$desc)

#Run ITS ----
#NB: Model can run with either annual or quarterly counts. Ensure to update model type text and to update folder and file names.
#Model options negbin/pois
for(k in dslist){
  print(k)
  ds <- input_data %>% 
    mutate(ipd = round(as.numeric(ipd)),
           pop = round(as.numeric(pop))) %>% 
    filter(desc == k) %>% 
    dplyr::select(-desc) 
  colnames(ds) <- c("date", "ipd", "pop")
  ipd.mod <- step_func(ds=ds, outcome_name='ipd', denom='pop', mod=distribution, other.covars=other_co_variates)
  assign(paste0("ipd.mod.",k), ipd.mod)
  ipd.mod.rr <- cbind(ipd.mod$rr.annual, ipd.mod$aic, ipd.mod$deviance, ipd.mod$df, ipd.mod$overdispersion)
  assign(paste0("ipd.mod.rr",k), ipd.mod.rr)
}

#Compile and format results ----

#For all models: Generate reductions and plot results
#List of model results
model_outputs <- mget(paste0("ipd.mod.rr", dslist))

#Combine results, Calculate rates
rr.all <- bind_rows((model_outputs)) %>% 
  mutate(ipd_model_rate = (preds.q.50 / denom)*100000,
         ipd_model_rate_2.5 = (preds.q.2.5 / denom)*100000,
         ipd_model_rate_97.5 = (preds.q.97.5 / denom)*100000,
         ipd_expected_rate = (preds.cf.q.50 / denom)*100000,
         ipd_expected_rate_2.5 = (preds.cf.q.2.5 / denom)*100000,
         ipd_expected_rate_97.5 = (preds.cf.q.97.5 / denom)*100000)

rr.all <- rr.all %>% 
  group_by(ds) %>% 
  mutate(pre_period = if_else(`date...7`<2009, 1, 0),
         ipd_true_19 = max(if_else(`date...7`==2019, ipd_true, NA), na.rm = TRUE),
         denom_19 = max(if_else(`date...7`==2019, denom, NA), na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(ds, pre_period) %>% 
  mutate(ipd_true_pre = if_else(pre_period==1, sum(ipd_true), NA),
         denom_pre = if_else(pre_period==1, sum(denom), NA)) %>% 
  ungroup() %>% 
  group_by(ds) %>% 
  mutate(ipd_true_pre = max(ipd_true_pre, na.rm=TRUE),
         denom_pre = max(denom_pre, na.rm=TRUE)) %>% 
  ungroup() %>% 
  dplyr::select(-pre_period)

rr.all <- cbind(rr.all,
                cbind(as.data.frame(epi.conf(as.matrix(rr.all %>% dplyr::select(ipd_true, denom)), ctype = "inc.rate", method = "exact", N = 1000, design = 1, 
                                       conf.level = 0.95) * 100000) %>% 
                  rename("ipd_true_rate" = "est",
                         "ipd_true_rate_2.5" = "lower",
                         "ipd_true_rate_97.5" = "upper")),
                cbind(as.data.frame(epi.conf(as.matrix(rr.all %>% dplyr::select(ipd_true_pre, denom_pre)), ctype = "inc.rate", method = "exact", N = 1000, design = 1, 
                                       conf.level = 0.95) * 100000) %>% 
                  rename("ipd_true_pre_rate" = "est",
                         "ipd_true_pre_rate_2.5" = "lower",
                         "ipd_true_pre_rate_97.5" = "upper"),
                  as.data.frame(epi.conf(as.matrix(rr.all %>% dplyr::select(ipd_true_19, denom_19)), ctype = "inc.rate", method = "exact", N = 1000, design = 1, 
                                        conf.level = 0.95) * 100000) %>% 
                    rename("ipd_true_19_rate" = "est",
                           "ipd_true_19_rate_2.5" = "lower",
                           "ipd_true_19_rate_97.5" = "upper")))

#IRRs for pre vs post
#Duplicate case numbers across all years
irr_calc <- rr.all  %>% 
  group_by(ds) %>% 
  summarise(ipd_true_19 = max(ipd_true_19, na.rm = TRUE),
            denom_19 = max(denom_19, na.rm = TRUE),
            ipd_true_pre = max(ipd_true_pre, na.rm = TRUE),
            denom_pre = max(denom_pre, na.rm = TRUE))

n_rows <- nrow(irr_calc)
res_out <- data.frame(IRR = numeric(n_rows), IRR_2.5 = numeric(n_rows), IRR_97.5 = numeric(n_rows))

for (i in 1:n_rows) {
  dat <- c(irr_calc$ipd_true_19[i], irr_calc$denom_19[i], irr_calc$ipd_true_pre[i], irr_calc$denom_pre[i])
  res <- epi.2by2(dat, method = "cohort.time", digits = 2,
                  conf.level = 0.95, units = 100000, interpret = FALSE, outcome = "as.columns")
  
  # Assign the IRR and its confidence intervals to the res_out dataframe
  res_out[i, "2005_2008 cases"] <- as.numeric(res$tab$`   Outcome +`[2])
  res_out[i, "2005_2008 population size"] <- as.numeric(res$tab$`   Time at risk`[2])
  res_out[i, "Risk ratio 2005_2008 vs 2019 (observed)"] <- res$massoc.summary$est[1]
  res_out[i, "Risk ratio 2005_2008 vs 2019 (observed) Lower 95%CI"] <- res$massoc.summary$lower[1]
  res_out[i, "Risk ratio 2005_2008 vs 2019 (observed) Upper 95%CI"] <- res$massoc.summary$upper[1]
  res_out[i, "Absolute risk difference 2005_2008 vs 2019 (observed)"] <- res$massoc.summary$est[2]
  res_out[i, "Absolute risk difference 2005_2008 vs 2019 (observed) Lower 95%CI"] <- res$massoc.summary$lower[2]
  res_out[i, "Absolute risk difference 2005_2008 vs 2019 (observed) Upper 95%CI"] <- res$massoc.summary$upper[2]
  
}
res_out <- cbind(irr_calc %>% dplyr::select(ds), res_out)

rr.all <- left_join(rr.all, res_out, by="ds")

#IRRs for predicted vs modelled
irr_calc <- rr.all  %>% 
  filter(`date...7`==2019) %>% 
  dplyr::select(ds, denom, preds.q.50, preds.cf.q.50)

res_out <- data.frame(IRR = numeric(n_rows), IRR_2.5 = numeric(n_rows), IRR_97.5 = numeric(n_rows))

for (i in 1:n_rows) {
  dat <- c(irr_calc$preds.q.50[i], irr_calc$denom[i], irr_calc$preds.cf.q.50[i], irr_calc$denom[i])
  res <- epi.2by2(dat, method = "cohort.time", digits = 2,
                  conf.level = 0.95, units = 100000, interpret = FALSE, outcome = "as.columns")
  
  # Assign the IRR and its confidence intervals to the res_out dataframe
  res_out[i, "Risk ratio expected vs modeled"] <- res$massoc.summary$est[1]
  res_out[i, "Risk ratio expected vs modeled Lower 95%CI"] <- res$massoc.summary$lower[1]
  res_out[i, "Risk ratio expected vs modeled Upper 95%CI"] <- res$massoc.summary$upper[1]
  res_out[i, "Absolute risk difference expected vs modeled"] <- res$massoc.summary$est[2]
  res_out[i, "Absolute risk difference expected vs modeled Lower 95%CI"] <- res$massoc.summary$upper[2]
  res_out[i, "Absolute risk difference expected vs modeled Upper 95%CI"] <- res$massoc.summary$lower[2]
  
}
res_out <- cbind(irr_calc %>% dplyr::select(ds), res_out)

rr.all <- left_join(rr.all, res_out, by="ds")

#Tidy up results
rr.all <- rr.all %>% 
  rename("Risk ratio (ITS)" = "rr.q.50",
         "Risk ratio (ITS) Lower 95%CI" = "rr.q.2.5",
         "Risk ratio (ITS) Upper 95%CI" = "rr.q.97.5",
         "Subgroup" = "ds",
         "Date" = "date...7",
         "Modeled cases (ITS)" = "preds.q.50",
         "Modeled cases (ITS) Lower 95%CI" = "preds.q.2.5",
         "Modeled cases (ITS) Upper 95%CI" = "preds.q.97.5",
         "Expected cases (ITS)" = "preds.cf.q.50",
         "Expected cases (ITS) Lower 95%CI" = "preds.cf.q.2.5",
         "Expected cases (ITS) Upper 95%CI" = "preds.cf.q.97.5",
         "Population size" = "denom",
         "Observed cases (Imputed)" = "ipd_true",
         "ITS model AIC" = "ipd.mod$aic",
         "ITS model deviance" = "ipd.mod$deviance",
         "ITS model overdispersion" = "ipd.mod$overdispersion",
         "Modeled incidence (ITS)" = "ipd_model_rate",
         "Modeled incidence (ITS) Lower 95%CI" = "ipd_model_rate_2.5",
         "Modeled incidence (ITS) Upper 95%CI" = "ipd_model_rate_97.5",
         "Expected incidence (ITS)" = "ipd_expected_rate",
         "Expected incidence (ITS) Lower 95%CI" = "ipd_expected_rate_2.5",
         "Expected incidence (ITS) Upper 95%CI" = "ipd_expected_rate_97.5",
         "Observed incidence" = "ipd_true_rate",
         "Observed incidence Lower 95%CI" = "ipd_true_rate_2.5",
         "Observed incidence Upper 95%CI" = "ipd_true_rate_97.5",
         "2005_2008 incidence" = "ipd_true_pre_rate",
         "2005_2008 incidence Lower 95%CI" = "ipd_true_pre_rate_2.5",
         "2005_2008 incidence Upper 95%CI" = "ipd_true_pre_rate_97.5",) %>% 
  mutate("Percentage risk difference 2005_2008 vs 2019 (observed)" = ((1-`Risk ratio 2005_2008 vs 2019 (observed)`)*100)*-1,
         "Percentage risk difference 2005_2008 vs 2019 (observed) Lower 95%CI" =((1-`Risk ratio 2005_2008 vs 2019 (observed) Lower 95%CI`)*100)*-1,
         "Percentage risk difference 2005_2008 vs 2019 (observed) Upper 95%CI" = ((1-`Risk ratio 2005_2008 vs 2019 (observed) Upper 95%CI`)*100)*-1,
         "Percentage risk difference expected vs modeled" = ((1-`Risk ratio expected vs modeled`)*100)*-1,
         "Percentage risk difference expected vs modeled Lower 95%CI" = ((1-`Risk ratio expected vs modeled Lower 95%CI`)*100)*-1,
         "Percentage risk difference expected vs modeled Upper 95%CI" = ((1-`Risk ratio expected vs modeled Upper 95%CI`)*100)*-1) %>% 
  dplyr::select("Date", "Subgroup", "Population size", 
         "Observed cases (Imputed)", 
         "Observed incidence", "Observed incidence Lower 95%CI", "Observed incidence Upper 95%CI",
         "2005_2008 cases", "2005_2008 population size",
         "2005_2008 incidence", "2005_2008 incidence Lower 95%CI",  "2005_2008 incidence Lower 95%CI",
         "Risk ratio 2005_2008 vs 2019 (observed)", "Risk ratio 2005_2008 vs 2019 (observed) Lower 95%CI", "Risk ratio 2005_2008 vs 2019 (observed) Upper 95%CI",
         "Absolute risk difference 2005_2008 vs 2019 (observed)", 
         "Absolute risk difference 2005_2008 vs 2019 (observed) Lower 95%CI", 
         "Absolute risk difference 2005_2008 vs 2019 (observed) Upper 95%CI",
         "Percentage risk difference 2005_2008 vs 2019 (observed)", 
         "Percentage risk difference 2005_2008 vs 2019 (observed) Lower 95%CI",
         "Percentage risk difference 2005_2008 vs 2019 (observed) Upper 95%CI",
         "Modeled cases (ITS)", "Modeled cases (ITS) Lower 95%CI", "Modeled cases (ITS) Upper 95%CI" , 
         "Modeled incidence (ITS)", "Modeled incidence (ITS) Lower 95%CI", "Modeled incidence (ITS) Upper 95%CI", 
         "Expected cases (ITS)" , "Expected cases (ITS) Lower 95%CI", "Expected cases (ITS) Upper 95%CI", 
         "Expected incidence (ITS)", "Expected incidence (ITS) Lower 95%CI" , "Expected incidence (ITS) Upper 95%CI",
         "Risk ratio expected vs modeled", "Risk ratio expected vs modeled Lower 95%CI", "Risk ratio expected vs modeled Upper 95%CI",
         "Absolute risk difference expected vs modeled", 
         "Absolute risk difference expected vs modeled Lower 95%CI", 
         "Absolute risk difference expected vs modeled Upper 95%CI",
         "Percentage risk difference expected vs modeled", 
         "Percentage risk difference expected vs modeled Lower 95%CI",
         "Percentage risk difference expected vs modeled Upper 95%CI",
         "Risk ratio expected vs modeled", 
         "Risk ratio expected vs modeled Lower 95%CI", "Risk ratio expected vs modeled Upper 95%CI", 
         "Absolute risk difference expected vs modeled", 
         "Absolute risk difference expected vs modeled Lower 95%CI", "Absolute risk difference expected vs modeled Upper 95%CI",
         "ITS model AIC", "ITS model deviance", "ITS model overdispersion") %>% 
  mutate("Model type" = model_to_run)

#Format category names
rr.all <- rr.all %>% 
  mutate(Category = case_when(Subgroup=="All" ~ "All",
                              Subgroup=="PCV7" ~ "PCV7",
                              Subgroup=="PCV13add" ~ "PCV13add", 
                              Subgroup=="nonPCV" ~ "nonPCV",
                              Subgroup=="<10w" ~ "< 10 weeks",
                              Subgroup=="<2y" ~ "< 2 years",
                              Subgroup=="2-4y" ~ "2-4 years",
                              Subgroup=="5-14y" ~ "5-14 years",
                              Subgroup=="15-24y" ~ "15-24 years",
                              Subgroup=="25-44y" ~ "25-44 years",
                              Subgroup=="45-64y" ~ "45-64 years",
                              Subgroup==">=65y" ~ "\U2265 65 years",
                              Subgroup=="<10w_PCV13add" ~ "< 10 weeks, PCV13add",
                              Subgroup=="<10w_PCV7" ~ "< 10 weeks, PCV7",
                              Subgroup=="<10w_nonPCV" ~ "< 10 weeks, non-PCV",
                              Subgroup=="15-24y_PCV13add" ~ "15-24 years, PCV13add",
                              Subgroup=="15-24y_PCV7" ~ "15-24 years, PCV7",
                              Subgroup=="15-24y_nonPCV" ~ "15-24 years, nonPCV",
                              Subgroup=="2-4y_PCV13add" ~ "2-4 years, PCV13add",
                              Subgroup=="2-4y_PCV7" ~ "2-4 years, PCV7",
                              Subgroup=="2-4y_nonPCV" ~ "2-4 years, nonPCV",
                              Subgroup=="25-44y_PCV13add" ~ "25-44 years, PCV13add",
                              Subgroup=="25-44y_PCV7" ~ "25-44 years, PCV7",
                              Subgroup=="25-44y_nonPCV" ~ "25-44 years, nonPCV",
                              Subgroup=="45-64y_PCV13add" ~ "45-64 years, PCV13add",
                              Subgroup=="45-64y_PCV7" ~ "45-64 years, PCV7",
                              Subgroup=="45-64y_nonPCV" ~ "45-64 years, nonPCV",
                              Subgroup=="5-14y_PCV13add" ~ "5-14 years, PCV13add",
                              Subgroup=="5-14y_PCV7" ~ "5-14 years, PCV7",
                              Subgroup=="5-14y_nonPCV" ~ "5-14 years, nonPCV",
                              Subgroup=="<2y_PCV13add" ~ "< 2 years, PCV13add",
                              Subgroup=="<2y_PCV7" ~ "< 2 years, PCV7",
                              Subgroup=="<2y_nonPCV" ~ "< 2 years, nonPCV",
                              Subgroup==">=65y_PCV13add" ~ "\U2265 65 years, PCV13add",
                              Subgroup==">=65y_PCV7" ~ "\U2265 65 years, PCV7",
                              Subgroup==">=65y_nonPCV" ~ "\U2265 65 years, nonPCV",
                              Subgroup=="pen_res" ~ "Penicillin-nonsusceptible",
                              Subgroup=="pen_sus" ~ "Penicillin-susceptible",
                              Subgroup=="multi_res" ~ "Multidrug-resistant",
                              Subgroup=="multi_sus"  ~ "Multidrug-susceptible")) %>% 
  mutate(Category = factor(Category, levels = c("All", "PCV7", "PCV13add", "nonPCV", "< 10 weeks", "< 10 weeks, PCV7", "< 10 weeks, PCV13add", "< 10 weeks, non-PCV", "< 2 years", "< 2 years, PCV7", "< 2 years, PCV13add", "< 2 years, nonPCV", "2-4 years", "2-4 years, PCV7", "2-4 years, PCV13add", "2-4 years, nonPCV", "5-14 years", "5-14 years, PCV7", "5-14 years, PCV13add", "5-14 years, nonPCV", "15-24 years", "15-24 years, PCV7", "15-24 years, PCV13add", "15-24 years, nonPCV", "25-44 years", "25-44 years, PCV7", "25-44 years, PCV13add", "25-44 years, nonPCV", "45-64 years", "45-64 years, PCV7", "45-64 years, PCV13add", "45-64 years, nonPCV", "\U2265 65 years", "\U2265 65 years, PCV7", "\U2265 65 years, PCV13add", "\U2265 65 years, nonPCV", "Penicillin-susceptible", "Penicillin-nonsusceptible", "Multidrug-susceptible", "Multidrug-resistant"))) %>% 
  mutate(fig_group = case_when(Category=="All" | Category=="PCV7" | Category=="PCV13add" | Category=="nonPCV" ~ "Overall",
                               grepl("<", Category) | grepl("2-4", Category) | grepl("5-14", Category) ~ "Children",
                               grepl("15-24", Category) | grepl("25-44", Category) | grepl("45-64", Category) | grepl("65", Category) ~ "Adults",
                               grepl("susce", Category) | grepl("res", Category) ~ "Abx"))

#Sort
rr.all %>% 
  arrange(Category, Date)

#Plot ITS Model output ----
plot_its <- function(grouping){
rr.all %>%
  filter(fig_group==grouping) %>% 
ggplot(aes(x = Date, y = `Observed incidence`)) +
  geom_point(aes(color = 'Observed incidence'), pch = 16, size = 1) +
  labs(y = "IPD cases/100,000 population", x = "Collection year") +
  geom_line(aes(x = Date, y = `Expected incidence (ITS)`, color = 'Expected incidence (no intervention)'), lty = 2, size = 0.5) +
  geom_line(aes(x = Date, y = `Modeled incidence (ITS)`, color = 'Modeled incidence (actual)'), lty = 1, size = 0.5) +
  geom_line(aes(x = Date, y = 0, color = 'PCV introduction'), lty = 4, size = 0.5) +
  geom_ribbon(aes(
    x = Date,
    ymin = `Modeled incidence (ITS) Lower 95%CI`,
    ymax = `Modeled incidence (ITS) Upper 95%CI`
  ), fill = 'black', alpha = 0.3) +
  geom_ribbon(aes(
    x = Date,
    ymin = ifelse(Date >= 2008, `Expected incidence (ITS) Lower 95%CI`, NA),
    ymax = ifelse(Date >= 2008, `Expected incidence (ITS) Upper 95%CI`, NA)
  ), fill = 'red', alpha = 0.3) +
  geom_vline(xintercept = 2009, col = 'gray', lty = 4, size = 0.8) +
  geom_vline(xintercept = 2011, col = 'gray', lty = 4, size = 0.8) +
  scale_y_continuous(expand = c(0, 0), breaks = my_breaks, limits = my_limits) +
  theme_classic() +
  theme(
    axis.line.x.bottom = element_line(size = 0.5),
    axis.line.y.left = element_line(size = 0.5),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x =element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.text = element_text(size = 8), 
    legend.position = "bottom",
    legend.margin = margin(-10, -5, -5, -5),
    legend.spacing.y= unit(-0.3, "cm")) +
  scale_color_manual(
    values = c('Expected incidence (no intervention)' = 'red', 
               'Modeled incidence (actual)' = 'black', 
               'Observed incidence' = '#606060',
               'PCV introduction' = 'gray'),
    guide = guide_legend(title = NULL,
                         override.aes = 
                           list(linetype = c(2, 1, NA, 4), 
                                pch = c(NA, NA, 16, NA), 
                                color = c('red', 'black', '#606060', 'gray')),
                         nrow = 1,
                         byrow = TRUE)) +
  #facet_wrap(~ Category, scales = "free_y") +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 8, color = "black"))+
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5)  
}

#Save plots
plot <- plot_its("Overall") +
  facet_wrap(~Category, scales = "free_y", nrow = 1, ncol = 4)
assign(paste0("plot_",model_desc,"_Overall"), plot)
ggsave(paste0("ITS results and plots/",model_desc,"/",model_desc,"_Overall.tiff"), plot=get(paste0("plot_",model_desc,"_Overall")), device="tiff", dpi=300, width = 7, height = 2.5)

plot <- plot_its("Children") +
  facet_wrap(~Category, scales = "free_y", nrow = 4, ncol = 4)
assign(paste0("plot_",model_desc,"_Children"), plot)
ggsave(paste0("ITS results and plots/",model_desc,"/",model_desc,"_Children.tiff"), plot=get(paste0("plot_",model_desc,"_Children")), device="tiff", dpi=300, width = 7, height = 7.5)

plot <- plot_its("Adults") +
  facet_wrap(~Category, scales = "free_y", nrow = 4, ncol = 4)
assign(paste0("plot_",model_desc,"_Adults"), plot)
ggsave(paste0("ITS results and plots/",model_desc,"/",model_desc,"_Adults.tiff"), plot=get(paste0("plot_",model_desc,"_Adults")), device="tiff", dpi=300, width = 7, height = 7.5)

plot <- plot_its("Abx") +
  facet_wrap(~Category, scales = "free_y", nrow = 2, ncol = 4)
assign(paste0("plot_",model_desc,"_Abx"), plot)
ggsave(paste0("ITS results and plots/",model_desc,"/",model_desc,"_Abx.tiff"), plot=get(paste0("plot_",model_desc,"_Abx")), device="tiff", dpi=300, width = 7, height = 2.5)


#Save output dataframe ----
rr.all <- rr.all %>%
  dplyr::select(`Model type`, everything())
assign(paste0("rr.all_",model_desc), rr.all)
write_xlsx(get(paste0("rr.all_",model_desc)), paste0("ITS results and plots/",model_desc,"/",model_desc,"_All results_",format(Sys.Date(), "%Y-%m-%d"),".xlsx"))

#Format main results  ----
main_res <- rr.all %>% 
  filter(Date==2019) %>% 
  mutate(`Modeled incidence` = paste0(round(`Modeled incidence (ITS)`, 1)," (", 
                                             round(`Modeled incidence (ITS) Lower 95%CI`, 1), " to ", 
                                             round(`Modeled incidence (ITS) Upper 95%CI`, 1), ")")) %>% 
  mutate(`Expected incidence`= paste0(round(`Expected incidence (ITS)`, 1)," (", 
                                       round(`Expected incidence (ITS) Lower 95%CI`, 1), " to ", 
                                       round(`Expected incidence (ITS) Upper 95%CI`, 1), ")")) %>% 
  mutate(`Absolute risk difference`= paste0(round(`Absolute risk difference expected vs modeled`, 1)," (", 
                                      round(`Absolute risk difference expected vs modeled Lower 95%CI`, 1), " to ", 
                                      round(`Absolute risk difference expected vs modeled Upper 95%CI`, 1), ")")) %>% 
  mutate(`Percentage risk difference`= paste0(round(`Percentage risk difference expected vs modeled`, 1)," (", 
                                            round(`Percentage risk difference expected vs modeled Lower 95%CI`, 1), " to ", 
                                            round(`Percentage risk difference expected vs modeled Upper 95%CI`, 1), ")")) %>% 
  dplyr::select(Category,`Modeled incidence`, `Expected incidence`, `Absolute risk difference`,`Percentage risk difference` )
write_xlsx(main_res, paste0("ITS results and plots/",model_desc,"/",model_desc,"_Main results_",format(Sys.Date(), "%Y-%m-%d"),".xlsx"))

#Comparison of AIC scores
#aic_scores <- bind_rows(filtered_data)  %>%
  #pivot_wider(id_cols = Subgroup, names_from = 'Model type', values_from = c(`ITS model AIC`, `ITS model overdispersion`, `ITS model deviance`))

#Export AIC deviance comparison
#write_xlsx(aic_scores , paste0("ITS results and plots/PCV Impact IPD 2005-2019 ITS Model Comparison ",format(Sys.Date(), "%Y-%m-%d"),".xlsx"))

#