#-----------------------------------#
#### PCV Impact on IPD 2005-2019 ####
#-----------------------------------#

#Created by: Jackie Kleynhans
#First created: 21 November 2023
#Purpose: Create incidence plots for IPD data
#Follow after PCV Impact PCV script

#Overall incidence plots ----
  
  ipd_pcv <- an_imp_ipd %>% 
    mutate(inc = (ipd/pop)*100000) %>% 
    pivot_wider(id_cols = "date", names_from = "desc", values_from = "inc") %>% 
    ggplot() +
    geom_line(aes(x = date, y = `PCV7`, color = 'PCV7'), lty = 1, size = 0.5) +
    geom_line(aes(x = date, y = `PCV13add`, color = 'PCV13add'), lty = 1, size = 0.5) +
    geom_line(aes(x = date, y = `nonPCV`, color = 'nonPCV'), lty = 1, size = 0.5) +
    geom_line(aes(x = date, y = `All`, color = 'All'), lty = 2, size = 0.5) +
    geom_line(aes(x = date, y = 0, color = 'PCV introduction'), lty = 4, size = 0.5) +
    ylab("IPD cases/100,000 population") +
    xlab("Collection year") +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_date(expand = c(0, 0), date_breaks = "3 years",  date_labels = "%Y") +
    geom_vline(xintercept = as.Date(c("2009-01-01")), col = 'gray', lty = 4, size = 0.8) +
    geom_vline(xintercept = as.Date(c("2011-01-01")), col = 'gray', lty = 4, size = 0.8) +
    theme(axis.line.x.bottom = element_line(size = 0.5), axis.line.y.left = element_line(size = 0.5)) + 
    theme(axis.text.x = element_text(size = 6), 
          axis.text.y = element_text(size = 6),
          axis.title.x = element_text(size = 8), 
          axis.title.y = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.position = "bottom",
          legend.margin = margin(-10, -5, -5, -5),
          legend.spacing.y= unit(-0.3, "cm")) +
    scale_color_manual(values = c("All" = "#000066", "PCV7" = '#0099cc', "PCV13add" = '#b3ecff', "nonPCV" = '#00cc66', 'PCV introduction' = '#606060'),
                       breaks = c("All", "PCV7", "PCV13add", "nonPCV", 'PCV introduction'),
                       guide = guide_legend(title = NULL,
                                            override.aes = 
                                              list(linetype = c(2, 1, 1, 1, 4), 
                                                   color = c("#000066", '#0099cc', '#b3ecff', '#00cc66', '#606060')),
                                            nrow = 2,
                                            byrow = TRUE)) +
    labs(color = NULL)
  ggsave("ITS results and plots/IPD Incidence Serotype.tiff", plot=ipd_pcv, device="tiff", dpi=300, width = 3.5, height = 2.8)
  
  
  ipd_age_all <- an_imp_ipd %>% 
    mutate(inc = (ipd/pop) * 100000) %>% 
    pivot_wider(id_cols = "date", names_from = "desc", values_from = "inc") %>% 
    ggplot(aes(x = date)) +
    geom_line(aes(y = `<2y`, color = "< 2 years"), lty = 1, size = 0.5) +
    geom_line(aes(y = `2-4y`, color = "2-4 years"), lty = 1, size = 0.5) +
    geom_line(aes(y = `5-14y`, color = "5-14 years"), lty = 1, size = 0.5) +
    geom_line(aes(y = `15-24y`, color = "15-24 years"), lty = 1, size = 0.5) +
    geom_line(aes(y = `25-44y`, color = "25-44 years"), lty = 1, size = 0.5) +
    geom_line(aes(y = `45-64y`, color = "45-64 years"), lty = 1, size = 0.5) +
    geom_line(aes(y = `>=65y`, color = "≥ 65 years"), lty = 1, size = 0.5) +
    geom_line(aes(y = 0, color = 'PCV introduction'), lty = 4, size = 0.5) +
    ylab("IPD cases/100,000 population") +
    xlab("Collection year") +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_date(expand = c(0, 0), date_breaks = "3 years", date_labels = "%Y") +
    geom_vline(xintercept = as.Date(c("2009-04-01")), col = '#606060', lty = 4, size = 0.8) +
    geom_vline(xintercept = as.Date(c("2011-01-01")), col = '#606060', lty = 4, size = 0.8) +
    theme(
      axis.line.x.bottom = element_line(size = 0.5),
      axis.line.y.left = element_line(size = 0.5),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.position = "bottom",
      legend.margin = margin(-10, -5, -5, -5),
      legend.spacing.y= unit(-0.3, "cm")) +
    scale_color_manual(
      values = c('< 2 years' = "#000066",'2-4 years' = '#0099cc', '5-14 years' = '#b3ecff', '15-24 years' = '#62deba',
                 '25-44 years' = '#99ff99', '45-64 years' = '#00cc66', '≥ 65 years' = '#008738', 'PCV introduction' = '#606060'),
      labels = c('< 2 years', '2-4 years', '5-14 years', '15-24 years', 
                 '25-44 years', '45-64 years', '\u2265 65 years', 'PCV introduction'),
      breaks = c('< 2 years', '2-4 years', '5-14 years', '15-24 years', 
                 '25-44 years', '45-64 years', '≥ 65 years', 'PCV introduction'),
      guide = guide_legend(title = NULL,
                           override.aes = 
                             list(linetype = c(1, 1, 1, 1, 1, 1, 1, 4), 
                                  color = c("#000066", '#0099cc', '#b3ecff', '#62deba', '#99ff99', 
                                            '#00cc66',  '#008738', '#606060')),
                           nrow= 2,
                           byrow = TRUE)) +
    labs(color = NULL)
  
  ipd_age_m2 <- an_imp_ipd %>% 
    mutate(inc = (ipd/pop)*100000) %>% 
    pivot_wider(id_cols = "date", names_from = "desc", values_from = "inc") %>% 
    ggplot() +
    geom_line(aes(x=date, y=`2-4y`, col='2-4y'), lty=1, size=0.5) +
    geom_line(aes(x=date, y=`5-14y`, col='5-14y'), lty=1, size=0.5) +
    geom_line(aes(x=date, y=`15-24y`, col='15-24y'), lty=1, size=0.5) +
    geom_line(aes(x=date, y=`25-44y`, col='25-44y'), lty=1, size=0.5) +
    geom_line(aes(x=date, y=`45-64y`, col='45-64y'), lty=1, size=0.5) +
    geom_line(aes(x=date, y=`>=65y`, col='>=65y'), lty=1, size=0.5) +
    geom_line(aes(x = date, y = 0, color = 'PCV introduction'), lty = 4, size = 0.5) +
    #ylab("IPD cases/100,000 population") +
    #xlab("Collection year") +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_date(expand = c(0, 0), date_breaks = "3 years",  date_labels = "%Y") +
    geom_vline(xintercept=as.Date(c( "2009-04-01") ), col='#606060', lty=4, size=0.8) +
    geom_vline(xintercept=as.Date(c( "2011-01-01") ), col='#606060', lty=4, size=0.8) +
    theme(axis.line.x.bottom=element_line(size=1.2), axis.line.y.left=element_line(size=1.2)) + 
    theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.line.x.bottom = element_line(size = 0.5),
          axis.line.y.left = element_line(size = 0.5)) +
    scale_color_manual(
      name = "Age",
      values = c('<2y' = "#000066",'2-4y' = '#0099cc',  '5-14y' = '#b3ecff',  '15-24y' = '#62deba' , '25-44y' = '#99ff99', '45-64y' = '#00cc66', '>=65y' = '#008738'),
      labels = c( `<2y` = "< 2 years", `2-4y` ="2-4 years", `5-14y` ="5-14 years", `15-24y` ="15-24 years", 
                  `25-44y` ="25-44 years", `45-64y` ="45-64 years",`>=65y` = "\u2265 65 years"),
      breaks = c( "< 2 years", "2-4 years", "5-14 years", "15-24 years", 
                  "25-44 years", "45-64 years", "\u2265 65 years"))
  
  ipd_age_comb <-
    ggdraw() +
    draw_plot(ipd_age_all) +
    draw_plot(ipd_age_m2, x = 0.5, y = .55, width = .5, height = .45)
  ggsave("ITS results and plots/IPD Incidence Age with insert.tiff", plot=ipd_age_comb, device="tiff", dpi=300, width = 7, height = 2.8)
  ggsave("ITS results and plots/IPD Incidence Age.tiff", plot=ipd_age_all, device="tiff", dpi=300, width = 7, height = 2.8)
  
#Imputation trends over time
an_imp_act_ipd <- read_csv("Annual imputed and actual IPD counts 2024-04-26.csv")

an_imp_act_ipd <- an_imp_act_ipd %>% 
  mutate(proportion_imputed = (((ipd_imputed-ipd_actual)/ipd_imputed)*100)) %>% 
  filter(desc=="All" | desc=="<10w" | desc=="PCV7" | 
           desc=="<10w_nonPCV" | desc=="<2y_nonPCV" | desc=="2-4y_nonPCV" | desc=="5-14y_nonPCV" | desc=="15-24y_nonPCV" | desc=="25-44y_nonPCV" | desc=="45-64y_nonPCV" | desc==">=65y_nonPCV" | desc=="pen_res") %>% 
  mutate(Description = factor(desc, 
                              levels = c("All", "<10w", "PCV7", "pen_res", "<10w_nonPCV", "<2y_nonPCV", "2-4y_nonPCV", 
                                         "5-14y_nonPCV", "15-24y_nonPCV", "25-44y_nonPCV", "45-64y_nonPCV", 
                                         ">=65y_nonPCV"), 
                              labels = c("All ages and serotypes", "Age", "Serotype", 
                                         "Antibiotic susceptability",
                                         "Serotypes in < 10 weeks", 
                                         "Serotypes in < 2 years", "Serotypes in  2-4 years", 
                                         "Serotypes in 5-14 years", "Serotypes in 15-24 years", 
                                         "Serotypes in 25-44 years", "Serotypes in 45-64 years", 
                                         "Serotypes in \u2265 65 years")))


imputation_plot <- ggplot(an_imp_act_ipd) +
  geom_line(aes(x=date, y=proportion_imputed)) +
  facet_wrap(~Description) +
    ylab("Percentage of cases imputed") +
    xlab("Collection year") +
    theme_classic() +
    theme(axis.line.x.bottom=element_line(size=1.2), axis.line.y.left=element_line(size=1.2)) + 
    theme(
      axis.line.x.bottom = element_line(size = 0.5),
      axis.line.y.left = element_line(size = 0.5),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10)) +
    theme(strip.background = element_blank(), 
          strip.text = element_text(size = 8, color = "black"))
ggsave("ITS results and plots/Percentage cases imputed.tiff", plot=imputation_plot, device="tiff", dpi=300, width = 7, height = 7.5)


  