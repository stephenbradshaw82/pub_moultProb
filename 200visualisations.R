#' Author: Stephen Bradshaw
#' Date: 22 Jan 2025
#' Title: Visualisation / Plotting for Moult Prob work
#' Details:
#'      - 20250123 Initial creation for public repo
#'      

#### LIBRARIES ####
rm(list=ls())
options(stringsAsFactors = FALSE)
options(scipen = 999)

#' Packages (not in environment) -------------------------------------------
list.of.packages <- c("magrittr", "tidyr", "plyr", "dplyr", "ggplot2"
                      , "stringr", "purrr", "rebus","ggpubr"
                      , "foreach", "doParallel"
                      , "scales"
                      , "gridExtra"
                      , "posterior"
                      , "rstan", "stringr", "StanHeaders"
                      , "bayesplot"
                      , 'gridtext', "grid"
                      , "scales"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#' Libraries ---------------------------------------------------------------
req_packages <- list.of.packages 
sapply(req_packages,require,character.only = TRUE, quietly=TRUE)
#####

#### FUNCTIONS ####
source("src/functions.R")
#####

#################################################
#### GENERATED DATA PLOT (SINGLE BASE MODEL) ####
#################################################
#' Original gen plot (single) as base model
#--> Get ggs plots for damage model ####
useDir <- "outputs/"

#--> SAVED OUTPUTS ####
## Specify the folder name
folder_name <- paste0(useDir,"plots")

## Check if the folder exists, and create it if it doesn't
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
  cat("Folder 'outputs' created.\n")
} else {
  cat("Folder 'outputs' already exists.\n")
}

dir(useDir) 

use_height <- 480
use_width <- 480

useTask <- c("gen30lib18grow52")

growth_gen <- readRDS(paste0(useDir, "out_moult_", useTask, "_FITMCMC.rds")) ##24 Jun 2024 --> 04 Oct 2024
dat_gen <- readRDS(paste0(useDir, "out_moult_", useTask, "_DATA.rds")) ##24 Jun 2024 --> 04 Oct 2024
fit_gen <- readRDS(paste0(useDir, "out_moult_", useTask, "_FIT.rds")) ##24 Jun 2024 --> 04 Oct 2024

fit_mcmc <- growth_gen

##FIT OBJECT
# fit <- rstan::extract(fit_mcmc)
##SUMMARY
# out_summ <- capture.output(print(fit_mcmc))

##PLOTS
samples <- ggmcmc::ggs(fit_mcmc)

#--> P PARAMS ####
samples$Parameter %>% unique()
samples_a <- samples[samples$Parameter %>% str_detect("p\\[1\\]|p\\[2\\]|p\\[3\\]|p\\[4\\]"),]
samples_b <- samples[samples$Parameter %>% str_detect("p\\[5\\]|p\\[6\\]|p\\[7\\]|p\\[8\\]"),]
samples_c <- samples[samples$Parameter %>% str_detect("p\\[9\\]|p\\[10\\]|p\\[11\\]|p\\[12\\]"),]

if ("hist"=="hist"){
  png(paste("outputs/plots/gendata_HISTparams_a.png", sep="")
      , width = use_width
      , height = use_height
      ,units="px")
  plot= ggmcmc::ggmcmc(
    D = samples_a
    , file = NULL
    , plot = 'ggs_histogram'
    , param_page=4
  ) 
  dev.off()
  
  png(paste("outputs/plots/gendata_HISTparams_b.png", sep="")
      , width = use_width
      , height = use_height
      ,units="px")
  plot= ggmcmc::ggmcmc(
    D = samples_b
    , file = NULL
    , plot = 'ggs_histogram'
    , param_page=4
  ) 
  dev.off()
  
  png(paste("outputs/plots/gendata_HISTparams_c.png", sep="")
      , width = use_width
      , height = use_height
      ,units="px")
  plot= ggmcmc::ggmcmc(
    D = samples_c
    , file = NULL
    , plot = 'ggs_histogram'
    , param_page=4
  ) 
  dev.off()
}

if ("trace"=="trace"){
  png(paste("outputs/plots/gendata_TRACEPLOTparams_a.png", sep="")
      , width = use_width
      , height = use_height
      ,units="px")
  plot= ggmcmc::ggmcmc(
    D = samples_a
    , file = NULL
    , plot = 'ggs_traceplot'
    , param_page=4
  ) 
  dev.off()
  
  png(paste("outputs/plots/gendata_TRACEPLOTparams_b.png", sep="")
      , width = use_width
      , height = use_height
      ,units="px")
  plot= ggmcmc::ggmcmc(
    D = samples_b
    , file = NULL
    , plot = 'ggs_traceplot'
    , param_page=4
  ) 
  dev.off()
  
  png(paste("outputs/plots/gendata_TRACEPLOTparams_c.png", sep="")
      , width = use_width
      , height = use_height
      ,units="px")
  plot= ggmcmc::ggmcmc(
    D = samples_c
    , file = NULL
    , plot = 'ggs_traceplot'
    , param_page=4
  ) 
  dev.off()
}

#--> OTHER PARAMS ####
samples$Parameter %>% unique()
samples_d <- samples[samples$Parameter %>% str_detect("growth_mu|growth_sigma|nogrowth_sigma"),]
samples_e <- samples[samples$Parameter %>% str_detect("damageIndicates_moulted|damageIndicates_nonmoulted|maturityIndicates_moulted"),]
samples_f <- samples[samples$Parameter %>% str_detect("maturityIndicates_nonmoulted|pleopodIndicates_moulted|pleopodIndicates_nonmoulted"),]

if ("hist"=="hist"){
  png(paste("outputs/plots/gendata_HISTparams_d.png", sep="")
      , width = use_width
      , height = use_height
      ,units="px")
  plot= ggmcmc::ggmcmc(
    D = samples_d
    , file = NULL
    , plot = 'ggs_histogram'
    , param_page=4
  ) 
  dev.off()
  
  png(paste("outputs/plots/gendata_HISTparams_e.png", sep="")
      , width = use_width
      , height = use_height
      ,units="px")
  plot= ggmcmc::ggmcmc(
    D = samples_e
    , file = NULL
    , plot = 'ggs_histogram'
    , param_page=4
  ) 
  dev.off()
  
  png(paste("outputs/plots/gendata_HISTparams_f.png", sep="")
      , width = use_width
      , height = use_height
      ,units="px")
  plot= ggmcmc::ggmcmc(
    D = samples_f
    , file = NULL
    , plot = 'ggs_histogram'
    , param_page=4
  ) 
  dev.off()
}

if ("trace"=="trace"){
  png(paste("outputs/plots/gendata_TRACEPLOTparams_d.png", sep="")
      , width = use_width
      , height = use_height
      ,units="px")
  plot= ggmcmc::ggmcmc(
    D = samples_d
    , file = NULL
    , plot = 'ggs_traceplot'
    , param_page=4
  ) 
  dev.off()
  
  png(paste("outputs/plots/gendata_TRACEPLOTparams_e.png", sep="")
      , width = use_width
      , height = use_height
      ,units="px")
  plot= ggmcmc::ggmcmc(
    D = samples_e
    , file = NULL
    , plot = 'ggs_traceplot'
    , param_page=4
  ) 
  dev.off()
  
  png(paste("outputs/plots/gendata_TRACEPLOTparams_f.png", sep="")
      , width = use_width
      , height = use_height
      ,units="px")
  plot= ggmcmc::ggmcmc(
    D = samples_f
    , file = NULL
    , plot = 'ggs_traceplot'
    , param_page=4
  ) 
  dev.off()
}
#--> SINGLE GEN PLOT ####
if(100==100){
  #--> Basic plot ####
  
  tmpGen <- data.frame(
    p = fit_gen$p %>% colMeans()
    , sd = fit_gen$p %>% apply(2, sd)
    , months = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    , gp = "Estimated (95% CI)"
  )
  tmpCoded <- data.frame(
    p = c(0.05, 0.05, 0.1, 0.4, 0.7, 0.4, 0.2, 0.15, 0.05, 0.02, 0.02, 0.01)
    , sd = 0
    , months = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    , gp = "Moult Probability"
  )
  
  tmp <- rbind(tmpGen, tmpCoded)
  tmp$months <- factor(tmp$months, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  pgen <- ggplot(tmp, aes(x = months, y = p, group = gp, color = gp, linetype = gp)) +
    geom_point(data = subset(tmp, gp == "Estimated (95% CI)"), position = position_jitter(width = 0), alpha = 0.5) +
    geom_line() +
    geom_errorbar(data = subset(tmp, gp == "Estimated (95% CI)"), aes(ymin = p - 2 * sd, ymax = p + 2 * sd), width = 0.2, position = position_dodge(0.75)) +
    scale_color_manual(values = c("Estimated (95% CI)" = "black", "Moult Probability" = "red")) +
    scale_linetype_manual(values = c("Estimated (95% CI)" = "solid", "Moult Probability" = "dashed")) +
    labs(x = "Calendar month", y = "Probability of moulting", color = "", linetype = "") +
    theme_minimal() +
    theme(
      legend.position = c(0.8, 0.95),  # Adjust the legend position to top-left corner
      legend.key.size = unit(0.75, "lines"),  # Adjust legend key size
      legend.text = element_text(size = 10),
    ) + ylim(-0.1, 1)  # Fix the vertical axis from 0 to 1.3
  
  
  #--> Add label ####
  use_label_1 <- c("Ancillary data: 30%")
  use_label_2 <- c("Liberty: 18 months")
  use_label_3 <- c("Growth: N(5,2)")
  use_label_4 <- c("Lines of evidence: All")
  
  p <- pgen + geom_text(aes(x = 1, y = 0.95), label = paste0(use_label_1), hjust = 0, vjust = 0, size = 3.5, color = "black",family = "Arial" )
  p <- p + geom_text(aes(x = 1, y = 0.87), label = paste0(use_label_2), hjust = 0, vjust = 0, size = 3.5, color = "black",family = "Arial" )
  p <- p + geom_text(aes(x = 1, y = 0.79), label = paste0(use_label_3), hjust = 0, vjust = 0, size = 3.5, color = "black",family = "Arial" )
  p <- p + geom_text(aes(x = 1, y = 0.71), label = paste0(use_label_4), hjust = 0, vjust = 0, size = 3.5, color = "black",family = "Arial" )
  p <- p +     theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

  ggsave("outputs/plots/plot_genplot.png", p, width = 4*480, height = 2*480, units = "px", dpi = 300)
  

}
#####

#############################
#### GENERATED DATA PLOT ####
#############################
#--> Get ggs plots for damage model ####
useDir <- "outputs/"

dir(useDir) 

use_height <- 480
use_width <- 480

tasks <- c("gen05lib18grow52", #"gen30lib18grow52", 
           "gen50lib18grow52", "gen80lib18grow52",
           
           "gen30lib24grow52", #"gen30lib18grow52", 
           "gen30lib12grow52", "gen30lib06grow52",
           
           "gen30lib18grow72", #"gen30lib18grow52", 
           "gen30lib18grow32", "gen30lib18grow11",
           
           "gen30lib18grow52_grDamage", #"gen30lib18grow52", 
           "gen30lib18grow52_grPleopod", "gen30lib18grow52_grMaturity"
           )

store_scenario1234 <- list()

for (nn in 1:length(tasks)){
  #--> Set scenario element ####
  useTask <- tasks[nn]
  
  #--> Load data ####
  # growth_gen <- readRDS("_outputs/moult_LLonly_21/out_moult_generated_FITMCMC.rds")
  growth_gen <- readRDS(paste0(useDir, "out_moult_", useTask, "_FITMCMC.rds")) ##24 Jun 2024 --> 04 Oct 2024
  dat_gen <- readRDS(paste0(useDir, "out_moult_", useTask, "_DATA.rds")) ##24 Jun 2024 --> 04 Oct 2024
  fit_gen <- readRDS(paste0(useDir, "out_moult_", useTask, "_FIT.rds")) ##24 Jun 2024 --> 04 Oct 2024
  # growth_gen <- readRDS("_outputs/moult/out_moult_generated_FITMCMC.rds")
  
  fit_mcmc <- growth_gen
  
  ##FIT OBJECT
  # fit <- rstan::extract(fit_mcmc)
  ##SUMMARY
  # out_summ <- capture.output(print(fit_mcmc))
  
  ##PLOTS
  samples <- ggmcmc::ggs(fit_mcmc)#, family="b")
  
  #--> Basic plot ####
  tmpGen <- data.frame(
    p = fit_gen$p %>% colMeans()
    , sd = fit_gen$p %>% apply(2, sd)
    , months = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    , gp = "Est. (95% CI)"
  )
  tmpCoded <- data.frame(
    p = c(0.05, 0.05, 0.1, 0.4, 0.7, 0.4, 0.2, 0.15, 0.05, 0.02, 0.02, 0.01)
    , sd = 0
    , months = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    , gp = "Moult Prob."
  )
  
  tmp <- rbind(tmpGen, tmpCoded)
  tmp$months <- factor(tmp$months, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  pgen <- ggplot(tmp, aes(x = months, y = p, group = gp, color = gp, linetype = gp)) +
    geom_point(data = subset(tmp, gp == "Est. (95% CI)"), position = position_jitter(width = 0), alpha = 0.5) +
    geom_line() +
    geom_errorbar(data = subset(tmp, gp == "Est. (95% CI)"), aes(ymin = p - 2 * sd, ymax = p + 2 * sd), width = 0.2, position = position_dodge(0.75)) +
    scale_color_manual(values = c("Est. (95% CI)" = "black", "Moult Prob." = "red")) +
    scale_linetype_manual(values = c("Est. (95% CI)" = "solid", "Moult Prob." = "dashed")) +
    labs(x = "Calendar month", y = "Probability of moulting", color = "", linetype = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  store_scenario1234[[nn]] <- pgen    
  
  ######
}

#--> Plot pgen ####
store_scenario1234[[1]]
store_scenario1234[[2]]
store_scenario1234[[3]]
store_scenario1234[[4]]

#--> Conditionally modify plots  ####
for (n in 1:length(store_scenario1234)){
  p <- store_scenario1234[[n]]
  
  #--> Update axis text size ####
  if (n %in% c(1:length(store_scenario1234))) {
    p <- p + theme(
      axis.text.y = element_text(size = 9),     # Optionally, increase axis text size as well
      axis.text.x = element_text(size = 9),     # Optionally, increase axis text size as well
    )
  }
  
  #--> Legend on top right plot only ####
  if (n %in% c(3)){
    p <- p + theme(
      legend.text = element_text(size = 8)  # Increase legend text size
      ,legend.position = c(0.95, 0.95)  # Adjust the legend position
      ,legend.justification = c(1, 1)  # Adjust justification for top-right
      ,legend.key.size = unit(0.75, "lines") # Adjust legend key size
    )
  } else {
    p <- p + theme(legend.position = "none")  # Completely remove the legend
  }
  
  #--> Add conditional labelling ####
  if (n %in% c(1:9)){
    p <- p + theme(
      axis.text.x = element_text(colour = "white")
    )
  }

  if (n %!in% c(1,4,7,10)){
    p <- p + theme(
      axis.text.y = element_text(colour = "white")
    )
  }
  
  #--> Add conditional titles ####
  if (n %in% c(1:length(store_scenario1234))) {
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  }
  
  #--> Adjust ylim for all ####
  p <- p + scale_y_continuous(breaks = c(0.0, 0.5, 1.0), labels = scales::number_format(accuracy = 0.1), expand = c(0, 0)) +
    coord_cartesian(ylim = c(-0.1, 1.1))

  #--> Add text label ####
  use_label_1 <- c(rep("Scenario 1", 3), rep("Scenario 2", 3), rep("Scenario 3", 3),rep("Scenario 4", 3))
  
  use_label_2 <- c("Ancillary data: 5%", "Ancillary data: 50%", "Ancillary data: 80%",
                   "Liberty: 24 months", "Liberty: 12 months", "Liberty: 6 months",
                   "Growth: N(7,2)", "Growth: N(3,2)", "Growth: N(1,1)",
                   "Lines of evidence: +Damage", "Lines of evidence: +Pleopod", "Lines of evidence: +Maturity")

  p <- p + geom_text(aes(x = 0.5, y = 0.925), label = paste0(use_label_1[n], " - ", use_label_2[n]), hjust = 0, vjust = 0, size = 3.5, color = "black",family = "Arial" )

  
  #--> Update list ####
  store_scenario1234[[n]] <- p
}

#--> Ordering ####
# Specify the order for grid.arrange
order <- c(1,2,3,4,5,6,7,8,9,10,11,12)  # Change the order as needed
order_use <- order + 0

## Arrange the plots using grid.arrange
arranged_plots <- arrangeGrob(grobs = store_scenario1234[order_use]
                              , ncol = 3
                              , order = order_use
                              , heights = rep(1, length(order_use) / 3)  # Make heights uniform, reduces space between rows
                              , widths = rep(1, 3)  # Make widths uniform, reduces space between columns
                              , padding = unit(0,"line")
                              , margin = unit(0,"line")
                              )

#--> Organise grob plots #####
yleft = richtext_grob(text="Probability of moulting"
                      , rot=90
                      , gp = gpar(col = "black", fontsize = 16)
)

bottom_month = richtext_grob(text = '<span style="color:black">Calendar month</span>'
                             ,gp = gpar(col = "black", fontsize = 16)
)

#--> Save Line Plots ####
png("outputs/plots/plot_out_scenario1234.png", width = 3*4*480, height = 4*2*480, units = "px", res = 600)
grid.arrange(arranged_plots
             , ncol=1
             , left = yleft
             , bottom = bottom_month
)
dev.off()
#####

##################################################################
##################################################################