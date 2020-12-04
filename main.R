# crosslink analysis for Marie (Kathrin Lang Lab)

# dependencies
library(tidyverse)
library(data.table)
library(ggthemes)

source("analyse_crosslinks.R")


# global theme for plots
theme_set(theme_bw())


# constants
output_folder <- "output"
decoy         <- "Decoy_"
rab           <- "sp\\|Q9H0U4\\|RAB1B"
drra          <- "sp\\|Q29ST3\\|DRRA"


# analyse
initial <- analyse_crosslinks(file = file.path("..", "initial_detection", "rerun", "20180207_Rab_DrrA.perc.inter.txt"),
                              decoy = decoy,
                              rab = rab,
                              drra = drra,
                              name = "initial_detection",
                              output_folder = file.path(output_folder, "initial_detection"))


confirm <- analyse_crosslinks(file = file.path("..", "confirmation", "rerun", "20180710_D82C_confirmation.perc.inter.txt"),
                              decoy = decoy,
                              rab = rab,
                              drra = drra,
                              name = "confirmation",
                              output_folder = file.path(output_folder, "confirmation"))

