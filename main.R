library(tidyverse)
library(data.table)
library(ggthemes)

source("analyse_crosslinks.R")

theme_set(theme_tufte())

name <- "confirmation"
output_folder <- "output"
decoy <- "Decoy_"
rab <- "sp\\|Q9H0U4\\|RAB1B"
drra <- "sp\\|Q29ST3\\|DRRA"


initial <- analyse_crosslinks(file = file.path("..", "initial_detection", "rerun", "20180207_MS_Rab1b_DrrA_complex_SEC_crosslinks.perc.inter.txt"),
                              decoy = decoy,
                              rab = rab,
                              drra = drra,
                              name = "initial_detection",
                              output_folder = file.path(output_folder, "initial_detection"))


confirm <- analyse_crosslinks(file = file.path("..", "confirmation", "rerun", "20180710_Marie.perc.inter.txt"),
                              decoy = decoy,
                              rab = rab,
                              drra = drra,
                              name = "confirmation",
                              output_folder = file.path(output_folder, "confirmation"))

