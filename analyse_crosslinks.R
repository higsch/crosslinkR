analyse_crosslinks <- function (file,
                                decoy,
                                rab,
                                drra,
                                name,
                                output_folder) {
  
  # create output folder
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
  # determine number of columns in input file
  column_number <- max(count.fields(file = file, sep = "\t"))
  
  # read data
  data <- read.table(file = file,
                     sep = "\t",
                     fill = TRUE,
                     header = FALSE,
                     col.names = c(1:column_number),
                     quote = "",
                     as.is = TRUE,
                     stringsAsFactors = FALSE)
  message("Data dimensions: ", dim(data))
  
  # repair headers
  headers <- data[1, ]
  headers[headers == ""] <- paste("Proteins_", 2:(sum(headers == "") + 1), sep = "")
  colnames(data) <- headers
  data <- data[-1, ]
  
  # analyse
  data_uni <- data %>%
    unite(Protein, starts_with("Protein"), sep = "__") %>%
    mutate(Protein = gsub("_*$", "", Protein)) %>% # aggregate protein names
    mutate_at(.var = c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13), .funs = as.numeric) %>% # make numerical columns numeric
    arrange(desc(Score)) %>% # sort by score
    mutate(order = 1:nrow(.)) %>% # introduce order numbers
    mutate(hasRab = grepl(rab, Protein)) %>% # check if Rab is present
    mutate(hasDrra = grepl(drra, Protein)) %>% # check if DrrA is present
    mutate(isCrosslink = hasRab & hasDrra) %>%  # check if both are present
    mutate(isDecoy = grepl(decoy, Protein)) %>% # check if it's a decoy
    mutate(qValue = NA) %>%
    mutate(isHit = isCrosslink & !isDecoy) # identify hit
  
  # calculate q-value
  for (i in 1:nrow(data_uni)) {
    isDecoy <- data_uni[1:i, "isDecoy"]
    data_uni[i, "qValue"] <- sum(isDecoy) / sum(!isDecoy)
  }
  
  # extract all hits
  hits <- data_uni %>%
    filter(isHit == TRUE) %>%
    separate(Protein, into = c("Protein_1", "Protein_2"), sep = "__") %>%
    mutate_at(.vars = c("Protein_1", "Protein_2"), .funs = function (x) if_else(grepl(rab, x), "Rab", if_else(grepl(drra, x), "DrrA", "unknown"))) %>%
    select(Protein_1, Protein_2, Peptide, Score, qValue)
  
  # write hits
  write.table(x = hits,
              file = file.path(output_folder, paste0("hits_", name, ".txt")),
              quote = FALSE,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  
  # plot score density with hits scores
  p <- ggplot(data = data_uni) +
    geom_density(aes(Score, group = isDecoy, color = isDecoy)) +
    scale_color_manual(values = c("#78c4c1", "#eb470c")) +
    geom_density(aes(Score), color = "#808d9c") +
    geom_point(data = hits, aes(x = Score, y = 0), color = "#247ee3", size = 3, alpha = .7) +
    ggtitle(label = "Score density and hits")
  
  pdf(file = file.path(output_folder, paste0("score_density_", name, ".pdf")))
  print(p)
  dev.off()
  
  # plot q-values by score with hit scores
  p <- ggplot(data = data_uni, aes(x = Score, y = qValue)) +
    geom_line(color = "#808d9c") +
    geom_point(data = hits, color = "#247ee3", size = 3, alpha = .7) +
    ylab("q-value") +
    ggtitle(label = "q-values vs. score and hits")
  
  pdf(file = file.path(output_folder, paste0("q_value_vs_score_", name, ".pdf")))
  print(p)
  dev.off()
  
  return(data_uni)
}