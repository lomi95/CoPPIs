## code to prepare `categories.interactome.hs` dataset goes here

ppi.int <- build_interactome.BETA(directory_interactome,9606,scores_threshold = c("experimental" = 150, "database" = 300))
categories.interactome.hs <- ppi.int$categories.interactome

usethis::use_data(categories.interactome.hs, overwrite = TRUE)
