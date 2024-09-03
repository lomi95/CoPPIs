## code to prepare `annotation.interactome.hs` dataset goes here
ppi.int <- build_interactome.BETA(directory_interactome,9606,scores_threshold = c("experimental" = 150, "database" = 300))
annotation.interactome.hs <- ppi.int$annotation.interactome

usethis::use_data(annotation.interactome.hs, overwrite = TRUE)
