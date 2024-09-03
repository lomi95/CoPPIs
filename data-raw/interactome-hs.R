## code to prepare `interactome.hs` dataset goes here

ppi.int <- build_interactome.BETA(directory_interactome = "C:/Users/WKS/Desktop/Coppi_new/Interactomes/9606.protein.links.detailed.v12.0.txt",
                                  9606,scores_threshold = c("experimental" = 150, "database" = 300))
interactome.hs <- ppi.int$interactome

usethis::use_data(interactome.hs, overwrite = TRUE, compress = "xz")
