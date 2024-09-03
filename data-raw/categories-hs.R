## code to prepare `categories-hs` vector goes here

categories.hs <- c("COMPARTMENTS","Component","DISEASES","Function","HPO",
                "InterPro","Keyword","NetworkNeighborAL","Pfam","Process",
                "RCTM","SMART","TISSUES","WikiPathways")
usethis::use_data(categories.hs, overwrite = TRUE, compress = "xz")
