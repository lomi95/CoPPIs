# Co-expressed Protein Protein Interactions (CoPPIs) algorithm

CoPPIs is an algorithm developed in R. Its purpose is to measure the level of correlation between physically/functionally interacting proteins in protein complexes, GO biological terms and pathways.
Moreover, in pairwise comparisons, provides a score indicating the condition/group where the complex/process/pathway is most correlated.
## 📦 Installation

To install this package, make sure you have RCy3 installed.
You can install it with:
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RCy3")
```

To install CoPPIs run the following command: 
```r
devtools::install_github("lomi95/CoPPIs")
```

## 🚀 Usage
Read your data
```r
dataset <- openxlsx::read.xlsx("yourDataTable.xslx")
nog <- c("Group_1", "Group_2", "Group_3") # GroupN
genes_id <- dataset$yourGeneColumn
```
Make sure that in `colnames(dataset)` are present "Group_1", "Group_2", "Group_N" and that `genes_id` has the same length of your dataset rows.
To obtain significant results, it is recommended to process a number of profiles per group of no less than 20.

If you want to visualize the results in Cytoscape, open the software before to run the CoPPIs code. 
Otherwise, add the parameter `Cytoscape = FALSE`
```r
CoPPI_results <- CoPPIs_pipeline(dataset = dataset,
                                   names_of_groups = nog,
                                   genes_id = genes_id,
                                   name_analysis = "CoPPIs_analysis")

```
The execution time depends on dataset size, and the resulting R object may occupy between 1 and 10 GB of memory.

Other parameters 
## 📂 Output

Results are stored in a directory specified by the `name_analysis` parameter. The output includes:
- A **Cytoscape file** in .cys format (if `Cytoscape = TRUE`, default).
- A set of **Excel files**, each corresponding to the selected categories (default: `categories = c("Component", "Process", "RCTM", "WikiPathways")`). Further categories may be inserted through `categories.hs` vector. However, some of them may be organism-dependent.
    The columns of Excel files are:
  - *description* and *term*: Respectively the name and the identifier of the significant term.
  - *preferredNames* and *number_of_genes*: Respectively the names and the number of the genes that are annotated with the term.
  -  *N.PPI* and *Nodes.PPI: Respectively the number of known Protein-Protein Interaction between *preferredNames* and the number of interacting proteins.
  -  *genes_found* and *number_of_genes_found*:  Respectively the names and the number of the genes in the analyzed dataset that are annotated with the term.
  -  *N.edges*: The number of known Protein-Protein Interaction between *genes_found*.
  -  *cor.GroupX*: The average correlation of the considered group.
  -  *Sign_edges.GroupX* and *Sign_nodes.GroupX*: Respectively the number of significant correlation and the number of proteins that significantly correlates.
  -  *Perc_edges.GroupX* and *Perc_nodes.GroupX*: Respectively the percentage of significant correlation on *N.edges* and the percentage of proteins that significantly correlates on *Nodes.PPI*.
  -  *p.adj GroupX vs GroupY* and *score GroupX vs GroupY*: Respectively the p value adjusted and the CoPPIs score (positive if significant for GroupX, negative if significant for GroupY)
- The Cytoscape file and Excel tables contain equivalent results.

### 🔍 Structure of `CoPPIs_results`
The `CoPPIs_results` object is structured as a hierarchical list containing multiple components:

- **`prot.annotation`**: A collection of category-specific objects, each containing a data frame with annotated biological information extracted from the dataset.
- **`gCOR.groups` and `gCOR.categories`**: Lists indexed by experimental groups, storing correlation metrics that can be reused in subsequent analyses to optimize computational efficiency (in CoPPIs_pipeline parameters **`gCOR.groups` and `gCOR.categories`**).
- **`parameters`**: A list detailing all parameters used in the execution.
- **`resultsCoPPI`**: A category-specific list containing:
  - *`table_summary`*: A data frame summarizing results, corresponding to the `node_table.terms` sheet in the generated Excel files.
  - *`gPath.All` and `gPath.Sign`*: Lists of igraph objects, representing the full set of annotated terms and the subset of statistically significant terms per group, respectively. Edges are weighted by correlation values.
  - *`sig_path`*: Data frames detailing significantly enriched terms for pairwise group comparisons.
  - *`tab_pathways_protein`*: A matrix where rows correspond to biological terms and columns indicate the number of associated proteins.
  - *`graph_similarity`*: A pairwise similarity matrix computed from significant protein correlations per term and used for the assessment of similarity among terms.

This structure facilitates further downstream analyses and integrative studies within the same experimental framework.


