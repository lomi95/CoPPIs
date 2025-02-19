# CoPPIs Co-expressed Protein Protein Interactions

CoPPIs is an algorithm developed to identify the most coordinated biological terms and pathways within protein-protein interaction (PPI) networks. 
By leveraging information derived from protein correlations and the structure of PPI interactions, 
CoPPIs enables the detection of highly connected and biologically relevant functional modules.

## 📦 Installation

To install this package, make sure you have RCy3 installed.
You can install it with:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RCy3")

devtools::install_github("lomi95/CoPPIs")
```

## 🚀 Usage
Read your data
```r
dataset <- openxlsx::read.xlsx("yourDataTable.xslx")
nog <- c("Group_1", "Group_2", "Group_3") # GroupN
genes_id <- dataset$GeneColumn
```
Make sure that in `colnames(dataset)` are present "Group_1", "Group_2", "Group_N"
and that `genes_id` has the same length of your dataset rows.

Open Cytoscape if you want to have your results in it, otherwise in your parameters add `Cytoscape = F`
```r
CoPPI_results <- CoPPIs_pipeline(dataset = dataset,
                                   names_of_groups = nog,
                                   genes_id = genes_id,
                                   name_analysis = "CoPPIs_analysis")

```
The execution time depends on dataset size, and the resulting R object may occupy between 1 and 10 GB of memory.

## 📂 Output

Results are stored in a directory specified by the `name_analysis` parameter. The output includes:
- A **Cytoscape file** (if `Cytoscape = TRUE`, default).
- A set of **Excel files**, each corresponding to a selected category (default: `categories = c("Component", "Process", "RCTM", "WikiPathways")`).
- The Cytoscape file and Excel tables contain equivalent data representations.

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
  - *`graph_similarity`*: A pairwise similarity matrix computed from significant protein correlations per term, used for pathway similarity assessment.

This structure facilitates further downstream analyses and integrative studies within the same experimental framework.


