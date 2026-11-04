# pertsig

Perturbation signatures from <https://alserglab.wustl.edu/files/pertsig/>
as gene sets for enrichment analysis with
[fgsea](https://bioconductor.org/packages/fgsea/).

## Installation

```r
# install.packages("remotes")
remotes::install_github("alserglab/pertsig")
```

## Usage

```r
library(pertsig)
library(fgsea)

signatures <- loadSignatures("mmu", geneId = "entrez")  # or "hsa", geneId = "symbol"

data(exampleRanks)
geneSets <- getGeneSets(signatures, topK = 200, universe = names(exampleRanks))
fgsea(geneSets, exampleRanks)
```

See the [tutorial](https://alserglab.github.io/pertsig/articles/pertsig-tutorial.html)
for a full example.
