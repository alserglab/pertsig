# Load perturbation signatures

Reads the perturbation signatures (per-gene logistic regression
coefficients) for human or mouse. The CSV file is downloaded from
<https://alserglab.wustl.edu/files/pertsig/> on the first call and read
from the cache afterwards.

## Usage

``` r
loadSignatures(
  species = c("hsa", "mmu"),
  geneId = c("symbol", "entrez"),
  cacheDir = tools::R_user_dir("pertsig", which = "cache"),
  force = FALSE
)
```

## Arguments

- species:

  Species code: `"hsa"` (human) or `"mmu"` (mouse).

- geneId:

  Gene identifiers to use as row names: `"symbol"` or `"entrez"`. If an
  identifier occurs more than once, only its first row is kept.

- cacheDir:

  Directory where downloaded files are cached.

- force:

  If `TRUE`, download the file again even if it is cached.

## Value

A numeric matrix of signature weights with genes in rows and signatures
in columns. Signatures that merge several perturbations are named after
the first one (`"IFNG + STAT1"` becomes `"IFNG"`).

## Examples

``` r
# \donttest{
signatures <- loadSignatures("mmu", geneId = "entrez")
signatures[1:3, 1:2]
#>        17_ESTRADIOL 2_3_7_8_TETRACHLORODIBENZO_P_DIOXIN
#> 18770     0.2714234                          -1.6146960
#> 12038    -0.3557221                          -0.9174224
#> 230603   -0.2913645                          -1.0929476
# }
```
