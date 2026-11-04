# Convert signatures to gene sets

Takes the `topK` genes with the largest (and/or smallest) weights in
each signature.

## Usage

``` r
getGeneSets(
  signatures,
  topK = 200,
  direction = c("up", "both", "down"),
  universe = rownames(signatures)
)
```

## Arguments

- signatures:

  Signature weight matrix (genes x signatures), as returned by
  [`loadSignatures()`](https://alserglab.github.io/pertsig/reference/loadSignatures.md).

- topK:

  Number of genes per gene set.

- direction:

  `"up"` (top weights) and `"down"` (bottom weights) make one set per
  signature, named as the signature; `"both"` makes two sets per
  signature, named with `_UP` and `_DN` suffixes.

- universe:

  Genes to choose from, for example the names of the statistic vector
  the gene sets will be tested against. Restricting to it keeps every
  gene set at `topK` genes after fgsea drops genes that are missing from
  the statistic.

## Value

A named list of character vectors.

## Examples

``` r
# \donttest{
signatures <- loadSignatures("mmu", geneId = "entrez")
geneSets <- getGeneSets(signatures, topK = 5)
geneSets[["TBX21"]]
#> [1] "57765"  "12766"  "16162"  "14945"  "235534"
# }
```
