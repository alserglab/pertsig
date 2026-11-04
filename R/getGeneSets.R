#' Convert signatures to gene sets
#'
#' Takes the `topK` genes with the largest (and/or smallest) weights in each
#' signature.
#'
#' @param signatures Signature weight matrix (genes x signatures), as
#'   returned by [loadSignatures()].
#' @param topK Number of genes per gene set.
#' @param direction `"up"` (top weights) and `"down"` (bottom weights) make
#'   one set per signature, named as the signature; `"both"` makes two sets
#'   per signature, named with `_UP` and `_DN` suffixes.
#' @param universe Genes to choose from, for example the names of the
#'   statistic vector the gene sets will be tested against. Restricting to
#'   it keeps every gene set at `topK` genes after fgsea drops genes that
#'   are missing from the statistic.
#'
#' @return A named list of character vectors.
#'
#' @examples
#' \donttest{
#' signatures <- loadSignatures("mmu", geneId = "entrez")
#' geneSets <- getGeneSets(signatures, topK = 5)
#' geneSets[["TBX21"]]
#' }
#' @export
getGeneSets <- function(signatures,
                        topK = 200,
                        direction = c("up", "both", "down"),
                        universe = rownames(signatures)) {
    direction <- match.arg(direction)
    signatures <- signatures[rownames(signatures) %in% universe, , drop = FALSE]
    topK <- min(topK, nrow(signatures))

    columns <- asplit(signatures, 2)
    topGenes <- function(weights, decreasing) {
        names(weights)[order(weights, decreasing = decreasing)[seq_len(topK)]]
    }
    upSets <- lapply(columns, topGenes, decreasing = TRUE)
    downSets <- lapply(columns, topGenes, decreasing = FALSE)

    if (direction == "up") {
        return(upSets)
    }
    if (direction == "down") {
        return(downSets)
    }
    names(upSets) <- paste0(names(upSets), "_UP")
    names(downSets) <- paste0(names(downSets), "_DN")
    c(upSets, downSets)
}
