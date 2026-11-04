signaturesBaseUrl <- "https://alserglab.wustl.edu/files/pertsig/"

#' Load perturbation signatures
#'
#' Reads the perturbation signatures (per-gene logistic regression
#' coefficients) for human or mouse. The CSV file is downloaded from
#' <https://alserglab.wustl.edu/files/pertsig/> on the first call and read
#' from the cache afterwards.
#'
#' @param species Species code: `"hsa"` (human) or `"mmu"` (mouse).
#' @param geneId Gene identifiers to use as row names: `"symbol"` or
#'   `"entrez"`. If an identifier occurs more than once, only its first row
#'   is kept.
#' @param cacheDir Directory where downloaded files are cached.
#' @param force If `TRUE`, download the file again even if it is cached.
#'
#' @return A numeric matrix of signature weights with genes in rows and
#'   signatures in columns. Signatures that merge several perturbations are
#'   named after the first one (`"IFNG + STAT1"` becomes `"IFNG"`).
#'
#' @examples
#' \donttest{
#' signatures <- loadSignatures("mmu", geneId = "entrez")
#' signatures[1:3, 1:2]
#' }
#' @export
loadSignatures <- function(species = c("hsa", "mmu"),
                           geneId = c("symbol", "entrez"),
                           cacheDir = tools::R_user_dir("pertsig", which = "cache"),
                           force = FALSE) {
    species <- match.arg(species)
    geneId <- match.arg(geneId)
    file <- downloadSignatures(species, cacheDir = cacheDir, force = force)

    table <- data.table::fread(file, check.names = FALSE)
    idColumns <- c(entrez = "gene", symbol = "gene_symbol")
    missingColumns <- setdiff(idColumns, names(table))
    if (length(missingColumns) > 0) {
        stop("Missing column(s) in ", file, ": ",
             paste(missingColumns, collapse = ", "))
    }

    geneIds <- as.character(table[[idColumns[[geneId]]]])
    weights <- as.matrix(table[, setdiff(names(table), idColumns), with = FALSE])
    keep <- !is.na(geneIds) & geneIds != "" & !duplicated(geneIds)
    weights <- weights[keep, , drop = FALSE]
    rownames(weights) <- geneIds[keep]
    colnames(weights) <- sub("\\s*\\+.*$", "", colnames(weights))
    weights
}

#' Download a signatures file into the cache
#'
#' @param species Species code, as in `loadSignatures()`.
#' @param cacheDir Cache directory.
#' @param force Download again even if the file is cached.
#' @return Path to the cached file.
#' @keywords internal
#' @noRd
downloadSignatures <- function(species, cacheDir, force = FALSE) {
    fileName <- sprintf("logreg_signatures_%s.csv", species)
    path <- file.path(cacheDir, fileName)
    if (file.exists(path) && !force) {
        return(path)
    }

    dir.create(cacheDir, recursive = TRUE, showWarnings = FALSE)
    oldOptions <- options(timeout = max(600, getOption("timeout")))
    on.exit(options(oldOptions), add = TRUE)

    partialPath <- paste0(path, ".part")
    utils::download.file(paste0(signaturesBaseUrl, fileName),
                         destfile = partialPath, mode = "wb")
    file.rename(partialPath, path)
    path
}
