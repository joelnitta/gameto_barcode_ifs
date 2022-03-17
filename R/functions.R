#' Calculate interspecific distances within a genus
#'
#' @param aln DNA alignment
#' @param genus_select Genus to calculate distances
#' @param g_s_tib Tibble with two columns: genus, species
#'
#' @return Tibble
#' 
calc_intersp_by_genus <- function(aln, genus_select, g_s_tib) {
  
  species_select <- g_s_tib %>%
    filter(genus == genus_select) %>%
    pull(species)
  
  aln <- aln[species_select, ] %>%
    ape::del.colgapsonly()
  
  dist <- ape::dist.dna(
    aln, model = "raw", as.matrix = FALSE, pairwise.deletion = TRUE) %>%
    as.list() %>%
    unlist()
  
  tibble(
    mean = mean(dist, na.rm = TRUE),
    sd = sd(dist, na.rm = TRUE),
    max = max(dist, na.rm = TRUE),
    min = min(dist, na.rm = TRUE),
    genus = genus_select
  )
}