#' Chelonoides_tree
#'
#' An ape 'phylo' object representing a time-scaled phylogenetic tree for the
#' 13 known species of Galapagos giant tortoides (Chelonoides spp.). In addition
#' to the branch lengths, the tree object also contains the root age ($root.age)
#' and a data.frame recording the node ages in millions-of-years ago and modern
#' endemic distributions (islands, coordinates) of each tip ($biogeography).
#'
#' @details Tortoise relationships were inferred and the tree timescales using
#' mitochondrial DNA, then origin islands inferred in BioGeoBEARS.
#'
#' Galapagos giant tortoises show highly endemic geographic distributions, even
#' where multiple species co-exist on one island. The longitude-latitude
#' coordinates correspond to the population epicentres for each species.
#' For the internal nodes of the tree, the ancestral island was inferred, but
#' not the precise coordinates, hence why the coordinates are NA for these
#' entries.
#'
#' Island abbreviations: Isa = Isabella, Santi = Santiago, SantaC = Santa Cruz,
#' SanC = San Cristobal, Pin = Pinta, Esp = Espagola, Flo = Floreana
#'
#' @source Time-scaled tree and node-origin islands from Poulakakis et al (2020).
#' Colonization history of Galapagos giant tortoises: Insights from mitogenomes
#' support the progression rule. Journal of Zoological Systematics and
#' Evolutionary Research, 58, 1262-1275
#' @source Modern species locations from Jensen et al (2022). The Galapagos giant
#' tortoise Chelonoidis phantasticus is not extinct. Communications Biology, 5,
#' 546
