
# Compare two taxonomy

tax_comp_circle(d_pq, "Genus__eukaryome_Glomero", "Genus")
tax_comp_circle(d_pq, "Family__eukaryome_Glomero", "Family")
tax_comp_circle <- function(
    physeq,
    tax_level1 = NULL,
    tax_level2 = NULL,
    suffix_1 = "_1",
    suffix_2 = "_2") {

  tab <- table(physeq@tax_table[, tax_level1], physeq@tax_table[, tax_level2])
  df_circle <- data.frame(
    from = rep(paste0(rownames(tab), suffix_1), times = ncol(tab)),
    to = rep(paste0(colnames(tab), suffix_2), each = nrow(tab)),
    value = as.vector(tab),
    stringsAsFactors = FALSE
  )

  suffix_name <- c(paste0(rownames(tab), suffix_1), paste0(colnames(tab), suffix_2))

  uniq_names <- unique(c(rownames(tab), colnames(tab)))

  grid.col <- circlize::rand_color(length(uniq_names))

  circlize::circos.par(gap.degree = 1)

  circlize::chordDiagram(df_circle)
  circlize::circos.clear()
}

tax_comp_sankey <- function(
    physeq,
    tax_level1 = NULL,
    tax_level2 = NULL){

  df_sank <- as.data.frame(unclass(physeq@tax_table)) %>%
    count(.data[[tax_level1]],.data[[tax_level2]])

  ggplot(df_sank, aes(axis1 = .data[[tax_level1]],
                      axis2 = .data[[tax_level2]],
                      y = n)
  ) +
    ggalluvial::geom_alluvium()
}




formula_taxo <- formula(paste0("~", paste0(colnames(d_pq@tax_table)[c(1:7)], collapse="/")))

data_taxo <- as.data.frame(d_pq@tax_table[,c(1:7)]) |>
  mutate_if(is.character, as.factor)

phy_tree <- ape::as.phylo.formula(formula_taxo, data = data_taxo, collapse=T)


p1 <- ggtree(phy_tree, layout='roundrect')
