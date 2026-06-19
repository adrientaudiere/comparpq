# Regenerate pre-computed vignettes from their .Rmd.orig sources.
# The output .Rmd files are committed to git so that build and build_site skip
# the heavy computations (DA methods, simulations) on every run.
#
# Requirements: ANCOMBC, ALDEx2, maaslin3, DESeq2, edgeR, limma, radEmu
#
# Run from the package root:
#   Rscript vignettes/precompile.R

library(knitr)

old_wd <- setwd("vignettes")
on.exit(setwd(old_wd), add = TRUE)

message("Precompiling compare-da-methods ...")
knitr::knit("compare-da-methods.Rmd.orig", output = "compare-da-methods.Rmd")

message("Precompiling benchmark-da-methods ...")
knitr::knit(
  "benchmark-da-methods.Rmd.orig",
  output = "benchmark-da-methods.Rmd"
)

message("Done. Commit the updated .Rmd files and vignettes/figure/ to git.")
