# ---------------------------------------------------------------------------
# Run this INSTEAD of devtools::check() when preparing a CRAN submission.
#
# Step 1 - build a source tarball:
pkg <- devtools::build()

# Step 2 - check THAT tarball (not the source directory) with --as-cran:
devtools::check_built(pkg, args = "--as-cran")

# Notes
# -----
# The package vignette is vignettes/rmsBMA.Rmd and builds to HTML
# (rmarkdown::html_vignette), so there is no PDF under inst/doc and nothing
# to compact with --compact-vignettes.
#
# The journal manuscript describing the package (rmsBMA.Rnw) lives in a
# SEPARATE project and must never be placed in vignettes/. It shares the
# base name "rmsBMA" with the vignette, so its .pdf/.tex output makes
# R CMD build fail with:
#   "Located more than one 'weave' output file ... for vignette with
#    name 'rmsBMA'"
# .Rbuildignore does not prevent this - keep the manuscript out of the
# package directory entirely.
# ---------------------------------------------------------------------------
