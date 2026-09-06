# =============================================================================
# rmsBMA 0.2.0 -- pre-submission checks
#
# Run these in order. Steps 1-5 are local and take a few minutes. Steps 6-8
# submit to remote check machines and return results by e-mail, so start them
# early and read the mail before step 9.
#
# IMPORTANT: start from a fresh R session (Session > Restart R, or Cmd+Shift+F10
# in RStudio). A stale lazy-load database from 0.1.2 will make a perfectly good
# package look broken.
# =============================================================================

pkg <- "~/Programming/rmsBMA"

install.packages(c("devtools", "rhub", "urlchecker", "spelling", "revdepcheck"))
# revdepcheck is not on CRAN:
# remotes::install_github("r-lib/revdepcheck")


# -----------------------------------------------------------------------------
# 1. Regenerate documentation, then check that nothing changed unexpectedly
# -----------------------------------------------------------------------------
# You have roxygen2 8.1.0; the man/ files in the repo were written with it.
# If document() produces a diff, inspect it before going further.

devtools::document(pkg)

# then, in a terminal:
#   cd ~/Programming/rmsBMA && git status --short && git diff -- man NAMESPACE


# -----------------------------------------------------------------------------
# 2. Install and run the test suite on its own
# -----------------------------------------------------------------------------
# Faster feedback than a full check, and it tells you whether a failure is in
# the package or in the check environment.

devtools::install(pkg, build_vignettes = TRUE, force = TRUE)
devtools::test(pkg)

# Sanity-check the two things that changed most:
library(rmsBMA)
data(modelSpace)
print(modelSpace)                       # new print method, 6-element object
b <- bma(modelSpace)
summary(b)
coef(b)
coef(b, prior = "random", conditional = TRUE)
model_sizes(b, type = "histogram")
model_pmp(b, top = 10, type = "histogram")


# -----------------------------------------------------------------------------
# 3. Full local check, as CRAN runs it
# -----------------------------------------------------------------------------
# This is the one that must come back clean. manual = TRUE builds the PDF
# reference manual, which is where the 0.1.2 trouble started.

devtools::check(
  pkg,
  document    = FALSE,   # already done in step 1
  manual      = TRUE,
  remote      = TRUE,    # enables the CRAN incoming checks (DOIs, URLs, version)
  incoming    = TRUE,
  run_dont_test = TRUE
)

# Equivalent from a terminal, if you prefer to see the raw log:
#   cd ~/Programming
#   R CMD build rmsBMA
#   R CMD check --as-cran rmsBMA_0.2.0.tar.gz
#   cat rmsBMA.Rcheck/00check.log


# -----------------------------------------------------------------------------
# 4. URLs and DOIs
# -----------------------------------------------------------------------------
# CRAN fetches every URL and DOI in DESCRIPTION, the Rd files and the README.
# A 403 or 404 here becomes a NOTE on submission.

urlchecker::url_check(pkg)

# The Madigan & York DOI (10.2307/1403615) resolves to JSTOR. JSTOR sometimes
# refuses automated requests, which shows up as "Status: 403". If that happens,
# it is a false positive -- say so in cran-comments.md rather than removing it.


# -----------------------------------------------------------------------------
# 5. Spelling
# -----------------------------------------------------------------------------
# inst/WORDLIST already holds the accepted technical terms. Add to it rather
# than rewording correct statistical vocabulary.

spelling::spell_check_package(pkg)
# spelling::update_wordlist(pkg)   # only after reading the list above


# -----------------------------------------------------------------------------
# 6. Windows check machines (results by e-mail, ~30 min)
# -----------------------------------------------------------------------------
# CRAN checks against R-devel first. Your Mac runs neither R-devel nor Windows,
# so this is the step that finds problems your local check cannot.

devtools::check_win_devel(pkg)
devtools::check_win_release(pkg)


# -----------------------------------------------------------------------------
# 7. macOS builder
# -----------------------------------------------------------------------------

devtools::check_mac_release(pkg)


# -----------------------------------------------------------------------------
# 8. R-hub: the flavours that broke 0.1.2
# -----------------------------------------------------------------------------
# R-hub v2 runs on GitHub Actions, so it needs the package in a GitHub repo.
# One-time setup, which commits a workflow file:
#
#   rhub::rhub_setup()
#   # commit and push the .github/workflows/rhub.yaml it creates
#
# Then:

rhub::rhub_platforms()        # list what is available

rhub::rhub_check(
  platforms = c(
    "linux",     # r-devel on Ubuntu
    "nold",      # no long double -- exercises the new numeric guard
    "atlas",     # alternative BLAS, same family of problem as BLIS
    "valgrind"   # slow, but catches memory issues; optional here
  )
)

# The BLIS failure itself: I reproduced the mechanism under BLIS in a Linux
# container and the guarded code returns Inf for a constant y on both BLIS and
# reference BLAS. "atlas" and "nold" are the closest R-hub equivalents.


# -----------------------------------------------------------------------------
# 9. Reverse dependencies
# -----------------------------------------------------------------------------
# rmsBMA has none, so this returns immediately. Run it anyway -- it is cheap
# and CRAN asks about it.

revdepcheck::revdep_check(pkg, num_workers = 4)
# or, without the package:
tools::package_dependencies("rmsBMA", reverse = TRUE,
                            db = available.packages())


# -----------------------------------------------------------------------------
# 10. Update cran-comments.md, then submit
# -----------------------------------------------------------------------------
# cran-comments.md currently asserts "0 errors | 0 warnings | 0 notes".
# Replace that with what the runs above actually produced, listing each
# platform. If the JSTOR DOI came back 403, say so there.

file.edit(file.path(pkg, "cran-comments.md"))

# Then:
devtools::release(pkg)

# release() re-runs the checks, asks the CRAN policy questions, and submits.
# You will get a confirmation e-mail that must be answered before the
# submission enters the queue.
