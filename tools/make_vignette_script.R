# ---------------------------------------------------------------------------
# Regenerate tools/vignette_script.R from vignettes/rmsBMA.Rmd.
#
# The script is a runnable companion to the vignette: every R chunk, in order,
# with the section headings kept as comments. Run this after editing the
# vignette so the two cannot drift apart.
#
#   Rscript tools/make_vignette_script.R
#
# Run from the package root.
# ---------------------------------------------------------------------------

rmd <- "vignettes/rmsBMA.Rmd"
out <- "tools/vignette_script.R"
stopifnot(file.exists(rmd))

txt  <- readLines(rmd, warn = FALSE)
yaml <- which(txt == "---")
if (length(yaml) >= 2) txt <- txt[(yaml[2] + 1):length(txt)]

header <- c(
'################################################################################',
'##',
'##  rmsBMA -- runnable companion script to the package vignette',
'##',
'##  Every R chunk from vignettes/rmsBMA.Rmd, in order, with the section',
'##  headings kept as comments so you can navigate and run it piece by piece.',
'##  Do not edit by hand: regenerate with tools/make_vignette_script.R.',
'##',
'##  HOW TO USE',
'##    - Run top to bottom, or step through a section at a time in RStudio.',
'##    - Fold sections with Alt+O (Windows) / Cmd+Alt+O (macOS).',
'##    - Takes roughly fifteen seconds with RUN_SLOW = FALSE.',
'##',
'##  RUN_SLOW = TRUE additionally builds the two large K = 17 model spaces by',
'##  exhaustive enumeration, and rebuilds the packaged `modelSpace` object from',
'##  scratch to check it still matches (about a minute in total). Those chunks',
'##  are eval = FALSE in the vignette, which is why its own build is quick.',
'##',
'################################################################################',
'',
'RUN_SLOW <- FALSE',
'',
'library(rmsBMA)',
'',
'section <- function(x) message("\\n=== ", x, " ===")',
'')

body <- character(0)
i <- 1L; n <- length(txt)
in_chunk <- FALSE; opts <- ""; buf <- character(0)

flush_chunk <- function(opts, code) {
  code <- code[!(seq_along(code) == length(code) & !nzchar(code))]
  if (grepl("include=FALSE", opts, fixed = TRUE)) return(character(0))
  eval_false <- grepl("eval\\s*=\\s*FALSE", opts)
  if (!eval_false) return(code)
  first <- trimws(paste(code, collapse = " "))
  if (grepl("^install\\.packages|^\\?", first)) {
    return(c("# (interactive / installation step, left commented)", paste("#", code)))
  }
  if (grepl("model_space\\(", first)) {
    extra <- character(0)
    if (grepl("modelSpace\\s*<-", first)) {
      extra <- c(
        '  # Same strings, but the packaged object carries a stale dim attribute on',
        '  # its dimnames, so compare values and names rather than attributes.',
        '  message("  rebuilt modelSpace matches the packaged object: ",',
        '          isTRUE(all.equal(unname(modelSpace[[2]]),',
        '                           unname(rmsBMA::modelSpace[[2]]))) &&',
        '          identical(as.character(colnames(modelSpace[[2]])),',
        '                    as.character(colnames(rmsBMA::modelSpace[[2]]))))')
    }
    return(c("if (RUN_SLOW) {", paste0("  ", code), extra, "}"))
  }
  paste("#", code)
}

while (i <= n) {
  line <- txt[i]
  if (!in_chunk && grepl("^```\\{r", line)) {
    in_chunk <- TRUE
    opts <- sub("^```\\{r\\s*", "", sub("\\}\\s*$", "", line))
    buf <- character(0)
  } else if (in_chunk && grepl("^```\\s*$", line)) {
    in_chunk <- FALSE
    got <- flush_chunk(opts, buf)
    if (length(got)) body <- c(body, got, "")
  } else if (in_chunk) {
    buf <- c(buf, line)
  } else if (grepl("^#{1,3} ", line)) {
    lvl   <- nchar(sub("^(#+).*$", "\\1", line))
    title <- trimws(sub("^#+ ", "", line))
    title <- gsub("[`$\\\\]", "", title)
    rule  <- if (lvl == 1L) strrep("#", 79) else NULL
    body <- c(body, "",
              if (lvl == 1L) c(rule, paste("##", toupper(title)), rule)
              else paste0("## ", strrep("-", if (lvl == 2L) 2 else 4), " ", title, " ",
                          strrep("-", max(4, 68 - nchar(title)))),
              sprintf('section("%s")', gsub('"', "'", title)), "")
  }
  i <- i + 1L
}

dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
writeLines(c(header, body, "", 'message("\\n=== script complete ===")', ""), out)
cat("wrote", out, "-", length(c(header, body)), "lines\n")
