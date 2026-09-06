# Names used inside ggplot2 aesthetics and tidyr calls that R CMD check would
# otherwise report as undefined global variables.
#
# These are declared here, in one place, rather than next to the functions that
# need them. Placing a utils::globalVariables() call between a roxygen block and
# the function it documents attaches the block to that call instead of to the
# function, so roxygen cannot derive a \usage section and R CMD check reports
#   Rd files without \usage: ...
#   \arguments should not be documented without \usage.
utils::globalVariables(c(".data", "ID", "Value", "Probability"))
