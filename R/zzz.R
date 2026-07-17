# =============================================================================
# zzz.R
#
# Package-level housekeeping. `networkY()` uses dplyr's non-standard
# evaluation (columns referenced unquoted inside mutate/filter/etc.), which
# R CMD check otherwise flags as "no visible binding for global variable".
# Declaring them here is the standard way to silence that check note - it
# has no effect on how the functions actually run.
# =============================================================================

utils::globalVariables(c(
  "kinship", "from", "to", "x", "y", "id",
  "fx", "fy", "tx", "ty", "midX", "midY",
  "toCenterX", "toCenterY", "ex", "ey", "cross"
))
