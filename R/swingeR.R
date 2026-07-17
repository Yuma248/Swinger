# =============================================================================
# swingeR.R
#
# R translation of a Perl script that selects breeding/family groups from a
# pairwise relatedness matrix while respecting relatedness thresholds.
#
# INPUTS ARE NOW TWO DATA FRAMES (as produced by, e.g., a kinship package):
#
#   kinship_df: long-format pairwise relatedness, one row per pair, e.g.
#       from        to        kinship
#       SPBR_002    SPBR_001  0.207414
#     Column names are configurable (defaults: "from", "to", "kinship").
#     Only needs to contain one direction of each pair; it will be
#     symmetrized automatically. Self-pairs (from == to) are ignored -
#     the diagonal of the relatedness matrix is always set to 0.
#
#   info_df: one row per individual, e.g.
#       ID          Sex   IR
#       SPBR_001    M     0.012010
#     Column names are configurable (defaults: "ID", "Sex", "IR").
#
# A file-based reader (read_relatedness_matrix) is still provided below for
# backward compatibility if you ever have data in the original flat-file
# format instead.
# =============================================================================

# ---- 1a. Build the internal data structure from two data frames -----------

#' Build the internal relatedness data structure from two data frames
#'
#' @param kinship_df Long-format data frame with columns \code{from},
#'   \code{to}, \code{kinship} (or custom names, see below).
#' @param info_df Data frame with one row per individual: \code{ID},
#'   \code{Sex}, \code{F} (or custom names, see below).
#' @param from_col,to_col,kin_col Column names in \code{kinship_df}.
#'   Default \code{"from"}, \code{"to"}, \code{"kinship"}.
#' @param id_col,sex_col,ir_col Column names in \code{info_df}. Default
#'   \code{"ID"}, \code{"Sex"}, \code{"F"}.
#'
#' @return A list with elements \code{sex}, \code{ir}, \code{id}, \code{rel}
#'   (a symmetric relatedness matrix with 0 diagonal), and \code{n}. Used
#'   internally by \code{\link{swingeR}}.
#' @export
build_relatedness_data <- function(kinship_df, info_df,
                                    from_col = "from", to_col = "to",
                                    kin_col = "kinship",
                                    id_col = "ID", sex_col = "Sex",
                                    ir_col = "F") {

  id  <- as.character(info_df[[id_col]])
  sex <- as.character(info_df[[sex_col]])
  ir  <- suppressWarnings(as.numeric(info_df[[ir_col]]))
  n   <- length(id)

  if (anyDuplicated(id)) {
    stop("info_df contains duplicate IDs")
  }
  if (!all(sex %in% c("M", "F"))) {
    stop("Sex column must contain only 'M' or 'F' (row(s): ",
         paste(which(!sex %in% c("M", "F")), collapse = ", "), ")")
  }
  if (any(is.na(ir) | ir < -1 | ir > 1)) {
    stop("Internal relatedness (F) values must be numeric and within [-1, 1]")
  }

  from <- as.character(kinship_df[[from_col]])
  to   <- as.character(kinship_df[[to_col]])
  kin  <- suppressWarnings(as.numeric(kinship_df[[kin_col]]))

  if (any(is.na(kin) | kin < -1 | kin > 1)) {
    stop("Kinship values must be numeric and within [-1, 1]")
  }

  idx_from <- match(from, id)
  idx_to   <- match(to, id)
  unmatched <- is.na(idx_from) | is.na(idx_to)
  if (any(unmatched)) {
    warning(sum(unmatched), " row(s) in kinship_df reference IDs not found ",
            "in info_df and will be ignored")
  }

  rel <- matrix(0, n, n, dimnames = list(id, id))
  ok <- !unmatched
  self_pair <- ok & (from == to)
  ok <- ok & !self_pair
  rel[cbind(idx_from[ok], idx_to[ok])] <- kin[ok]
  rel[cbind(idx_to[ok], idx_from[ok])] <- kin[ok]
  diag(rel) <- 0

  list(sex = sex, ir = ir, id = id, rel = rel, n = n)
}

# ---- 1b. (Legacy) Read & validate the relatedness matrix from a flat file --

#' Read a relatedness matrix from a flat file (legacy format)
#'
#' @param file Path to a tab-delimited file with columns Sex, F, ID,
#'   followed by one relatedness column per individual (see package README
#'   for the exact expected layout).
#' @return A list with elements \code{sex}, \code{ir}, \code{id}, \code{rel},
#'   \code{n} — same structure as \code{\link{build_relatedness_data}}.
#' @export
read_relatedness_matrix <- function(file) {
  raw <- read.delim(file, header = TRUE, stringsAsFactors = FALSE,
                     check.names = FALSE)

  sex <- as.character(raw[[1]])
  ir  <- suppressWarnings(as.numeric(raw[[2]]))
  id  <- as.character(raw[[3]])
  rel <- as.matrix(raw[, -(1:3), drop = FALSE])
  storage.mode(rel) <- "double"

  if (!all(sex %in% c("M", "F"))) {
    stop("Sex column must contain only 'M' or 'F' (row(s): ",
         paste(which(!sex %in% c("M", "F")), collapse = ", "), ")")
  }
  if (any(is.na(ir) | ir < -1 | ir > 1)) {
    stop("Internal relatedness (F) values must be numeric and within [-1, 1]")
  }
  if (any(is.na(rel) | rel < -1 | rel > 1)) {
    stop("Pairwise relatedness values must be numeric and within [-1, 1]")
  }
  if (nrow(rel) != ncol(rel)) {
    stop("Relatedness block is not square: check that the file has exactly ",
         "one relatedness column per individual")
  }

  dimnames(rel) <- list(id, id)
  list(sex = sex, ir = ir, id = id, rel = rel, n = length(id))
}

# ---- 2. Helper: check every pairwise value in a combo is <= a threshold ----

.combo_within_threshold <- function(idx, rel, threshold) {
  if (length(idx) < 2) return(TRUE)
  pairs <- utils::combn(idx, 2)
  vals  <- rel[cbind(pairs[1, ], pairs[2, ])]
  all(vals <= threshold)
}

# ---- 3. Generate valid same-sex combinations of a given size --------------

.valid_same_sex_combos <- function(pool, size, rel, max_pair_rel) {
  if (length(pool) < size) return(list())
  combo_mat <- utils::combn(pool, size)
  keep <- apply(combo_mat, 2, .combo_within_threshold,
                rel = rel, threshold = max_pair_rel)
  lapply(seq_len(ncol(combo_mat))[keep], function(j) combo_mat[, j])
}

# ---- 4. Build candidate family groups (nF females + nM males) -------------

.build_family_candidates <- function(female_combos, male_combos, rel,
                                      max_pair_rel, max_group_avg_rel,
                                      max_group_couple_avg_rel) {
  candidates <- list()
  for (fc in female_combos) {
    for (mc in male_combos) {
      members <- c(fc, mc)

      # cross-sex pairwise relatedness check (THI)
      cross_pairs <- expand.grid(f = fc, m = mc)
      cross_vals  <- rel[cbind(cross_pairs$f, cross_pairs$m)]
      if (any(cross_vals > max_pair_rel)) next

      # whole-group average relatedness (all pairs)
      all_pairs <- utils::combn(members, 2)
      group_avg <- mean(rel[cbind(all_pairs[1, ], all_pairs[2, ])])

      # cross-sex ("couples") average relatedness
      couple_avg <- mean(cross_vals)

      if (group_avg <= max_group_avg_rel &&
          couple_avg <= max_group_couple_avg_rel) {
        key <- paste(sort(members), collapse = "_")
        candidates[[key]] <- list(members = members,
                                   group_avg = group_avg,
                                   couple_avg = couple_avg)
      }
    }
  }
  candidates
}

# ---- 5. Search for the best set of n_groups non-overlapping candidates ----

.best_group_combinations <- function(candidates, n_groups, max_combo_avg_rel,
                                      max_results) {
  if (length(candidates) == 0) return(list())

  if (n_groups == 1) {
    ord <- order(vapply(candidates, function(x) x$group_avg, numeric(1)))
    picks <- candidates[ord]
    picks <- picks[seq_len(min(max_results, length(picks)))]
    return(lapply(picks, function(x) list(groups = list(x),
                                           combo_avg = x$group_avg)))
  }

  idx_combos <- utils::combn(length(candidates), n_groups, simplify = FALSE)
  results <- list()

  for (ic in idx_combos) {
    grp_set <- candidates[ic]
    all_members <- unlist(lapply(grp_set, function(x) x$members))
    if (anyDuplicated(all_members)) next  # groups must not share individuals

    combo_avg <- mean(vapply(grp_set, function(x) x$group_avg, numeric(1)))
    if (combo_avg <= max_combo_avg_rel) {
      results[[length(results) + 1]] <- list(groups = grp_set,
                                              combo_avg = combo_avg)
    }
  }

  if (length(results) == 0) return(list())
  ord <- order(vapply(results, function(x) x$combo_avg, numeric(1)))
  results <- results[ord]
  results[seq_len(min(max_results, length(results)))]
}

# ---- 6. Main user-facing function ------------------------------------------

#' Select breeding/family groups from a pairwise relatedness matrix
#'
#' The main function of the \pkg{swingeR} package. Given a pairwise
#' kinship/relatedness data frame and per-individual sex + internal
#' relatedness data, selects sets of non-overlapping breeding/family groups
#' that satisfy the specified relatedness thresholds.
#'
#' @param kinship_df Long-format data frame with columns `from`, `to`,
#'   `kinship` (column names configurable via from_col/to_col/kin_col).
#' @param info_df Data frame with one row per individual: `ID`, `Sex`, `F`
#'   (column names configurable via id_col/sex_col/ir_col).
#' @param from_col,to_col,kin_col Column names in kinship_df.
#' @param id_col,sex_col,ir_col Column names in info_df.
#' @param n_groups Number of family groups to select (must not share members).
#' @param n_females,n_males Number of females / males per family group.
#' @param max_ir_f,max_ir_m Max internal relatedness allowed for females /
#'   males before they're excluded from consideration.
#' @param max_pair_rel_f,max_pair_rel_m Max pairwise relatedness allowed
#'   within same-sex combos (female-female, male-male).
#' @param max_pair_rel Max pairwise relatedness allowed between the sexes
#'   (i.e. within a breeding "couple").
#' @param max_group_avg_rel Max average relatedness across all pairs within
#'   a family group.
#' @param max_group_couple_avg_rel Max average relatedness across cross-sex
#'   pairs only, within a family group.
#' @param max_combo_avg_rel Max average relatedness across the n_groups
#'   selected together.
#' @param auto_tune If TRUE, relax thresholds by +0.01 (up to 2 times) when
#'   no valid solution is found, mirroring the original script's behaviour.
#' @param max_results How many top (lowest relatedness) combinations to
#'   return.
#' @param verbose Print progress messages.
#'
#' @return A list of result sets. Each result set contains a list of
#'   family groups (member IDs + relatedness stats) and the overall
#'   combination average relatedness.
#' @export
#'
#' @examples
#' \dontrun{
#' results <- swingeR(kinship_df, info_df, n_groups = 2, n_females = 2, n_males = 1)
#' print_breeding_groups(results)
#' }
swingeR <- function(kinship_df, info_df,
                     from_col = "from", to_col = "to",
                     kin_col = "kinship",
                     id_col = "ID", sex_col = "Sex",
                     ir_col = "F",
                     n_groups, n_females, n_males,
                     max_ir_f = 1, max_ir_m = 1,
                     max_pair_rel_f = 1, max_pair_rel_m = 1,
                     max_pair_rel = 1,
                     max_group_avg_rel = 1,
                     max_group_couple_avg_rel = 1,
                     max_combo_avg_rel = 1,
                     auto_tune = FALSE,
                     max_results = 3,
                     verbose = TRUE) {

  data <- build_relatedness_data(kinship_df, info_df,
                                  from_col, to_col, kin_col,
                                  id_col, sex_col, ir_col)
  id  <- data$id
  sex <- data$sex
  ir  <- data$ir
  rel <- data$rel

  attempt <- 0
  repeat {
    attempt <- attempt + 1

    females <- which(sex == "F" & ir <= max_ir_f)
    males   <- which(sex == "M" & ir <= max_ir_m)

    if (verbose) {
      message(sprintf(
        "[attempt %d] %d females / %d males pass internal-relatedness filter",
        attempt, length(females), length(males)))
    }

    female_combos <- .valid_same_sex_combos(females, n_females, rel,
                                             max_pair_rel_f)
    male_combos   <- .valid_same_sex_combos(males, n_males, rel,
                                             max_pair_rel_m)

    if (verbose) {
      message(sprintf("  %d valid female combos, %d valid male combos",
                       length(female_combos), length(male_combos)))
    }

    candidates <- .build_family_candidates(female_combos, male_combos, rel,
                                            max_pair_rel, max_group_avg_rel,
                                            max_group_couple_avg_rel)

    if (verbose) {
      message(sprintf("  %d candidate family groups pass group-level filters",
                       length(candidates)))
    }

    if (length(candidates) >= n_groups) {
      results <- .best_group_combinations(candidates, n_groups,
                                           max_combo_avg_rel, max_results)
      if (length(results) > 0) break
    } else {
      results <- list()
    }

    if (!auto_tune || attempt >= 3) break

    # relax thresholds, mirroring the Perl script's auto-tune behaviour
    max_ir_f <- max_ir_f + 0.01
    max_ir_m <- max_ir_m + 0.01
    max_pair_rel_f <- max_pair_rel_f + 0.01
    max_pair_rel_m <- max_pair_rel_m + 0.01
    max_pair_rel   <- max_pair_rel + 0.01
    max_group_avg_rel <- max_group_avg_rel + 0.01
    max_group_couple_avg_rel <- max_group_couple_avg_rel + 0.01
    max_combo_avg_rel <- max_combo_avg_rel + 0.01
    if (verbose) message("  No solution found - relaxing thresholds by 0.01")
  }

  if (length(results) == 0) {
    warning("No valid combination of ", n_groups,
            " family group(s) found. Consider relaxing thresholds or ",
            "setting auto_tune = TRUE.")
    return(invisible(NULL))
  }

  # attach readable IDs
  for (r in seq_along(results)) {
    for (g in seq_along(results[[r]]$groups)) {
      results[[r]]$groups[[g]]$ids <- id[results[[r]]$groups[[g]]$members]
    }
  }

  results
}

# ---- 7. Pretty-print helper ------------------------------------------------

#' Print the results of swingeR() in a readable format
#'
#' @param results Output of \code{\link{swingeR}}.
#' @return Invisibly \code{NULL}; called for its printed output.
#' @export
print_breeding_groups <- function(results) {
  for (r in seq_along(results)) {
    cat(sprintf("\n=== Combination %d (avg relatedness = %.4f) ===\n",
                r, results[[r]]$combo_avg))
    for (g in seq_along(results[[r]]$groups)) {
      grp <- results[[r]]$groups[[g]]
      cat(sprintf("  Group %d: %s\n", g, paste(grp$ids, collapse = ", ")))
      cat(sprintf("    group avg rel = %.4f | couples avg rel = %.4f\n",
                   grp$group_avg, grp$couple_avg))
    }
  }
}

# ---- 8. Example usage -------------------------------------------------------
# kinship_df:
#        from        to     kinship
#   1  SPBR_002  SPBR_001  0.207414
#   2  SPBR_003  SPBR_001  0.010021
#   ...
#
# info_df:
#        ID        Sex   IR
#   1  SPBR_001    M    0.012010
#   2  SPBR_002    F    0.030500
#   ...
#
# results <- swingeR(
#   kinship_df = kinship_df, info_df = info_df,
#   n_groups = 2, n_females = 2, n_males = 1,
#   max_ir_f = 0.1, max_ir_m = 0.1,
#   max_pair_rel_f = 0.2, max_pair_rel_m = 0.2, max_pair_rel = 0.1,
#   max_group_avg_rel = 0.15, max_group_couple_avg_rel = 0.1,
#   max_combo_avg_rel = 0.15,
#   auto_tune = TRUE, max_results = 3
# )
# print_breeding_groups(results)
