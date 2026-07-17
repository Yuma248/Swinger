# =============================================================================
# vcf_tools.R
#
# Helper functions for summarizing SNP- and individual-level metrics from
# VCF files (via vcfR) and for converting a VCF to GDS format (via SNPRelate)
# for downstream population genetics analyses (e.g. with select_breeding_groups()).
# =============================================================================

#' Safely coerce a vector to numeric
#'
#' Wraps \code{as.numeric()} and silences the coercion warning that appears
#' when non-numeric strings are turned into \code{NA}. Used internally to
#' coerce VCF fields like QUAL, which are read in as character.
#'
#' @param x A vector (typically character) to coerce.
#' @return A numeric vector, with non-numeric entries as \code{NA}.
#' @keywords internal
safe_num <- function(x) suppressWarnings(as.numeric(x))

#' Summarize per-SNP metrics from a VCF
#'
#' Builds a data frame with one row per variant, combining the VCF's fixed
#' fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) with derived
#' per-locus metrics: numeric QUAL, total read depth across samples, number
#' of non-missing genotype calls, and percent missingness.
#'
#' @param vcf A \code{vcfR} object, as returned by \code{vcfR::read.vcfR()}.
#'
#' @return A data frame with one row per SNP/variant and columns:
#'   \describe{
#'     \item{CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO}{the VCF fixed
#'       fields, as character (from \code{vcf@fix})}
#'     \item{QUAL_num}{QUAL coerced to numeric}
#'     \item{DP_raw}{summed read depth across all samples at that locus}
#'     \item{n_genotypes}{number of samples with a non-missing depth/genotype
#'       call at that locus}
#'     \item{missing}{percent of samples with a missing genotype
#'       (\code{NA}, \code{"./."}, or \code{".|."}) at that locus}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' library(vcfR)
#' vcf <- read.vcfR("SPBRsimulated.vcf.gz")
#' snp_df <- vcf2df(vcf)
#' head(snp_df)
#' }
vcf2df <- function(vcf) {
  ourSNPs <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)

  ourSNPs$QUAL_num <- safe_num(ourSNPs$QUAL)

  dp_mat <- vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE)
  ourSNPs$DP_raw      <- rowSums(dp_mat, na.rm = TRUE)
  ourSNPs$n_genotypes <- rowSums(!is.na(dp_mat))

  gt_mat <- vcfR::extract.gt(vcf, element = "GT")
  ourSNPs$missing <- rowMeans(is.na(gt_mat) | gt_mat %in% c("./.", ".|.")) * 100

  return(ourSNPs)
}

#' Summarize per-individual metrics from a VCF
#'
#' Builds a data frame with one row per sample/individual, giving percent
#' missingness and observed heterozygosity for each.
#'
#' @param vcf A \code{vcfR} object, as returned by \code{vcfR::read.vcfR()}.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{ID}{sample identifier (from the VCF genotype matrix column names)}
#'     \item{miss_pct}{percent of loci with a missing genotype call for that
#'       individual}
#'     \item{Ho}{observed heterozygosity: proportion of that individual's
#'       non-missing genotypes that are heterozygous, from
#'       \code{\link{calc_ind_het}}}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' library(vcfR)
#' vcf <- read.vcfR("SPBRsimulated.vcf.gz")
#' ind_df <- vcf2inddf(vcf)
#' head(ind_df)
#' }
vcf2inddf <- function(vcf) {
  gt_mat <- vcfR::extract.gt(vcf, element = "GT")

  miss_pct <- colMeans(is.na(gt_mat) | gt_mat %in% c("./.", ".|.")) * 100

  ourInd <- data.frame(ID = colnames(gt_mat), miss_pct = miss_pct,
                        row.names = NULL)

  ind_het <- calc_ind_het(vcf, plot = FALSE)
  ourInd$Ho <- as.numeric(ind_het[ourInd$ID])

  return(ourInd)
}

#' Convert a VCF file to GDS format and open it
#'
#' Thin wrapper around \code{SNPRelate::snpgdsVCF2GDS()} +
#' \code{SNPRelate::snpgdsOpen()}: converts a (bgzipped) VCF to a GDS file
#' on disk, keeping biallelic SNPs only, and returns the open GDS connection.
#'
#' @param vcffile Path to the input VCF (or VCF.gz) file.
#' @param gdsfile Path to write the output GDS file to.
#'
#' @return An open \code{gds.class}/\code{gdsn.class} connection object (as
#'   returned by \code{SNPRelate::snpgdsOpen()}), ready for use with
#'   \code{SNPRelate}/\code{gdsfmt} functions.
#'
#' @note The returned connection is an external pointer and will NOT survive
#'   saving/reloading your R workspace (e.g. via \code{.RData} or restarting
#'   RStudio with "restore workspace" on). If you get
#'   \code{Error in print.gdsn.class... Invalid GDS node object}, just
#'   reopen it: \code{genofile <- SNPRelate::snpgdsOpen(gdsfile)}. Close it
#'   explicitly with \code{SNPRelate::snpgdsClose(genofile)} when done, and
#'   avoid opening the same GDS file twice without closing the first
#'   connection.
#' @export
#'
#' @examples
#' \dontrun{
#' genofile <- vcf2gds("SPBRsimulated.vcf.gz", "simulated.gds")
#' SNPRelate::snpgdsSummary("simulated.gds")
#' SNPRelate::snpgdsClose(genofile)
#' }
vcf2gds <- function(vcffile, gdsfile) {
  SNPRelate::snpgdsVCF2GDS(vcffile, gdsfile, method = "biallelic.only")
  genofile <- SNPRelate::snpgdsOpen(gdsfile)
  return(genofile)
}

#' Calculate observed heterozygosity per individual from a VCF
#'
#' For each sample, computes the proportion of its non-missing genotype
#' calls that are heterozygous (i.e. the two alleles differ).
#'
#' @param vcf A \code{vcfR} object, as returned by \code{vcfR::read.vcfR()}.
#' @param plot Currently unused (reserved for a future plotting option).
#'   Default \code{TRUE}.
#' @param sep Regular expression used to split genotype strings into
#'   alleles. Default \code{"[/|]"} handles both unphased (\code{"0/1"}) and
#'   phased (\code{"0|1"}) genotype notation.
#'
#' @return A named numeric vector (invisibly), one value per sample, giving
#'   the proportion of that sample's non-missing loci that are
#'   heterozygous. Also prints a dimension summary and a
#'   \code{summary()} of the heterozygosity values as a side effect.
#' @export
#'
#' @examples
#' \dontrun{
#' library(vcfR)
#' vcf <- read.vcfR("SPBRsimulated.vcf.gz")
#' het <- calc_ind_het(vcf)
#' }
calc_ind_het <- function(vcf, plot = TRUE, sep = "[/|]") {

  # Step 1: extract genotypes
  gt_raw <- vcfR::extract.gt(vcf, element = "GT")
  cat("Genotype matrix dimensions:\n")
  cat("  Loci   :", nrow(gt_raw), "\n")
  cat("  Samples:", ncol(gt_raw), "\n\n")

  # Step 2: for each locus x sample, check if the two alleles differ
  # A genotype like "0|1" or "0/1" is heterozygous  (alleles differ)
  # A genotype like "0|0" or "1/1" is homozygous    (alleles same)
  # A genotype like "./." or ".|." is missing        (return NA)
  is_het <- apply(gt_raw, 2, function(sample_gts) {
    alleles <- strsplit(sample_gts, sep)
    sapply(alleles, function(a) {
      if (any(is.na(a)) || length(a) < 2) return(NA)
      if (any(a == "."))                  return(NA)
      a[1] != a[2]
    })
  })

  # Step 3: mean proportion of heterozygous loci per individual
  ind_het <- colMeans(is_het, na.rm = TRUE)
  names(ind_het) <- colnames(gt_raw)

  # Step 4: print summary
  cat("Individual heterozygosity summary:\n")
  print(round(summary(ind_het), 4))
  cat("\n")

  # Step 5: optional plot (not yet implemented; reserved via the `plot` arg)

  # Step 6: return the vector invisibly so it can be saved
  return(invisible(ind_het))
}


#' Calculate pairwise relatedness from a VCF
#'
#' For each pair, computes the relatedness Queller & Goodnigth (1989) estimator.
#'
#' @param vcf A \code{vcfR} object, as returned by \code{vcfR::read.vcfR()}.
#' @param plot Currently unused (reserved for a future plotting option).
#'   Default \code{TRUE}.
#' @param sep Regular expression used to split genotype strings into
#'   alleles. Default \code{"[/|]"} handles both unphased (\code{"0/1"}) and
#'   phased (\code{"0|1"}) genotype notation.
#'
#' @return A square matrix of pairwise relatedness.
#'   \code{summary()} of the heterozygosity values as a side effect.
#' @export
#'
#' @examples
#' \dontrun{
#' library(vcfR)
#' vcf <- read.vcfR("SPBRsimulated.vcf.gz")
#' QG<-QGrelt(vcf)
#' }


QGrel <- function(vcf_obj) {
  # 1. Extract the genotype matrix (GT element) from the vcfR object
  message("Extracting genotype matrix...")
  gt_matrix <- vcfR::extract.gt(vcf_obj, element = "GT")
  
  # 2. Convert text genotypes (e.g., "0/0", "0/1", "1/1") into numeric allele dosages (0, 1, 2)
  # Missing data "." or "./." becomes NA
  message("Converting genotypes to allele counts...")
  numeric_gt <- matrix(NA, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix))
  colnames(numeric_gt) <- colnames(gt_matrix)
  rownames(numeric_gt) <- rownames(gt_matrix)
  
  numeric_gt[gt_matrix == "0/0" | gt_matrix == "0|0"] <- 0
  numeric_gt[gt_matrix == "0/1" | gt_matrix == "0|1" | gt_matrix == "1/0" | gt_matrix == "1|0"] <- 1
  numeric_gt[gt_matrix == "1/1" | gt_matrix == "1|1"] <- 2
  
  # 3. Calculate reference population allele frequency (p) for each locus
  # p is the frequency of the alternate allele '1'
  message("Calculating allele frequencies...")
  p <- rowMeans(numeric_gt, na.rm = TRUE) / 2
  
  # Filter out invariant/fixed SNPs where p is 0 or 1 to prevent division by zero
  valid_loci <- which(p > 0 & p < 1)
  if (length(valid_loci) == 0) {
    stop("Error: No polymorphic SNPs found in your dataset.")
  }
  
  numeric_gt <- numeric_gt[valid_loci, ]
  p <- p[valid_loci]
  
  num_ind <- ncol(numeric_gt)
  ind_names <- colnames(numeric_gt)
  
  # Initialize the final symmetric relationship matrix
  qg_matrix <- matrix(NA, nrow = num_ind, ncol = num_ind)
  colnames(qg_matrix) <- ind_names
  rownames(qg_matrix) <- ind_names
  
  message("Calculating pairwise Queller & Goodnight estimates...")
  
  # 4. Pairwise loop to calculate relatedness
  for (i in 1:num_ind) {
    for (j in i:num_ind) {
      if (i == j) {
        qg_matrix[i, j] <- 1.0 # Self-comparison anchor
        next
      }
      
      # Extract genotypes for individuals X and Y
      X <- numeric_gt[, i]
      Y <- numeric_gt[, j]
      
      # Retain only loci where both individuals have data (shared non-missing loci)
      complete_cases <- !is.na(X) & !is.na(Y)
      
      if (sum(complete_cases) == 0) {
        qg_matrix[i, j] <- NA
        qg_matrix[j, i] <- NA
        next
      }
      
      X_c <- X[complete_cases]
      Y_c <- Y[complete_cases]
      p_c <- p[complete_cases]
      
      # Numerator and Denominators for directional evaluations
      # Using 2*p because the genotype data is scaled to 0, 1, 2 alleles
      num_xy <- (X_c - 2 * p_c) * (Y_c - 2 * p_c)
      den_x  <- 2 * p_c * (1 - p_c)
      
      # Directional calculations (X as reference vs Y as reference)
      r_xy <- sum(num_xy) / sum(den_x)
      
      # Average both directions to yield the standard symmetrical QG estimate
      qg_matrix[i, j] <- r_xy
      qg_matrix[j, i] <- r_xy
    }
  }
  
  message("Done!")
  return(qg_matrix)
}

