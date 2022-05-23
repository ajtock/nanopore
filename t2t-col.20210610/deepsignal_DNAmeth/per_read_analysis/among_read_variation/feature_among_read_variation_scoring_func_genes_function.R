#!/applications/R/R-4.0.0/bin/Rscript

# For each feat:
# 1. calculate a measure of among-read agreement in methylation state (e.g., Fleiss' kappa)

# Get reads that overlap each featName
fOverlapsStrand <- function(chr_featGR, chr_tabGR_str) {
  ## Note: findOverlaps() approach does not work where a window does not overlap
  ##       any positions in chr_tabGR, which can occur with smaller genomeBinSize
  # Identify overlapping windows and midpoint coordinates
  fOverlaps_str <- findOverlaps(query = chr_featGR,
                                subject = chr_tabGR_str,
                                type = "any",
                                select = "all",
                                ignore.strand = T)
  fOverlaps_str
}

# Function to calculate among-read agreement for a given feature x
makeDFx_strand <- function(fOverlaps_str, chr_tabGR_str, chr_featGR, featNum) {

  chr_tabGR_str_x <- chr_tabGR_str[subjectHits(fOverlaps_str[queryHits(fOverlaps_str) == featNum])]

  if(length(chr_tabGR_str_x) > 0) {

    chr_tabGR_str_x <- sortSeqlevels(chr_tabGR_str_x)
    chr_tabGR_str_x <- sort(chr_tabGR_str_x, by = ~ read + start)

    df_str_x <- data.frame(pos = start(chr_tabGR_str_x),
                           read = chr_tabGR_str_x$read,
                           call = chr_tabGR_str_x$call)

    pwider_str_x <- as.data.frame(tidyr::pivot_wider(data = df_str_x,
                                                     names_from = read,
#                                                     names_prefix = "read_",
                                                     values_from = call))
    pwider_str_x <- pwider_str_x[ with(data = pwider_str_x, expr = order(pos)), ]
    rownames(pwider_str_x) <- pwider_str_x[,1]
    pwider_str_x <- pwider_str_x[ , -1, drop = F]

    # kappam.fleiss() uses only rows (cytosines) with complete information
    # across all columns (reads)
    # Therefore, remove columns (reads) containing >= NAmax proportion NAs to
    # to retain more cytosines in the data.frame for kappa calculation
    mask_cols <- apply(pwider_str_x, MARGIN = 2, FUN = function(col) sum(is.na(col)) >= nrow(pwider_str_x) * NAmax)    
    # Report proportion of columns (reads) to be retained
    prop_reads_retained_str_x <- sum(!(mask_cols)) / ncol(pwider_str_x)
    # Conditionally remove columns (reads) containing >= NAmax proportion NAs
    if(sum(mask_cols) > 0) {
      pwider_str_x <- pwider_str_x[ , !(mask_cols), drop = F]
    }
    # Report number of columns (reads) to be retained for kappa and other calculations
    num_reads_retained_str_x <- ncol(pwider_str_x) 

    # Identify rows (cytosines) containing any NAs across the retained columns (reads),
    # as these will not be used by kappam.fleiss()
    mask_rows <- apply(pwider_str_x, MARGIN = 1, FUN = function(row) sum(is.na(row)) > 0)
    # Report proportion of rows (cytosines) to be retained
    prop_Cs_retained_str_x <- sum(!(mask_rows)) / nrow(pwider_str_x) 
    # Keep rows (cytosines) with NAs for other calculations
    stocha_pwider_str_x <- pwider_str_x
    # Conditionally remove rows (cytosines) containing any NAs
    if(sum(mask_rows) > 0) {
      pwider_str_x <- pwider_str_x[ !(mask_rows), , drop = F]
    }
    # Report number of rows (cytosines) to be retained for kappa and other calculations
    num_Cs_retained_str_x <- nrow(pwider_str_x)
    stocha_num_Cs_retained_str_x <- nrow(stocha_pwider_str_x)

    # Calculate mean methylation for region
    mean_mC_pwider_str_x <- mean(as.matrix(stocha_pwider_str_x), na.rm = T)
    
    # Calculate mean and standard deviation of per-read mean methylation
    mean_read_mC_pwider_str_x <- mean(colMeans(stocha_pwider_str_x, na.rm = T), na.rm = T)
    sd_read_mC_pwider_str_x <- sd(colMeans(stocha_pwider_str_x, na.rm = T), na.rm = T)

    # Calculate Fleiss' kappa
    if(nrow(pwider_str_x) >= min_Cs && nrow(pwider_str_x) <= max_Cs &&
       ncol(pwider_str_x) >= min_reads && ncol(pwider_str_x) <= max_reads) {

      # Calculate Fleiss' kappa
      fkappa_pwider_str_x <- irr::kappam.fleiss(pwider_str_x, detail = F)

      # Sanity checks
      stopifnot(fkappa_pwider_str_x$raters == num_reads_retained_str_x)
      stopifnot(fkappa_pwider_str_x$subjects == num_Cs_retained_str_x)

      fkappa_pwider_str_x_kappa <- fkappa_pwider_str_x$value
      fkappa_pwider_str_x_pval <- fkappa_pwider_str_x$p.value
      fkappa_pwider_str_x_zstat <- fkappa_pwider_str_x$statistic
      fkappa_pwider_str_x_reads <- fkappa_pwider_str_x$raters
      fkappa_pwider_str_x_Cs <- fkappa_pwider_str_x$subjects

    } else {

      fkappa_pwider_str_x_kappa <- NaN
      fkappa_pwider_str_x_pval <- NaN
      fkappa_pwider_str_x_zstat <- NaN
      fkappa_pwider_str_x_reads <- NaN
      fkappa_pwider_str_x_Cs <- NaN

    }

    # Calculate Krippendorff's alpha, an inter-rater reliability or agreement metric that can handle incomplete data,
    # and for which "computed reliabilities are comparable across any numbers of coders [raters], values, ... and unequal sample sizes.";
    # https://en.wikipedia.org/wiki/Krippendorff%27s_alpha
    # Also calculate site-to-site stochasticity
    if(nrow(stocha_pwider_str_x) >= min_Cs && nrow(stocha_pwider_str_x) <= max_Cs &&
       ncol(stocha_pwider_str_x) >= min_reads && ncol(stocha_pwider_str_x) <= max_reads) {

      # Calculate Krippendorff's alpha
      kalpha_pwider_str_x <- irr::kripp.alpha(t(stocha_pwider_str_x), method = "nominal")

      # Sanity checks
      stopifnot(kalpha_pwider_str_x$raters == num_reads_retained_str_x)
      stopifnot(kalpha_pwider_str_x$subjects == stocha_num_Cs_retained_str_x)

      kalpha_pwider_str_x_alpha <- kalpha_pwider_str_x$value
      kalpha_pwider_str_x_nmatchval <- kalpha_pwider_str_x$nmatchval


      # Calculate absolute differences between methylation statuses of neighbouring Cs within each read 
      absdiff_pwider_str_x <- abs(diff(as.matrix(stocha_pwider_str_x)))
      # Calculate the mean absolute difference for each read
      colMeans_absdiff_pwider_str_x <- colMeans(absdiff_pwider_str_x, na.rm = T)
      # Across all reads overlapping a given feature, calculate the mean and median of mean absolute differences
      mean_stocha_pwider_str_x <- mean(colMeans_absdiff_pwider_str_x, na.rm = T)
      median_stocha_pwider_str_x <- median(colMeans_absdiff_pwider_str_x, na.rm = T)
      # Across all reads overlapping a given feature, calculate the sd of mean absolute differences
      sd_stocha_pwider_str_x <- sd(colMeans_absdiff_pwider_str_x, na.rm = T)

      # Report number of rows (cytosines) retained for other calculations
      stocha_pwider_str_x_Cs <- nrow(stocha_pwider_str_x)


      # Calculate autocorrelations between methylation statuses of neighbouring Cs within each read
      acf_pwider_str_x_list <- apply(stocha_pwider_str_x, MARGIN = 2,
                                     FUN = function(col) acf(col, lag.max = 10, plot = F, na.action = na.pass))
      mean_min_acf_pwider_str_x <- mean(sapply(seq_along(acf_pwider_str_x_list), function(col) {
        if(sum(acf_pwider_str_x_list[[col]]$acf, na.rm = T) != 0) {
          min(as.vector(acf_pwider_str_x_list[[col]]$acf)[-1], na.rm = T)
        } else {
          NA
        }
      }), na.rm = T)
      mean_max_acf_pwider_str_x <- mean(sapply(seq_along(acf_pwider_str_x_list), function(col) {
        if(sum(acf_pwider_str_x_list[[col]]$acf, na.rm = T) != 0) {
          max(as.vector(acf_pwider_str_x_list[[col]]$acf)[-1], na.rm = T)
        } else {
          NA
        }
      }), na.rm = T)
      mean_mean_acf_pwider_str_x <- mean(sapply(seq_along(acf_pwider_str_x_list), function(col) {
        if(sum(acf_pwider_str_x_list[[col]]$acf, na.rm = T) != 0) {
          mean(as.vector(acf_pwider_str_x_list[[col]]$acf)[-1], na.rm = T)
        } else {
          NA
        }
      }), na.rm = T)
      mean_median_acf_pwider_str_x <- mean(sapply(seq_along(acf_pwider_str_x_list), function(col) {
        if(sum(acf_pwider_str_x_list[[col]]$acf, na.rm = T) != 0) {
          median(as.vector(acf_pwider_str_x_list[[col]]$acf)[-1], na.rm = T)
        } else {
          NA
        }
      }), na.rm = T)

    } else {

      kalpha_pwider_str_x_alpha <- NaN
      kalpha_pwider_str_x_nmatchval <- NaN
      mean_stocha_pwider_str_x <- NaN
      median_stocha_pwider_str_x <- NaN
      sd_stocha_pwider_str_x <- NaN
      stocha_pwider_str_x_Cs <- NaN
      mean_min_acf_pwider_str_x <- NaN
      mean_max_acf_pwider_str_x <- NaN
      mean_mean_acf_pwider_str_x <- NaN
      mean_median_acf_pwider_str_x <- NaN

    }

    fk_df_str_win_x <- data.frame(chr = seqnames(chr_featGR[featNum]),
                                  start = start(chr_featGR[featNum]),
                                  end = end(chr_featGR[featNum]),
                                  midpoint = round((start(chr_featGR[featNum])+end(chr_featGR[featNum]))/2),
                                  strand = strand(chr_featGR[featNum]),
                                  name = chr_featGR[featNum]$name,
                                  score = chr_featGR[featNum]$score,
                                  feature_width = chr_featGR[featNum]$feature_width,
                                  exons_width = chr_featGR[featNum]$exons_width,
                                  exons_width_prop = chr_featGR[featNum]$exons_width_prop,
                                  exons_count = chr_featGR[featNum]$exons_count,
                                  exons_count_per_kb = chr_featGR[featNum]$exons_count_per_kb,
                                  introns_width = chr_featGR[featNum]$introns_width,
                                  introns_width_prop = chr_featGR[featNum]$introns_width_prop,
                                  introns_count = chr_featGR[featNum]$introns_count,
                                  introns_count_per_kb = chr_featGR[featNum]$introns_count_per_kb,

                                  mean_mC_str = mean_mC_pwider_str_x,
                                  mean_read_mC_str = mean_read_mC_pwider_str_x,
                                  sd_read_mC_str = sd_read_mC_pwider_str_x,

                                  fk_kappa_str = fkappa_pwider_str_x_kappa,
                                  fk_pval_str = fkappa_pwider_str_x_pval,
                                  fk_zstat_str = fkappa_pwider_str_x_zstat,
                                  fk_reads_str = fkappa_pwider_str_x_reads,
                                  fk_Cs_str = fkappa_pwider_str_x_Cs,

                                  ka_alpha_str = kalpha_pwider_str_x_alpha,
                                  ka_nmatchval_str = kalpha_pwider_str_x_nmatchval,
                                  mean_stocha_str = mean_stocha_pwider_str_x,
                                  median_stocha_str = median_stocha_pwider_str_x,
                                  sd_stocha_str = sd_stocha_pwider_str_x,
                                  stocha_Cs_str = stocha_pwider_str_x_Cs,
                                  mean_min_acf_str = mean_min_acf_pwider_str_x,
                                  mean_max_acf_str = mean_max_acf_pwider_str_x,
                                  mean_mean_acf_str = mean_mean_acf_pwider_str_x,
                                  mean_median_acf_str = mean_median_acf_pwider_str_x
                                 ) 

  } else {

    fk_df_str_win_x <- data.frame(chr = seqnames(chr_featGR[featNum]),
                                  start = start(chr_featGR[featNum]),
                                  end = end(chr_featGR[featNum]),
                                  midpoint = round((start(chr_featGR[featNum])+end(chr_featGR[featNum]))/2),
                                  strand = strand(chr_featGR[featNum]),
                                  name = chr_featGR[featNum]$name,
                                  score = chr_featGR[featNum]$score,
                                  feature_width = chr_featGR[featNum]$feature_width,
                                  exons_width = chr_featGR[featNum]$exons_width,
                                  exons_width_prop = chr_featGR[featNum]$exons_width_prop,
                                  exons_count = chr_featGR[featNum]$exons_count,
                                  exons_count_per_kb = chr_featGR[featNum]$exons_count_per_kb,
                                  introns_width = chr_featGR[featNum]$introns_width,
                                  introns_width_prop = chr_featGR[featNum]$introns_width_prop,
                                  introns_count = chr_featGR[featNum]$introns_count,
                                  introns_count_per_kb = chr_featGR[featNum]$introns_count_per_kb,

                                  mean_mC_str = NaN,
                                  mean_read_mC_str = NaN,
                                  sd_read_mC_str = NaN,

                                  fk_kappa_str = NaN,
                                  fk_pval_str = NaN,
                                  fk_zstat_str = NaN,
                                  fk_reads_str = NaN,
                                  fk_Cs_str = NaN,

                                  ka_alpha_str = NaN,
                                  ka_nmatchval_str = NaN,
                                  mean_stocha_str = NaN,
                                  median_stocha_str = NaN,
                                  sd_stocha_str = NaN,
                                  stocha_Cs_str = NaN,
                                  mean_min_acf_str = NaN,
                                  mean_max_acf_str = NaN,
                                  mean_mean_acf_str = NaN,
                                  mean_median_acf_str = NaN
                                 ) 

  }

fk_df_str_win_x

}

