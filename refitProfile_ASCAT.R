merge_consecutive_positions <- function(segments) {
  ## merge consecutive rounded segments that share chr/nMajor/nMinor and are contiguous
  segments_df <- data.frame(segments)
  n <- nrow(segments_df)
  
  if (n == 0) return(segments)
  
  new_segment <- c(
    TRUE,
    !(segments_df$chr[-1]     == segments_df$chr[-n] &
        segments_df$nMajor[-1]  == segments_df$nMajor[-n] &
        segments_df$nMinor[-1]  == segments_df$nMinor[-n] &
        segments_df$startpos[-1] == segments_df$endpos[-n] + 1)
  )
  
  segments_df$group <- cumsum(new_segment)
  
  merged_df <- segments_df %>%
    group_by(group) %>%
    summarise(
      startpos = first(startpos),
      endpos   = last(endpos),
      chr      = first(chr),
      nMajor   = first(nMajor),
      nMinor   = first(nMinor),
      .groups  = "drop"
    )
  
  merged_list <- list(
    startpos = merged_df$startpos,
    endpos   = merged_df$endpos,
    chr      = merged_df$chr,
    nMajor   = merged_df$nMajor,
    nMinor   = merged_df$nMinor
  )
  
  merged_list$sample <- rep(segments$sample[1], length(merged_list$startpos))
  merged_list <- merged_list[c("sample", "startpos", "endpos", "chr", "nMajor", "nMinor")]
  
  return(merged_list)
}

refitProfileASCAT <- function(ASCATobj, refMinor, refMajor, chr, gamma = 0.55, wm = NULL) {
  
  segments <- ASCATobj$segments_raw
  
  org_purity <- ASCATobj$purity[[1]]
  org_nA     <- segments[, "nAraw"]
  org_nB     <- segments[, "nBraw"]
  org_psi    <- ASCATobj$psi[[1]]     # sample-level normalizing ploidy: psi = rho*ploidy + 2*(1-rho)
  # org_ploidy (tumor-only ploidy, ASCATobj$ploidy[[1]]) is NOT used here â€” LogR is normalized by psi, not ploidy
  
   ### LogR = gamma * log2( (rho*(nA+nB) + (1-rho)*2) / psi )
  LogRref <- gamma * log2((org_purity * (org_nA + org_nB) + (1 - org_purity) * 2) / org_psi)
  
  ### BAF = (1 - rho + rho*nB) / (2 - 2*rho + rho*(nA+nB))
  refBAF <- (1 - org_purity + org_purity * org_nB) /
    (2 - 2 * org_purity + org_purity * (org_nA + org_nB))
  
  segments$logR <- LogRref
  segments$BAF  <- refBAF
  
  cn. <- segments[segments$"chr" == chr, ]
  sizes <- as.numeric(as.character(cn.$endpos)) - as.numeric(as.character(cn.$startpos))
  
  ### if index of segment on chromosome is not specified, pick the largest segment
  if (is.null(wm)) { wm <- which.max(sizes)[1] }
  
  logR_seg <- cn.[wm, "logR"]
  BAF_seg  <- cn.[wm, "BAF"]
  
  ### --- Solve BAF = (1-rho+rho*nB) / (2-2rho+rho*(refMajor+refMinor)) for rho ---
  ### rho = (2*BAF-1) / ((2*BAF-1) - BAF*(refMajor+refMinor) + refMinor)
  rho <- (2 * BAF_seg - 1) /
    ((2 * BAF_seg - 1) - BAF_seg * (refMajor + refMinor) + refMinor)
  
  ### --- Solve LogR = gamma*log2((rho*nt+2(1-rho))/psi) for the new psi (sample ploidy) ---
  psi  <- (2 * (1 - rho) + rho * (refMajor + refMinor)) / (2^(logR_seg / gamma))
  psit <- (psi - 2 * (1 - rho)) / rho                      # tumor-only ploidy
  
  if (!is.finite(psit) || psit <= 0 || !is.finite(rho) || rho <= 0 || rho > 1) {
    stop("Not possible: ploidy<0 or purity \u2209 [0,1]")
  }
  
  ### recompute allele-specific copy number for every segment under the new rho/psit,
  ### holding each segment's (logR, BAF) fixed ---
  newA <- (rho - 1 - (refBAF - 1) * 2^(LogRref / gamma) * ((1 - rho) * 2 + rho * psit)) / rho
  newB <- (rho - 1 +  refBAF       * 2^(LogRref / gamma) * ((1 - rho) * 2 + rho * psit)) / rho
  
  ASCATobj$segments_raw$nAraw  <- newA
  ASCATobj$segments_raw$nBraw  <- newB
  ASCATobj$segments_raw$nMajor <- round(newA, 0)
  ASCATobj$segments_raw$nMinor <- round(newB, 0)
  
  ### merge rounded segments with the same value
  ASCATobj$segments <- merge_consecutive_positions(ASCATobj$segments_raw)
  
  ASCATobj$purity[[1]] <- rho
  ASCATobj$ploidy[[1]] <- psit
  ASCATobj$psi[[1]]    <- psi
  
  return(ASCATobj)
}