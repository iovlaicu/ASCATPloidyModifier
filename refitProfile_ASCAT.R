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
  
  ## same layout as ASCAT's own $segments data frame
  merged <- data.frame(
    sample   = rep(segments_df$sample[1], nrow(merged_df)),
    chr      = merged_df$chr,
    startpos = merged_df$startpos,
    endpos   = merged_df$endpos,
    nMajor   = merged_df$nMajor,
    nMinor   = merged_df$nMinor,
    stringsAsFactors = FALSE
  )

  return(merged)
}

refitProfileASCAT <- function(ASCATobj, refMinor, refMajor, chr, gamma = 0.55, wm = NULL) {
  
  if (refMajor == refMinor) {
    stop("A balanced state (major = minor) carries no allelic imbalance, so it cannot be used to refit purity; choose an unbalanced state")
  }

  segments <- ASCATobj$segments_raw

  org_purity <- ASCATobj$purity[[1]]
  org_nA     <- segments[, "nAraw"]
  org_nB     <- segments[, "nBraw"]
  org_psi    <- ASCATobj$psi[[1]]     # ASCAT psi

  ### ASCAT model: LogR = gamma * log2( (2*(1-rho) + rho*(nA+nB)) / (2*(1-rho) + rho*psi) )
  LogRref <- gamma * log2((org_purity * (org_nA + org_nB) + (1 - org_purity) * 2) /
                            (2 * (1 - org_purity) + org_purity * org_psi))
  
  ### BAF = (1 - rho + rho*nB) / (2 - 2*rho + rho*(nA+nB))
  refBAF <- (1 - org_purity + org_purity * org_nB) /
    (2 - 2 * org_purity + org_purity * (org_nA + org_nB))
  
  segments$logR <- LogRref
  segments$BAF  <- refBAF
  
  cn. <- segments[gsub("^chr", "", segments$"chr") == gsub("^chr", "", as.character(chr)), ]
  if (nrow(cn.) == 0) stop("No segments found on chromosome ", chr)
  sizes <- as.numeric(as.character(cn.$endpos)) - as.numeric(as.character(cn.$startpos))
  
  ### if index of segment on chromosome is not specified, pick the largest segment
  if (is.null(wm)) { wm <- which.max(sizes)[1] }
  
  logR_seg <- cn.[wm, "logR"]
  BAF_seg  <- cn.[wm, "BAF"]
  
  ### --- Solve BAF = (1-rho+rho*nB) / (2-2rho+rho*(refMajor+refMinor)) for rho ---
  ### rho = (2*BAF-1) / ((2*BAF-1) - BAF*(refMajor+refMinor) + refMinor)
  rho <- (2 * BAF_seg - 1) /
    ((2 * BAF_seg - 1) - BAF_seg * (refMajor + refMinor) + refMinor)
  
  ### --- Solve the LogR equation above for the new psi ---
  normaliser <- (2 * (1 - rho) + rho * (refMajor + refMinor)) / (2^(logR_seg / gamma))   # = 2*(1-rho) + rho*psi
  psi <- (normaliser - 2 * (1 - rho)) / rho

  if (!is.finite(psi) || psi <= 0 || !is.finite(rho) || rho <= 0 || rho > 1) {
    stop("Not possible: ploidy<0 or purity \u2209 [0,1]")
  }

  ### recompute allele-specific copy number for every segment under the new rho/psi,
  ### holding each segment's (logR, BAF) fixed (same formulas as ASCAT) ---
  newA <- (rho - 1 - (refBAF - 1) * 2^(LogRref / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
  newB <- (rho - 1 +  refBAF       * 2^(LogRref / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
  
  ASCATobj$segments_raw$nAraw  <- newA
  ASCATobj$segments_raw$nBraw  <- newB
  ASCATobj$segments_raw$nMajor <- round(newA, 0)
  ASCATobj$segments_raw$nMinor <- round(newB, 0)
  
  ### merge rounded segments with the same value
  ASCATobj$segments <- merge_consecutive_positions(ASCATobj$segments_raw)
  
  ASCATobj$purity[[1]] <- rho
  ASCATobj$psi[[1]]    <- psi
  ### tumour ploidy as ASCAT reports it (length-weighted mean total copy number); the app
  ### replaces it with ASCAT's own probe-weighted value via updatePloidyGoF()
  len <- as.numeric(ASCATobj$segments_raw$endpos) - as.numeric(ASCATobj$segments_raw$startpos)
  ASCATobj$ploidy[[1]] <- sum((newA + newB) * len) / sum(len)
  
  return(ASCATobj)
}