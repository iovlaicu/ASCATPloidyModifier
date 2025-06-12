merge_consecutive_positions <- function(segments) {
  
  ## merge small rounded segments if same value
  segments_df <- data.frame(segments)
  

  segments_df <- segments_df %>%
    mutate(
      group = cumsum(
        c(TRUE, 
          diff(startpos) != (endpos[-length(endpos)] - startpos[-1]) |
            nMajor != lag(nMajor) |
            nMinor != lag(nMinor) |
            chr != lag(chr)
        )
      )
    )
  

  merged_df <- segments_df %>%
    group_by(group) %>%
    summarise(
      startpos = first(startpos),
      endpos = last(endpos),
      chr = first(chr),
      nMajor = first(nMajor),
      nMinor = first(nMinor),
      .groups = "drop"
    )
  

  merged_list <- list(
    startpos = merged_df$startpos,
    endpos = merged_df$endpos,
    chr = merged_df$chr,
    nMajor = merged_df$nMajor,
    nMinor = merged_df$nMinor
  )
  
  merged_list$sample <- rep(segments$sample[1], length(merged_list$startpos))
  
  merged_list <- merged_list[c("sample", "startpos", "endpos", "chr", "nMajor", "nMinor")]
  
  
  return(merged_list)
}


refitProfileASCAT <- function(ASCATobj, refMinor, refMajor, chr, gamma=0.55, wm=NULL) {
  
  segments <- ASCATobj$segments_raw
  
  org_purity <- ASCATobj$purity[[1]]
  org_nA <- segments[, "nAraw"]
  org_nB <- segments[, "nBraw"]
  org_psi <- ASCATobj$psi[[1]]
  org_ploidy <- ASCATobj$ploidy[[1]]
  
  ### calculate logR using original ploidy, purity, and unrounded nMajor & nMinor
  
  LogRref<- gamma * log((org_purity *(org_nA + org_nB) + (1 - org_purity)*2) / org_ploidy)
  
  ### calculate BAF using original purity, and unrounded nMajor & nMinor
  
  refBAF <- (1 - org_purity + org_purity * org_nB) / (2 - 2 * org_purity + org_purity*(org_nA + org_nB) )
  
  segments$logR <- LogRref
  segments$BAF <- refBAF
  
  cn. <- segments[segments$"chr"==chr, ]
  sizes <- as.numeric(as.character(cn.$endpos))-as.numeric(as.character(cn.$startpos))
  
  ### if index of segment on chromosome is not specified, pick the largest segment
  if(is.null(wm)) {wm <- which.max(sizes)[1]}
 
  ### get logR and BAF for the segment
  logR_seg <- cn.[wm, "logR"]
  BAF_seg <- cn.[wm, "BAF"]
  
  ### caluclate new purity based on BAF and new nMajor & nMinor
  rho=(2*BAF_seg-1) / ((2*BAF_seg-1) + BAF_seg*(refMinor+refMajor) + refMinor)
  #psi = (rho*(refMajor+refMinor)+2-2*rho)/(2^(logR_seg/gamma))
  
  ### caluclate new ploidy based on logR and new nMajor & nMinor
  
  psi=(2*(1-rho)+rho*(refMajor+refMinor))/(2^(logR_seg/gamma))
  psit = (psi-2*(1-rho))/rho
  
  if(psit == Inf | psit == -Inf | rho == 0 | rho>1 ) stop("Not possible: ploidy<0 or purity ∉ [0,1]")
  
  ### recalculate nMajor & nMinor values for all segments
  newA = (rho-1 - (refBAF-1)*2^(LogRref/gamma) * ((1-rho)*2+rho*psit))/rho
  newB = (rho-1 + refBAF*2^(LogRref/gamma) * ((1-rho)*2+rho*psit))/rho
  
  ASCATobj$segments_raw$nAraw <- newA
  ASCATobj$segments_raw$nBraw <- newB
  ASCATobj$segments_raw$nMajor <- round(newA, 0)
  ASCATobj$segments_raw$nMinor <- round(newB, 0)
  
  ### merge rounded segments with the same value
  ASCATobj$segments <- merge_consecutive_positions(ASCATobj$segments_raw)
  
  ASCATobj$purity[[1]] <- rho
  ASCATobj$ploidy[[1]] <- psit
  ASCATobj$psi[[1]] <- psi
  

  return (ASCATobj)
  
}




