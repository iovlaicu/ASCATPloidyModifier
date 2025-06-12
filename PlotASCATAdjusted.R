library(RColorBrewer)

ascat.plotAdjustedAscatProfile=function(ASCAT_output_object, REF, SAMPLE, y_limit=5, plot_unrounded=FALSE, hideSeg=FALSE, highlight=NULL) {
  
  nonrounded=ASCAT_output_object$segments_raw[, c(1:4, 7:8)]
  colnames(nonrounded)[5:6]=c("nMajor", "nMinor")
  nonrounded$nMajor=nonrounded$nMajor+nonrounded$nMinor
  colourAN = "lightpink" 
  colourBN = "lightskyblue" 
  
  nonrounded$nMajor=ifelse(nonrounded$nMajor>y_limit, y_limit+0.1, nonrounded$nMajor)
  nonrounded$nMinor=ifelse(nonrounded$nMinor>y_limit, y_limit+0.1, nonrounded$nMinor)
  
  
  SEGMENTS=ASCAT_output_object$segments
  SEGMENTS$nMajor=SEGMENTS$nMajor-0.1
  SEGMENTS$nMinor=SEGMENTS$nMinor+0.1
  colourA = "#E03546" # red
  colourB = "#3557E0" # blue
  
  colourH = "chartreuse"
  
  SEGMENTS$nMajor=ifelse(SEGMENTS$nMajor>y_limit, y_limit+0.1, SEGMENTS$nMajor)
  SEGMENTS$nMinor=ifelse(SEGMENTS$nMinor>y_limit, y_limit+0.1, SEGMENTS$nMinor)
  
  if (REF=="hg19") {
    REF=data.frame(chrom=c(1:22, "X"),
                   start=rep(1, 23),
                   end=c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
                         135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                         59128983, 63025520, 48129895, 51304566, 155270560))
  } else if (REF=="hg38") {
    REF=data.frame(chrom=c(1:22, "X"),
                   start=rep(1, 23),
                   end=c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
                         133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                         58617616, 64444167, 46709983, 50818468, 156040895))
  } else if (REF=="CHM13") {
    REF=data.frame(chrom=c(1:22, "X"),
                   start=rep(1, 23),
                   end=c(248387328, 242696752, 201105948, 193574945, 182045439, 172126628, 160567428, 146259331, 150617247,
                         134758134, 135127769, 133324548, 113566686, 101161492, 99753195, 96330374, 84276897, 80542538,
                         61707364, 66210255, 45090682, 51324926, 154259566))
  } else {
    stopifnot(is.data.frame(REF))
    stopifnot(identical(colnames(REF), c("chrom", "start", "end")))
  }
  
  SEGMENTS$chr=gsub("^chr", "", SEGMENTS$chr)
  nonrounded$chr=gsub("^chr", "", nonrounded$chr)
  
  #stopifnot(all(ASCAT_output_object$segments$chr %in% REF$chrom))
  REF$size=REF$end-REF$start+1
  REF$middle=0
  for (i in 1:nrow(REF)) {
    if (i==1) {
      REF$middle[i]=REF$size[i]/2
    } else {
      REF$middle[i]=sum(as.numeric(REF$size[1:(i-1)]))+REF$size[i]/2
    }
  }; rm(i)
  REF$cumul=cumsum(as.numeric(REF$size))
  REF$add=cumsum(as.numeric(c(0, REF$size[1:(nrow(REF)-1)])))
  
  SEGMENTS$startpos_adjusted=SEGMENTS$startpos
  SEGMENTS$endpos_adjusted=SEGMENTS$endpos
  
  nonrounded$startpos_adjusted=nonrounded$startpos
  nonrounded$endpos_adjusted=nonrounded$endpos
  
  for (CHR in unique(REF$chrom)) {
    INDEX=which(SEGMENTS$chr==CHR)
    if (length(INDEX)>0) {
      SEGMENTS$startpos_adjusted[INDEX]=SEGMENTS$startpos_adjusted[INDEX]+REF$add[which(REF$chrom==CHR)]
      SEGMENTS$endpos_adjusted[INDEX]=SEGMENTS$endpos_adjusted[INDEX]+REF$add[which(REF$chrom==CHR)]
      
      
    }
    rm(INDEX)
  }; rm(CHR)
  
  for (CHR in unique(REF$chrom)) {
    INDEX=which(nonrounded$chr==CHR)
    if (length(INDEX)>0) {
      
      nonrounded$startpos_adjusted[INDEX]=nonrounded$startpos_adjusted[INDEX]+REF$add[which(REF$chrom==CHR)]
      nonrounded$endpos_adjusted[INDEX]=nonrounded$endpos_adjusted[INDEX]+REF$add[which(REF$chrom==CHR)]
    }
    rm(INDEX)
  }; rm(CHR)
  
  
  SEGS=SEGMENTS#[which(SEGMENTS$sample==SAMPLE), ]
  SEGS_nonrounded=nonrounded #[which(nonrounded$sample==SAMPLE), ]
  #if (nrow(SEGS)==0) warning(paste0("No segments for sample: ", SAMPLE))
  maintitle = paste("Ploidy: ", sprintf("%1.2f", ASCAT_output_object$ploidy[1]), ", purity: ", sprintf("%2.0f", ASCAT_output_object$purity[1]*100), "%, goodness of fit: ", sprintf("%2.1f", ASCAT_output_object$goodnessOfFit[1]), "%", ifelse(isTRUE(ASCAT_output_object$nonaberrantarrays[1]), ", non-aberrant", ""), sep="")
  #png(filename = paste0(png_prefix, SAMPLE, ".adjusted", ifelse(plot_unrounded, "rawprofile", "ASCATprofile"), ".png"), width = 2000, height = (y_limit*100), res = 200)
  par(mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
  ticks=seq(0, y_limit, 1)
  plot(c(1, REF$cumul[nrow(REF)]), c(0, y_limit), type = "n", xaxt = "n", yaxt="n", main = maintitle, xlab = "", ylab = "")
  axis(side = 2, at = ticks)
  abline(h=ticks, col="lightgrey", lty=1)
  
  
  #segments(SEGS_nonrounded$startpos_adjusted, (SEGS_nonrounded$nMajor-0.07), SEGS_nonrounded$endpos_adjusted, (SEGS_nonrounded$nMajor+0.07), col=colourAN)
  #segments(SEGS_nonrounded$startpos_adjusted, (SEGS_nonrounded$nMinor-0.07), SEGS_nonrounded$endpos_adjusted, (SEGS_nonrounded$nMinor+0.07), col=colourBN)
  if(plot_unrounded){
    
    rect(SEGS_nonrounded$startpos_adjusted, (SEGS_nonrounded$nMajor-0.07), SEGS_nonrounded$endpos_adjusted, (SEGS_nonrounded$nMajor+0.07), col=ifelse(SEGS_nonrounded$nMajor>=y_limit, adjustcolor(colourAN, red.f=0.75, green.f=0.75, blue.f=0.75, red.f=0.75, green.f=0.75, blue.f=0.75), colourAN), border=ifelse(SEGS_nonrounded$nMajor>=y_limit, adjustcolor(colourAN, red.f=0.75, green.f=0.75, blue.f=0.75), colourAN))
    rect(SEGS_nonrounded$startpos_adjusted, (SEGS_nonrounded$nMinor-0.07), SEGS_nonrounded$endpos_adjusted, (SEGS_nonrounded$nMinor+0.07), col=ifelse(SEGS_nonrounded$nMinor>=y_limit, adjustcolor(colourBN, red.f=0.75, green.f=0.75, blue.f=0.75), colourBN), border=ifelse(SEGS_nonrounded$nMinor>=y_limit, adjustcolor(colourBN, red.f=0.75, green.f=0.75, blue.f=0.75), colourBN))
    
  }
  if(!hideSeg){
    rect(SEGS$startpos_adjusted, (SEGS$nMajor-0.07), SEGS$endpos_adjusted, (SEGS$nMajor+0.07), col=ifelse(SEGS$nMajor>=y_limit, adjustcolor(colourA, red.f=0.75, green.f=0.75, blue.f=0.75), colourA), border=ifelse(SEGS$nMajor>=y_limit, adjustcolor(colourA, red.f=0.75, green.f=0.75, blue.f=0.75), colourA))
    rect(SEGS$startpos_adjusted, (SEGS$nMinor-0.07), SEGS$endpos_adjusted, (SEGS$nMinor+0.07), col=ifelse(SEGS$nMinor>=y_limit, adjustcolor(colourB, red.f=0.75, green.f=0.75, blue.f=0.75), colourB), border=ifelse(SEGS$nMinor>=y_limit, adjustcolor(colourB, red.f=0.75, green.f=0.75, blue.f=0.75), colourB))
    
  }
  if(!is.null(highlight)){
    
  
    if(length(highlight$x)==2){
     
      SEGS=SEGMENTS[which(SEGMENTS$startpos_adjusted <= highlight$x[1] & SEGMENTS$endpos_adjusted >= highlight$x[1]), ]
      
      rect(SEGS$startpos_adjusted, (SEGS$nMajor-0.07), SEGS$endpos_adjusted, (SEGS$nMajor+0.07), col=ifelse(SEGS$nMajor>=y_limit, adjustcolor(colourH, red.f=0.75, green.f=0.75, blue.f=0.75), colourH), border=ifelse(SEGS$nMajor>=y_limit, adjustcolor(colourH, red.f=0.75, green.f=0.75, blue.f=0.75), colourH))
      rect(SEGS$startpos_adjusted, (SEGS$nMinor-0.07), SEGS$endpos_adjusted, (SEGS$nMinor+0.07), col=ifelse(SEGS$nMinor>=y_limit, adjustcolor(colourH, red.f=0.75, green.f=0.75, blue.f=0.75), colourH), border=ifelse(SEGS$nMinor>=y_limit, adjustcolor(colourH, red.f=0.75, green.f=0.75, blue.f=0.75), colourH))
      
    }
    
    else {
      
      SEGS=SEGMENTS[which(SEGMENTS$startpos_adjusted <= highlight$x & SEGMENTS$endpos_adjusted >= highlight$x), ]
      
      if(SEGS$nMajor == highlight$y) {
        rect(SEGS$startpos_adjusted, (SEGS$nMajor-0.07), SEGS$endpos_adjusted, (SEGS$nMajor+0.07), col=ifelse(SEGS$nMajor>=y_limit, adjustcolor(colourH, red.f=0.75, green.f=0.75, blue.f=0.75), colourH), border=ifelse(SEGS$nMajor>=y_limit, adjustcolor(colourH, red.f=0.75, green.f=0.75, blue.f=0.75), colourH))
        
      }
      else {
        
        rect(SEGS$startpos_adjusted, (SEGS$nMinor-0.07), SEGS$endpos_adjusted, (SEGS$nMinor+0.07), col=ifelse(SEGS$nMinor>=y_limit, adjustcolor(colourH, red.f=0.75, green.f=0.75, blue.f=0.75), colourH), border=ifelse(SEGS$nMinor>=y_limit, adjustcolor(colourH, red.f=0.75, green.f=0.75, blue.f=0.75), colourH))
        
      }
      
    
    }
   
  }
   abline(v=c(1, REF$cumul), lty=1, col="lightgrey")
  text(REF$middle, y_limit, REF$chrom, pos = 1, cex = 2)
  #dev.off()
  #rm(SEGS, ticks, maintitle)
  #; rm(SAMPLE)
}

#ASCAT_output_object <- 
# ASCAT_output_object <- ascat.output
# REF <- "hg38"
# y_limit=5
# plot_unrounded=FALSE
# # png_prefix="test"
# # 
# ascat.plotAdjustedAscatProfile(res_ASCAT, REF = "hg38")  

