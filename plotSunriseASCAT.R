library(RColorBrewer)
source("findMinima.R")
ascat.plotSunrise<-function(d, psi_opt1, rho_opt1, minim=TRUE, optima=FALSE) {
  
  par(mar = c(5, 5, 0.5, 0.5), cex=0.75, cex.lab=2, cex.axis=2)
  tryCatch({
    
    if (minim) {
      hmcol = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
    } else {
      hmcol = colorRampPalette(brewer.pal(10, "RdBu"))(256)
    }
    image(log(d), col = hmcol, axes = FALSE, xlab = "Ploidy", ylab = "Purity")
    
    ploidy_min<-as.numeric(rownames(d)[1])
    ploidy_max<-as.numeric(rownames(d)[nrow(d)])
    purity_min<-as.numeric(colnames(d)[1])
    purity_max<-as.numeric(colnames(d)[ncol(d)])
    
    if(optima) {
      
        bao1 <- findLocalMinima(d)
        ao <- bao1$bao
        
        text((ao[,1] - 0.5) / nrow(d),
             (ao[,2] - 0.5) / ncol(d),
             label = 1:nrow(ao),
             col = "chartreuse", pch = 19)
      
     
    }
    
    PLOIDY_LABELS=seq(ploidy_min, ploidy_max, length.out=max(c((ploidy_max-ploidy_min)+1, 2)))
    axis(1, at = (PLOIDY_LABELS-ploidy_min) / (ploidy_max-ploidy_min), labels = round(PLOIDY_LABELS, 2))
    if (purity_min==0.1 && purity_max==1.05) {
      PURITY_LABELS=round(seq(purity_min, purity_max, by = 0.3),2) # default Y-axis
    } else {
      PURITY_LABELS=round(seq(purity_min, purity_max, length.out=4),2)
    }
    axis(2, at = round((PURITY_LABELS-purity_min) / (purity_max-purity_min),2), labels = PURITY_LABELS)
    
    if (psi_opt1>0 && rho_opt1>0) {
  
        psi_index <- which.min(abs(as.numeric(rownames(d)) - psi_opt1))
        rho_index <- which.min(abs(as.numeric(colnames(d)) - rho_opt1))
        points((psi_index - 0.5) / nrow(d), 
               (rho_index - 0.5) / ncol(d), 
               col = "green", pch = 4, cex = 2)
      
      
      }
    
    if (optima) ao else invisible(NULL)
  
  },
  error=function(e) {
    print(e)
  },
  warning=function(w){
    print(w)
  })

}
