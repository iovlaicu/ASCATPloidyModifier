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
      #print(ao)
      #msg = ao
      
      #cmd = paste0("echo ", msg)
     
      new_pl <- as.numeric(rownames(d)[ao[,1]])
      new_pr <- as.numeric(colnames(d)[ao[,2]])
      
      # cmd = paste0("echo ", new_pl)
      # system(cmd)
      
      text((new_pl-ploidy_min) / (ploidy_max-ploidy_min),
           (new_pr-purity_min) / (purity_max-purity_min), label=1:nrow(ao),
           col = "chartreuse", pch = 19)
      
      
      #ao <- bao1$ao
      # if(plotClust){
      #   text(ao[,2]/ncol(errs),
      #        1-ao[,1]/nrow(errs), label=1:nrow(ao),
      #        col = RColorBrewer::brewer.pal(12,"Paired")[bao1$clusts], cex=.6)
      #   
      # }
      
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
      points((psi_opt1-ploidy_min) / (ploidy_max-ploidy_min), (rho_opt1-purity_min) / (purity_max-purity_min), col="green", pch=4, cex = 2)
    }
    
    ao
  
  },
  error=function(e) {
    print(e)
  },
  warning=function(w){
    print(w)
  })

}

#distance_marix <- ascat.output$distance_matrix[[1]]

#ascat.plotSunrise(distance_marix, psi_opt1 = res_ASCAT$allSolutions[[1]]$ploidy, rho_opt1 = res_ASCAT$allSolutions[[1]]$purity)
