
myascat.plotSegmentedData = function(ASCATobj, plotLogR=TRUE, logr.y_values=c(-2, 2)) {
  for (arraynr in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    Select_nonNAs = rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]])
    AllIDs = 1:dim(ASCATobj$Tumor_LogR)[1]
    names(AllIDs) = rownames(ASCATobj$Tumor_LogR)
    HetIDs = AllIDs[Select_nonNAs]
    
    #png(filename = file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr], ".ASPCF.png", sep="")), width = 2000, height = 1000, res = 200)
    
    par(mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main=3, cex.axis = 2)
    
    r = ASCATobj$Tumor_LogR_segmented[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), arraynr]
    beta = ASCATobj$Tumor_BAF_segmented[[arraynr]][, , drop=FALSE]
    
    if(plotLogR) {
      
      plot(c(1, length(r)), logr.y_values, type = "n", xaxt = "n", main = paste(colnames(ASCATobj$Tumor_BAF)[arraynr], ", LogR", sep=""), xlab = "", ylab = "")
      points(ASCATobj$Tumor_LogR[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), arraynr], col = "red", pch=ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000, ".", 20))
      points(ASCATobj$Tumor_LogR[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), arraynr], col = "#77000011", pch=ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000, ".", 20))
      points(r, col="#1b38ae", pch=16)
      abline(v=0.5, lty=1, col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ASCATobj$ch)) {
        chrk = intersect(ASCATobj$ch[[j]], HetIDs)
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2
        text(tpos, logr.y_values[2], ASCATobj$chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5, lty=1, col="lightgrey")
      }
      
    }
    else {
      
      plot(c(1, length(beta)), c(0, 1), type = "n", xaxt = "n", main = paste(colnames(ASCATobj$Tumor_BAF)[arraynr], ", BAF", sep=""), xlab = "", ylab = "")
      points(ASCATobj$Tumor_BAF[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), arraynr], col = "red", pch=ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000, ".", 20))
      points(ASCATobj$Tumor_BAF[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), arraynr], col = "#77000011", pch=ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000, ".", 20))

      points(beta, col = "#1b38ae", pch=16)
      points(1-beta, col = "#1b38ae", pch=16)
      abline(v=0.5, lty=1, col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ASCATobj$ch)) {
        chrk = intersect(ASCATobj$ch[[j]], HetIDs)
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2
        text(tpos, 1, ASCATobj$chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5, lty=1, col="lightgrey")
      }
      
    }

    #dev.off()
  }
}
