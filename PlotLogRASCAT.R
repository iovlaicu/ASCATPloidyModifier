myascat.plotRawData = function(ASCATobj, logr.y_values=c(-2, 2)) {
 
  for (i in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    #png(filename = file.path(img.dir, paste(img.prefix, ASCATobj$samples[i], ".tumour.png", sep="")), width = 2000, height = 1000, res = 200)
    #par(mar = c(0.5, 5, 5, 0.5), mfrow = c(2, 1), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000, ".", 20))
    p <- plot(c(1, dim(ASCATobj$Tumor_LogR)[1]), logr.y_values, type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", LogR", sep = ""), xlab = "", ylab = "")
    p <- p + points(ASCATobj$Tumor_LogR[, i], col="red")
    p <- p + points(ASCATobj$Tumor_LogR[, i], col="#77000011")
    p <- p + abline(v=0.5, lty=1, col="lightgrey")
    chrk_tot_len = 0
    
    for (j in 1:length(ASCATobj$ch)) {
      chrk = ASCATobj$ch[[j]]
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2
      p <- p + text(tpos, logr.y_values[2], ASCATobj$chrs[j], pos = 1, cex = 2)
      p <- p + abline(v=vpos+0.5, lty=1, col="lightgrey")
    }

  }
  
  return(p)
  
}