#' Compare mutation load against pcawg cohorts

#' @param col color vector for input pcawg cohort. Default gray70 and black.

#' @param bg_col background color. Default'#E0E0E0', '#FFFFFF'

#' @param medianCol color for median line. Default red.


bg_col = c('#E0E0E0', '#FFFFFF')

medianCol = 'red'

pcawg.cohort = data.table::fread("PCAWGid_TRAperSample_tumor.txt", sep = '\t', stringsAsFactors = FALSE)

#class(pcawg.cohort)

pcawg.cohort = pcawg.cohort[,.(Tumor_Sample_Barcode, total, cohort)]

pcawg.cohort$total = as.numeric(as.character(pcawg.cohort$total))


pcawg.cohort.med = pcawg.cohort[,.(.N, median(total)),cohort][order(V2, decreasing = TRUE)]

pcawg.cohort$cohort = factor(x = pcawg.cohort$cohort,levels = pcawg.cohort.med$cohort)

colnames(pcawg.cohort.med) = c('Cohort', 'Cohort_Size', 'Median_Mutations')

pcawg.cohort$pcawg = 'pcawg'

#return(pcawg.cohort)

pcawg.cohort = split(pcawg.cohort, as.factor(pcawg.cohort$cohort))

plot.dat = lapply(seq_len(length(pcawg.cohort)), function(i){
  
  x = pcawg.cohort[[i]]
  
  x = data.table::data.table(rev(seq(i-1, i, length.out = nrow(x))),
                             
                             x[order(total, decreasing = T), total],
                             
                             x[,pcawg])
  
  x
  
})

names(plot.dat) = names(pcawg.cohort)

y_lims = range(log10(data.table::rbindlist(l = plot.dat)[,V2]))

#y_lims = range(log10(unlist(lapply(plot.dat, function(x) max(x[,V2], na.rm = TRUE)))))

y_max = ceiling(max(y_lims))

y_min = floor(min(y_lims))

y_lims = c(y_min, y_max)

y_at = pretty(y_lims)





par(mar = c(4, 2, 1, 0))

plot(NA, NA, xlim = c(0, length(plot.dat)), ylim = y_lims, axes = FALSE, xlab = NA, ylab = NA)

rect(xleft = seq(0, length(plot.dat)-1, 1), ybottom = min(y_lims), xright = seq(1, length(plot.dat), 1),
     
     ytop = y_max, col = grDevices::adjustcolor(col = bg_col, alpha.f = 0.2),
     
     border = NA)

abline(h = pretty(y_lims), lty = 2, col = "gray70")

#abline(v = seq(1, length(plot.dat)), lty = 1, col = "gray70")


lapply(seq_len(length(plot.dat)), function(i){
  
  x = plot.dat[[i]]
  
  points(x$V1, log10(x$V2), pch = 16, cex = 0.3, col = 'black')
  
})


axis(side = 2, at = y_at, las = 2, line = -1, tick = FALSE, cex.axis = 0.9)

axis(side = 1, at = seq(0.5, length(plot.dat)-0.5, 1), labels = names(plot.dat),
     
     las = 2, tick = FALSE, line = -1, cex.axis = 0.9)

mtext(text = "log10 ( TRA per sample )", side = 2, line = 1.2, cex.lab = 0.9)



pcawg.cohort.med[, Median_Mutations_log10 := log10(Median_Mutations)]

lapply(seq_len(nrow(pcawg.cohort.med)), function(i){
  
  segments(x0 = i-0.8, x1 = i-0.2, y0 = pcawg.cohort.med[i, Median_Mutations_log10],
           
           y1 = pcawg.cohort.med[i, Median_Mutations_log10], col = medianCol)
  
})
