colfunc <- colorRampPalette(c("gray87", "gray24"))
colors <- colfunc(99)
h_series <- seq(0.005,0.495,by=0.005)

these <- c(1,3,4,5,6,11) #selected statistics

for (g in 1:length(genes)){
  pdf(paste("../sample_output/",toupper(genes[g]),"_ecdf",".pdf",sep = ""),paper = "USr",width = 11,height = 8.5)
  
  FPR_neu <- read.csv(paste("../sample_output/",genes[g],"_h+hs_FPR_thresh-neu.csv",sep=""),stringsAsFactors = F)

  par(mfrow=c(2,3))
  
  for(i in these){
    test <- FPR_neu[FPR_neu$h==0,]
    test_pI.ordered = sort(test[,i+2])
    plot(test_pI.ordered, (1:99)/99, type = 's', ylim = c(0, 1), xlab = 'FPR', xlim=c(0,max(FPR_neu[,i+2])),
         ylab = '', main = paste(paste(toupper(genes[g]),' - ',colnames(test)[i+2],sep = ""),col="blue",lwd=3))
         
         
         for (h in 1:length(h_series)){
           this <- FPR_neu[FPR_neu$h==h_series[h],]
           this_pI.ordered = sort(this[,i+2])
           lines(this_pI.ordered, (1:99)/99,col=colors[h])
         }
         
         
         this <- FPR_neu[FPR_neu$h==0.5,]
         this_pI.ordered = sort(this[,i+2])
         lines(this_pI.ordered, (1:99)/99,col="black",lwd=3)
         
         
         this <- FPR_neu[FPR_neu$h=="hs",]
         this_pI.ordered = sort(this[,i+2])
         lines(this_pI.ordered, (1:99)/99,col="green",lwd=3)
         #legend("bottomright",legend = c("h=0.5","h=0","hs"),col=c("black","blue","green"),lwd=3,lty=1)
  }
  dev.off()
}
