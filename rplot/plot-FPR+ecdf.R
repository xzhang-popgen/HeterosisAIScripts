###1. get 5% critical values from neutral simulations 

genes <- c("hyal2") #in practice, the list should consists of 26 AI candidate region names

for(i in 1:length(genes)){
  setwd(paste("/Users/xinjunzhang/Desktop/HeterosisAIScripts/stats/",genes[i],"/",sep=""))
  
  neu <- read.table(paste("neu/",genes[i],"-0.0_human_neu_windows.txt",sep=""))
  neg <- read.table(paste(genes[i],"-0.0_human_neg_windows.txt",sep=""))
  
  title <- c("n","h","meanpI","start","end","freq_before","freq_after","pI","Dstat","fD","Het","divratioavg","Q_1_100_q95","Q_1_100_q90","Q_1_100_max","U_1_0_100","U_1_20_100","U_1_50_100","U_1_80_100")
  
  colnames(neu) <- title
  colnames(neg) <- title
  
  
  neg$window <- c(1:99)
  neg$pos <- round((neg$end+neg$start)/2)
  
  neu$window <- c(1:99)
  neu$pos <- round((neu$end+neu$start)/2)
  
  
  ptable_all <- data.frame(matrix(ncol = 5, nrow = 0)) #gene,stat,thresh,popx12
  colnames(ptable_all) <- c("gene","h","stat","thresh_neu","thresh_neg")
  
  
  for (g in c(0)){ #h=0 only because here we are only getting the deleterious model's critical values, in contrast to the neutral model
    neg_g <- neg[neg$h==g,]
    neu_g <- neu
    
    p_neu <- c()
    p_neg <- c()
    
    ptable <- data.frame(matrix(ncol = 5, nrow = 12)) #gene,stat,thresh,popx12
    colnames(ptable) <- c("gene","h","stat","thresh_neu","thresh_neg")
    for (p in c(8:19)){
      index <- p
      
      #############
      if(p==12){
        thresh_neu = quantile(neu_g[,index][complete.cases(neu_g[,index])],probs = 0.05)
      }
      if (p!=12){
        thresh_neu = quantile(neu_g[,index][complete.cases(neu_g[,index])],probs = 0.95)
      }
      
      if(p==12){
        thresh_neg = quantile(neg_g[,index][complete.cases(neg_g[,index])],probs = 0.05)
      }
      if (p!=12){
        thresh_neg = quantile(neg_g[,index][complete.cases(neg_g[,index])],probs = 0.95)
      }  
      
      p_neu <- c(p_neu,thresh_neu)
      p_neg <- c(p_neg,thresh_neg)
      
    }
    
    ptable$thresh_neu <- p_neu
    ptable$thresh_neg <- p_neg
    ptable$gene <- "hyal2"
    ptable$h <- g
    ptable$stat <- title[8:19]
    
    ptable_all <- rbind(ptable_all,ptable)
    
  }
  
  write.csv(ptable_all,paste("../sample_output/",genes[i],"_thresh_neu+neg.csv",sep=""),row.names = F)
  
}

###########################################################################

###2. get FPRs

for(g in 1:length(genes)){
  setwd(paste("/Users/xinjunzhang/Desktop/HeterosisAIScripts/stats/",genes[g],"/",sep=""))
  
  #part1: h in [0,0.5]
  thresh <- read.csv(paste("../sample_output/",genes[g],"_thresh_neu+neg.csv",sep=""),stringsAsFactors = F)

  allfiles <- list.files(pattern=".txt")
  
  title <- c("n","h","meanpI","start","end","freq_before","freq_after","pI","Dstat","fD","Het","divratioavg","Q_1_100_q95","Q_1_100_q90","Q_1_100_max","U_1_0_100","U_1_20_100","U_1_50_100","U_1_80_100")
  
  hvalues <- seq(0.0,0.5,by=0.005)
  
  
  FPR_neg <- as.data.frame(t(c("h","window","pI","D","fD","Het","RD","Q95","Q90","Qmax","U0","U20","U50","U80")))
  colnames(FPR_neg) <- c("h","window","pI","D","fD","Het","RD","Q95","Q90","Qmax","U0","U20","U50","U80")
  FPR_neg<-FPR_neg[-1,] 
  
  FPR_neu <- as.data.frame(t(c("h","window","pI","D","fD","Het","RD","Q95","Q90","Qmax","U0","U20","U50","U80")))
  colnames(FPR_neu) <- c("h","window","pI","D","fD","Het","RD","Q95","Q90","Qmax","U0","U20","U50","U80")
  FPR_neu<-FPR_neu[-1,] 
  
  for (i in 1:length(hvalues)){
    h <- hvalues[i]
    
    if (i ==1){file <- read.table(paste(genes[g],"-","0.0","_human_neg_windows.txt",sep = ""),stringsAsFactors = F)}
    else{file <- read.table(paste(genes[g],"-",h,"_human_neg_windows.txt",sep = ""),stringsAsFactors = F)}
    
    colnames(file) <- title
    
    file$window <- c(1:99)
    
    for (w in 1:99){
      file_w <- file[file$window==w,]
      
      fpr_neg_w <- c()
      fpr_neu_w <- c()
      
      for (s in 1:12){
        stat_neu <- thresh$thresh_neu[s]
        stat_neg <- thresh$thresh_neg[s]
        
        stat_sims <- file_w[,s+7]
        total <- length(stat_sims)
        stat_sims <- stat_sims[complete.cases(stat_sims)]
        
        if (s!=5){
          fpr_neu_s <- length(stat_sims[stat_sims>=stat_neu])/total
          fpr_neg_s <- length(stat_sims[stat_sims>=stat_neg])/total
        }
        else{
          fpr_neu_s <- length(stat_sims[stat_sims<=stat_neu])/total
          fpr_neg_s <- length(stat_sims[stat_sims<=stat_neg])/total
        }
        
        fpr_neg_w <- c(fpr_neg_w,fpr_neg_s)
        fpr_neu_w <- c(fpr_neu_w,fpr_neu_s)
        
      }
      
      neg_row <- c(h,w,fpr_neg_w)
      neu_row <- c(h,w,fpr_neu_w)
      
      neg_frame <- as.data.frame(t(neg_row))
      neu_frame <- as.data.frame(t(neu_row))
      
      colnames(neg_frame) <- c("h","window","pI","D","fD","Het","RD","Q95","Q90","Qmax","U0","U20","U50","U80")
      colnames(neu_frame) <- c("h","window","pI","D","fD","Het","RD","Q95","Q90","Qmax","U0","U20","U50","U80")
      
      FPR_neg <- rbind(FPR_neg ,neg_frame)
      FPR_neu <- rbind(FPR_neu ,neu_frame)
      
    }
    
  }
  
  write.csv(FPR_neg,paste("../sample_output/",genes[g],"_FPR_thresh-neg.csv",sep=""),row.names = F)
  write.csv(FPR_neu,paste("../sample_output/",genes[g],"_FPR_thresh-neu.csv",sep=""),row.names = F)
  
  
  #part2: add hs to the table
  
  FPR_neghs <- as.data.frame(t(c("h","window","pI","D","fD","Het","RD","Q95","Q90","Qmax","U0","U20","U50","U80")))
  colnames(FPR_neghs) <- c("h","window","pI","D","fD","Het","RD","Q95","Q90","Qmax","U0","U20","U50","U80")
  FPR_neghs<-FPR_neghs[-1,] 
  
  FPR_neuhs <- as.data.frame(t(c("h","window","pI","D","fD","Het","RD","Q95","Q90","Qmax","U0","U20","U50","U80")))
  colnames(FPR_neuhs) <- c("h","window","pI","D","fD","Het","RD","Q95","Q90","Qmax","U0","U20","U50","U80")
  FPR_neuhs<-FPR_neuhs[-1,] 
  
  for (i in 1){
    h <- "hs"
    
    if (i ==1){file <- read.table(paste("../hs/",genes[g],"-hs_human_neg_windows.txt",sep = ""),stringsAsFactors = F)}
    
    
    colnames(file) <- title
    
    file$window <- c(1:99)

    for (w in 1:99){
      file_w <- file[file$window==w,]
      
      fpr_neg_w <- c()
      fpr_neu_w <- c()
      
      for (s in 1:12){
        stat_neu <- thresh$thresh_neu[s]
        stat_neg <- thresh$thresh_neg[s]
        
        stat_sims <- file_w[,s+7]
        total <- length(stat_sims)
        stat_sims <- stat_sims[complete.cases(stat_sims)]
        
        if (s!=5){
          fpr_neu_s <- length(stat_sims[stat_sims>=stat_neu])/total
          fpr_neg_s <- length(stat_sims[stat_sims>=stat_neg])/total
        }
        else{
          fpr_neu_s <- length(stat_sims[stat_sims<=stat_neu])/total
          fpr_neg_s <- length(stat_sims[stat_sims<=stat_neg])/total
        }
        
        fpr_neg_w <- c(fpr_neg_w,fpr_neg_s)
        fpr_neu_w <- c(fpr_neu_w,fpr_neu_s)
        
      }
      
      neg_row <- c(h,w,fpr_neg_w)
      neu_row <- c(h,w,fpr_neu_w)
      
      neg_frame <- as.data.frame(t(neg_row))
      neu_frame <- as.data.frame(t(neu_row))
      
      colnames(neg_frame) <- c("h","window","pI","D","fD","Het","RD","Q95","Q90","Qmax","U0","U20","U50","U80")
      colnames(neu_frame) <- c("h","window","pI","D","fD","Het","RD","Q95","Q90","Qmax","U0","U20","U50","U80")
      
      FPR_neghs <- rbind(FPR_neghs ,neg_frame)
      FPR_neuhs <- rbind(FPR_neuhs ,neu_frame)
      
    }
    
  }
  
  write.csv(FPR_neghs,paste("../sample_output/",genes[g],"_hs_FPR_thresh-neg.csv",sep = ""),row.names = F)
  write.csv(FPR_neuhs,paste("../sample_output/",genes[g],"_hs_FPR_thresh-neu.csv",sep = ""),row.names = F)
  
  
  
  FPR_neg_all <- rbind(FPR_neg,FPR_neghs)
  FPR_neu_all <- rbind(FPR_neu,FPR_neuhs)
  
  write.csv(FPR_neg_all,paste("../sample_output/",genes[g],"_h+hs_FPR_thresh-neg.csv",sep = ""),row.names = F)
  write.csv(FPR_neu_all,paste("../sample_output/",genes[g],"_h+hs_FPR_thresh-neu.csv",sep = ""),row.names = F)
  
}


###########################################################################

###3. plot the FPRs
color <- rainbow(12)

#FPR plot
for (g in 1:length(genes)){
  pdf(paste("../sample_output/",toupper(genes[g]),"_FPR",".pdf",sep = ""),paper = "USr",width = 11,height = 8.5)
  
  par(mfrow=c(3,2))
  FPR_neu <- read.csv(paste("../sample_output/",genes[g],"_h+hs_FPR_thresh-neu.csv",sep=""),stringsAsFactors = F)
  
  these <- c(1,3,4,5,6,10) #plot these statistics (by column index)
  #these <- c(1:12) #plot all statistics
  for(i in these){
    boxplot(FPR_neu[,i+2]~FPR_neu[,1],data=FPR_neu,col=c(rep(color[i],101),"gray"),range=0,xaxt = "n",xlab="Dominance - Deleterious mutations",ylab="FPR",
            cex.lab=1.4,cex.axis=1.2,main=paste("HYAL2 - ",colnames(FPR_neu)[i+2],sep = ""),cex.main=1.5)
    abline(h=0.05,lwd=3,col="gray60",lty=1)
    axis(side=1, at=c(seq(from=1, to=110,by=10 ),102), labels=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,"hs"), cex.axis=1,las=2)
  }
  
  dev.off()
}

#ecdf plot

colfunc <- colorRampPalette(c("gray87", "gray24"))
colors <- colfunc(99)
h_series <- seq(0.005,0.495,by=0.005)

these <- c(1,3,4,5,6,10) #selected statistics

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




