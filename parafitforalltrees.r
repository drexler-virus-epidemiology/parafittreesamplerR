library (ape)
library(vegan)


trees<-read.nexus("XXXXXXXXXXXXXXXXXX.t") #change 

#nexus trees extract and save in folder
sample.trees<-function(trees, burnin, final.number){           
  
  #NUMBER OF TREES IN ORIGINAL FILE
  original.number<-as.numeric(length(trees))                           
  
  #THIS CREATES THE POST BURNIN PORTION
  post.burnin.trees<-trees[(burnin*original.number):original.number]   
  
  #THIS DOWNSAMPLES THE COLLECTION OF TREES
  random.trees<-sample(post.burnin.trees, final.number)                 
  
  for (i in 1:length(random.trees)) { # create a file name for the current tree
    
    # write the current tree to the file
    write.tree(random.trees[[i]], file = paste0("tree", i, ".newick"))
  }

}

#run the function, you get 10 trees that will be in your folder
sample.trees(trees, .25, 5) # 25% burnin and 10 trees

#read all trees in the folder 

folder <- ("C:/XXXXXXXXXXXXXXXXXXXXXXX") #change

tree_files <- list.files(folder, pattern = ".newick", full.names = TRUE)

#remove the outgroups change
for (file in tree_files) {
  
  tree <- read.tree(file)
  
  tree_pruned <- drop.tip(tree, "XXXXXXXXXXXX") #change accordingly 

  newick_file <- gsub(".newick", ".newick", file) # modify the original file name if you want
  write.tree(tree_pruned, file = newick_file)
  
}

#use the tree for the parafit functions in a loop 

tree_fileswoog <- list.files(folder, pattern = ".newick", full.names = TRUE)

#run the parafit function in all files

for (file in tree_fileswoog) {
  
  
  links = "C:/XXXXXXXXXXXXXXXXXXXXX.csv" #change
  tree_host ="C:/XXXXXXXXXXXXXXXXXXXXXXXXx.nwk" #change
  tree_parasite = (file)
  permutations = "permut"
  N.perm = 1000
  output_file <- "output.txt"
  
  HP <- as.matrix(read.csv(links, header=TRUE, sep= ";"))
  NLinks = sum(HP)
  TreeH <- read.tree(tree_host)
  TreeP <- read.tree(tree_parasite)
  
  host.D <- cophenetic(TreeH)
  para.D <- cophenetic(TreeP)
  host.D <- host.D[rownames(HP),rownames(HP)]
  para.D <- para.D[colnames(HP),colnames(HP)]
  
  
  PACo <- function (H.dist, P.dist, HP.bin){
    HP.bin <- which(HP.bin > 0, arr.in=TRUE)
    H.PCo <- pcoa(H.dist, correction="cailliez")$vectors
    P.PCo <- pcoa(P.dist, correction="cailliez")$vectors
    H.PCo <- H.PCo[HP.bin[,1],]
    P.PCo <- P.PCo[HP.bin[,2],]
    list (H.PCo = H.PCo, P.PCo = P.PCo)
  }
  
  PACo.fit <- PACo(host.D, para.D, HP)
  HP.proc <- procrustes(PACo.fit$H.PCo, PACo.fit$P.PCo)
  
  HostX <- HP.proc$X
  ParY <- HP.proc$Yrot
  plot(HostX, asp=1, pch=46)
  points(ParY, pch=1)
  arrows(ParY[,1], ParY[,2], HostX[,1], HostX[,2], length=0.12, angle=15,xpd=FALSE)
  HostX <- unique(HP.proc$X)
  ParY <- unique(HP.proc$Yrot)
  #identify(ParY[,1], ParY[,2], rownames(ParY), offset=0.3, xpd=FALSE, cex=0.8)
  #identify(HostX[,1], HostX[,2], rownames(HostX),offset=0.3, xpd=TRUE, cex= 0.8)
  
  
  m2.obs <- HP.proc$ss
  
  P.value = 0
  set.seed(5)
  
  
  for (n in c(1:N.perm))
  {if (NLinks <= nrow(HP) | NLinks <= ncol(HP))
  {flag2 <- TRUE
  while (flag2 == TRUE)
  {HP.perm <- t(apply(HP,1,sample))
  if(any(colSums(HP.perm) == NLinks)) flag2 <- TRUE
  else flag2 <- FALSE}
  } else { HP.perm <- t(apply(HP,1,sample))}
    PACo.perm <- PACo(host.D, para.D, HP.perm)
    m2.perm <- procrustes(PACo.perm$H.PCo, PACo.perm$P.PCo)$ss
    write(m2.perm, file=permutations, sep="\t", append=TRUE)
    if (m2.perm <= m2.obs){
      cat("olala\n")
      P.value = P.value + 1
    }
  }
  
  
  P.value <- P.value/N.perm
  
  cat("for parasite tree" , file, "The observed m2 is ", m2.obs, "\np-value is", P.value, " based on ", N.perm,"permutations.", file = "output.txt", append = TRUE)
  
  
  HP.ones <- which(HP > 0, arr.in=TRUE)
  SQres.jackn <- matrix(rep(NA, NLinks**2), NLinks)
  colnames (SQres.jackn) <- paste(rownames(HostX),rownames(ParY), sep="-")
  t.critical = qt(0.975,NLinks-1)
  for(i in c(1:NLinks))
  {HP.ind <- HP
  HP.ind[HP.ones[i,1],HP.ones[i,2]]=0
  PACo.ind <- PACo(host.D, para.D, HP.ind)
  Proc.ind <- procrustes(PACo.ind$H.PCo, PACo.ind$P.PCo)
  res.Proc.ind <- c(residuals(Proc.ind))
  res.Proc.ind <- append (res.Proc.ind, NA, after= i-1)
  SQres.jackn [i, ] <- res.Proc.ind}
  SQres.jackn <- SQres.jackn**2
  SQres <- (residuals (HP.proc)**2)
  SQres.jackn <- SQres.jackn*(-(NLinks-1))
  SQres <- SQres*NLinks
  SQres.jackn <- t(apply(SQres.jackn, 1, "+", SQres))
  phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE)
  phi.UCI <- apply(SQres.jackn, 2, sd, na.rm = TRUE)
  phi.UCI <- phi.mean + t.critical * phi.UCI/sqrt(NLinks)
  pat.bar <- barplot(phi.mean, names.arg = " ", space = 0.25, col="white", xlab=
                       "Host-parasite link", ylab= "Squared residuals", ylim=c(0, max(phi.UCI)),
                     cex.lab=1.2)
  text(pat.bar, par("usr")[3] - 0.001, srt = 330, adj = 0, labels =
         colnames(SQres.jackn), xpd = TRUE, font = 1, cex=0.6)
  arrows(pat.bar, phi.mean, pat.bar, phi.UCI, length= 0.05, angle=90)
  abline(a=median(phi.mean), b=0, lty=2)
  
  
  res <- parafit(host.D, para.D, HP, nperm=N.perm, test.links=TRUE, correction="lingoes")
  
  write.csv (res$link.table, file=paste0(sub("(.*/|\\.csv)", "", file), ".csv"))
  
  sink (output_file, append = TRUE)
  print(res)
  
  #listofhosts <- cat("list of hosts\n")
  listofhosts <-  for (i in 1:length(rownames(host.D))){
    cat(i, rownames(host.D)[i],"\n")
  }

  print (listofhosts)
  #listofparasites <- cat("list of parasites\n")
  listofparasites <- for (i in 1:length(rownames(para.D))){
    cat(i, rownames(para.D)[i],"\n")
  }
print (listofparasites)
sink()
}
