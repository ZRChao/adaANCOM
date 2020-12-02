#' @title adaANCOM
#'
#' @description The main function for testing on the internal and leaf nodes. It first test on each internal node by
#' building log-ratios transformation and then test for each leaf node by build log-ratios transformation adaptively based on the first step.
#'
#' @param otu a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the taxa.
#' @param Y interested outcome.
#' @param tree the phylogenetic tree.
#' @param tfun the testing procedure, such as t.test or wilcox rank sum test for two-group comparison.
#' @param smooth the strategy for zero smooth. Can be choosen from a constant (e.g 0.5/1), or DM, ZIDM, MIX which corresponding the posterior mean or after likelihood ratio test to choose the best fitted model's.
#' @param alpha the significant level both for internal and leaf nodes.
#'
#'
#' @return Returns a list:
#' \itemize{
#' \item V testing results for internal nodes V .
#' \item L testing results for leaf nodes L.
#' }
#'
#' @examples
#' library(phyloseq)
#' data(COMBO)
#' otu.tab <- otu_table(COMBO)
#' tree <- phy_tree(COMBO)
#' group <- sample_data(COMBO)$bmicat1norm2ow3ob
#' group[group==3] <- 2
#' fit <- adaANCOM(otu=t(data.frame(otu.tab)), Y=group, tree=tree, tfun = t.test, smooth=0.5)
#' summary(fit)
#' @export



adaANCOM = function(otu=NULL, Y = NULL, tree, tfun = t.test, smooth=0.5, alpha=0.05){
  # otu can save as one matrix (N*p) for leaves
  # matrix one with colnames match OTU names;
  edge  <- tree$edge
  Nnod  <- tree$Nnode
  root  <- min(edge[,1])
  a <- alpha

  ## check the OTU and sample size,
  otu_names <- tree$tip.label # matching 1:p
  G1    <- unique(Y)[1] # two group analysis
  N1    <- length(G1)
  p     <- length(otu_names)
  tab_names <- colnames(otu)
  if(length(intersect(otu_names, tab_names)) != p) stop('Mismatch OTU between table and tree.')

  if(is.matrix(otu) | is.data.frame(otu)) {
    otu1 <- otu[Y==G1, ]; res1 <- Ytreefun(otu1, tree)
    otu2 <- otu[Y!=G1, ]; res2 <- Ytreefun(otu2, tree)
    Yotu1 <- res1;
    Yotu2 <- res2;
    # cat('# Start analysis for', p, ' OTUs of one table .\n')

    ## pseudo on the tree to smooth the count data
    P1 <- Smooth_Tree(Yotu1, tree, type=smooth, rel = F, outlier=T)
    P2 <- Smooth_Tree(Yotu2, tree, type=smooth, rel = F, outlier=T)
    Potu1 <- P1$Ytree; Altb1 <- P1$Altb
    Potu2 <- P2$Ytree; Altb2 <- P2$Altb
  }

  Tfun_check <- function(Group, tfun) {
    tf <- sum(all.equal(tfun, t.test)==T, all.equal(tfun, wilcox.test)==T)
    if(length(unique(Group)) == 1) {
      stop('Only one group of this  data !')
    }
    if(length(unique(Group)) == 2) {
      l1 <- unique(Group)[1]
      if(sum(Group==l1)<2 | sum(Group !=l1)<2){
        stop('Not enough sample of any groups.')
      }
      if(tf == 0)
        stop('Two group test method mismatched !')
    }else{
      if(tf >0 & length(unique(Group))>2)
        stop('Multiple group test method mismatch !') # kruskal.test
    }
    tfun
  }

  logratio_test <- function(x, y, tfun=tf) {
    f <- formula("x ~ G")
    xlr <- data.frame(x=c(x, y), G=c(rep('A', length(x)), rep('B', length(y))))
    xlr <- xlr[!is.infinite(xlr$x), ]
    xlr <- xlr[!is.na(xlr$x), ]
    if(length(unique(xlr$x))==1) {
      pv <- 1
    }else{
      t1 <- try(tfun(f, data=xlr)$p.value, silent = T)
      t2 <- try(wilcox.test(f, data=xlr)$p.value, silent = T)
      t2 <- ifelse(class(t2)=='try-error'|is.na(t2), 1, t2)
      pv <- ifelse(class(t1)=='try-error', t2, t1)
    }
    pv
  }

  tf <- Tfun_check(Y, tfun)

  ## test for internal nodes
  taxa_info <- data.frame(taxa=as.character((p+1):(p+Nnod)), p.value=1, p.adj=1)
  taxa.pv <- c()
  for(i in 1:Nnod) {
    x <- log(Potu1[[i]][, 1]) - log(Potu1[[i]][, 2])
    y <- log(Potu2[[i]][, 1]) - log(Potu2[[i]][, 2])
    taxa.pv[i] <- logratio_test(x, y, tfun = tf)
  }
  taxa_info$p.value <- taxa.pv
  taxa_info$p.adj   <- taxa.pvj <- p.adjust(taxa.pv, method = 'fdr')

  ## test for terminal nodes
  potu <- rep(1, p); deg <- rep(F, p)
  for(i in 1:p) {
    b0 <- find_ascent(edge, i) # delete itself
    b <- b0[-1]
    (pa<- b[which(taxa.pvj[b-p]<=a)])
    (na<- length(pa))
    if(na==0) {
      potu[i] <- taxa.pv[b[1]-p]
    }else{
      op <- edge[edge[, 1]==pa[na], 2]
      (op <- as.character(setdiff(op, b0)))
      x <- log(Altb1[, i]) - log(Altb1[, op]);
      y <- log(Altb2[, i]) - log(Altb2[, op])
      potu[i] <- logratio_test(x, y, tfun=tf)
      if(na>1 & potu[i]>a) deg[i] <- T
      # while(na>1 & deg[i]) {
      #   na <- na - 1; deg[i] <- F
      #   op <- as.character(setdiff(edge[edge[, 1]==pa[na], 2] , b0))
      #   x <- log(Altb1[, i]) - log(Altb1[, op])
      #   y <- log(Altb2[, i]) - log(Altb2[, op])
      #   potu[i] <- logratio_test(x, y, tfun=tf)
      #   if(potu[i]>a) deg[i] <- T
      # }
    }
  }
  otu_info <- data.frame(otu=otu_names, p.value= potu, deg=deg)

  result <- list(V=taxa_info, L=otu_info)
  result
}


