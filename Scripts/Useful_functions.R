
fun.p <- function(i, j,fact,resp) {
  fact2 <- droplevels(fact[as.numeric(fact) %in% c(i, 
                                                   j)])
  resp2 <- as.matrix(resp)
  rows <- which(fact %in% levels(fact2))
  resp2 <- as.dist(resp2[rows, rows])
  vegan::adonis(resp2 ~ fact2, permutations = nperm)$aov.tab[1, 
                                                             "Pr(>F)"]
}
fun.p1 <- function(i, j,fact,resp) {
  fact2 <- droplevels(fact[as.numeric(fact) %in% c(i, 
                                                   j)])
  resp2 <- as.matrix(resp)
  rows <- which(fact %in% levels(fact2))
  resp2 <- as.dist(resp2[rows, rows])
  vegan::adonis(resp2 ~ fact2, permutations = nperm)$aov.tab[1, 
                                                             "F.Model"]
}

My_pairwise.perm.manova=function (resp, fact, nperm = 999, p.method = "fdr") 
{
  fact <- factor(fact)
  Pairs=data.frame(t(combn(levels(fact),m=2)))
  Pairs$P_value=NA
  Pairs$F_value=NA
  Pairs$X1=factor(Pairs$X1,levels=levels(fact))
  Pairs$X2=factor(Pairs$X2,levels=levels(fact))
  
  for (k in 1:dim(Pairs)[1])
  {
    Pairs$P_value[k]=fun.p(i=Pairs$X1[k],j=Pairs$X2[k],fact=fact,resp=resp)
    Pairs$F_value[k]=fun.p1(Pairs$X1[k],Pairs$X2[k],fact=fact,resp=resp)
  }
  
  p.adjust(Pairs$P_value, method = "fdr", n = length(Pairs$P_value))
  return(Pairs)
}




pairwise.tableF=function (compare.levels, level.names) 
{
  ix <- setNames(seq_along(level.names), level.names)
  pp <- outer(ix[-1L], ix[-length(ix)], function(ivec, jvec) sapply(seq_along(ivec), 
                                                                    function(k) {
                                                                      i <- ivec[k]
                                                                      j <- jvec[k]
                                                                      if (i > j) 
                                                                        compare.levels(i, j)
                                                                      else NA
                                                                    }))
  pp
}

pairwise.perm.manovaF=function (resp, fact, test = c("Pillai", "Wilks", "Hotelling-Lawley", 
                               "Roy", "Spherical"), nperm = 999, progress = TRUE, p.method = "fdr") 
{
  call <- match.call()
  dname <- paste0(deparse(call$resp), " by ", deparse(substitute(fact)), 
                  "\n", nperm, " permutations")
  fact <- factor(fact)
  if ("dist" %in% class(resp)) {
    fun.p <- function(i, j) {
      fact2 <- droplevels(fact[as.numeric(fact) %in% c(i, 
                                                       j)])
      resp2 <- as.matrix(resp)
      rows <- which(fact %in% levels(fact2))
      resp2 <- as.dist(resp2[rows, rows])
      vegan::adonis(resp2 ~ fact2, permutations = nperm)$aov.tab[1, 
                                                                 "F.Model"]
    }
    multcomp <- pairwise.tableF(fun.p, levels(fact))
    method <- "permutation MANOVAs on a distance matrix"
  }
  else {
    if (nrow(resp) != length(fact)) {
      stop(paste("'", deparse(substitute(resp)), "' and '", 
                 deparse(substitute(fact)), "' lengths differ", 
                 sep = ""))
    }
    test <- match.arg(test)
    if (!is.matrix(resp)) {
      resp <- as.matrix(resp)
    }
    if (!is.factor(fact)) {
      fact <- factor(fact)
    }
    fun.p <- function(i, j) {
      resp2 <- resp[as.numeric(fact) %in% c(i, j), ]
      fact2 <- droplevels(fact[as.numeric(fact) %in% c(i, 
                                                       j)])
      perm.manova(resp2, fact2, test = test, nperm = nperm, 
                  progress)
    }
    multcomp <- pairwise.table(fun.p, levels(fact), p.adjust.method = p.method)
    method <- paste0("permutation MANOVAs (test: ", test, 
                     ")")
  }
  result <- list(method = method, data.name = dname, F = multcomp)
  return(result)
}

#Function to format aov table 

format_aov_table=function(adonis2table,data_type)
{
  adonis2table=data.frame(adonis2table)
  colnames(adonis2table)[5]="P_value"
  row.names(adonis2table)=paste(data_type,row.names(adonis2table),sep="_")
  for (i in 1:4){adonis2table[,i]=round(adonis2table[,i],2)}
  adonis2table[,5]=round(adonis2table[,5],4)
  return(adonis2table)
}




# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.

# This function has been copied from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}