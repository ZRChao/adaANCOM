\name{adaANCOM-package}
\alias{adaANCOM-package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Differential Abundance Analysis of Microbiome Data Incorporating Phylogeny
}
\description{
\strong{adaANCOM} simplified of adaptive ANCOM which do differential abundance analysis based on the tree both for internal and leaf nodes. The package here provides function for tree structure exploration, several parameter estimation functions for ZIDM or ZIDTM, zero smooth and outliers detection by using likelihood ratio test and so on.

    The complete list of functions can be displayed with library(help = adaANCOM).

    More information on adaANCOM can be found :
    \url{https://github.com/ZRChao/adaANCOM}
}
%\usage{
%adaANCOM(x)
%}
%- maybe also 'usage' for other objects documented here.
%\arguments{
%  \item{x}{
%%     ~~Describe \code{x} here~~
%}
%}
\details{
      Package: adaANCOM

      Type: R package

      Version: 0.1.0

      Date: 2020-12-2

      License: GPL-2
}

\author{
Chao Zhou

Maintainer : Chao ZHOU <Supdream8@sjtu.edu.cn>
}


\references{
Chao Zhou, Tao Wang. (2020) Differential Abundance Analysis of Microbiome Data Incorporating Phylogeny (Under review).

Mandal, Siddhartha, et al. "Analysis of composition of microbiomes: a novel method for studying microbial composition." Microbial ecology in health and disease 26.1 (2015): 27663.
}



%% ~Make other sections like Warning with \section{Warning }{....} ~



\examples{
library(phyloseq)
data(COMBO)
otu.tab <- otu_table(COMBO)
tree <- phy_tree(COMBO)
group <- sample_data(COMBO)$bmicat1norm2ow3ob
group[group==3] <- 2
fit <- adaANCOM(otu=t(data.frame(otu.tab)), Y=group, tree=tree, tfun = t.test, smooth=0.5)
summary(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.




