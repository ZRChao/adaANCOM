## adaANCOM

We introduced adaANCOM for microbiome differential abundance analysis by incorporating phylengy.

## Installation and loading the package

```R
devtools::install_github("ZRChao/adaANCOM")  
library(adaANCOM)
```
## Usage

* Data generation from zero-inflated Dirichlet multinomial (ZIDM) distribution.

```R
N <- 100
D1 <- D2 <- round(runif(N, 1, 100))
Pi <- rzidirichlet(N=N, alpha =  c(3, 6), pi=0.1) #  ZID
which(Pi[,1]==0) # true zero
X1 <- matrix(NA, N, 2)
for(i in 1:N) X1[i, ] <- rmultinom(1, D1[i], prob=Pi[i, ])
```

* Parameter estimation

```R
est0 <- mom.dm(X1) # moment methods for DM
est0

est1 <- est.dm.NR(X1) # Newton-Raphson for DM
est1

est2 <- est.zidm.EM(X1) # expectation-maximization (EM) for ZIDM
est2
```


* Model selection with two steps likelihood ratio test, return the p-value for testing `MN < DM < ZIDM`, and the loglikelihood of this three models, 

```R
LikRTest(X1)
```

* Posterior mean transformation, return a matrix with the transformed relative data

```R
prop <- Smooth_X(X1, type='MIX', rel=TRUE)
head(prop)
```

* Outlier detection, return the index of the outliers and values of the log-ratio for the transformed data

```R
res <- Outliers(X1, prop)
res
```


* Next, we illustrate the differential abundance testing procedure based on the example data

```R
library(phyloseq)
data(example)
data <- example_data$data
fit <- adaANCOM(otu=data.frame(otu_table(data)), 
                Y=sample_data(data)$group, 
                tree=phy_tree(data), tfun = t.test, smooth=0.5)
fit$L # results for the leaves
fit$V # results for the internal nodes
```

More examples could be found each functions in the package. 

Any suggestions or problem, please contact _Chao Zhou（Supdream8@sjtu.edu.cn)_ .
