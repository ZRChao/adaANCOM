## adaANCOM

We introduced adaANCOM for microbiome differential abundance analysis by incorporating phylengy.

## Installation

```R
devtools::install_github("ZRChao/adaANCOM")  
```
## Usage

* Following, we use the `COMBO` data as example.

```{R}
library(phyloseq)
library(adaANCOM)
data(COMBO)
COMBO  # data information

otu.tab <- otu_table(COMBO)
tree <- phy_tree(COMBO)
group <- sample_data(COMBO)$bmicat1norm2ow3ob
group[group==3] <- 2
```

* Data transformation

```R

```

```R
fit <- adaANCOM(otu=t(data.frame(otu.tab)), Y=group, tree=tree, tfun = t.test, smooth=0.5)
summary(fit)
```

More examples could be found in the functions ```adaANCOM```. 

Any suggestions or problem, please contact _Chao Zhou（Supdream8@sjtu.edu.cn)_ .
