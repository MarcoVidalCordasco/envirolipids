# brGDGT-NPP

R package to predict Net Primary Productivity (NPP) from the relative abundance
of individual brGDGTs, using a centred log-ratio (CLR) transformation and a
pre-trained random forest model. It includes diagnostic tools to assess
whether new samples fall within the brGDGT composition space of the training
data. When mean annual temperature and precipitation are available, the package
also computes potential NPP with the Miami and Thornthwaite Memorial models,
allowing comparison between climate-based and brGDGT-based NPP estimates to
identify deviations from climatic expectations.

## Installation

```r
install.packages("devtools")
devtools::install_github("MarcoVidalCordasco/brGDGT-NPP",
                         build_vignettes = TRUE,
                         dependencies = TRUE)
```

## Usage

```r
library(envirolipids)
data(example_data)

# Do the samples fall within the training data?
check <- npp_brgdgt_applicable(example_data)
check
autoplot(check, type = "pca")    # also "mahalanobis", "ternary", "map"

# NPP from brGDGTs
npp <- npp_brgdgt(example_data, sample_type = "Lacustrine Sediment")
autoplot(npp, x = example_data$Age)

# Potential NPP from climate, and the deficit
climate <- npp_miami(mat = my_mat, map = my_map)   # or npp_thornthwaite()
deficit <- hanpp(npp_empirical = npp, npp_model = climate)
autoplot(deficit, x = example_data$Age)
```

See `vignette("envirolipids-demo")` for the full workflow.
