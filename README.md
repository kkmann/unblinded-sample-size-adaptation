![run](https://github.com/kkmann/unplanned-sample-size-adaptation/workflows/run/badge.svg)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kkmann/unblinded-sample-size-adaptation/master?urlpath=rstudio)
[![DOI](https://zenodo.org/badge/274669694.svg)](https://zenodo.org/badge/latestdoi/274669694)
[![arXiv](https://img.shields.io/badge/arXiv-2010.06567-b31b1b.svg)](https://arxiv.org/abs/2010.06567)

# Sample Size Recalculation, When and How?

The preprint is available at https://arxiv.org/abs/2010.06567

Adapting the final sample size of a trial to the evidence accruing during the trial is a natural way to address planning uncertainty. Designs with adaptive sample size need to account for their optional stopping to guarantee strict type-I error-rate control. A variety of different methods to maintain type-I error-rate control after unplanned changes of the initial sample size have been proposed in the literature. This makes interim analyses for the purpose of sample size recalculation feasible in a regulatory context. Since the sample size is usually determined via an argument based on the power of the trial, an interim analysis raises the question of how the final sample size should be determined conditional on the accrued information. Conditional power is a concept often put forward in this context. Since it depends on the unknown effect size, we take a strict estimation perspective and compare assumed conditional power, observed conditional power, and predictive power with respect to their properties as estimators of the unknown conditional power. We then demonstrate that pre-planning an interim analysis using methodology for unplanned interim analyses is ineffective and naturally leads to the concept of optimal two-stage designs. We conclude that unplanned design adaptations should only be conducted as reaction to trial-external new evidence, operational needs to violate the originally chosen design, or post hoc changes in the objective criterion. Finally, we show that commonly discussed sample size recalculation rules can lead to paradoxical outcomes and propose two alternative ways of reacting to newly emerging trial-external evidence. 

Feel free to explore the repository interactively using the binder badge powered by https://mybinder.org or clone and explore locally.

`R` script files for reproducing the results (roughy ordered by section) are
available in the `R` folder.
The required dependencies can be installed using the
[`renv`](https://rstudio.github.io/renv/articles/renv.html) 
package and the provided lockfile `renv.lock`.
A bash script `run` to execute all code is provided. 
On Linux/Unix, execute 
```
./run
```
(might take a while).
Plots are stored in `output/figures`.
