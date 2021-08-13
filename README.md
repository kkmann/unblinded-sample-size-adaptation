![run](https://github.com/kkmann/unplanned-sample-size-adaptation/workflows/run/badge.svg)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kkmann/unblinded-sample-size-adaptation/master?urlpath=rstudio)
[![DOI](https://zenodo.org/badge/274669694.svg)](https://zenodo.org/badge/latestdoi/274669694)
[![arXiv](https://img.shields.io/badge/arXiv-2010.06567-b31b1b.svg)](https://arxiv.org/abs/2010.06567)

# Sample Size Recalculation, When and How?

The preprint is available at https://arxiv.org/abs/2010.06567

Adapting the final sample size of a trial to the evidence accruing
during the trial is a natural way to address planning uncertainty.
Since the sample size is usually determined by an
argument based on the power of the trial,
an interim analysis raises the question of
how the final sample size should be determined conditional on the accrued information.
To this end, we first review and compare common approaches to estimating conditional power
which is often used in heuristic sample size recalculation rules.
We then discuss the connection of heuristic sample size recalculation and
optimal two-stage designs demonstrating that the latter is the superior approach in a fully pre-planned setting.
Hence, unplanned design adaptations should only be conducted
as reaction to trial-external new evidence, operational needs to violate the originally chosen design,
or \textit{post~hoc} changes in the optimality criterion but not as a reaction to trial-internal data.
We are able to show that commonly discussed sample size recalculation rules lead to paradoxical adaptations
where an initially planned optimal design is not invariant under the adaptation rule
even if the planning assumptions do not change.
Finally, we propose two alternative ways of reacting to
newly emerging trial-external evidence in ways that are consistent with the originally
planned design to avoid such inconsistencies.

Feel free to explore the repository interactively using the binder badge powered by https://mybinder.org or clone and explore locally.

`R` script files for reproducing the results are available in the `R` folder.
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
