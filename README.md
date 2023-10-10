# StableMate
StableMate [(Deng et al.2023)](https://www.biorxiv.org/content/10.1101/2023.09.26.559658v1) is a regression and variable selection framework in a multi-environment regression setting. Given a response, StableMate makes distincitions between stable and unstable predictors that have consistent or changing functional dependencies (e.g, fitted regression models) on the response across environments, respectively. Stable predictors entile strong causal implication and can be selected to enforce generalizability of regression models. Unstable predictors (which we refered to as environment-specific predictors in our manuscript) reveals the impact of environment pertubation on the studied system of predictors and response.

StableMate is built based on the theoretical foundation of **Stabilized Regression** [(Pfister et al. 2021)](https://arxiv.org/abs/1911.01850), but implements a different algorithm, **Stochastic Stepwise Variable Selection (ST2)** [(Xin et al. 2012)](https://www.tandfonline.com/doi/abs/10.1080/10618600.2012.679223), to discern stability, since the original algorithm is slow and inacurrate when applied to large scale data. We made further modification on the ST2 algorithm to improve computational efficiency, and proposed a new statistical concept called **pseudo-predictor** that adds on intuitive and accurate benchmark of ST2 selections. 

![haha](./figures/intro.png)

## StableMate analysis
Go to the [analysis](./analysis) folder for the code and data used in [Deng et al. (2023)](https://www.biorxiv.org/content/10.1101/2023.09.26.559658v1). 

## StableMate package
A proper package of StableMate is under development. However, we provide the source code of all SateleMate functions, along with a detailed vigenette, in the [pkg](./pkg) folder. What we have now only comes short in documentatiobn and environment management compared to a proper R package. Other than that, users can easily run StableMate by sourceing, [stablemate_pkg.R](./pkg/stablemate_pkg.R), the R script we provided 

## Reference
Deng, Y., Mao, J., Choi, J., & Le Cao, K. A. (2023). StableMate: a new statistical method to select stable predictors in omics data. bioRxiv, 2023-09.

Pfister, N., Williams, E. G., Peters, J., Aebersold, R., & Bühlmann, P. (2021). Stabilizing variable selection and regression. The Annals of Applied Statistics, 15(3), 1220-1246.

Xin, L., & Zhu, M. (2012). Stochastic stepwise ensembles for variable selection. Journal of Computational and Graphical Statistics, 21(2), 275-294.
