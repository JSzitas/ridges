
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ridges

<!-- badges: start -->
<!-- badges: end -->

## About

**Work in Progress**

A tiny ridge regression library for multiple targets, with optional
automatic selection of the ridge parameter (Using [the estimator ‘k36’
from this paper (see
references)](article.ajtas.org/pdf/10.11648.j.ajtas.20211002.13.pdf),
for all, not only Poisson regression cases. )

Obviously this is highly experimental and has no mathematical
guarantees… but I try to have a curated set of benchmarks to showcase
whether it works in practice.

The goal here is absolutely **NOT** to provide an interface like **glm**
or **glmnet**, but rather to have a lean, easy to understand set of
solvers, with models capable of producing predictions. These solvers
should be easy to plug-in where-ever you need to solve a ridge
regression problem, and not much else.

It also serves as a simple educational tool - since the implementation
is primarily concerned with simplicity and does no unnecessary
calculations, it can be used to reasonably showcase ridge regression in
a generalized linear model context. While elastic-net could be
considered a better approach in many contexts, I believe most
implementations are opaque and hard to understand.

To the goal of simplicity, this package currently supports no
scaling/re-scaling, no addition of intercepts, no formulas, and
absolutely no offsets. The scaling and intercepts **might** be added at
some point, but I would prefer to keep this contained as a very
bare-bones implementation.

## Installation

Currently only from github.

## References

-   [Etaga, Harrison & Florence, Aforka & Abidemi, Awopeju & Njideka,
    Etaga. (2021). “Poisson Ridge Regression Estimators: A Performance
    Test.” American Journal of Theoretical and Applied Statistics. 10.
    111-121.
    10.11648/j.ajtas.20211002.13.](article.ajtas.org/pdf/10.11648.j.ajtas.20211002.13.pdf)
