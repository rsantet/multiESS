# multiESS

Compute the multivariate effective sample size (*mESS*) of a Markov chain, using the multivariate dependence structure of the process.

This is a Julia implementation based on the [Python implementation][1] by Gabriel Perren, itself based on the [MATLAB implementation][2] by Luigi Acerbi, of the *mESS* estimation method described in [Vats et al. (2015)][3]. Both the MATLAB and Python codes (and consequently this one) have some minor tweaks for the choice of batch size *b* for the computation of the Monte Carlo covariance matrix.

See the R package [mcmcse][4] for a separate implementation.

*Disclaimer:* This is a stripped down version of the original MATLAB code. Most notably it does neither accept nor return the Monte Carlo covariance matrix and does not return the batch size. It also can process only *one* chain at a time.

## Details

The effective sample size of a Markov chain is the size of an i.i.d. sample with the same covariance structure as the current chain.

*mESS* is given by

    mESS = n |Λ|^{1/p}/ |Σ|^{1/p}

where *n* is the current sample size, Λ is the sample covariance matrix, *p* is the number of parameters and Σ is an estimate of the Monte Carlo covariance matrix for the Markov chain (here obtained by batch estimation).

### Reference

[Vats, D., Flegal, J. M., & Jones, G. L. "Multivariate Output Analysis for Markov chain Monte Carlo", arXiv preprint arXiv:1512.07713 (2015)][3]

[1]: https://github.com/Gabriel-p/multiESS
[2]: https://github.com/lacerbi/multiESS
[3]: http://arxiv.org/abs/1512.07713
[4]: https://cran.r-project.org/web/packages/mcmcse/index.html