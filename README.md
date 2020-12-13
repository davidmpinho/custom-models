# Custom models

These are some of the models I implemented in [Stan](https://mc-stan.org/),
partly for easy reference, and partly to share models that I have developed or
extended.

## Models currently implemented 

* **Autoregressive Conditional Poisson/Negative-Binomial (ARCP/ARCNG)**,
 which works similarly to GARCH models, but for counts instead of squared
observations).  The advantage of this model over structural time series is that
it only needs 2-4 parameters to model the 'residual' component of a time series
-- the component that is (ideally) stationary and cannot be explained away by
other covariates.  This model was originally developed by 
[Heinen's (2003)](http://dx.doi.org/10.2139/ssrn.1117187) article.  I added the
following changes and extensions:
    - Allow for complex multiplicative seasonality and external predictors
      (i.e., exogenous terms). 
    - Reparameterization of many variables to make the model easier to
      estimate, which makes the use of variational inference more viable.
      (The reparameterizations can also be applied to GARCH models.)

* **Bayesian structural time series**. I base this on the code of 
[Phalen](https://peterphalen.github.io/ceasefire/bsts). I
optimized it and simplified it. The model struggles with larger samples (>>2000),
but it tends to be better when doing inference.

## Future goals

* Add generated quantities block to all models (to do posterior predictive
checks and generate predictions). 

* For the ARCP/ARCNG model:
    - Let the unconditional mean ("omega") vary through time to allow for
      trends. Currently, the sum of alpha and beta parameters is always too high
      when modeling non-stationary time series. The plan is to use splines to
      model omega, although they cannot be too flexible. Something similar to the
      [C-GARCH](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5848) can also 
      work, although that is bad for extrapolating -- the process always reverts 
      to a long-run mean. 
    - Experiment with "time before/after event" factors. These are only useful for 
      events with some irregularity (like sales before/after Easter) instead of 
      events (like Christmas). I need a mixture of the functions "additive_trend" 
      and "additive_factor". 
    - Allow for overdispersion in the seasonality. Right now I need to specify
      holiday dates for the yearly seasonality work well. But even there,
      the seasonality can still struggle in the days before those dates. The solution
      is to add seasonality innovations (the "epsilon" variables) with student-t(nu,
      0, sigma), where nu is somewhere between 3 and 7 (or it could be estimated),
      and sigma is then chosen so that the standard deviation of the innovations is
      one.  
    - Support for multiple time series. To be complete, this requires two
      things: allow for a multivariate residual component, and hierarchical
      relationships between seasonality and trends. This last one is done in the
      following way (for yearly seasonality as an example): instead of estimating one
      yearly trend for every time series, first estimate the average trend, and then
      allow for hierarchically-structured deviations from that mean; only after that
      do you use indexing in the code. This could be applied to splines in a similar
      way, I think, but it is less straightforward in that case.  
    - Allow missing data. This is simple: in the loop, if y is not
      missing proceed as normal; if it is missing use the expectation (which is
      easy to compute, there is no need to add another set of parameters to estimate
      missing variables).

* For time series more generally, apply what I learned with ARC model to ARMA
models. The motivation is to have something like 
[Prophet](https://facebook.github.io/prophet/), but with better performance and
without all the tacky choices that they have -- like having to specify bounds,
no support for count variables or hierarchical/multivariate relationships,
modeling seasonality with Fourier series, or having two 'trends'.  So the
objective is to have a middle term between conventional ARMA-type models (which
work well enough in most cases but are very limiting) and black-box models like
neural networks or tree-based models (which do not even work that well most of
the time).

* Gaussian random field approximation with splines (i.e., for spatial data).
Like before, the motivation is to have a middle term between the
classical approach (estimate a full covariance matrix) and the black-box
algorithm approach. 

* "Default" splines for regression. Currently, that I know of, there is no easy
way to use p-splines (penalized splines) in typical regression problems where
there are many variables and interactions that are plausibly non-linear.
[Kharratzadeh](https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html)
has shown how to do this to some extent, but the indexes and knot locations
have to be specified outside of the code itself in a very manual way. I would
like something that is about as easy to fit as polynomial regression, while
also allowing for regularization (and not having the drawbacks of polynomial
regression in general). 

* Hierarchical splines. The idea is to take the simple "default" splines and
allow deviations from them with a hierarchical structure. 

 
