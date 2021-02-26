# Custom models

These are some of the models I implemented in [Stan](https://mc-stan.org/),
partly for easy reference, and partly to share models that I have developed or
extended.

Feel free to ask me questions or make pull requests.

## Models currently implemented 

### Regression

* [Fast piecewise linear function (similar to first-degree splines)](https://github.com/davidmpinho/custom-models/blob/main/regression/fast_1d_splines.stan).
The objective with these is to specify a large number of knots (>30) to approximate a smooth function. When using many knots, it's faster
than trying to do it with [Kharratzadeh's](https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html) spline implementation because the
coefficients at the knots are estimated in an identical way Bayesian structural time series models, with points between parameters being interpolated 
in a similar way to first-degree splines. This scales better and is particularly advantageous when setting priors: there is a more intuitive
relationship between the local and global flexibilities of the function because, a priori, it behaves like a random walk. 
In my particular implementation, I've made the priors 
(mostly) invariant to changes in the number of knots, such that, e.g., sigma=0.5 is equal to a change of 0.5 in Y for every 1 standard deviation change in X. 
The disadvantage of this method is that interpolation to values outside the range is not really possible (but that's a tricky subject for every model anyways). 
It also doesn't work very well with exponential relationships because the piecewise slopes have (on average) a constant 'rate of change', and that causes 
overfitting in the ranges of values where Y is relatively closer to zero. You can fix that by just exponentiating phi if Y > 0. I tried to find a more 
general solution, one that would work when Y has negative values. Because I like this piecewise linear function idea, I put a piecewise linear function inside 
my piecewise linear function: the flexibility of the spline is allowed to change over the range of X in an exponential way
([fast_1d_splines_exp.stan](https://github.com/davidmpinho/custom-models/blob/main/regression/fast_1d_splines.stan)). Seems 
to take ~50% longer to fit on toy data sets. But I haven't really tested this, so proceed with caution (always!).

### Time series

* [Autoregressive Conditional Poisson/Negative-Binomial (ARCP/ARCNG)](https://github.com/davidmpinho/custom-models/blob/main/time_series/autoregressive_conditional.stan),
 which works similarly to GARCH models, but for counts instead of squared
observations).  The advantage of this model over structural time series is that
it only needs 2-4 parameters to model the 'residual' component of a time series
-- the component that is (ideally) stationary and cannot be explained away by
other covariates.  This model was originally developed by 
[Heinen's (2003)](http://dx.doi.org/10.2139/ssrn.1117187) article. I added the
following changes and extensions:
    - Allow for complex multiplicative seasonality and external predictors
      (i.e., exogenous terms). 
    - Reparameterization of many variables to make the model easier to
      estimate, which makes the use of variational inference more viable.
      (The reparameterizations can also be applied to GARCH models.)

* [Bayesian structural time series](https://github.com/davidmpinho/custom-models/blob/main/time_series/bayes_state_space.stan). I base this on the code of 
[Phalen](https://peterphalen.github.io/ceasefire/bsts). I
optimized it and simplified it. The model struggles with larger samples (>>2000) when using
HMC, but tends to be more robust when doing inference in general problems -- there are no strong assumptions
about how the latent state should evolve.

## Reparameterizations

* **Student-t**: this is different from the reparameterization from the 
[Stan User's guide](Reparameterization). In my case, one parameter sets 
the standard deviation of the whole distribution, while the other one 
sets (what would be) the standard deviation of a standard student-t, 
which maps directly to the parameter nu.
With this, the problematic posterior is (mostly) gone, though there is  
still see so. 
* **Beta and alpha parameters in GARCH/ARCP models**: this is the one I
mentioned in the previous section. Here I eliminate all the unrelated code. 


## Future goals

* "Default" splines for regression. Currently, that I know of, there is
no easy way to use penalized splines in a principled way. 
[Kharratzadeh](https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html)
has shown how to do this to some extent, but I see some issues:
  * indexes and knot locations need to be specified outside of the
  code itself in a very manual way;
  * there is no principled way of setting knots; we could use 
  many of them to get an overparameterized model and then regularize, but...
  * ...using too many knots is slow, and the priors have to be adapted to the
  number of knots. 

* Hierarchical splines. The idea is to take the simple "default" splines and
allow deviations from them with a hierarchical structure. Check the work 
of Dr. Gavin Simpson on this. 


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
      
* For time series, model events (e.g., days before/after Christmas).
This is mostly 'feature engineering' -- indexes 1, 2, ..., n/2-1 are 
the days before an event, day n/2 is the event, and n/2+1, ..., n are the 
days after the event. I think I only need to have separate scales for the 
days before and after. That said, if there are issues of identifiability, 
either restrict the number of days that are being modelled (probably necessary
anyways because of the factor model assumptions), or add a positive 
(and negative) 'slope' to the random walk process in the days before 
(and after). 

* Factor model (yes, again) for survival analysis. Just restrict the 
differences to be negative and model them on the log scale. Mind the 
censoring, though.  

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

* Add generated quantities block to all models (to do posterior predictive
checks and generate predictions). 
 
