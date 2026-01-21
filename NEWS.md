# simsurv 1.0.1.9000
 
#### New features
   * Current development version

---

# simsurv 1.0.1 (2026-01-21)
 
#### Bug fixes
   * Fixed NEWS.md format to comply with CRAN requirements (ISO 8601 dates and proper heading levels)

---

# simsurv 1.0.0 (2021-01-09)
 
#### New features
   * Version number bumped to 1.0.0 to correspond to publication of JSS paper.

---

# simsurv 0.2.5 (2019-02-22)
 
#### Bug fixes
   * Fix bug where infinite survival time were incorrectly returned as NaN (due to: Inf * d = NaN when d = 0 is for a censored individual)

---

# simsurv 0.2.4 (2019-02-06)

#### New features
   * Allow negative shape parameters for the Gompertz distribution (note that this can lead to infinite survival times, in which case the survival time is set to 'Inf' and a warning is printed).

#### Bug fixes
   * The simsurv() function documentation and the technical vignette have been updated to better clarify the parameterisations for the parametric distributions.
	 
---

# simsurv 0.2.3 (2019-02-01)

#### New features
   * The `rootfun` argument has been added. This allows the user to apply any transformation to each side of the root finding equation. The default is `rootfun = log` which corresponds to `uniroot` solving -H(t) - log(u) = 0.
   * The `rootsolver` argument has been added. This allows the user to choose between stats::uniroot or BB::dfsane for the root finding.
	 
---

# simsurv 0.2.2 (2018-05-18)

#### New features
   * The `uniroot` function (called internally by `simsurv`) now solves - H(t) - log(u) = 0 instead of exp(-H(t)) - u = 0, i.e. it is on log scale.
   
---

# simsurv 0.2.1 (2018-05-16)

#### Bug fixes
   * The `interval` argument of `simsurv` supports a lower limit of zero (thanks to Alessandro Gasparini, @ellessenne).
   
---

# simsurv 0.2.0 (2018-03-06) 
 
#### New features    
   * Added technical and example usage vignettes
   * Added user-specified cumulative hazard or log cumulative hazard
   * Added analytical forms for the inverted survival function when generating survival times from standard distributions (instead of using numerical root finding). This has lead to about a 5-fold increase in speed when simulating event times from standard parametric distributions.
 
---

# simsurv 0.1.0 (2017-07-27)
  
Initial CRAN release    
