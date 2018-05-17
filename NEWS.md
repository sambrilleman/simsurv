### simsurv 0.2.2-9000 (X/X/2018) 
 
Current development version
 
---

### simsurv 0.2.2 (18/5/2018)

#### New features
   * The `uniroot` function (called internally by `simsurv`) now solves - H(t) - log(u) = 0 instead of exp(-H(t)) - u = 0, i.e. it is on log scale.
   
---

### simsurv 0.2.1 (16/5/2018)

#### Bug fixes
   * The `interval` argument of `simsurv` supports a lower limit of zero (thanks to Alessandro Gasparini, @ellessenne).
   
---

### simsurv 0.2.0 (6/3/2018) 
 
#### New features    
   * Added technical and example usage vignettes
   * Added user-specified cumulative hazard or log cumulative hazard
   * Added analytical forms for the inverted survival function when generating survival times from standard distributions (instead of using numerical root finding). This has lead to about a 5-fold increase in speed when simulating event times from standard parametric distributions.
 
---

### simsurv 0.1.0 (27/7/2017)
  
Initial CRAN release

---

NEWS Template 

### Version (Date) 
 
#### Serious bug fixes   
   * List here if any 
   
#### Bug fixes    
   * List here if any 
   
#### New features    
   * List here if any    
