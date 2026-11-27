# Rceattle: Fits the Multispecies Assessment Model (CEATTLE) Using TMB

Implements the CEATTLE model using Template Model Builder ('TMB';
Kristensen et al. 2016
[doi:10.18637/jss.v070.i05](https://doi.org/10.18637/jss.v070.i05) ),
which can be installed following
<https://github.com/kaskr/adcomp/wiki/Download>. Structured similar to
the original manuscript in terms of modularization. Separate functions
estimate retrospective temperature- and size-specific predator rations,
prey preference, and weight-at-age. These are then used as inputs to the
CEATTLE model to evaluate how predation mortality, recruitment, and
survival of three target species change under historical climate
conditions and harvest rates.

## References

- Adams, G.D., Holsman, K.K., Barbeaux, S.J., Dorn, M.W., Ianelli, J.N.,
  Spies, I., Stewart, I.J., Punt, A.E. (2022). An ensemble approach to
  understand predation mortality for groundfish in the Gulf of Alaska.
  Fisheries Research, 251, 106303. doi:10.1016/j.fishres.2022.106303

- Holsman, K.K., Ianelli, J., Aydin, K., Punt, A.E., Moffitt, E.A.
  (2016). A comparison of fisheries biological reference points
  estimated from temperature-specific multi-species and single-species
  climate-enhanced stock assessment models. Deep Sea Research Part II:
  Topical Studies in Oceanography, 134, 360-378.

- Wassermann, S.N., Adams, G.D., Haltuch, M.A., Kaplan, I.C., Marshall,
  K.N., Punt, A.E. (2025). Even low levels of cannibalism can bias
  population estimates for Pacific hake. ICES Journal of Marine Science,
  82(1), fsae064.

- Kristensen, K., Nielsen, A., Berg, C.W., Skaug, H., Bell, B.M. (2016).
  TMB: automatic differentiation and Laplace approximation. Journal of
  Statistical Software, 70(5). doi:10.18637/jss.v070.i05.

## See also

Useful links:

- <https://afsc-assessments.github.io/Rceattle/>

- <https://github.com/afsc-assessments/Rceattle>

- Report bugs at <https://github.com/afsc-assessments/Rceattle/issues>

## Author

**Maintainer**: Grant Adams <grant.adams@noaa.gov>
([ORCID](https://orcid.org/0000-0003-0297-8347))

Authors:

- Grant Adams <grant.adams@noaa.gov>
  ([ORCID](https://orcid.org/0000-0003-0297-8347))

Other contributors:

- Kirstin K. Holsman (original ADMB CEATTLE model; Eastern Bering Sea
  data) \[contributor\]

- James N. Ianelli (original ADMB CEATTLE model; AMAK and ebswp
  conventions) \[contributor\]

- Kerim Aydin (original ADMB CEATTLE model) \[contributor\]

- Cole C. Monnahan (bias-adjustment controls in fit_control())
  \[contributor\]

- Raquel Ruiz Diaz (predicted diet averaged over predator ages, weighted
  by numbers) \[contributor\]

- Tim Miller ([ORCID](https://orcid.org/0000-0003-1411-1206)) (OSA
  composition code, via OSA_multivariate_dists and WHAM) \[contributor\]

- Brian Stock ([ORCID](https://orcid.org/0000-0002-2393-6747)) (WHAM
  co-author (Stock and Miller 2021)) \[contributor\]

- Anders Nielsen (OSA composition code, incl. osa_order())
  \[contributor\]

- Christoffer Moesgaard Albertsen (OSA composition code, from
  OSA_multivariate_dists) \[contributor\]

- Kasper Kristensen (OSA composition code, from OSA_multivariate_dists)
  \[contributor\]

- Chris M. Legault (OSA composition code, from OSA_multivariate_dists)
  \[contributor\]

- Giancarlo M. Correa (growth curves, from WHAM's growth branch; OFI
  2023 examples) \[contributor\]

- Jane Sullivan (whamGrowthData, from the OFI 2023 workshop's WHAM
  examples) \[contributor\]

- James Thorson ([ORCID](https://orcid.org/0000-0001-7415-1010))
  (TMBAIC(), from TMBhelper) \[contributor\]

- Matthew L. H. Cheng (construct_Q(), from GMRF_WAA) \[contributor\]

- Vanessa Trijoulet (OSA composition code; predator-prey vulnerability
  and suitability) \[contributor\]

- Gavin Fay (TMBphase()) \[contributor\]

- Arni Magnusson (rich.colors(), from gplots) \[contributor\]

- Ian G. Taylor (rich.colors.short(), from r4ss) \[contributor\]
