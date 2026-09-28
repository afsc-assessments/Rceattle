# Rceattle: Fits the Multispecies Assessment Model (CEATTLE) Using TMB

Implements the CEATTLE model using Template Model Builder ('TMB';
Kristensen et al. 2015), which can be installed following
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

- Kristensen, K., Nielsen, A., Berg, C.W., Skaug, H., Bell, B. (2015).
  TMB: automatic differentiation and Laplace approximation. arXiv
  preprint arXiv:1509.00660.

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
