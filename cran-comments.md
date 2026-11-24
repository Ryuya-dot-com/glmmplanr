## R CMD check results

0 errors | 0 warnings | 1 note

* This is the first submission of this package to CRAN.
* `cmdstanr` is listed in `Suggests` (used only when the user explicitly selects
  `engine = "stan_hmc"` with CmdStan) and its repository is declared in
  `Additional_repositories: https://stan-dev.r-universe.dev/`.
* Possibly misspelled words flagged by `win-builder` in `DESCRIPTION`
  (`Boeck`, `De`, `GLMM`, `GLMMs`, `Rasch`, `Westfall`, `al`, `et`,
  `psycholinguistics`, `testlet`) are author surnames (`Westfall`, `De Boeck`),
  the conventional citation abbreviation `et al.`, or standard
  domain/statistical terminology.

## Test environments

* Local macOS (aarch64-apple-darwin25.4.0), R 4.6.1
  - `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran --no-manual`: 0 errors, 0 warnings, 1 note
  - `_R_CHECK_DEPENDS_ONLY_=true R CMD check --no-manual` (optional `Suggests` hidden): 0 errors, 0 warnings
* GitHub Actions (ubuntu-latest, macOS-latest, windows-latest; R-release): 0 errors, 0 warnings
* win-builder (R Under development (unstable), 2026-11-17 r90926 ucrt): 0 errors, 0 warnings, 1 note (total check time 304s; tests 81s; examples 15s; vignettes 10s)
* win-builder (Windows Server 2022, R 4.6.2, 2026-10-29 ucrt): 0 errors, 0 warnings, 1 note (total check time 158s; tests 55s)
