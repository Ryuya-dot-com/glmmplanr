## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission (`glmmplanr` 0.1.0).
* `cmdstanr` is an optional dependency in `Suggests` (used only when a user
  explicitly selects `engine = "stan_hmc"` with CmdStan) and is available from
  `Additional_repositories: https://stan-dev.r-universe.dev/`. All examples,
  vignettes, and tests run without `cmdstanr` or `instantiate` installed
  (verified with `_R_CHECK_DEPENDS_ONLY_=true`).
* Words flagged by spell checking in `DESCRIPTION` (`Boeck`, `De`, `GLMM`,
  `GLMMs`, `Rasch`, `Westfall`, `al`, `et`, `psycholinguistics`, `testlet`) are
  author surnames (`Westfall`, `De Boeck`), the standard abbreviation `et al.`,
  or established statistical/psychometric terms.

## Test environments

* Local macOS (`aarch64-apple-darwin23`), R 4.6.1:
  * `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran --no-manual`: 0 errors, 0 warnings, 1 note
  * `_R_CHECK_DEPENDS_ONLY_=true R CMD check --no-manual`: 0 errors, 0 warnings, 0 notes
* GitHub Actions (`ubuntu-latest` R-oldrel / R-release / R-devel, `macos-latest` R-release, `windows-latest` R-release):
  * 0 errors, 0 warnings across all 5 matrix configurations
* Win-builder (Windows Server 2022 x64, GCC 14.3.0):
  * R Under development (unstable) (`2026-11-17 r90926 ucrt`): 0 errors, 0 warnings, 1 note
  * R 4.6.2 (`2026-10-29 ucrt`): 0 errors, 0 warnings, 1 note
