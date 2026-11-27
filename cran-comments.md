## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission (`glmmplanr` 0.1.0; `1,319,579` bytes,
  SHA-256 `2ac8769f78cc41c8313214018d1971733a8e26fb43c4b8b2a6cccdb6533db3e6`).
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
  * Full `NOT_CRAN=true` test suite: 1,328 assertions passed, 0 failures, 0 warnings
  * `tools:::check_url_db()` on built source archive (90 unique URLs across `DESCRIPTION`, `README.md`, `NEWS.md`, `man/*.Rd`, and `inst/doc/*.html`): 0 invalid URLs
* Pre-submission external checks to record on the final release commit before uploading to CRAN:
  * GitHub Actions (`ubuntu-latest`, `macos-latest`, `windows-latest`; R-oldrel, R-release, R-devel): pending push of final commit
  * Win-builder (`R-release` and `R-devel`, Windows Server x86_64): pending upload of final `glmmplanr_0.1.0.tar.gz`
