## Submission status

Preparation for the first CRAN submission. The current source has not yet
been checked by win-builder or the three-platform GitHub Actions workflow.
The historical results below must not be reported as checks of this source.
The post-correction default-bootstrap validation is in progress.

## Current local candidate

Source archive: `glmmplanr_0.1.0.tar.gz` (1,252,586 bytes).
SHA-256: `a170a0d4cd7a70ae4beda98d4c471fc882bd10bead3efab0ab8d109256158400`.

0 errors | 0 warnings | 1 note

* macOS, R 4.6.1, `aarch64-apple-darwin23`.
* The archive check used `--as-cran --no-manual --no-tests`, with
  `_R_CHECK_FORCE_SUGGESTS_=false`, `_R_CHECK_DEPENDS_ONLY_=true`,
  `_R_CHECK_CRAN_INCOMING_REMOTE_=false`, and `_R_CHECK_SYSTEM_CLOCK_=false`.
* The NOTE concerns the Apple developer-tool temporary cache `xcrun_db`.
  Package installation, examples, and vignette reconstruction passed.
* A separate regression run passed 713 assertions with zero failures or
  warnings; one opt-in Stan test was skipped. The archive's 37 R/C++/test
  files match that tested source. This is combined stage evidence, not a
  full single-command CRAN check of the archive.
* The optional `instantiate` package was unavailable. The archive's 12
  README-relative links resolve within the archive.
* Logs and source/archive identity are retained in
  `validation/l2_uncertainty_comparison/` (excluded from the CRAN tarball).

## Package notes

* `cmdstanr` is listed in `Suggests` (used only when the user explicitly selects
  `engine = "stan_hmc"` with CmdStan) and its repository is declared in
  `Additional_repositories: https://stan-dev.r-universe.dev/`.
* Possibly misspelled words flagged by `win-builder` in `DESCRIPTION`
  (`Boeck`, `De`, `GLMM`, `GLMMs`, `Rasch`, `Westfall`, `al`, `et`,
  `psycholinguistics`, `testlet`) are author surnames (`Westfall`, `De Boeck`),
  the conventional citation abbreviation `et al.`, or standard
  domain/statistical terminology.

## Historical environments (earlier source, not current-candidate evidence)

* Local macOS (aarch64-apple-darwin25.4.0), R 4.6.1
  - `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran --no-manual`: 0 errors, 0 warnings, 1 note
  - `_R_CHECK_DEPENDS_ONLY_=true R CMD check --no-manual` (optional `Suggests` hidden): 0 errors, 0 warnings
* GitHub Actions (ubuntu-latest, macOS-latest, windows-latest; R-release): 0 errors, 0 warnings
* win-builder (R Under development (unstable), 2026-11-17 r90926 ucrt): 0 errors, 0 warnings, 1 note (total check time 304s; tests 81s; examples 15s; vignettes 10s)
* win-builder (Windows Server 2022, R 4.6.2, 2026-10-29 ucrt): 0 errors, 0 warnings, 1 note (total check time 158s; tests 55s)
