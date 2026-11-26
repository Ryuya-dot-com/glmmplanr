## Submission status

Preparation for the first CRAN submission; no formal CRAN submission has
been made. Win-builder R-release and R-devel uploads were acknowledged;
results are pending. The three-platform GitHub Actions check is being
prepared on the release branch. Historical checks are not current evidence.

The post-correction default-bootstrap confirmation completed 3,072 outer
datasets and 361,728 retained inner fits, all valid. The validation guide
retains indeterminate planning/calibration checks and an adverse nominal-null
flag. Successful fitting is not presented as uniform inferential calibration.

## Current local candidate

Source archive: `glmmplanr_0.1.0.tar.gz` (1,254,503 bytes).
SHA-256: `b9b025ca54039867b32795661b63c20f02f9f50631b72c98e19002eb07b31087`.

0 errors | 0 warnings | 1 note

* macOS, R 4.6.1, `aarch64-apple-darwin23`.
* The final difference check used `--as-cran --no-manual --no-tests
  --no-examples --no-vignettes`, with `_R_CHECK_FORCE_SUGGESTS_=false`,
  `_R_CHECK_DEPENDS_ONLY_=true`, `_R_CHECK_CRAN_INCOMING_REMOTE_=true`, and
  `_R_CHECK_SYSTEM_CLOCK_=false`.
* The NOTE reports `New submission` and `cmdstanr` in an additional
  repository; its availability was confirmed. Remote URL checking passed
  after correcting two canonical CRAN links and two unpublished figure links.
* The preceding archive (SHA-256
  `219ebc6c952997324b19c1804b381b901846bece80c6ac35f28ccc94af94c0c6`)
  passed installation, examples and vignette reconstruction. The final
  archive differs only in DESCRIPTION's packaging metadata, the validation
  guide and its generated HTML, and three printed environment addresses in
  the otherwise identical planning-workflow HTML. The final build also
  reconstructed the vignettes successfully.
* A separate regression run passed 713 assertions with zero failures or
  warnings; one opt-in Stan test was skipped. The archive's 37 R/C++/test
  files match that tested source. These are combined stage results, not a
  full single-command check of the final archive.
* The optional `instantiate` package was unavailable. All 12 README-relative
  links resolve within the archive. Compiler/compiled-code checks passed.
* Logs and archive/source identity are retained in
  `validation/release_candidate/`; the full regression log is in
  `validation/l2_uncertainty_comparison/` (both excluded from the CRAN tarball).

## Current external check requests

* Win-builder R-release: acknowledged upload of the pre-link-fix archive
  `219ebc6c...`; results pending. Its R/C++/test files match the current source.
* Win-builder R-devel: acknowledged upload of the final archive `b9b025ca...`;
  results pending. Upload acknowledgment is not a passing check.
* GitHub Actions: current-source macOS, Windows and Linux checks pending.

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
