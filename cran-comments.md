## Submission status

Preparation for the first CRAN submission; no formal CRAN submission has
been made. Win-builder R-release and R-devel uploads were acknowledged;
results are pending. The preceding source passed three-platform GitHub
Actions. The current candidate adds plot styling and annotation controls;
its local checks are recorded separately below.

The post-correction default-bootstrap confirmation completed 3,072 outer
datasets and 361,728 retained inner fits, all valid. The validation guide
retains indeterminate planning/calibration checks and an adverse nominal-null
flag. Successful fitting is not presented as uniform inferential calibration.

## Current local candidate

Source archive: `glmmplanr_0.1.0.tar.gz` (1,240,738 bytes).
SHA-256: `5a80fc0b2858a5479ef166f90428e218dacb5cade833218324e58989cce63c5f`.

0 errors | 0 warnings | 1 note

* macOS, R 4.6.1, `aarch64-apple-darwin23`.
* The plot-update check used `--as-cran --no-manual --no-tests`,
  with `_R_CHECK_FORCE_SUGGESTS_=false`,
  `_R_CHECK_DEPENDS_ONLY_=true`, `_R_CHECK_CRAN_INCOMING_REMOTE_=true`, and
  `_R_CHECK_SYSTEM_CLOCK_=false`.
* The NOTE reports `New submission` and `cmdstanr` in an additional
  repository; its availability was confirmed. Installation, remote URL
  checks, examples, help files and vignette reconstruction passed.
* The checked archive has SHA-256
  `649b428ab2129e642e3e3310de18191bfe042e63e8f90659a54e7a90fff1f3bc`.
  The final archive corrects one existing section anchor's placement in the
  validation guide and its HTML. Other differences are DESCRIPTION's build
  metadata and three printed environment addresses in the planning-workflow
  HTML. No code, tests, help, examples or vignette source changed after the
  check; the final build reconstructed the vignettes successfully.
* The targeted plot/plan-method tests passed 79 assertions with zero
  failures, warnings or skips. Ten static plot variants were rendered;
  annotation toggles preserve all built layer data. Plotly layout and trace
  data were checked, without a browser interaction test. Saved inputs and
  the regenerated default-bootstrap figure's numerical summaries are unchanged.
* The earlier full regression passed 713 assertions with zero failures or
  warnings and one opt-in Stan skip. Source comparison confirms that native
  kernels, estimators, optimization and uncertainty calculations are unchanged
  by this display update. The full regression was not repeated.
* The final archive's 37 R/C++/test files match the current source. These are
  combined stage results, not a full single-command check of the final archive.
  The optional `instantiate` package was unavailable. All 13 README-relative
  links resolve within the archive. Compiler/compiled-code checks passed.
* Plot-update logs, figures and archive/source identity are retained in
  `validation/plot_defaults/`. Earlier checks are in
  `validation/release_candidate/` and `validation/l2_uncertainty_comparison/`
  (all excluded from the CRAN tarball).

## Current external check requests

* Win-builder R-release: acknowledged upload of the pre-link-fix archive
  `219ebc6c...`; results pending.
* Win-builder R-devel: acknowledged upload of the preceding archive `b9b025ca...`;
  results pending. Upload acknowledgment is not a passing check.
  Neither upload includes the new plot defaults; no duplicate request was sent.
* GitHub Actions: preceding-source macOS, Windows and Linux, R 4.6.2, all
  `Status: OK`. Each platform reports 706 passed assertions, zero failures
  or warnings, and one opt-in Stan skip. Examples and vignette reconstruction
  passed on all three. CI uses `NOT_CRAN=true` and disables incoming
  feasibility checks; the separate local archive check enables them.
  Candidate source: `29aadf4a99016ed21fa1ffca6bc710369368ed6f`.
  Run: https://github.com/Ryuya-dot-com/glmmplanr/actions/runs/45795084138
  This is evidence for the unchanged numerical implementation, not a
  three-platform check of the new plotting code.

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
