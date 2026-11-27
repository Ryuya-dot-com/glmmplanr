#!/bin/sh
# Usage:
#   sh validation/run-isolated.sh SCRIPT_NAME [SCRIPT_ARGS...]
# Set VALIDATION_COPY_ALL=1 to copy the full archived validation/ tree (~213 MB)
# instead of staging only validation/*.R and script-specific input dependencies.
set -eu
if [ "$#" -lt 1 ]; then
  echo "Usage: sh validation/run-isolated.sh SCRIPT_NAME [SCRIPT_ARGS...]" >&2
  exit 2
fi
name=$1
shift
case "$name" in ''|*[!a-zA-Z0-9_]*) echo "Invalid validation script name" >&2; exit 2 ;; esac
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test -f "$root/validation/$name.R" || { echo "Unknown validation script: $name" >&2; exit 2; }
run=$(mktemp -d "${TMPDIR:-/tmp}/glmmplanr-validation.XXXXXXXX")
echo "Validation directory: $run"
mkdir -p "$run/source/validation" "$run/library"
cd "$root"
cp -R DESCRIPTION NAMESPACE LICENSE R src man inst tests "$run/source/"

stage_rel() {
  for rel in "$@"; do
    if [ -e "$root/validation/$rel" ]; then
      parent=$(dirname "$run/source/validation/$rel")
      mkdir -p "$parent"
      cp -R "$root/validation/$rel" "$parent/"
      echo "validation/$rel" >> "$run/staged_inputs.txt"
    fi
  done
}

: > "$run/staged_inputs.txt"
if [ "${VALIDATION_COPY_ALL:-0}" = "1" ]; then
  cp -R validation "$run/source/"
  echo "validation/ (full copy via VALIDATION_COPY_ALL=1)" >> "$run/staged_inputs.txt"
else
  cp validation/*.R "$run/source/validation/"
  echo "validation/*.R" >> "$run/staged_inputs.txt"
  dep_list=$(Rscript --vanilla - "$root" "$name" <<'EOF'
args <- commandArgs(TRUE)
root <- args[1L]
name <- args[2L]
val_dir <- file.path(root, "validation")
existing <- setdiff(list.files(val_dir, all.files = FALSE), c("README.md"))
existing_non_r <- existing[!grepl("\\.(R|sh|py)$", existing)]
safe_strs <- function(pd) {
  raw <- pd$text[pd$token == "STR_CONST"]
  unlist(lapply(raw, function(s) {
    tryCatch(as.character(eval(parse(text = s))), error = function(e) character(0))
  }))
}
entry <- file.path(val_dir, paste0(name, ".R"))
visited <- character(0)
queue <- entry
all_strs <- character(0)
raw_tokens <- character(0)
exclude <- character(0)
is_reader_script <- grepl("^(summarize|plot|check|compare|confirm|inspect)_", name)
while (length(queue) > 0L) {
  f <- queue[1L]; queue <- queue[-1L]
  if (f %in% visited || !file.exists(f)) next
  visited <- c(visited, f)
  pd <- getParseData(parse(f, keep.source = TRUE))
  strs <- safe_strs(pd)
  raw_tokens <- c(raw_tokens, strs)
  queue <- c(queue, file.path(root, grep("^validation/.*\\.R$", strs, value = TRUE)))
  lines <- readLines(f, warn = FALSE)
  creates_out <- any(grepl("(dir\\.create|checkpoint_dir)\\s*\\(", lines))
  if (identical(f, entry) && creates_out && !is_reader_script) {
    for (ol in grep("^\\s*(out|out_dir)\\s*<-", lines, value = TRUE)) {
      opd <- tryCatch(getParseData(parse(text = ol, keep.source = TRUE)), error = function(e) NULL)
      if (!is.null(opd)) exclude <- c(exclude, grep("^validation/", safe_strs(opd), value = TRUE))
    }
  }
  for (wl in grep("(write\\.csv|saveRDS|ggsave|png|pdf)\\s*\\(", lines, value = TRUE)) {
    wpd <- tryCatch(getParseData(parse(text = wl, keep.source = TRUE)), error = function(e) NULL)
    if (!is.null(wpd)) exclude <- c(exclude, grep("^validation/", safe_strs(wpd), value = TRUE))
  }
  all_strs <- c(all_strs, grep("^validation/", strs, value = TRUE))
}
# Match constructed basenames/prefixes in validation/ (e.g. file.path(out, paste0("stem_", stage)))
matched_rel <- character(0)
toks <- unique(sub("^validation/+", "", raw_tokens[nzchar(raw_tokens)]))
toks <- toks[!toks %in% c("validation", "pilot", "large", "medium", "smoke", "summary.csv", "fits.csv", "fits.csv.gz", "diagnostics.csv", "manifest.rds")]
for (tk in toks) {
  if (tk %in% existing_non_r) matched_rel <- c(matched_rel, tk)
  if (nchar(tk) >= 6L && grepl("_", tk)) {
    hits <- existing_non_r[startsWith(existing_non_r, tk) | startsWith(existing_non_r, paste0(tk, "_"))]
    matched_rel <- c(matched_rel, hits)
  }
}
if (is_reader_script) {
  stem <- sub("^(summarize|plot|check|compare|confirm|inspect)_", "", name)
  if (nchar(stem) >= 4L) {
    hits <- existing_non_r[existing_non_r == stem | startsWith(existing_non_r, paste0(stem, "_"))]
    hits <- hits[!grepl("\\.(png|pdf)$", hits)]
    matched_rel <- c(matched_rel, hits)
  }
}
norm_p <- function(x) {
  if (!length(x)) return(character(0))
  v <- vapply(unique(x), function(s) if (grepl("%", s)) dirname(sub("%.*$", "", s)) else s, character(1L), USE.NAMES = FALSE)
  unique(sub("/+$", "", v))
}
cands <- unique(c(setdiff(norm_p(all_strs), norm_p(exclude)), paste0("validation/", matched_rel)))
cands <- setdiff(cands, norm_p(exclude))
cands <- cands[!grepl("\\.R$", cands) & cands != "validation"]
rels <- unique(sub("^validation/", "", cands))
rels <- rels[file.exists(file.path(val_dir, rels))]
for (r in rels) cat(r, "\n", sep = "")
EOF
  )
  if [ -n "$dep_list" ]; then
    OLD_IFS=$IFS
    IFS='
'
    for rel in $dep_list; do
      [ -n "$rel" ] && stage_rel "$rel"
    done
    IFS=$OLD_IFS
  fi
fi

cd "$run/source"
export R_LIBS="$run/library${R_LIBS:+:$R_LIBS}"
R CMD INSTALL --preclean --clean --library="$run/library" . > "$run/install.log" 2>&1 ||
  { cat "$run/install.log" >&2; exit 1; }
Rscript --vanilla -e '
  args <- commandArgs(TRUE)
  stopifnot(normalizePath(find.package("glmmplanr")) == normalizePath(file.path(args[1], "glmmplanr")))
  print(sessionInfo())
  for (p in c("glmmplanr", "lme4", "TAM", "rstan", "cmdstanr", "simr", "PowRPriori"))
    if (nzchar(system.file(package=p))) cat(p, as.character(packageVersion(p)), "\n")
' "$run/library" > "$run/environment.txt"

if [ "$#" -eq 0 ]; then
  if [ "$name" = coverage_recheck ]; then
    set -- "$run/source/validation"
  elif [ "$name" = covariance_mode_replay ]; then
    set -- --fresh "$run/replay_out"
  fi
fi

Rscript --vanilla "validation/$name.R" "$@" > "$run/run.log" 2>&1 ||
  { tail -n 30 "$run/run.log" >&2; exit 1; }
echo "Results: $run/source/validation"
if [ -d "$run/replay_out" ]; then echo "Fresh replay output: $run/replay_out"; fi
echo "Staged inputs: $run/staged_inputs.txt"
echo "Console: $run/run.log"
