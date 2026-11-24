#!/bin/sh
# Usage: sh validation/run-isolated.sh rasch_fixed_bank
set -eu
if [ "$#" -ne 1 ]; then
  echo "Usage: sh validation/run-isolated.sh SCRIPT_NAME" >&2
  exit 2
fi
name=$1
case "$name" in ''|*[!a-zA-Z0-9_]*) echo "Invalid validation script name" >&2; exit 2 ;; esac
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test -f "$root/validation/$name.R" || { echo "Unknown validation script: $name" >&2; exit 2; }
run=$(mktemp -d "${TMPDIR:-/tmp}/glmmplanr-validation.XXXXXXXX")
echo "Validation directory: $run"
mkdir "$run/source" "$run/library"
cd "$root"
cp -R DESCRIPTION NAMESPACE LICENSE R src man inst tests validation "$run/source/"
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
set --
if [ "$name" = coverage_recheck ]; then set -- "$run/source/validation"; fi
Rscript --vanilla "validation/$name.R" "$@" > "$run/run.log" 2>&1 ||
  { tail -n 30 "$run/run.log" >&2; exit 1; }
echo "Results: $run/source/validation"
echo "Console: $run/run.log"
