"""Independent arithmetic audit; R is used only to decode the frozen RDS draws."""
import csv
import gzip
import json
import math
import statistics as st
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

root = Path(__file__).resolve().parents[1]
out = Path(sys.argv[1]) if len(sys.argv) > 1 else root / "validation/interval_mechanism_audit"
normal = st.NormalDist()
z = normal.inv_cdf(.975)


def read(name):
    path = out / name
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return list(csv.DictReader(f))


def key(row):
    return row["case"], row["condition"], row["replicate"]


def close(a, b, tol=1e-10):
    assert abs(float(a) - float(b)) < tol, (a, b)


def wilson(k, n):
    p = k / n
    den = 1 + z*z/n
    mid = (p + z*z/(2*n)) / den
    half = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return mid-half, mid+half


rows = read("intervals.csv.gz")
index = {(*key(v), v["method"]): v for v in rows}
assert len(rows) == len(index) == 12864
groups = defaultdict(list)
for v in rows:
    valid = v["valid"] == "TRUE"
    lo, hi, truth = map(float, (v["lower"], v["upper"], v["truth"]))
    assert valid and math.isfinite(lo) and math.isfinite(hi) and lo < hi
    cover = lo <= truth <= hi
    detect = lo > 0 or hi < 0
    width = hi-lo <= 10
    for name, value in [("covered",cover),("detected",detect),("width_ok",width),("joint",detect and width)]:
        assert (v[name] == "TRUE") == value
    groups[(v["case"],v["condition"],v["method"])].append(v)
for s in read("interval_summary.csv"):
    data = groups[(s["case"],s["condition"],s["method"])]
    assert len(data) == int(s["n"]) == int(s["n_valid"])
    for metric, name in [("coverage","covered"),("detection","detected"),("joint","joint")]:
        count = sum(v[name] == "TRUE" for v in data)
        assert count == int(s["n_" + name])
        close(count/len(data),s[metric])
        lo, hi = wilson(count,len(data))
        close(lo,s[metric+"_lo"]); close(hi,s[metric+"_hi"])
    expected = "pass" if float(s["joint_lo"]) >= .8 else "fail" if float(s["joint_hi"]) < .8 else "indeterminate"
    assert expected == s["joint_decision"]
for s in read("paired.csv"):
    a = groups[(s["case"],s["condition"],s["a"])]
    b = {key(v):v for v in groups[(s["case"],s["condition"],s["b"])]}
    pair = [(v,b[key(v)]) for v in a if key(v) in b]
    assert len(pair) == int(s["n"])
    for metric, name in [("joint","joint"),("coverage","covered"),("detection","detected")]:
        delta = [int(y[name]=="TRUE")-int(x[name]=="TRUE") for x,y in pair]
        close(st.mean(delta),s[metric+"_change"])
        close(st.stdev(delta)/math.sqrt(len(delta)),s[metric+"_mcse"])

# Quantile interpolation independently follows the documented normal-score rule.
def quantile(draws, p):
    ordered = sorted(draws)
    rank = (len(draws)+1)*p
    k = math.floor(rank)
    if k == rank:
        return ordered[k-1]
    weight = (normal.inv_cdf(p)-normal.inv_cdf(k/(len(draws)+1))) / (
        normal.inv_cdf((k+1)/(len(draws)+1))-normal.inv_cdf(k/(len(draws)+1)))
    return ordered[k-1] + weight*(ordered[k]-ordered[k-1])


program = '''
p <- "validation/l2_default_bootstrap_confirm"
rows <- list()
for(f in list.files(p,"^job-.*[.]rds$",full.names=TRUE)) {
 r<-readRDS(f)
 for(id in names(r$bootstrap_draws)) rows[[length(rows)+1L]] <- data.frame(
 case=r$trials$case[1],condition=r$trials$condition[1],replicate=as.integer(id),
 value=100*r$bootstrap_draws[[id]])
}
write.csv(do.call(rbind,rows),stdout(),row.names=FALSE)
'''
with tempfile.TemporaryFile(mode="w+t") as decoded:
    subprocess.run(["Rscript","--vanilla","-e",program],cwd=root,stdout=decoded,check=True)
    decoded.seek(0)
    draws = defaultdict(list)
    for v in csv.DictReader(decoded):
        draws[key(v)].append(float(v["value"]))
assert len(draws) == 3072 and sum(map(len,draws.values())) == 361728
maximum_error = 0.0
for k, b in draws.items():
    theta = float(index[(*k,"joint_delta")]["estimate"])
    for B in ([99,399] if len(b)==399 else [99]):
        x = b[:B]
        qlo, qhi = quantile(x,.025), quantile(x,.975)
        se = st.stdev(x)
        endpoints = {f"wald_{B}":(theta-z*se,theta+z*se),
                     f"percentile_{B}":(qlo,qhi),f"basic_{B}":(2*theta-qhi,2*theta-qlo)}
        for method, (lo, hi) in endpoints.items():
            actual = index[(*k,method)]
            maximum_error = max(maximum_error, abs(lo-float(actual["lower"])),abs(hi-float(actual["upper"])))
            close(lo,actual["lower"]); close(hi,actual["upper"])
for d in read("error_diagnostics.csv"):
    data = groups[(d["case"],d["condition"],d["method"])]
    errors = [float(v["estimate"])-float(v["truth"]) for v in data]
    ses = [float(v["se"]) for v in data]
    standardized = [e/s for e,s in zip(errors,ses)]
    close(st.mean(errors),d["bias"])
    close(st.stdev(errors),d["empirical_sd"])
    close(math.sqrt(st.mean(s*s for s in ses)),d["rms_se"])
    close(st.stdev(standardized),d["z_sd"])
    assert sum(v < -z for v in standardized) == int(d["z_left"])
    assert sum(v > z for v in standardized) == int(d["z_right"])
report = dict(interval_rows=len(rows),outer_datasets=len(draws),saved_draws=sum(map(len,draws.values())),
              max_endpoint_difference_pp=maximum_error,endpoint_tolerance_pp=1e-10,
              summary_rows=len(groups),paired_rows=len(read("paired.csv")),new_fits=0,
              checked="draw-based endpoints, counts, Wilson bounds, paired MCSEs, error diagnostics")
(out/"verification.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2))
