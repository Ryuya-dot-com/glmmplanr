"""Compare this saved-response replay to the preceding curvature contract repair."""
import csv
import hashlib
import json
from pathlib import Path

out = Path(__file__).resolve().parent
prior = out.parent / 'curvature_contract_replay'
results = {}
for name, keys, columns in [
    ('comparison.csv', ('case', 'replicate'),
     ('estimate', 'se', 'valid', 'early_valid', 'covered', 'detected', 'width_ok', 'joint', 'response_hash')),
    ('reference_fits.csv', ('case', 'replicate', 'method'),
     ('estimate', 'se', 'valid', 'covered', 'detected', 'width_ok', 'joint', 'response_hash')),
    ('reference_summary.csv', ('case', 'replicate'),
     ('native_valid', 'profile_valid', 'early_valid', 'evaluations', 'failed_evaluations')),
]:
    def rows(path):
        with path.open() as f:
            return list(csv.DictReader(f))
    a, b = rows(prior / name), rows(out / name)
    aa = {tuple(r[k] for k in keys): r for r in a}
    bb = {tuple(r[k] for k in keys): r for r in b}
    assert len(aa) == len(a) == len(b) == len(bb)
    assert aa.keys() == bb.keys()
    changes = {c: sum(aa[k][c] != bb[k][c] for k in aa) for c in columns}
    assert not any(changes.values()), (name, changes)
    results[name] = {'rows': len(a), 'changed': changes}
    if 'estimate' in columns:
        results[name]['max_point_difference_pp'] = max(abs(float(bb[k]['estimate']) - float(aa[k]['estimate'])) for k in aa)
        results[name]['max_se_difference_pp'] = max(abs(float(bb[k]['se']) - float(aa[k]['se'])) for k in aa)
    if name == 'reference_summary.csv':
        results[name]['failed_evaluations'] = sum(int(r['failed_evaluations']) for r in b)
results['source_sha256'] = hashlib.sha256((out / 'source_glmm_laplace.cpp').read_bytes()).hexdigest()
results['comparison'] = 'Current curvature guard versus preceding curvature contract repair; no new responses.'
(out / 'verification.json').write_text(json.dumps(results, indent=2) + '\n')
print(json.dumps(results, indent=2))
