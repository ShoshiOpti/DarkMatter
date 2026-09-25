"""Development audit of every native C# bootstrap value; not used by the application."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
for ancestor in HERE.parents:
    packages = ancestor / "tmp/python_packages"
    if packages.is_dir():
        sys.path.insert(0, str(packages))
        break
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--actual", type=Path, default=PROJECT / "output/pointwise_comparison/results/bootstrap_summaries.npz")
parser.add_argument("--reference", type=Path, default=PROJECT.parent / "Holonomy_Cubic_TwoField_Investigation_2026-09-24/numerics/pointwise_comparison_2026_09_25/results/bootstrap_summaries.npz")
parser.add_argument("--output", type=Path, default=HERE / "pointwise_bootstrap_array_check.json")
args = parser.parse_args()
atol, rtol = 1e-8, 1e-10
expected_shapes = {"means": (99999, 24), "tstats": (99999, 24), "pooled_means": (99999, 4)}
checks = {}
passed = True
with np.load(args.reference, allow_pickle=False) as reference, np.load(args.actual, allow_pickle=False) as actual:
    keys_match = set(reference.files) == set(actual.files) == set(expected_shapes)
    passed &= keys_match
    for key, shape in expected_shapes.items():
        if key not in reference.files or key not in actual.files:
            checks[key] = {"passed": False, "reason": "missing array"}
            passed = False
            continue
        left, right = actual[key], reference[key]
        shapes_match = left.shape == right.shape == shape
        finite = bool(np.isfinite(left).all() and np.isfinite(right).all())
        record = {"native_shape": list(left.shape), "reference_shape": list(right.shape),
                  "expected_shape": list(shape), "native_dtype": str(left.dtype), "reference_dtype": str(right.dtype),
                  "shapes_match": shapes_match, "all_values_finite": finite}
        if shapes_match and finite:
            error = np.abs(left - right)
            limits = atol + rtol * np.abs(right)
            mismatches = int(np.count_nonzero(error > limits))
            index = np.unravel_index(np.argmax(error), error.shape)
            record.update(compared_values=int(left.size), mismatched_values=mismatches,
                          max_absolute_error=float(error[index]),
                          max_fraction_of_tolerance=float(np.max(error / limits)),
                          max_absolute_error_index=[int(i) for i in index],
                          native_at_max_error=float(left[index]), reference_at_max_error=float(right[index]),
                          passed=mismatches == 0)
        else:
            record["passed"] = False
        checks[key] = record
        passed &= record["passed"]

report = {"status": "pass" if passed else "fail", "all_expected_array_keys_present": keys_match,
          "absolute_tolerance": atol, "relative_tolerance": rtol,
          "tolerance_source": "Original independent_check/verify_pointwise.py close() defaults",
          "compared_values": sum(c.get("compared_values", 0) for c in checks.values()),
          "arrays": checks, "numpy_version": np.__version__,
          "native_path": str(args.actual.resolve()), "reference_path": str(args.reference.resolve()),
          "native_sha256": hashlib.sha256(args.actual.read_bytes()).hexdigest(),
          "reference_sha256": hashlib.sha256(args.reference.read_bytes()).hexdigest(),
          "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          "scope": "All values generated natively by C# compared to original stored NumPy bootstrap arrays; Python used only for this independent development audit."}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
raise SystemExit(0 if passed else 1)
