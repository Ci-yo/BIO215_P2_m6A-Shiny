from pathlib import Path
import csv
import re


ROOT = Path(__file__).resolve().parents[1]
REQUIRED = [
    "README.md",
    "app.R",
    "m6A_imbalanced.csv",
    "m6A_input_example.csv",
    "rf_fit.rds",
    "assets/app-preview.png",
    "www/Mech_figure2.png",
]

missing = [name for name in REQUIRED if not (ROOT / name).is_file()]
assert not missing, f"Missing runtime files: {missing}"

with (ROOT / "m6A_imbalanced.csv").open(encoding="utf-8-sig", newline="") as handle:
    reader = csv.reader(handle)
    header = next(reader)
    row_count = sum(1 for _ in reader)

expected = {
    "gc_content",
    "RNA_type",
    "RNA_region",
    "exon_length",
    "distance_to_junction",
    "evolutionary_conservation",
    "DNA_5mer",
}
assert expected.issubset(header), f"Dataset schema is missing: {sorted(expected - set(header))}"
assert row_count == 48055, f"Expected 48,055 sites, found {row_count:,}"
assert (ROOT / "rf_fit.rds").stat().st_size > 1_000_000, "Model file is unexpectedly small"

app_text = (ROOT / "app.R").read_text(encoding="utf-8")
for filename in REQUIRED[2:5]:
    assert filename in app_text, f"app.R does not reference {filename}"

secret_pattern = re.compile(r"(?i)(api[_-]?key|token|password)\s*[:=]\s*['\"][^'\"]+['\"]")
assert not secret_pattern.search(app_text), "Possible hard-coded secret in app.R"
for forbidden in (".Renviron", ".Rhistory", ".RData", "rsconnect"):
    assert not (ROOT / forbidden).exists(), f"Remove local-only item: {forbidden}"

print(f"OK: Shiny runtime complete; {row_count:,} sites and model asset validated")
