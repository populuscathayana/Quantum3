"""Offline integrity and scientific arithmetic checks; standard library only."""
from pathlib import Path
import hashlib, json, re
root = Path(__file__).resolve().parents[1]
manifest = json.loads((root/'data/publication-20261023.json').read_text())
for name, record in manifest['files'].items():
    p = root/name
    assert p.stat().st_size == record['bytes'], name
    assert hashlib.sha256(p.read_bytes()).hexdigest() == record['sha256'], name
tex = (root/'experiment-report/cycloaddition_mechanism.tex').read_text()
labels = re.findall(r'\\label\{([^}]+)\}', tex)
assert len(labels) == len(set(labels)), 'Duplicate LaTeX labels'
for image in re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}', tex):
    assert (root/'experiment-report'/image).is_file(), image
html = (root/'experiment-report/maleic-media.html').read_text()
for link in re.findall(r'(?:src|href)="(figures/[^"]+)"', html):
    assert (root/'experiment-report'/link).is_file(), link
energy = json.loads((root/'analysis/maleic-endo-exo-20261022/endo-exo-energy-comparison.json').read_text())
for program, values in energy['programs'].items():
    for symbol in ['E', 'G']:
        endo, exo = (values[f'{symbol}_{a}_barrier_kcal_mol'] for a in ['endo','exo'])
        assert abs(exo-endo-values[f'{symbol}_exo_minus_endo_kcal_mol']) < 1e-9
        print(f'{program} {symbol}: endo={endo:.2f}, exo={exo:.2f}, difference={exo-endo:.2f} kcal/mol')
print(f"PASS: {len(manifest['files'])} file hashes, figure links, unique labels, energy differences")
