#!/usr/bin/env python3
"""Offline scientific checks and saved-data figure regeneration."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--output', type=Path, default=ROOT/'build/figures')
    args = parser.parse_args()
    subprocess.run([sys.executable, str(ROOT/'scripts/verify_results.py')], check=True, cwd=ROOT)
    if args.check_only:
        return
    output = args.output.resolve()
    protected = [ROOT/'sources', ROOT/'data', ROOT/'experiment-report']
    if any(output == p or p in output.parents for p in protected):
        parser.error('Choose a derived-output directory outside sources, data and experiment-report')
    output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, QUANTUM3_OUTPUT=str(output))
    for name in ['plot_energy_ladder.py', 'plot_paths.py', 'plot_gibbs.py', 'plot_irc.py']:
        subprocess.run([sys.executable, str(ROOT/'scripts'/name), '--replace-render'],
                       check=True, cwd=ROOT, env=env)
    print(f'Regenerated figures: {output}')


if __name__ == '__main__':
    main()
