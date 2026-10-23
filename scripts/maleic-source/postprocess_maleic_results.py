"""Source-backed thermochemistry, complete Gaussian IRC and stopped ORCA audit.

Offline only: never launches an electronic-structure calculation or alters runs.
"""
import json
import math
import re
import subprocess
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.data import atomic_numbers
from ase.io import read
from ase.units import Bohr
from audit_stage1 import fchk_fields
from audit_thermochemistry import automorphisms, KCAL, R_J, KB, H, C, AMU_KG
from maleic_common import graph, section
from maleic_orca_backend import sha
from resume_maleic_irc import rows, trajectory, join_forward

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT/'analysis/maleic-postprocessing-20261021'
T = 298.15


def load(path):
    return json.loads((ROOT/path).read_text())


def one(pattern, text):
    values = re.findall(pattern, text, re.M)
    assert len(values) == 1, (pattern, values)
    return float(values[0])


def symmetry(x, numbers, masses):
    a = np.zeros((len(x), len(x)), dtype=bool)
    for i, j in graph(Atoms(numbers=numbers, positions=x)):
        a[i, j] = a[j, i] = True
    centered = x-np.average(x, axis=0, weights=masses)
    fits = []
    for perm in automorphisms(a, numbers, masses):
        target = centered[perm]
        u, _, vt = np.linalg.svd(centered.T@target)
        rotation = u@np.diag([1, 1, np.linalg.det(u@vt)])@vt
        assert np.linalg.det(rotation) > .999999
        fits.append(dict(permutation=perm,
                         max_error_A=float(np.linalg.norm(centered@rotation-target, axis=1).max())))
    counts = {str(tol): sum(p['max_error_A'] < tol for p in fits)
              for tol in [1e-4, 3e-4, 1e-3, 3e-3]}
    # Report tolerance sensitivity explicitly; the 0.001 A group must be closed.
    accepted = [p for p in fits if p['max_error_A'] < 1e-3]
    perms = {tuple(p['permutation']) for p in accepted}
    assert perms and all(tuple(p[q[i]] for i in range(len(x))) in perms for p in perms for q in perms)
    assert counts['0.001'] == counts['0.003'], 'Unstable symmetry assignment'
    return dict(sigma=len(accepted), tolerance_A=.001, tolerance_sweep=counts,
                proper_operations=accepted, graph_automorphism_count=len(fits))


def entropy(x, masses, frequencies, sigma):
    r = (x-np.average(x, axis=0, weights=masses))*1e-10
    inertia = sum(m*AMU_KG*(np.dot(v, v)*np.eye(3)-np.outer(v, v)) for m, v in zip(masses, r))
    b = H/(8*math.pi**2*np.linalg.eigvalsh(inertia))
    sr = R_J*(math.log(math.sqrt(math.pi)/sigma*(KB*T/H)**1.5/np.sqrt(np.prod(b)))+1.5)
    st = R_J*(math.log((2*math.pi*sum(masses)*AMU_KG*KB*T/H**2)**1.5*KB*T/101325.)+2.5)
    positive = np.array([v for v in frequencies if v > 0])
    z = H*C*100*positive/(KB*T)
    modes = R_J*(z/np.expm1(z)-np.log(-np.expm1(-z)))
    return dict(rotational_J_mol_K=sr, translational_J_mol_K=st,
                vibrational_J_mol_K=float(sum(modes)), total_J_mol_K=float(sr+st+sum(modes)),
                low_modes=[dict(frequency_cm_1=float(v), entropy_J_mol_K=float(s))
                           for v, s in zip(positive, modes) if v < 100],
                rotational_constants_Hz=b.tolist())


def record(program, name, directory, accepted):
    d = ROOT/directory
    out = (d/'calc.out').read_text()
    assert accepted['passed'] and abs(accepted['gradient_max_au']) < 1e-5
    if program == 'Gaussian':
        f = fchk_fields(d/'calc.fchk')
        x = np.array(f['Current cartesian coordinates']).reshape(-1, 3)*Bohr
        masses = np.array(f['Real atomic weights']); numbers = f['Atomic numbers']
        frequencies = np.array(f['Vib-E2'][:f['Number of Normal Modes']])
        assert 'Normal termination' in out[-2000:]
        sigma0 = int(one(r'Rotational symmetry number\s+(\d+)\.', out))
        temp = one(r'Temperature\s+([\d.]+) Kelvin', out)
        pressure = one(r'Pressure\s+([\d.]+) Atm', out)
        hraw = one(r'Sum of electronic and thermal Enthalpies=\s*([-\d.]+)', out)
        graw = one(r'Sum of electronic and thermal Free Energies=\s*([-\d.]+)', out)
        sraw = one(r'^\s*Total\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s*$', out)*4.184
        provenance = dict(output_sha256=sha(d/'calc.out'), hessian_sha256=sha(d/'calc.fchk'))
    else:
        hs = (d/'calc.hess').read_text()
        lines = [line.split() for line in section(hs, 'atoms')[1:]]
        x = np.array([list(map(float, line[2:])) for line in lines])*Bohr
        masses = np.array([float(line[1]) for line in lines]); numbers = [atomic_numbers[line[0]] for line in lines]
        frequencies = np.array([float(line.split()[1]) for line in section(hs, 'vibrational_frequencies')[1:]])[6:]
        assert 'ORCA TERMINATED NORMALLY' in out[-3000:] and 'Quasi RRHO          ...    False' in out
        sigma0 = int(one(r'Symmetry Number:\s+(\d+)', out))
        temp = one(r'^Temperature\s+\.\.\.\s+([\d.]+) K', out)
        pressure = one(r'^Pressure\s+\.\.\.\s+([\d.]+) atm', out)
        hraw = one(r'Total Enthalpy\s+\.\.\.\s+([-\d.]+) Eh', out)
        graw = one(r'Final Gibbs free energy\s+\.\.\.\s+([-\d.]+) Eh', out)
        sraw = one(r'Final entropy term\s+\.\.\.\s+([-\d.]+) Eh', out)*KCAL*4184/T
        # Explicitly check that all real frequencies entered the RRHO calculation.
        used = np.array(list(map(float, re.findall(r'^freq\.\s+([\d.]+)\s+E\(vib\)', out, re.M))))
        positive = frequencies[frequencies > 0]
        assert len(used) == len(positive) and max(abs(used-positive)) < .006
        provenance = dict(output_sha256=sha(d/'calc.out'), hessian_sha256=sha(d/'calc.hess'))
    assert (temp, pressure) == (T, 1.)
    assert len(frequencies) == 3*len(x)-6
    assert sum(frequencies < 0) == (1 if name == 'TS' else 0)
    assert max(abs(frequencies-np.array(accepted['frequency_cm_1']))) < .001
    sym = symmetry(x, numbers, masses)
    ent = entropy(x, masses, frequencies, sigma0)
    assert abs(ent['total_J_mol_K']-sraw) < .015, (program, name, ent['total_J_mol_K'], sraw)
    reconstruction_error = (hraw-T*ent['total_J_mol_K']/(KCAL*4184)-graw)*KCAL
    assert abs(reconstruction_error) < .002
    correction = R_J*T*math.log(sym['sigma']/sigma0)/4184
    return dict(program=program, species=name, source_directory=directory, **provenance,
                E_hartree=accepted['energy_hartree'], H_raw_hartree=hraw, G_raw_hartree=graw,
                G_symmetry_corrected_hartree=graw+correction/KCAL,
                printed_sigma=sigma0, symmetry=sym, symmetry_correction_kcal_mol=correction,
                entropy_reconstruction=ent, printed_entropy_J_mol_K=sraw,
                G_reconstruction_error_kcal_mol=reconstruction_error,
                all_positive_frequencies_in_RRHO=True, frequencies_cm_1=frequencies.tolist(),
                positions_A=x.tolist(), atomic_numbers=list(map(int, numbers)))


def terminal_audit():
    root = ROOT/'Gau3-local-orca-maleic-irc-resume1-20261021'
    old = ROOT/'Gau3-local-orca-maleic-irc-20261021'
    for folder in [root, old]:
        for rel, digest in json.loads((folder/'manifest.json').read_text()).items():
            assert sha(folder/rel) == digest, ('Immutable file changed', folder, rel)
    launch = json.loads((root/'launch.json').read_text())
    process = json.loads((root/'irc-forward-down/process.json').read_text())
    ps = subprocess.run(['ps', '-axo', 'pid=,stat=,command='], capture_output=True, text=True, check=True).stdout
    for pid in [launch['pid'], process['pid']]:
        assert not re.search(r'^\s*'+str(pid)+r'\s', ps, re.M), 'Expected stopped process is alive'
    assert not any('software/orca-6.1.1/orca' in line for line in ps.splitlines())
    out = (root/'irc-forward-down/calc.out').read_text()
    failure = (root/'failure.txt').read_text()
    assert 'subprocess.TimeoutExpired:' in failure and '7200.0 seconds' in failure
    assert 'ORCA TERMINATED NORMALLY' not in out and 'THE IRC HAS CONVERGED' not in out
    interrupted = json.loads((root/'irc-forward-down/interrupted.json').read_text())
    assert interrupted['returncode'] == -15 and process['timeout_seconds'] == 7200
    rec = rows(out)
    a, ae = trajectory(root/'source-records/forward-before_trj.xyz')
    b, be = trajectory(root/'irc-forward-down/calc_IRC_D_trj.xyz')
    assert len(b) == len(rec) and len(a) == 35
    highprec = read(root/'irc-forward-down/calc_IRC_D.xyz')
    assert max(abs(b[-1].positions-highprec.positions).ravel()) < 6e-7
    assert all(abs(r['energy']-e) < 5.1e-7 for r, e in zip(rec, be))
    joined = join_forward(a, ae, b, be, read(root/'source-records/resume-point.xyz'))
    assert not any((root/name).exists() for name in ['irc-backward', 'irc-forward-endpoint', 'irc-backward-endpoint', 'result.ready'])
    return dict(state='stopped_at_second_7200s_call_limit_not_converged',
                source_directory=str(root.relative_to(ROOT)), unchanged_manifests=True,
                ended_unix=interrupted['ended_unix'], driver_and_ORCA_exited=True,
                normal_termination=False, IRC_converged=False,
                new_segment_steps=len(b)-1, cumulative_steps=joined['total_steps'],
                last_printed_gradient_au=rec[-1], final_forming_distances_A=joined['forming_distances_A'][-1],
                final_energy_hartree=be[-1], join_energy_difference_hartree=joined['join_energy_difference_hartree'],
                no_backward_endpoints_or_NEB_started=True,
                policy='No further quantum execution without new explicit authorization.', combined_partial_path=joined)


def analyze():
    initial = load('analysis/maleic-stage1-review-20261021/audit.json')
    ts = load('analysis/maleic-ts-review-20261021/audit.json')
    ep = load('analysis/maleic-endpoints-review-20261021/audit.json')['Gaussian']
    paths = load('analysis/maleic-irc-review-20261021/audit.json')['gaussian']
    records = {'Gaussian': {}, 'ORCA': {}}
    names = {'butadiene-opt':'diene', 'maleic-anhydride-opt':'maleic', 'cis-adduct-opt':'independent_product'}
    for job, name in names.items():
        records['Gaussian'][name] = record('Gaussian', name, 'received/maleic-stage1-20261021T0531/raw/jobs/'+job, initial['gaussian'][job])
    records['Gaussian']['TS'] = record('Gaussian', 'TS', ts['source_directory']+'/jobs/ts-refine', ts)
    for job, name in [('reactant-complex-opt','RC'), ('cis-adduct-endpoint-opt','product')]:
        records['Gaussian'][name] = record('Gaussian', name, ep[job]['source_directory'], ep[job])
    for job, name in [('butadiene','diene'), ('maleic-anhydride','maleic'), ('cis-adduct','product'), ('ts','TS')]:
        r = initial['orca'][job]
        records['ORCA'][name] = record('ORCA', name, r['source_directory'], r)
    relative = {}
    for prog, species in records.items():
        keys = ['E_hartree', 'G_raw_hartree', 'G_symmetry_corrected_hartree']
        base = {key: sum(species[n][key] for n in ['diene', 'maleic']) for key in keys}
        relative[prog] = {name:{key.replace('_hartree', '_kcal_mol'):(r[key]-base[key])*KCAL for key in keys}
                          for name, r in species.items() if name not in ['diene', 'maleic']}
    eref = sum(records['Gaussian'][n]['E_hartree'] for n in ['diene', 'maleic'])
    curve = []
    for direction, sign in [('forward', -1), ('reverse', 1)]:
        assert paths[direction]['passed'] and paths[direction]['all_correctors_converged']
        for point in paths[direction]['points']:
            curve.append(dict(point, coordinate_sqrt_amu_bohr=sign*point['coordinate_magnitude_sqrt_amu_bohr'],
                              relative_energy_kcal_mol=(point['energy_hartree']-eref)*KCAL, direction=direction))
    curve.append(dict(coordinate_sqrt_amu_bohr=0., relative_energy_kcal_mol=(ts['energy_hartree']-eref)*KCAL,
                      coordinates_A=ts['positions_A'], forming_distances_A=ts['forming_distances_A'], direction='TS', point=0))
    curve.sort(key=lambda p:p['coordinate_sqrt_amu_bohr'])
    assert len(curve) == 53 and max(curve, key=lambda p:p['relative_energy_kcal_mol'])['direction'] == 'TS'
    return dict(temperature_K=T, pressure_atm=1., energy_unit='kcal/mol',
                reference='Fresh separated twisted near-s-cis diene + maleic anhydride, separately within each program.',
                conventions='Unscaled RRHO, each positive vibration included; proper rotational symmetry only. No conformer population model.',
                sources=['https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/thermochemistry.html'],
                species=records, relative_to_separated_reactants=relative, Gaussian_IRC=curve,
                TS_comparison=ts['comparison_to_ORCA'], ORCA_IRC_terminal=terminal_audit(),
                NEB_CI_status='not_computed_for_maleic_reaction')


if __name__ == '__main__':
    result = analyze()
    DEST.mkdir(exist_ok=True)
    with (DEST/'data.json').open('x') as handle:
        json.dump(result, handle, indent=2, ensure_ascii=False); handle.write('\n')
    print(json.dumps(dict(relative=result['relative_to_separated_reactants'],
                         sigma={p:{n:r['symmetry']['sigma'] for n,r in v.items()} for p,v in result['species'].items()},
                         stopped_at=result['ORCA_IRC_terminal']['cumulative_steps']), indent=2))
