"""Geometry/Hessian gates for the explicitly authorized A3011/B3303 branch.

No electronic calculations are started by this module. Hessian projection follows
the independently checked CI400 audit, generalized to minima and ring torsion.
"""
import hashlib
import itertools
import json
import re
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.units import Hartree, Bohr

EDGES = {(0,1),(1,2),(2,3),(4,5),(0,4),(3,5),(0,6),(0,7),(1,8),
         (2,9),(3,10),(3,11),(4,12),(4,13),(5,14),(5,15)}
RING = [0,1,2,3,5,4]
H_GROUPS = [[6,7],[8],[9],[10,11],[12,13],[14,15]]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fit(x, y):
    xc, yc = x-x.mean(0), y-y.mean(0)
    u, _, vt = np.linalg.svd(xc.T@yc)
    rot = u@np.diag([1.,1.,np.linalg.det(u@vt)])@vt
    aligned = xc@rot+y.mean(0)
    error = np.linalg.norm(aligned-y, axis=1)
    assert np.linalg.det(rot) > .999999
    return float(np.sqrt(np.mean(error**2))), float(error.max()), aligned


def best_match(x, y):
    choices = itertools.product(*(itertools.permutations(g) for g in H_GROUPS))
    candidates = []
    for groups in choices:
        perm = list(range(6))+[i for group in groups for i in group]
        rms, maximum, aligned = fit(x[perm], y)
        candidates.append((rms, maximum, perm, aligned))
    return min(candidates, key=lambda v: v[0])


def geometry(atoms, reference=None):
    assert list(atoms.numbers) == [6]*6+[1]*10 and not atoms.constraints
    x = atoms.positions
    assert np.isfinite(x).all()
    d = atoms.get_all_distances()
    found = {(i,j) for i in range(16) for j in range(i+1,16)
             if d[i,j] < (1.85 if j<6 else 1.25 if i<6 else .8)}
    assert found == EDGES, 'Cyclohexene connectivity changed'
    for i,j in EDGES:
        assert (1.1 if j<6 else .9) < d[i,j] < (1.85 if j<6 else 1.25)
    assert 1.25 < d[1,2] < 1.45, 'C2=C3 double bond changed'
    if reference is not None:
        for c in [0,3,4,5]:
            ns = sorted(j if i==c else i for i,j in EDGES if c in (i,j))
            def volume(z):
                p = z[ns]
                return np.linalg.det(p[:3]-p[3])
            assert volume(x)*volume(reference) > 0, 'Mapped tetrahedral orientation inverted'


def torsions(x):
    angles = []
    for k in range(6):
        p = x[[RING[(k+j)%6] for j in range(4)]]
        b0, b1, b2 = p[0]-p[1], p[2]-p[1], p[3]-p[2]
        b1 /= np.linalg.norm(b1)
        v, w = b0-np.dot(b0,b1)*b1, b2-np.dot(b2,b1)*b1
        angles.append(np.degrees(np.arctan2(np.dot(np.cross(b1,v),w), np.dot(v,w))))
    return np.array(angles)


def mode_metrics(x, mode):
    mode = mode/np.linalg.norm(mode)
    delta = torsions(x+.001*mode)-torsions(x-.001*mode)
    delta = ((delta+180)%360-180)/.002
    stretches = {f'{i+1}-{j+1}': float(np.dot((x[j]-x[i])/np.linalg.norm(x[j]-x[i]), mode[j]-mode[i]))
                 for i,j in sorted(EDGES)}
    geometry(Atoms('C6H10', positions=x+.1*mode), x)
    geometry(Atoms('C6H10', positions=x-.1*mode), x)
    return dict(ring_torsion_derivatives_deg_per_A=delta.tolist(), bond_derivatives=stretches,
                max_bond_derivative=max(map(abs,stretches.values())),
                max_ring_torsion_derivative_deg_per_A=float(abs(delta).max()))


def section(text, key):
    return text.split('$'+key+'\n',1)[1].split('$',1)[0].strip().splitlines()


def matrix(text, key):
    lines = section(text,key)
    dims = list(map(int,lines[0].split()))
    n,m = dims if len(dims)==2 else [dims[0]]*2
    a = np.full((n,m),np.nan)
    i = 1
    while not np.isfinite(a).all():
        cols = list(map(int,lines[i].split())); i += 1
        for _ in range(n):
            vals = lines[i].split(); i += 1
            a[int(vals[0]),cols] = list(map(float,vals[1:]))
    return a


def project(h, x, mass):
    mv = np.repeat(mass,3)
    hw = h/np.sqrt(np.outer(mv,mv))
    centered = x-np.average(x,axis=0,weights=mass)
    rigid = [(np.sqrt(mass)[:,None]*np.broadcast_to(a,x.shape)).ravel() for a in np.eye(3)]
    rigid += [(np.sqrt(mass)[:,None]*np.cross(a,centered)).ravel() for a in np.eye(3)]
    u,s,_ = np.linalg.svd(np.array(rigid).T,full_matrices=True)
    assert sum(s>1e-8)==6
    internal = u[:,6:]
    eig,vec = np.linalg.eigh(internal.T@hw@internal)
    mode = (internal@vec[:,0])/np.sqrt(mv)
    mode /= np.linalg.norm(mode)
    factor = np.sqrt(4.3597447222071e-18/(5.29177210903e-11**2*1.66053906660e-27))/(2*np.pi*299792458*100)
    return eig, np.sign(eig)*np.sqrt(abs(eig))*factor, mode.reshape(16,3)


def frequency_gate(directory, atoms, imaginary_count):
    directory = Path(directory)
    inp = (directory/'calc.inp').read_text()
    assert 'NumFreq' in inp and 'AnFreq' not in inp and 'Restart false' in inp
    assert 'nprocs 8 nprocs_group 1' in inp and 'QuasiRRHO false' in inp
    out = (directory/'calc.out').read_text()
    assert 'ORCA TERMINATED NORMALLY' in out[-3000:]
    assert re.search(r'Number of Electrons\s+NEL\s+\.+\s+46\b',out)
    steps = list(map(int,re.findall(r'Calculating gradient on displaced geometry\s+(\d+)\s+\(of\s+90\)',out)))
    assert sorted(steps)==list(range(1,91))
    assert 'Central differences            ... used' in out
    assert 'Numerical increment            ... 5.000e-03 bohr' in out
    hs = (directory/'calc.hess').read_text()
    al = section(hs,'atoms'); assert int(al[0])==16
    rows = [v.split() for v in al[1:17]]
    assert [v[0] for v in rows]==['C']*6+['H']*10
    mass = np.array([float(v[1]) for v in rows])
    x = np.array([list(map(float,v[2:])) for v in rows])*Bohr
    assert np.max(abs(x-atoms.positions))<1e-7
    h, modes = matrix(hs,'hessian'), matrix(hs,'normal_modes')
    assert h.shape==modes.shape==(48,48) and np.max(abs(h-h.T))<1e-10
    fr = section(hs,'vibrational_frequencies'); assert int(fr[0])==48
    freq = np.array([float(v.split()[1]) for v in fr[1:]])
    assert len(freq)==48 and max(abs(freq[:6]))<1e-8
    vib = freq[6:]
    eig, rebuilt, mode = project(h,x,mass)
    assert sum(vib<0)==sum(eig<0)==imaginary_count, 'Wrong vibrational/negative-curvature count'
    assert max(abs(rebuilt-vib))<.03
    en = [s.strip() for s in (directory/'calc.engrad').read_text().splitlines() if s.strip() and not s.startswith('#')]
    assert int(en[0])==16
    gradient = np.array(list(map(float,en[2:50]))).reshape(16,3)
    ex = np.array([list(map(float,s.split()[1:])) for s in en[50:66]])*Bohr
    assert np.max(abs(ex-x))<4e-8
    result = json.loads((directory/'result.json').read_text())
    assert abs(float(en[1])-result['energy_hartree'])<1e-10
    force_max = float(np.linalg.norm(gradient*Hartree/Bohr,axis=1).max())
    assert force_max<.001 and np.max(abs(gradient))<1e-5, 'Residual gradient is not a stationary point'
    record = dict(accepted=True,frequencies_cm_1=vib.tolist(),negative_curvature_count=int(sum(eig<0)),
                  max_real_force_eV_A=force_max,gradient_max_au=float(np.max(abs(gradient))),
                  energy_hartree=float(en[1]),hessian_sha256=sha(directory/'calc.hess'),
                  input_sha256=sha(directory/'calc.inp'),positions_A=atoms.positions.tolist(),
                  Hessian_coordinate_error_A=float(np.max(abs(x-atoms.positions))),
                  rebuilt_frequency_max_difference_cm_1=float(max(abs(rebuilt-vib))),
                  numerical_displacements=90,rigid_projected_modes=6,vibrational_modes=42)
    if imaginary_count==1:
        printed = modes[:,6].reshape(16,3); printed /= np.linalg.norm(printed)
        assert abs(float(np.sum(printed*mode)))>.9999
        metrics = mode_metrics(x,mode)
        assert metrics['max_ring_torsion_derivative_deg_per_A']>5, 'Mode is not appreciably ring-torsional'
        assert metrics['max_bond_derivative']<.1, 'Mode has substantial covalent stretching'
        record.update(mode_metrics=metrics,imaginary_mode_cartesian=mode.tolist(),
                      mode_assignment='Ring torsion/puckering; connectivity still needs two-sided IRC')
    return record


def irc_converged(text):
    assert 'THE IRC HAS CONVERGED' in text, 'IRC did not report convergence'
    assert 'MAXIMUM NUMBER OF ITERATIONS REACHED - STOPPING IRC RUN' not in text, 'IRC iteration bound reached'
    assert 'ORCA TERMINATED NORMALLY' in text[-3000:]
