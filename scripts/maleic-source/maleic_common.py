"""Geometry and independent numerical-Hessian checks for the new reaction."""
from io import StringIO
import json
from pathlib import Path
import re
import numpy as np
from ase import Atoms
from ase.io import read
from ase.units import Hartree, Bohr
from maleic_orca_backend import sha, validate_atoms

# Fresh Open Babel --gen3d output (2026-10-21), not a previous QC geometry.
MALEIC_SMILES = 'O=C1OC(=O)C=C1'
ADDUCT_SMILES = 'O=C1OC(=O)[C@@H]2CC=CC[C@H]12'
MALEIC_RAW = '''9
Open Babel 3.1.1 --gen3d, O=C1OC(=O)C=C1
O 2.25010 -0.45042 -0.04927
C 1.09294 -0.07818 -0.02628
O -0.00810 -0.85570 -0.00144
C -1.10334 -0.07009 0.02076
O -2.26322 -0.43363 0.04448
C -0.66312 1.31637 0.01032
C 0.66298 1.31138 -0.02003
H -1.33563 2.15511 0.02446
H 1.34167 2.14512 -0.03704
'''
ADDUCT_RAW = '''19
Open Babel 3.1.1 --gen3d, O=C1OC(=O)[C@@H]2CC=CC[C@H]12
O 0.55696 2.55262 -0.04449
C 0.55950 1.35257 0.17866
O 1.53047 0.51273 -0.23635
C 1.27372 -0.74933 0.16362
O 2.00717 -1.69893 -0.06629
C -0.00841 -0.82298 0.93781
C -0.95557 -1.83327 0.28707
C -1.73179 -1.25030 -0.85476
C -2.16146 0.01702 -0.81258
C -1.87210 0.88276 0.37575
C -0.48629 0.62212 0.96821
H 0.24788 -1.15111 1.95262
H -1.67473 -2.19483 1.03288
H -0.40536 -2.71346 -0.06521
H -1.97225 -1.88846 -1.70199
H -2.75648 0.42772 -1.62486
H -2.63994 0.68491 1.13449
H -1.98134 1.93490 0.08927
H -0.43921 1.00868 1.99374
'''
MALEIC_MAP = [5,6,3,1,4,0,2,7,8]
# Exchange the chemically equivalent two ends of the generated cis adduct.
# This is an atom permutation, not a spatial reflection: mapped tetrahedral
# orientations then agree with addition from the selected reactant approach.
ADDUCT_MAP = [6,7,8,9,5,10,3,1,4,0,2,12,13,14,15,16,17,11,18]
NUMBERS = [6]*8+[8]*3+[1]*8
FORMING = {(0,4),(3,5)}
DIENE_EDGES = {(0,1),(1,2),(2,3),(0,4),(0,5),(1,6),(2,7),(3,8),(3,9)}
MALEIC_EDGES = {(0,1),(0,2),(1,3),(2,4),(3,5),(2,6),(3,6),(0,7),(1,8)}
BASE_EDGES = {(0,1),(1,2),(2,3),(4,5),(4,6),(5,7),(6,8),(7,9),(6,10),(7,10),
              (0,11),(0,12),(1,13),(2,14),(3,15),(3,16),(4,17),(5,18)}


def diene(twisted=True):
    x = [[-1.41,1.16,0],[-.74,0,0],[.74,0,0],[1.41,1.16,0]]
    for parent, dx, dy in [(0,.5,np.sqrt(3)/2),(0,-1,0),(1,-.5,-np.sqrt(3)/2),
                           (2,.5,-np.sqrt(3)/2),(3,-.5,np.sqrt(3)/2),(3,1,0)]:
        x.append((np.array(x[parent])+[1.09*dx,1.09*dy,0]).tolist())
    a = Atoms('C4H6', positions=x)
    if twisted:
        angle = np.pi/6
        for i in (3,7,8,9):
            xx,yy,zz = a.positions[i]
            a.positions[i] = [xx, yy*np.cos(angle)-zz*np.sin(angle), yy*np.sin(angle)+zz*np.cos(angle)]
    return a


def seeds():
    m = read(StringIO(MALEIC_RAW), format='xyz')[MALEIC_MAP]
    p = read(StringIO(ADDUCT_RAW), format='xyz')[ADDUCT_MAP]
    # Proper rigid orientation: alkene above diene termini, anhydride below its centre.
    center = m.positions[:2].mean(0)
    ex = m.positions[1]-m.positions[0]; ex /= np.linalg.norm(ex)
    ey = center-m.positions[2:4].mean(0)
    ey -= np.dot(ey,ex)*ex; ey /= np.linalg.norm(ey)
    ez = np.cross(ex,ey)
    aligned = (m.positions-center)@np.column_stack([ex,ey,ez])
    # Remove only the negligible MM out-of-plane noise in this planar reactant.
    aligned[:,2] = 0.
    m.positions = aligned
    planar = diene(False)
    # Equal target distances; the two terminal x offsets are equal by construction.
    alkene_half = (aligned[1,0]-aligned[0,0])/2
    aligned[:2,0] = [-alkene_half,alkene_half]
    target = 2.30
    height = np.sqrt(target**2-(1.41-alkene_half)**2)
    shifted = aligned + [0,1.16,height]
    combined = Atoms(numbers=NUMBERS, positions=np.vstack([
        planar.positions[:4], shifted[:7], planar.positions[4:], shifted[7:]]))
    return {'butadiene': diene(), 'maleic-anhydride': m, 'cis-adduct': p,
            'ts-distance-seed': combined}


def graph(atoms):
    d = atoms.get_all_distances()
    assert np.min(d+np.eye(len(atoms))*100) > .7, 'Atomic overlap'
    found = set()
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            pair = frozenset((atoms.numbers[i],atoms.numbers[j]))
            cutoff = 1.25 if 1 in pair and len(pair)>1 else .8 if pair=={1} else 1.85 if pair=={6} else 1.65 if pair=={6,8} else 1.55
            if d[i,j] < cutoff:
                found.add((i,j))
    return found


def volume(x, center):
    neighbors = [0,5,6,17] if center==4 else [3,4,7,18]
    p = x[neighbors]
    return np.linalg.det(p[:3]-p[3])


def geometry(atoms, kind, reference=None):
    validate_atoms(atoms)
    found = graph(atoms)
    if kind=='butadiene':
        assert list(atoms.numbers)==[6]*4+[1]*6 and found==DIENE_EDGES
    elif kind=='maleic-anhydride':
        assert list(atoms.numbers)==[6]*4+[8]*3+[1]*2 and found==MALEIC_EDGES
    else:
        assert list(atoms.numbers)==NUMBERS
        if kind=='cis-adduct':
            assert found==BASE_EDGES|FORMING, f'Incorrect adduct graph: {found^(BASE_EDGES|FORMING)}'
            assert 1.25<atoms.get_distance(1,2)<1.45, 'Product C2=C3 changed'
        else:
            assert kind in ('ts-distance-seed','ts')
            assert found-FORMING == BASE_EDGES, 'Unintended bond changes'
            assert all(1.85<atoms.get_distance(*b)<2.80 for b in FORMING)
        if reference is not None:
            assert all(volume(atoms.positions,c)*volume(reference,c)>0 for c in (4,5)), 'Mapped stereochemistry changed'
    return dict(formula=atoms.get_chemical_formula(), atoms=len(atoms), electrons=int(sum(atoms.numbers)),
                vibrations=3*len(atoms)-6, edges=[list(v) for v in sorted(found)],
                forming_distances_A=[float(atoms.get_distance(*v)) for v in sorted(FORMING)] if len(atoms)==19 else None)


def section(text,key):
    return text.split('$'+key+'\n',1)[1].split('$',1)[0].strip().splitlines()


def matrix(text,key):
    lines=section(text,key); dims=list(map(int,lines[0].split()))
    n,m=dims if len(dims)==2 else [dims[0]]*2
    a=np.full((n,m),np.nan); i=1
    while not np.isfinite(a).all():
        cols=list(map(int,lines[i].split())); i+=1
        for _ in range(n):
            vals=lines[i].split(); i+=1
            a[int(vals[0]),cols]=list(map(float,vals[1:]))
    return a


def project(h,x,mass):
    mv=np.repeat(mass,3); hw=h/np.sqrt(np.outer(mv,mv))
    centered=x-np.average(x,axis=0,weights=mass)
    rigid=[(np.sqrt(mass)[:,None]*np.broadcast_to(v,x.shape)).ravel() for v in np.eye(3)]
    rigid += [(np.sqrt(mass)[:,None]*np.cross(v,centered)).ravel() for v in np.eye(3)]
    u,s,_=np.linalg.svd(np.array(rigid).T,full_matrices=True)
    assert sum(s>1e-8)==6
    internal=u[:,6:]; eig,vec=np.linalg.eigh(internal.T@hw@internal)
    mode=(internal@vec[:,0])/np.sqrt(mv); mode/=np.linalg.norm(mode)
    factor=np.sqrt(4.3597447222071e-18/(5.29177210903e-11**2*1.66053906660e-27))/(2*np.pi*299792458*100)
    return eig,np.sign(eig)*np.sqrt(abs(eig))*factor,mode.reshape(len(x),3)


def frequency_gate(directory,atoms,imaginary_count):
    directory=Path(directory); n=len(atoms); ndof=3*n
    inp=(directory/'calc.inp').read_text(); out=(directory/'calc.out').read_text()
    assert 'NumFreq' in inp and 'AnFreq' not in inp and 'Restart false' in inp
    assert 'nprocs 8 nprocs_group 1' in inp and 'QuasiRRHO false' in inp
    assert 'ORCA TERMINATED NORMALLY' in out[-3000:]
    assert re.search(r'Number of Electrons\s+NEL\s+\.+\s+'+str(sum(atoms.numbers))+r'\b',out)
    displacements=re.findall(r'Calculating gradient on displaced geometry\s+(\d+)\s+\(of\s+(\d+)\)',out)
    totals={int(t) for _,t in displacements}
    assert len(totals)==1 and totals.issubset({6*n,6*n-6})
    total=totals.pop()
    assert sorted(int(k) for k,_ in displacements)==list(range(1,total+1)), 'Incomplete central differences'
    assert 'Central differences            ... used' in out
    assert 'Numerical increment            ... 5.000e-03 bohr' in out
    hs=(directory/'calc.hess').read_text(); al=section(hs,'atoms'); assert int(al[0])==n
    rows=[v.split() for v in al[1:]]; assert len(rows)==n
    assert [v[0] for v in rows]==list(atoms.get_chemical_symbols())
    mass=np.array([float(v[1]) for v in rows]); x=np.array([list(map(float,v[2:])) for v in rows])*Bohr
    assert np.max(abs(x-atoms.positions))<1e-7
    h,modes=matrix(hs,'hessian'),matrix(hs,'normal_modes')
    assert h.shape==modes.shape==(ndof,ndof) and np.max(abs(h-h.T))<1e-10
    fr=section(hs,'vibrational_frequencies'); assert int(fr[0])==ndof
    freq=np.array([float(v.split()[1]) for v in fr[1:]])
    assert len(freq)==ndof and max(abs(freq[:6]))<1e-8
    vib=freq[6:]; eig,rebuilt,mode=project(h,x,mass)
    assert sum(vib<0)==sum(eig<0)==imaginary_count, 'Wrong negative-curvature count'
    assert max(abs(rebuilt-vib))<.03
    en=[v.strip() for v in (directory/'calc.engrad').read_text().splitlines() if v.strip() and not v.startswith('#')]
    assert int(en[0])==n
    gradient=np.array(list(map(float,en[2:2+ndof]))).reshape(n,3)
    ex=np.array([list(map(float,v.split()[1:])) for v in en[2+ndof:2+ndof+n]])*Bohr
    assert np.max(abs(ex-x))<4e-8
    result=json.loads((directory/'result.json').read_text())
    assert abs(float(en[1])-result['energy_hartree'])<1e-10
    force=float(np.linalg.norm(gradient*Hartree/Bohr,axis=1).max())
    assert force<.001 and np.max(abs(gradient))<1e-5, 'Stationary-point gradient not accepted'
    record=dict(accepted=True, frequencies_cm_1=vib.tolist(), negative_curvature_count=int(sum(eig<0)),
                atoms=n,electrons=int(sum(atoms.numbers)),vibrational_modes=ndof-6,rigid_projected_modes=6,
                numerical_displacements=total,max_real_force_eV_A=force,
                gradient_max_au=float(np.max(abs(gradient))),energy_hartree=float(en[1]),
                hessian_sha256=sha(directory/'calc.hess'),input_sha256=sha(directory/'calc.inp'),
                positions_A=atoms.positions.tolist(),hessian_coordinate_error_A=float(np.max(abs(x-atoms.positions))),
                rebuilt_frequency_max_difference_cm_1=float(max(abs(rebuilt-vib))))
    if imaginary_count:
        assert n==19 and imaginary_count==1
        printed=modes[:,6].reshape(n,3); printed/=np.linalg.norm(printed)
        assert abs(float(np.sum(printed*mode)))>.9999
        derivatives=[float(np.dot((x[j]-x[i])/np.linalg.norm(x[j]-x[i]),mode[j]-mode[i])) for i,j in sorted(FORMING)]
        assert derivatives[0]*derivatives[1]>0 and min(map(abs,derivatives))>.15, 'Not a two-bond addition mode'
        record.update(imaginary_mode_cartesian=mode.tolist(),forming_bond_mode_derivatives=derivatives,
                      mode_assignment='Cooperative formation/breaking of C1-C5 and C4-C6; IRC pending')
    return record
