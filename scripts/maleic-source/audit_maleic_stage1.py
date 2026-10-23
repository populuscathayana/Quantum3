"""Read-only scientific checks for the new Gaussian and ORCA reaction stages."""
import argparse
import json
from pathlib import Path
import re
import numpy as np
from scipy.linalg import null_space
from ase import Atoms
from ase.io import read
from ase.units import Bohr, Hartree
from audit_stage1 import fchk_fields, output_geometry, dihedral
from gaussian_convergence import optimization_completion
from maleic_common import geometry, matrix, section, FORMING
from maleic_orca_backend import sha, save_json

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT/'Gau3-local-orca-maleic-da-20261021'
INPUT = ROOT/'Gau3-maleic-stage1-20261021'
JOBIDS = {'butadiene-opt':4546,'maleic-anhydride-opt':4547,'cis-adduct-opt':4548,'ts-constrained':4549}


def independent_curvatures(h, x, mass):
    # Independently construct the null space of the six mass-weighted rigid motions.
    centered = x - np.average(x,axis=0,weights=mass)
    rigid = []
    for axis in np.eye(3):
        rigid.append((np.tile(axis,(len(x),1))*np.sqrt(mass)[:,None]).ravel())
        rigid.append((np.cross(centered,axis)*np.sqrt(mass)[:,None]).ravel())
    internal = null_space(np.array(rigid))
    assert internal.shape == (3*len(x),3*len(x)-6)
    rootmass = np.repeat(np.sqrt(mass),3)
    projected = internal.T@(h/rootmass[:,None]/rootmass[None,:])@internal
    vals,vec = np.linalg.eigh(projected)
    factor = np.sqrt(4.3597447222071e-18/(5.29177210903e-11**2*1.66053906660e-27))/(2*np.pi*29979245800)
    mode = ((internal@vec[:,0])/rootmass).reshape(-1,3)
    mode /= np.linalg.norm(mode)
    return np.sign(vals)*np.sqrt(abs(vals))*factor, mode


def orca():
    assert (LOCAL/'stage1.ready').is_file()
    assert json.loads((LOCAL/'status.json').read_text())['phase']=='stage1_complete_pending_path_preparation'
    for rel,digest in json.loads((LOCAL/'manifest.json').read_text()).items():
        assert sha(LOCAL/rel)==digest
    records={}
    for kind in ('butadiene','maleic-anhydride','cis-adduct','ts'):
        opt = LOCAL/('ts-opt' if kind=='ts' else kind+'/opt')
        freq = LOCAL/('ts-frequency' if kind=='ts' else kind+'/frequency')
        atoms=read(opt/'calc.xyz'); n=len(atoms); ndim=3*n
        expected=1 if kind=='ts' else 0
        seed=read(LOCAL/'seeds'/('ts-distance-seed.xyz' if kind=='ts' else kind+'.xyz'))
        g=geometry(atoms,kind,seed.positions if n==19 else None)
        for folder in (opt,freq):
            out=(folder/'calc.out').read_text()
            result=json.loads((folder/'result.json').read_text())
            assert json.loads((folder/'exit.json').read_text())['returncode']==0
            assert 'ORCA TERMINATED NORMALLY' in out[-3000:]
            assert sha(folder/'calc.inp')==result['input_sha256'] and sha(folder/'calc.out')==result['output_sha256']
            assert f'NEL' in out and re.search(r'Number of Electrons\s+NEL\s+\.+\s+'+str(g['electrons'])+r'\b',out)
        assert 'THE OPTIMIZATION HAS CONVERGED' in (opt/'calc.out').read_text()
        inp=(freq/'calc.inp').read_text(); out=(freq/'calc.out').read_text()
        assert 'NumFreq' in inp and 'Restart false' in inp and 'CentralDiff true' in inp
        assert 'Increment 0.005' in inp and 'QuasiRRHO false' in inp and 'nprocs 8 nprocs_group 1' in inp
        coords=inp.split('* xyz 0 1\n')[1].split('*')[0].strip().splitlines()
        assert np.max(abs(np.array([list(map(float,line.split()[1:])) for line in coords])-atoms.positions))<1e-9
        displacements=re.findall(r'Calculating gradient on displaced geometry\s+(\d+)\s+\(of\s+(\d+)\)',out)
        total=int(displacements[0][1]); assert total in (6*n,6*n-6)
        assert sorted(int(a) for a,b in displacements)==list(range(1,total+1)) and all(int(b)==total for a,b in displacements)
        hs=(freq/'calc.hess').read_text(); rows=[line.split() for line in section(hs,'atoms')[1:]]
        mass=np.array([float(r[1]) for r in rows]); x=np.array([list(map(float,r[2:])) for r in rows])*Bohr
        assert [r[0] for r in rows]==list(atoms.get_chemical_symbols()) and np.max(abs(x-atoms.positions))<1e-7
        h=matrix(hs,'hessian'); assert h.shape==(ndim,ndim) and np.max(abs(h-h.T))<1e-10
        frequencies=np.array([float(line.split()[1]) for line in section(hs,'vibrational_frequencies')[1:]])
        assert len(frequencies)==ndim and max(abs(frequencies[:6]))<1e-8
        rebuilt,mode=independent_curvatures(h,x,mass)
        assert sum(rebuilt<0)==sum(frequencies[6:]<0)==expected
        err=float(max(abs(rebuilt-frequencies[6:]))); assert err<.03
        en=[line.strip() for line in (freq/'calc.engrad').read_text().splitlines() if line.strip() and not line.startswith('#')]
        assert int(en[0])==n
        gradient=np.array(list(map(float,en[2:2+ndim]))).reshape(n,3)
        ex=np.array([list(map(float,line.split()[1:])) for line in en[2+ndim:2+ndim+n]])*Bohr
        assert np.max(abs(ex-x))<4e-8
        maxforce=float(np.linalg.norm(gradient*Hartree/Bohr,axis=1).max())
        assert np.max(abs(gradient))<1e-5 and maxforce<.001
        energy=float(en[1]); assert abs(energy-float(re.findall(r'FINAL SINGLE POINT ENERGY\s+([-\d.]+)',out)[-1]))<1e-10
        record=dict(g,passed=True,energy_hartree=energy,frequency_cm_1=frequencies[6:].tolist(),negative_curvatures=expected,
                    rebuilt_frequency_error_cm_1=err,max_real_force_eV_A=maxforce,gradient_max_au=float(abs(gradient).max()),
                    hessian_sha256=sha(freq/'calc.hess'),opt_geometry_sha256=sha(opt/'calc.xyz'),
                    numerical_displacements=total,source_directory=str(freq.relative_to(ROOT)))
        if expected:
            derivative=[float(np.dot((x[j]-x[i])/np.linalg.norm(x[j]-x[i]),mode[j]-mode[i])) for i,j in sorted(FORMING)]
            assert derivative[0]*derivative[1]>0 and min(map(abs,derivative))>.15
            record.update(forming_mode_derivatives=derivative,imaginary_mode_cartesian=mode.tolist())
        records[kind]=record
    return records


def gaussian(raw):
    manifest=json.loads((INPUT/'manifest.json').read_text())
    assert (raw/'manifest.json').read_bytes()==(INPUT/'manifest.json').read_bytes()
    records={}
    for name,jobid in JOBIDS.items():
        spec=manifest['jobs'][name]; folder=raw/'jobs'/name
        assert sha(folder/'calc.gjf')==spec['input_sha256'] and sha(folder/'gaussian.slurm')==spec['slurm_sha256']
        assert (folder/'submitted.jobid').read_text().strip()==str(jobid)==(folder/'.submission-lock/sbatch.stdout').read_text().strip()
        assert f'job_id={jobid}' in (folder/'result.ready').read_text() and 'Read checkpoint file' in (folder/'formchk.log').read_text()
        out=(folder/'calc.out').read_text(); completion=optimization_completion(out)
        assert completion['step']<=spec['optimization_max_steps']
        f=fchk_fields(folder/'calc.fchk'); n=spec['atoms']; ndim=3*n
        assert f['Number of atoms']==n and f['Number of electrons']==spec['electrons']
        assert f['Charge']==0 and f['Multiplicity']==1 and f['Number of alpha electrons']==f['Number of beta electrons']==spec['electrons']//2
        atoms=Atoms(numbers=f['Atomic numbers'],positions=np.array(f['Current cartesian coordinates']).reshape(n,3)*Bohr)
        seed=read(INPUT/'structures'/f"{spec['kind']}.xyz")
        g=geometry(atoms,spec['kind'],seed.positions if n==19 else None)
        nums,xyz=output_geometry(out); assert nums==f['Atomic numbers']
        distance_error=float(np.max(abs(atoms.get_all_distances()-np.linalg.norm(xyz[:,None]-xyz[None,:],axis=2))))
        assert distance_error<3e-6
        energy=f['Total Energy']; assert abs(energy-float(re.findall(r'SCF Done:\s+E\([^)]+\)\s*=\s*([-\d.]+)',out)[-1]))<1e-7
        fr=[float(v) for row in re.findall(r'Frequencies --\s+([^\n]+)',out) for v in row.split()]
        assert len(fr)==spec['expected_frequencies']
        record=dict(g,passed=True,jobid=jobid,energy_hartree=energy,optimization_completion=completion,
                    frequency_cm_1=fr,positions_A=atoms.positions.tolist(),atomic_numbers=f['Atomic numbers'],
                    fchk_sha256=sha(folder/'calc.fchk'),output_fchk_distance_error_A=distance_error)
        if name=='ts-constrained':
            assert all(abs(v-2.30)<1e-6 for v in g['forming_distances_A'])
            record['accepted_for']='unconstrained Gaussian TS search only'
        else:
            assert f['Number of Normal Modes']==ndim-6 and min(fr)>0
            assert len(f['Vib-Modes'])==(ndim-6)*ndim
            assert max(abs(np.array(fr)-f['Vib-E2'][:ndim-6]))<1e-3
            h=np.zeros((ndim,ndim)); h[np.tril_indices(ndim)]=f['Cartesian Force Constants']; h=h+h.T-np.diag(h.diagonal())
            rebuilt,_=independent_curvatures(h,atoms.positions,np.array(f['Real atomic weights']))
            assert min(rebuilt)>0 and max(abs(rebuilt-fr))<.2
            grad=np.array(f['Cartesian Gradient']); assert max(abs(grad))<1e-4
            conditions=re.findall(r'Temperature\s+([\d.]+) Kelvin\.\s+Pressure\s+([\d.]+) Atm\.',out)
            assert tuple(map(float,conditions[-1]))==(298.15,1.)
            record.update(G_raw_hartree=float(re.findall(r'Sum of electronic and thermal Free Energies=\s*([-\d.]+)',out)[-1]),
                          negative_curvatures=0,rebuilt_frequency_error_cm_1=float(max(abs(rebuilt-fr))),gradient_max_au=float(max(abs(grad))),
                          rotational_symmetry_number=int(re.findall(r'Rotational symmetry number\s+(\d+)',out)[-1]))
        if name in ('butadiene-opt','ts-constrained'):
            record['diene_torsion_deg']=dihedral(atoms.positions.copy()); assert abs(record['diene_torsion_deg'])<90
        records[name]=record
    return records


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--raw',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); assert not a.output.exists()
    report=dict(passed=True,gaussian=gaussian(a.raw.resolve()),orca=orca(),
                scope='Independent raw-output/geometry/curvature/mode/gradient checks, not IRC connectivity; thermochemical symmetry correction pending')
    a.output.parent.mkdir(parents=True,exist_ok=True); save_json(a.output,report)
    print(json.dumps({engine:{name:{k:r[k] for k in ('energy_hartree','frequency_cm_1','forming_distances_A') if k in r}
                     for name,r in report[engine].items()} for engine in ('gaussian','orca')},indent=2))
