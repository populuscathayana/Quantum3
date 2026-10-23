"""Independent final-output/Hessian audit for the fresh maleic Gaussian TS."""
import argparse
import itertools
import json
from pathlib import Path
import re
import numpy as np
from ase import Atoms
from ase.io import read
from ase.units import Bohr, Hartree
from audit_stage1 import fchk_fields, output_geometry
from audit_maleic_stage1 import independent_curvatures
from gaussian_convergence import optimization_completion
from maleic_common import geometry, FORMING, BASE_EDGES
from maleic_orca_backend import sha, save_json

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT/'received/maleic-ts-4551-20261021T0625'
INPUT = ROOT/'Gau3-maleic-ts-20261021'
REMOTE_CHK_SHA = '095c6bae2a7455775ba9d18748cc2080482b122f92562072de086b6a20a7f499'
REMOTE_FCHK_SHA = 'e3426dd54159b9d50b5e9a061ed8d66e62249d30d839e52f73708811315710c6'


def proper_fit(x, y):
    x0=x-x.mean(axis=0); y0=y-y.mean(axis=0)
    u,_,vt=np.linalg.svd(x0.T@y0)
    r=u@np.diag([1.,1.,np.linalg.det(u@vt)])@vt
    assert np.linalg.det(r)>.999999
    return float(np.sqrt(np.mean(np.sum((x0@r-y0)**2,axis=1)))),r


def audit(raw=RAW, jobid=4551, input_root=INPUT, checkpoint_sha=REMOTE_CHK_SHA,
          fchk_sha=REMOTE_FCHK_SHA, local_root=None, local_record=None):
    folder=raw/'jobs/ts-refine'; out=(folder/'calc.out').read_text()
    manifest=json.loads((input_root/'manifest.json').read_text())
    assert (raw/'manifest.json').read_bytes()==(input_root/'manifest.json').read_bytes()
    for line in (raw/'SHA256SUMS').read_text().splitlines():
        digest,rel=line.split('  ',1); assert sha(raw/rel)==digest
    assert sha(folder/'calc.gjf')==manifest['input_sha256']
    assert sha(folder/'calc.fchk')==fchk_sha
    assert (folder/'submitted.jobid').read_text().strip()==(folder/'.submission-lock/sbatch.stdout').read_text().strip()==str(jobid)
    assert f'job_id={jobid}\n' in (folder/'result.ready').read_text()
    assert 'Read checkpoint file' in (folder/'formchk.log').read_text()
    assert 'Normal termination of Gaussian' in out[-1500:]
    completion=optimization_completion(out); assert completion['step']<=200
    f=fchk_fields(folder/'calc.fchk'); n=19; dim=57; modes=51
    assert [f[k] for k in ('Number of atoms','Number of electrons','Charge','Multiplicity','Number of alpha electrons','Number of beta electrons')]==[19,80,0,1,40,40]
    atoms=Atoms(numbers=f['Atomic numbers'],positions=np.array(f['Current cartesian coordinates']).reshape(n,3)*Bohr)
    seed=read(input_root/'structures/ts-start.xyz'); g=geometry(atoms,'ts',seed.positions)
    nums,ox=output_geometry(out); assert nums==f['Atomic numbers']
    distance_error=float(np.max(abs(atoms.get_all_distances()-np.linalg.norm(ox[:,None]-ox[None,:],axis=2))))
    assert distance_error<3e-6
    freqs=np.array([float(v) for row in re.findall(r'Frequencies --\s+([^\n]+)',out) for v in row.split()])
    assert len(freqs)==f['Number of Normal Modes']==modes and sum(freqs<0)==1
    assert len(f['Vib-Modes'])==dim*modes
    assert max(abs(freqs-f['Vib-E2'][:modes]))<.001
    h=np.zeros((dim,dim)); h[np.tril_indices(dim)]=f['Cartesian Force Constants']; h=h+h.T-np.diag(h.diagonal())
    rebuilt,mode=independent_curvatures(h,atoms.positions,np.array(f['Real atomic weights']))
    assert sum(rebuilt<0)==1 and max(abs(rebuilt-freqs))<.2
    printed=np.array(f['Vib-Modes'][:dim]).reshape(n,3); printed/=np.linalg.norm(printed)
    cosine=float(abs(np.sum(printed*mode))); assert cosine>.9999
    deriv=[float(np.dot((atoms.positions[j]-atoms.positions[i])/atoms.get_distance(i,j),mode[j]-mode[i])) for i,j in sorted(FORMING)]
    assert deriv[0]*deriv[1]>0 and min(map(abs,deriv))>.15
    grad=np.array(f['Cartesian Gradient']).reshape(n,3)
    maxforce=float(np.linalg.norm(grad*Hartree/Bohr,axis=1).max())
    assert abs(grad).max()<1e-5 and maxforce<.001
    energy=f['Total Energy']; assert abs(energy-float(re.findall(r'SCF Done:\s+E\([^)]+\)\s*=\s*([-\d.]+)',out)[-1]))<1e-7
    assert tuple(map(float,re.findall(r'Temperature\s+([\d.]+) Kelvin\.\s+Pressure\s+([\d.]+) Atm\.',out)[-1]))==(298.15,1.)
    prior=json.loads((ROOT/'analysis/maleic-stage1-review-20261021/audit.json').read_text())
    local=read((local_root or ROOT/'Gau3-local-orca-maleic-da-20261021')/'ts-opt/calc.xyz')
    assert list(local.numbers)==list(atoms.numbers)
    # Preserve carbon/oxygen identities; enumerate only H attached to the same carbon.
    groups=[]
    for carbon in range(8):
        hs=sorted(j if i==carbon else i for i,j in BASE_EDGES if carbon in (i,j) and max(i,j)>=11)
        if len(hs)>1: groups.append(hs)
    choices=[]
    for perms in itertools.product(*(list(itertools.permutations(group)) for group in groups)):
        p=np.arange(n)
        for group,perm in zip(groups,perms): p[group]=perm
        rms,r=proper_fit(atoms.positions[p],local.positions); choices.append((rms,p,r))
    rms,p,r=min(choices,key=lambda row:row[0]); assert rms<.02
    local_ts=local_record or prior['orca']['ts']
    localmode=np.array(local_ts['imaginary_mode_cartesian'])
    overlap=float(abs(np.sum((mode[p]@r)*localmode))); assert overlap>.99
    ge=prior['gaussian']; oe=prior['orca']
    barrier=(energy-ge['butadiene-opt']['energy_hartree']-ge['maleic-anhydride-opt']['energy_hartree'])*627.5094740631
    obarrier=(local_ts['energy_hartree']-oe['butadiene']['energy_hartree']-oe['maleic-anhydride']['energy_hartree'])*627.5094740631
    return dict(g,passed=True,jobid=jobid,accepted_for='stationary-point comparison and postprocessing',
        source_directory=str(raw.relative_to(ROOT)),source_checkpoint_sha256=checkpoint_sha,
        fchk_sha256=sha(folder/'calc.fchk'),output_sha256=sha(folder/'calc.out'),
        input_sha256=sha(folder/'calc.gjf'),optimization_completion=completion,
        energy_hartree=energy,frequency_cm_1=freqs.tolist(),negative_curvatures=1,
        positions_A=atoms.positions.tolist(),atomic_numbers=nums,imaginary_mode_cartesian=mode.tolist(),
        mode_hessian_cosine_abs=cosine,forming_mode_derivatives=deriv,
        rebuilt_frequency_error_cm_1=float(max(abs(rebuilt-freqs))),
        max_real_force_eV_A=maxforce,gradient_max_au=float(abs(grad).max()),
        output_fchk_distance_error_A=distance_error,
        G_raw_hartree=float(re.findall(r'Sum of electronic and thermal Free Energies=\s*([-\d.]+)',out)[-1]),
        rotational_symmetry_number_reported=int(re.findall(r'Rotational symmetry number\s+(\d+)',out)[-1]),
        symmetry_correction='pending separate thermochemical audit',
        comparison_to_ORCA=dict(proper_rotation_determinant=float(np.linalg.det(r)),
            permutation_zero_based=p.tolist(),RMSD_A=rms,imaginary_mode_cosine_abs=overlap,
            Gaussian_electronic_barrier_kcal_mol=barrier,ORCA_electronic_barrier_kcal_mol=obarrier,
            reference='Each program uses its own fresh separated reactants; no mixed absolute energies.'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=audit();a.output.parent.mkdir(parents=True,exist_ok=True);save_json(a.output,result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('positions_A','imaginary_mode_cartesian','edges','frequency_cm_1')},indent=2))
