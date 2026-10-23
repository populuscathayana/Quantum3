"""Independent raw Gaussian IRC endpoint geometry/frequency/Hessian review."""
import argparse
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
from maleic_common import NUMBERS, BASE_EDGES, graph, geometry
from maleic_orca_backend import sha, save_json
from run_maleic_irc import path_geometry, endpoint_kind, product_match

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT/'Gau3-maleic-endpoints-20261021'
JOBS = {'reactant-complex-opt':4579, 'cis-adduct-endpoint-opt':4580}


def atoms_from(fields):
    return Atoms(numbers=fields['Atomic numbers'],
                 positions=np.array(fields['Current cartesian coordinates']).reshape(-1,3)*Bohr)


def audit(raw):
    manifest = json.loads((PACKAGE/'manifest.json').read_text())
    assert (raw/'manifest.json').read_bytes()==(PACKAGE/'manifest.json').read_bytes()
    reference = atoms_from(fchk_fields(ROOT/'received/maleic-stage1-20261021T0531/raw/jobs/cis-adduct-opt/calc.fchk'))
    records = {}
    for name,jobid in JOBS.items():
        spec=manifest['jobs'][name]; folder=raw/'jobs'/name
        assert sha(folder/'calc.gjf')==spec['input_sha256']
        assert sha(folder/'gaussian.slurm')==spec['slurm_sha256']
        assert (folder/'submitted.jobid').read_text().strip()==str(jobid)==(folder/'.submission-lock/sbatch.stdout').read_text().strip()
        assert f'job_id={jobid}' in (folder/'result.ready').read_text()
        assert 'Read checkpoint file' in (folder/'formchk.log').read_text()
        source=raw/'source-records'/spec['source_direction']/'calc.fchk'
        assert sha(source)==spec['source_fchk_sha256']
        initial=read(raw/'structures'/f'{name}.xyz')
        assert np.max(abs(initial.positions-atoms_from(fchk_fields(source)).positions))<1e-9
        out=(folder/'calc.out').read_text(); completion=optimization_completion(out)
        # Freq repeats a two-step stationarity check after the main optimization.
        # Gaussian's printed internal limit can be below the requested MaxCycles.
        optimization_steps=re.findall(r'Step number\s+(\d+) out of a maximum of\s+(\d+)',out)
        assert optimization_steps and all(int(s)<=int(limit)<=200 for s,limit in optimization_steps)
        assert 'MaxCycles=200' in (folder/'calc.gjf').read_text()
        f=fchk_fields(folder/'calc.fchk'); atoms=atoms_from(f)
        assert f['Number of atoms']==19 and f['Atomic numbers']==NUMBERS
        assert f['Number of electrons']==80 and f['Number of alpha electrons']==f['Number of beta electrons']==40
        assert f['Charge']==0 and f['Multiplicity']==1
        distances=path_geometry(atoms)
        kind=endpoint_kind(atoms)
        if spec['kind']=='reactant-complex':
            assert kind=='reactant-complex' and graph(atoms)==BASE_EDGES
        else:
            assert kind=='product'
            geometry(atoms,'cis-adduct',initial.positions)
        nums,coords=output_geometry(out); assert nums==NUMBERS
        error=float(np.max(abs(atoms.get_all_distances()-np.linalg.norm(coords[:,None]-coords[None,:],axis=2))))
        assert error<3e-6
        energy=f['Total Energy']
        assert abs(energy-float(re.findall(r'SCF Done:\s+E\([^)]+\)\s*=\s*([-\d.]+)',out)[-1]))<1e-7
        fr=np.array([float(x) for row in re.findall(r'Frequencies --\s+([^\n]+)',out) for x in row.split()])
        assert len(fr)==f['Number of Normal Modes']==51 and min(fr)>0
        assert len(f['Vib-Modes'])==51*57
        assert max(abs(fr-f['Vib-E2'][:51]))<.001
        h=np.zeros((57,57));h[np.tril_indices(57)]=f['Cartesian Force Constants']
        h=h+h.T-np.diag(h.diagonal())
        rebuilt,_=independent_curvatures(h,atoms.positions,np.array(f['Real atomic weights']))
        assert min(rebuilt)>0 and max(abs(rebuilt-fr))<.2
        gradient=np.array(f['Cartesian Gradient']).reshape(19,3)
        assert abs(gradient).max()<1e-4
        conditions=re.findall(r'Temperature\s+([\d.]+) Kelvin\.\s+Pressure\s+([\d.]+) Atm\.',out)
        assert tuple(map(float,conditions[-1]))==(298.15,1.)
        records[name]=dict(passed=True,jobid=jobid,source_jobid=spec['source_jobid'],
            formula=atoms.get_chemical_formula(),atoms=19,electrons=80,vibrations=51,negative_curvatures=0,
            endpoint_kind=kind,forming_distances_A=distances,energy_hartree=energy,
            G_raw_hartree=float(re.findall(r'Sum of electronic and thermal Free Energies=\s*([-\d.]+)',out)[-1]),
            printed_rotational_symmetry_number=int(re.findall(r'Rotational symmetry number\s+(\d+)',out)[-1]),
            frequency_cm_1=fr.tolist(),rebuilt_frequency_error_cm_1=float(max(abs(rebuilt-fr))),
            gradient_max_au=float(abs(gradient).max()),max_real_force_eV_A=float(np.linalg.norm(gradient*Hartree/Bohr,axis=1).max()),
            positions_A=atoms.positions.tolist(),output_fchk_distance_error_A=error,
            optimization_completion=completion,fchk_sha256=sha(folder/'calc.fchk'),
            source_fchk_sha256=sha(source),source_directory=str(folder.relative_to(ROOT)))
        if kind=='product':
            records[name]['comparison_with_4548']=product_match(atoms.positions,reference.positions)
            ref_energy=fchk_fields(ROOT/'received/maleic-stage1-20261021T0531/raw/jobs/cis-adduct-opt/calc.fchk')['Total Energy']
            records[name]['energy_difference_from_4548_kcal_mol']=(energy-ref_energy)*627.5094740631
    return dict(passed=True,Gaussian=records,scope='IRC terminal minima independently confirmed from new Opt/Freq and projected raw Hessians. Rotational symmetry correction and RRHO comparison remain separate postprocessing.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();assert not args.output.exists()
    report=audit(args.raw.resolve());args.output.parent.mkdir(parents=True,exist_ok=True)
    save_json(args.output,report)
    print(json.dumps({name:{k:r[k] for k in ('energy_hartree','forming_distances_A','gradient_max_au','negative_curvatures','comparison_with_4548') if k in r} for name,r in report['Gaussian'].items()},indent=2))
