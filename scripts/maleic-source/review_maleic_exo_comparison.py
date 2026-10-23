"""Offline exo stationary-point and same-program RRHO comparison; no jobs."""
import json
from ase import Atoms
from audit_maleic_gaussian_ts import ROOT, audit
from maleic_approach import approach
from maleic_orca_backend import sha, save_json
from postprocess_maleic_results import record, KCAL

DEST=ROOT/'analysis/maleic-endo-exo-20261022'
RAW=ROOT/'received/maleic-exo-ts-4640-20261022T0350'

def main():
    archive=ROOT/'received/maleic-exo-ts-4640-20261022T0350.tar.gz'
    assert sha(archive)=='28288efe70e20f63fc3ff9a54a91a836d47e400b858b909c26e1d01a35b820ed'
    local=json.loads((DEST/'exo-orca-independent-audit.json').read_text())
    result=audit(RAW,jobid=4640,input_root=ROOT/'Gau3-maleic-exo-ts-20261022',
        checkpoint_sha='d37b2352bd65880609e49c3fe15a6cb754f520e037eb675d3c802aa03c2e4d65',
        fchk_sha='5fd7d0ba8f1cb8ae07508d5a477af565b5cfb10d9f07c3fb1e5816abe53cf946',
        local_root=ROOT/'Gau3-local-orca-maleic-exo-20261022',local_record=local)
    result['approach']=approach(Atoms(numbers=result['atomic_numbers'],positions=result['positions_A']))
    assert result['approach']['approach']=='exo'
    if not (DEST/'exo-gaussian-independent-audit.json').exists():
        save_json(DEST/'exo-gaussian-independent-audit.json',result)
    records={
      'Gaussian':record('Gaussian','TS',str((RAW/'jobs/ts-refine').relative_to(ROOT)),result),
      'ORCA':record('ORCA','TS','Gau3-local-orca-maleic-exo-20261022/ts-frequency',local)}
    if (DEST/'exo-thermochemistry.json').exists():
        assert json.loads((DEST/'exo-thermochemistry.json').read_text())==records
    else:
        save_json(DEST/'exo-thermochemistry.json',records)
    prior=json.loads((ROOT/'analysis/maleic-postprocessing-20261021/data.json').read_text())
    comparison={}
    for program, exo in records.items():
        species=prior['species'][program]
        values={}
        for field, label in [('E_hartree','E'),('G_symmetry_corrected_hartree','G')]:
            reference=sum(species[name][field] for name in ('diene','maleic'))
            values[label+'_exo_barrier_kcal_mol']=(exo[field]-reference)*KCAL
            values[label+'_endo_barrier_kcal_mol']=(species['TS'][field]-reference)*KCAL
            values[label+'_exo_minus_endo_kcal_mol']=(exo[field]-species['TS'][field])*KCAL
        comparison[program]=values
    save_json(DEST/'endo-exo-energy-comparison.json',dict(
        temperature_K=298.15,pressure_atm=1,reference='same-program separated reactants',
        free_energy='RRHO with independently audited rotational symmetry',programs=comparison))
    print(json.dumps(comparison,indent=2))
    print(json.dumps(dict(Gaussian=result['comparison_to_ORCA'],thermochemistry={
        k:{'G':v['G_symmetry_corrected_hartree'],'sigma':v['symmetry']['sigma']} for k,v in records.items()}),indent=2))

if __name__=='__main__':main()
