"""Independent exo raw-file audit adapted from the validated stage1 review.
Read-only quantum sources; no path calculations.
"""
from audit_maleic_stage1 import *
from maleic_approach import approach
LOCAL = ROOT/'Gau3-local-orca-maleic-exo-20261022'

def orca():
    assert (LOCAL/'stage1.ready').is_file()
    assert json.loads((LOCAL/'status.json').read_text())['phase']=='exo_TS_complete_pending_review'
    for rel,digest in json.loads((LOCAL/'manifest.json').read_text()).items():
        assert sha(LOCAL/rel)==digest
    records={}
    for kind in ('ts',):
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


if __name__ == '__main__':
    ts=orca()['ts']
    descriptor=approach(read(LOCAL/'ts-opt/calc.xyz'))
    assert descriptor['approach']=='exo'
    base=json.loads((ROOT/'analysis/maleic-stage1-review-20261021/audit.json').read_text())['orca']
    reference=base['butadiene']['energy_hartree']+base['maleic-anhydride']['energy_hartree']
    ts.update(approach=descriptor,electronic_barrier_kcal_mol=(ts['energy_hartree']-reference)*627.509474,
              exo_minus_endo_TS_kcal_mol=(ts['energy_hartree']-base['ts']['energy_hartree'])*627.509474,
              scope='Stationary-point audit only; no IRC/NEB followup under latest user instruction')
    print(json.dumps(ts,indent=2))
