"""Verify the public snapshot against raw scientific files; no quantum calls."""
import hashlib
import json
import math
from pathlib import Path
import re
import unittest
import numpy as np
from ase.units import Hartree, Bohr
from pypdf import PdfReader
from read_fchk import fchk_fields
from conformer_path_common import section, matrix, project

ROOT = Path(__file__).resolve().parents[1]
KCAL = 627.5094740631


def load(name):
    return json.loads((ROOT/'data'/name).read_text())


class ScientificChecks(unittest.TestCase):
    def test_manifest(self):
        for name, record in json.loads((ROOT/'MANIFEST.json').read_text())['files'].items():
            p = ROOT/name
            self.assertEqual(p.stat().st_size, record['bytes'], name)
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), record['sha256'], name)

    def test_gaussian_raw_stationary_points(self):
        for name, row in load('thermochemistry.json')['jobs'].items():
            directory = ROOT/row['source_directory']
            f = fchk_fields(directory/'calc.fchk')
            out = (directory/'calc.out').read_text()
            self.assertIn('Normal termination of Gaussian', out[-3000:])
            self.assertAlmostEqual(f['Total Energy'], row['raw_Gaussian_hartree']['E'], places=10)
            n = 3*f['Number of atoms']-6
            negative = sum(v < 0 for v in f['Vib-E2'][:n])
            self.assertEqual(negative, 1 if name in ['ts-opt', 'cyclohexene-conf-ts'] else 0)
            g = float(re.findall(r'Sum of electronic and thermal Free Energies=\s*([-\d.]+)',out)[-1])
            self.assertEqual(g, row['raw_Gaussian_hartree']['G'])

    def test_symmetry_and_thermochemistry(self):
        d = load('thermochemistry.json')
        expected = {'ethylene-opt':4, 'butadiene-opt':2, 'cyclohexene-opt':2}
        for name, row in d['jobs'].items():
            sigma = row['rotational_symmetry']['sigma']
            self.assertEqual(sigma, expected.get(name, 1))
            correction = 8.31446261815324*298.15/4184*math.log(sigma/row['printed_sigma'])
            self.assertAlmostEqual(correction,row['symmetry_G_correction_kcal_mol'],places=7)
        g = d['relative_to_separated_reactants_kcal_mol']['ts-opt']['G_symmetry_corrected']
        self.assertAlmostEqual(g,26.469316558565676,places=7)

    def test_orca_hessians_and_gradients(self):
        a = load('conformer-audit.json')['stationary_points']
        for label, imaginary in [('A',0),('B',0),('TS',1),('cycloaddition-CI400',1)]:
            dirname = f'conformer-{label}' if label in a else label
            directory = ROOT/'sources/orca'/dirname
            text = (directory/'calc.hess').read_text()
            rows = [v.split() for v in section(text,'atoms')[1:17]]
            self.assertEqual([v[0] for v in rows],['C']*6+['H']*10)
            mass = np.array([float(v[1]) for v in rows])
            xyz = np.array([list(map(float,v[2:])) for v in rows])*Bohr
            eig, rebuilt, _ = project(matrix(text,'hessian'),xyz,mass)
            frequencies = np.array([float(v.split()[1]) for v in section(text,'vibrational_frequencies')[1:]])
            self.assertEqual(len(frequencies),48)
            self.assertEqual(sum(eig<0),imaginary)
            self.assertLess(max(abs(rebuilt-frequencies[6:])),.03)
            out = (directory/'calc.out').read_text()
            self.assertIn('ORCA TERMINATED NORMALLY',out[-3000:])
            steps = re.findall(r'Calculating gradient on displaced geometry\s+(\d+)\s+\(of\s+90\)',out)
            self.assertEqual(sorted(map(int,steps)),list(range(1,91)))
            en = [s.strip() for s in (directory/'calc.engrad').read_text().splitlines()
                  if s.strip() and not s.startswith('#')]
            gradient = np.array(list(map(float,en[2:50]))).reshape(16,3)
            self.assertLess(np.linalg.norm(gradient*Hartree/Bohr,axis=1).max(),.001)
            if label in a:
                self.assertAlmostEqual(float(en[1]),a[label]['energy_hartree'],places=10)
                self.assertLess(abs(xyz-np.array(a[label]['positions_A'])).max(),1e-7)

    def test_ladder_references(self):
        d, t = load('energy-ladder.json'), load('thermochemistry.json')['jobs']
        ref = t['ethylene-opt']['raw_Gaussian_hartree']['E']+t['butadiene-opt']['raw_Gaussian_hartree']['E']
        mapping = {'RC':'reactant-endpoint-opt','TS_DA':'ts-opt','B':'product-endpoint-opt',
                   'TS_conf_candidate':'cyclohexene-conf-ts','A':'cyclohexene-opt'}
        for row in d['Gaussian_global_states'][1:]:
            expected = (t[mapping[row['key']]]['raw_Gaussian_hartree']['E']-ref)*KCAL
            self.assertAlmostEqual(row['E_kcal_mol'],expected,places=8)
        st = load('conformer-audit.json')['stationary_points']
        for side in ['A','B']:
            height = (st['TS']['energy_hartree']-st[side]['energy_hartree'])*KCAL
            self.assertAlmostEqual(height,d['conformer_rows'][2][side],places=9)
        self.assertTrue(d['A_B_CI_NEB_converged'])
        self.assertFalse(d['A_B_IRC_confirmed'])

    def test_conformer_sampled_peaks(self):
        d, bands = load('energy-ladder.json'),load('conformer-bands.json')
        for index, phase in enumerate(['NEB','CI-NEB']):
            energy = max(p['energy_eV'] for p in bands[phase]['images'])/Hartree
            for side in ['A','B']:
                ref = d['ORCA_endpoint_reference_E_hartree'][side]
                self.assertAlmostEqual(d['conformer_rows'][index][side],(energy-ref)*KCAL,places=8)
        self.assertLess(load('conformer-audit.json')['bands']['CI-NEB']['whole_band_fmax_eV_A'],.01)

    def test_full_irc(self):
        d = load('irc.json')
        self.assertEqual(len(d['points']),55)
        self.assertEqual(d['points'][34]['branch'],'TS')
        self.assertEqual(sum(p['branch']=='reactant' for p in d['points']),34)
        self.assertEqual(sum(p['branch']=='product' for p in d['points']),20)
        self.assertTrue(all(a['s']<b['s'] for a,b in zip(d['points'],d['points'][1:])))
        for p in d['points']:
            x = np.array(p['coordinates_A'])
            self.assertTrue(np.allclose([np.linalg.norm(x[0]-x[4]),np.linalg.norm(x[3]-x[5])],p['forming_bonds_A'],atol=1e-7))
            self.assertAlmostEqual(p['relative_electronic_energy_kcal_mol'],(p['energy_hartree']-d['ts_energy_hartree'])*KCAL,places=8)

    def test_cycloaddition_barriers(self):
        d = load('path-comparison.json')
        for key, expected in [('NEB',15.736025496514529),('CI_NEB',15.738138404366314),('Gaussian_TS',15.740742207185244)]:
            self.assertAlmostEqual(d[key]['forward_complex_barrier_kcal_mol'],expected,places=8)
        self.assertLess(d['CI_NEB']['whole_band_effective_fmax_eV_A'],.02)

    def test_report_and_pdf(self):
        directory = ROOT/'experiment-report'
        tex = (directory/'cycloaddition_mechanism.tex').read_text()
        figures = re.findall(r'\\includegraphics\[[^]]*\]\{([^}]+)\}',tex)
        self.assertEqual(len(figures),11)
        for name in figures:
            self.assertTrue((directory/name).is_file())
        pdf = PdfReader(directory/'cycloaddition_mechanism.pdf')
        self.assertEqual(len(pdf.pages),15)
        url = 'https://github.com/populuscathayana/Quantum3'
        self.assertIn(url,tex)
        self.assertNotIn(url,pdf.pages[0].extract_text())
        self.assertIn(url,pdf.pages[-1].extract_text())
        self.assertGreater(tex.index(url),tex.index('\\end{thebibliography}'))
        self.assertTrue(any(x.get_object().get('/A',{}).get('/URI')==url for x in pdf.pages[-1].get('/Annots',[])))
        for value in ['5.5907','0.0272','26.47']:
            self.assertIn(value,'\n'.join(p.extract_text() for p in pdf.pages))


if __name__ == '__main__':
    unittest.main(verbosity=2)
