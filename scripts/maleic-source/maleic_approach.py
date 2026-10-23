"""Source-backed endo/exo approach geometry; no electronic calculations."""
from pathlib import Path
import json
import numpy as np
from ase.io import read
from maleic_common import NUMBERS, BASE_EDGES, FORMING, graph, geometry
from maleic_orca_backend import ROOT, sha

SOURCE = ROOT/'Gau3-local-orca-maleic-da-20261021'
MALEIC_INDICES = [4,5,6,7,8,9,10,17,18]
DIENE_INDICES = [0,1,2,3,11,12,13,14,15,16]


def approach(atoms):
    """Signed anhydride projection toward the diene interior (TS/complex only)."""
    assert list(atoms.numbers)==NUMBERS
    x=atoms.positions
    axis=x[3]-x[0];axis/=np.linalg.norm(axis)
    inside=x[[1,2]].mean(0)-x[[0,3]].mean(0)
    inside-=inside@axis*axis;inside/=np.linalg.norm(inside)
    substituent=x[[6,7,10]].mean(0)-x[[4,5]].mean(0)
    value=float(substituent@inside)
    return dict(projection_toward_diene_interior_A=value,
                approach='endo' if value>.35 else 'exo' if value<-.35 else 'intermediate',
                definition='C7/C8/O11 centroid relative to C5/C6, projected toward C2/C3 from C1/C4, perpendicular to terminal axis; not a product stereoisomer label')


def exo_seed():
    manifest=json.loads((SOURCE/'manifest.json').read_text())
    rel='seeds/ts-distance-seed.xyz'
    assert sha(SOURCE/rel)==manifest[rel]
    old=read(SOURCE/rel);geometry(old,'ts-distance-seed')
    assert approach(old)['approach']=='endo'
    new=old.copy();x=old.positions
    center=x[[4,5]].mean(0);axis=x[5]-x[4];axis/=np.linalg.norm(axis)
    # Rodrigues rotation at pi: determinant +1, never reflection.
    rotation=2*np.outer(axis,axis)-np.eye(3)
    new.positions[MALEIC_INDICES]=(x[MALEIC_INDICES]-center)@rotation.T+center
    assert np.linalg.det(rotation)>.999999999
    assert np.max(abs(new.positions[DIENE_INDICES]-x[DIENE_INDICES]))==0
    for indices in [MALEIC_INDICES,DIENE_INDICES]:
        assert np.max(abs(new[indices].get_all_distances()-old[indices].get_all_distances()))<1e-12
    assert graph(new)==BASE_EDGES
    assert all(abs(new.get_distance(*b)-2.30)<1e-8 for b in FORMING)
    assert approach(new)['approach']=='exo'
    geometry(new,'ts-distance-seed')
    record=dict(source_directory=str(SOURCE),source_file=rel,source_sha256=sha(SOURCE/rel),
                source_is_unoptimized_seed=True,old=approach(old),new=approach(new),
                rotation=rotation.tolist(),rotation_determinant=float(np.linalg.det(rotation)),
                rotated_atoms_one_based=[i+1 for i in MALEIC_INDICES],
                same_fragment_internal_geometry=True,forming_distances_A=[new.get_distance(*b) for b in sorted(FORMING)],
                source_optimized_TS_or_IRC_used=False,geometry=geometry(new,'ts-distance-seed'))
    return new,record
