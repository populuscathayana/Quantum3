"""Front/side measured geometries and separately labelled endo/exo seed views."""
import json
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read
from PIL import ImageDraw
from render_maleic_results import molecule,fit,DATA,DEST
from postprocess_static import canvas,font,INK,MUTED,GRID
from maleic_approach import ROOT,NUMBERS,approach,exo_seed,SOURCE
from maleic_orca_backend import sha,save_json


def rx(deg):
    c,s=np.cos(np.radians(deg)),np.sin(np.radians(deg))
    return np.array([[1,0,0],[0,c,-s],[0,s,c]])


def ry(deg):
    c,s=np.cos(np.radians(deg)),np.sin(np.radians(deg))
    return np.array([[c,0,s],[0,1,0],[-s,0,c]])


def reference_frame(x):
    ex=x[3]-x[0];ex/=np.linalg.norm(ex)
    ey=x[[1,2]].mean(0)-x[[0,3]].mean(0);ey-=ey@ex*ex;ey/=np.linalg.norm(ey)
    frame=np.column_stack([ex,ey,np.cross(ex,ey)])
    assert np.linalg.det(frame)>.999999
    return frame


def render():
    data=json.loads(DATA.read_text());g=data['species']['Gaussian']
    target=np.array(g['TS']['positions_A']);basis=reference_frame(target)
    audit=json.loads((ROOT/'analysis/maleic-ts-review-20261021/audit.json').read_text())
    mode=np.array(audit['imaginary_mode_cartesian'])@basis
    cameras=[rx(35)@ry(12),rx(12)@ry(85)]
    im=canvas('Maleic anhydride addition: endo pathway',
              'Gaussian optimized geometries | Two views of each structure | TS: 442.956i cm−1',size=(2400,1650))
    d=ImageDraw.Draw(im);checks={}
    for cx,name,label in [(430,'RC','Reactant complex'),(1220,'TS','Endo transition state'),(2010,'product','Cis adduct')]:
        original=np.array(g[name]['positions_A']);x=fit(original,target)@basis
        d.text((cx,240),label,font=font(34,True),fill=INK,anchor='mm')
        for row,(cy,camera) in enumerate(zip([555,1130],cameras)):
            molecule(im,x,NUMBERS,(cx,cy),scale=100,mode=mode if name=='TS' else None,camera=camera)
        distances=[float(np.linalg.norm(original[j]-original[i])) for i,j in [(0,4),(3,5)]]
        d.text((cx,1460),f'{distances[0]:.4f} / {distances[1]:.4f} Å',font=font(32),fill=INK,anchor='mm')
        checks[name]=dict(original_positions_A=original.tolist(),display_positions_A=x.tolist(),forming_distances_A=distances,
                         pairwise_distance_max_error_A=float(np.max(abs(Atoms(numbers=NUMBERS,positions=x).get_all_distances()-Atoms(numbers=NUMBERS,positions=original).get_all_distances()))))
        assert checks[name]['pairwise_distance_max_error_A']<1e-12
    d.text((70,310),'Front',font=font(25,True),fill=MUTED)
    d.line((70,840,2330,840),fill=GRID,width=2)
    d.text((70,905),'Side',font=font(25,True),fill=MUTED)
    d.text((100,1540),'C: blue   O: red   H: pale gray | Distances: C1–C5 / C4–C6',font=font(28),fill=MUTED)
    d.text((100,1590),'Orange arrows: the TS imaginary-mode displacement; both views use the same coordinates.',font=font(26),fill=MUTED)
    image=DEST/'maleic-structures-multiview.png';im.save(image,dpi=(300,300))
    # A separate preview makes orientation explicit without presenting guesses as results.
    exo,seed_audit=exo_seed();endo=read(SOURCE/'seeds/ts-distance-seed.xyz')
    seed_basis=reference_frame(endo.positions)
    preview=canvas('Endo / exo approach: initial geometries',
        'Unoptimized distance-guided seeds | Same fragments; two forming distances fixed at 2.30 Å',size=(2000,1500))
    pd=ImageDraw.Draw(preview)
    for cx,a,label in [(520,endo,'ENDO: anhydride inward'),(1510,exo,'EXO: anhydride outward')]:
        pd.text((cx,250),label,font=font(33,True),fill=INK,anchor='mm')
        for cy,camera in zip([560,1130],cameras):molecule(preview,a.positions@seed_basis,NUMBERS,(cx,cy),100,camera=camera)
    pd.text((45,310),'Front',font=font(24,True),fill=MUTED);pd.text((45,890),'Side',font=font(24,True),fill=MUTED)
    pd.line((70,825,1930,825),fill=GRID,width=2)
    pd.text((90,1420),'These are starting geometries, not optimized TS structures or calculated barriers.',font=font(26),fill=MUTED)
    destination=ROOT/'analysis/maleic-endo-exo-20261022/approach-seeds.png';preview.save(destination,dpi=(300,300))
    record=dict(source=str(DATA.relative_to(ROOT)),source_sha256=sha(DATA),Gaussian_TS_classification=approach(Atoms(numbers=NUMBERS,positions=target)),
        measured_structure_figure=str(image.relative_to(ROOT)),measured_structure_figure_sha256=sha(image),
        initial_guess_figure=str(destination.relative_to(ROOT)),seed_audit=seed_audit,cameras=[c.tolist() for c in cameras],
        proper_cameras=True,checks=checks,comparison_result_figure='Pending new exo TS and its verified endpoints; no fabricated structures/energies')
    save_json(ROOT/'analysis/maleic-endo-exo-20261022/figure-audit.json',record)
    print(json.dumps({'measured_structure_figure':str(image),'initial_guess_figure':str(destination)},indent=2))


if __name__=='__main__':render()
