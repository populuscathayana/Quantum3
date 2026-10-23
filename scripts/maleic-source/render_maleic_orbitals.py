"""Render and validate true frontier orbital phases from the new ORCA GBWs."""
import json
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read
from ase.units import Bohr
from PIL import ImageDraw
from maleic_common import graph
from maleic_orca_backend import sha, save_json
from render_frontier import read_cube, surface, interpolate, fit, COLORS
from postprocess_static import canvas, font, view_matrix, INK, MUTED, GRID

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/'analysis/maleic-orbitals-20261021'


def render(im,item,center,scale=150):
    d=ImageDraw.Draw(im);view=view_matrix(35);cube=item['cube']
    shift=-cube['xyz'].mean(0)*Bohr
    atoms=(cube['xyz']*Bohr+shift)@view.T
    def xy(x):return x[...,:2]*[scale,-scale]+center
    primitives=[]
    for phase,triangles in enumerate(item['meshes']):
        for t in (triangles+shift)@view.T:
            n=np.cross(t[1]-t[0],t[2]-t[0]);n/=max(np.linalg.norm(n),1e-20)
            light=.60+.40*abs(n@np.array([.2,-.3,.9327379]))
            color=tuple(np.clip(COLORS[phase]*light+255*(1-light)*.15,0,255).astype(int))
            primitives.append((t[:,2].mean(),'tri',xy(t),color))
    for i,j in item['edges']:
        a,b=atoms[i],atoms[j]
        for t in np.linspace(.05,.95,20):
            p=a+(b-a)*t;q=a+(b-a)*min(t+.045,.96)
            primitives.append(((p[2]+q[2])/2,'bond',xy(np.array([p,q])),None))
    for i,p in enumerate(atoms):primitives.append((p[2],'atom',xy(p),cube['numbers'][i]))
    for _,kind,p,color in sorted(primitives,key=lambda q:q[0]):
        if kind=='tri':d.polygon([tuple(v) for v in p],fill=color)
        elif kind=='bond':d.line([tuple(v) for v in p],fill='#637682',width=8)
        else:
            radius=12 if color==1 else 18
            fill='#f0f2f4' if color==1 else '#ae514b' if color==8 else '#5b6977'
            d.ellipse((p[0]-radius,p[1]-radius,p[0]+radius,p[1]+radius),fill=fill,outline=INK,width=2)
    # Label only the reacting carbon sites, through a light skeleton overlay.
    for i,label in item['labels'].items():
        p=xy(atoms[i]);d.ellipse((p[0]-22,p[1]-22,p[0]+22,p[1]+22),fill='white',outline=INK,width=2)
        d.text(tuple(p),label,font=font(18,True),fill=INK,anchor='mm')


def main():
    exported=json.loads((DIR/'export.json').read_text());items=[];audit={}
    for kind,labels in [('butadiene',{0:'C1',3:'C4'}),('maleic-anhydride',{0:'C5',1:'C6'})]:
        r=exported[kind];cube=read_cube(ROOT/r['cube']);source=read(ROOT/r['source_directory']/'calc.xyz')
        assert sha(ROOT/r['cube'])==r['cube_sha256'] and cube['mo']==r['orbital_index_zero_based']
        assert cube['numbers']==source.numbers.tolist()
        error=float(np.max(abs(cube['xyz']*Bohr-source.positions)));assert error<6e-7
        norm=float(np.sum(cube['grid']**2)*abs(np.linalg.det(cube['steps'])))
        boundary=max(float(np.max(abs(np.take(cube['grid'],[0,-1],axis=k)))) for k in range(3))
        assert abs(norm-1)<.005 and boundary<.003
        item=dict(name=kind,cube=cube,edges=graph(source),labels=labels,meshes=[surface(cube,v) for v in [.03,-.03]])
        items.append(item)
        audit[kind]=dict(r,norm_integral=norm,boundary_max_amplitude=boundary,coordinate_max_error_A=error,
                         grid_shape=list(cube['grid'].shape),grid_steps_bohr=cube['steps'].tolist(),isovalue_bohr_m32=.03)
    im=canvas('Butadiene + maleic anhydride: frontier orbital phases',
              'ORCA B3LYP/G-D3(BJ)/def2-SVP | Isolated optimized reactants | Isovalue ±0.030 bohr−3/2',size=(2400,1500))
    d=ImageDraw.Draw(im)
    for i,item in enumerate(items):
        x=90+i*1190;kind=item['name'];r=exported[kind]
        title='1,3-Butadiene HOMO' if i==0 else 'Maleic anhydride LUMO'
        d.text((x,240),title,font=font(38,True),fill=INK)
        d.text((x,305),f'Energy: {r["energy_eV"]:+.4f} eV',font=font(30),fill=MUTED)
        render(im,item,(x+500,800))
    d.line((1200,235,1200,1240),fill=GRID,width=2)
    d.text((90,1320),'Blue / orange: opposite orbital phases. The overall sign of each orbital is arbitrary.',font=font(29),fill=MUTED)
    d.text((90,1380),f'HOMO(diene) → LUMO(anhydride) orbital-energy gap: {exported["maleic-anhydride"]["energy_eV"]-exported["butadiene"]["energy_eV"]:.4f} eV',font=font(31,True),fill=INK)
    im.save(ROOT/'experiment-report/figures/maleic-frontier-phases.png',dpi=(300,300))
    # Evaluate facing lobes after a rigid, proper fit of each isolated fragment to the TS.
    ts=read(ROOT/'Gau3-local-orca-maleic-da-20261021/ts-opt/calc.xyz').positions
    maps=[[0,1,2,3,11,12,13,14,15,16],[4,5,6,7,8,9,10,17,18]]
    transforms=[fit(item['cube']['xyz']*Bohr,ts[mapping]) for item,mapping in zip(items,maps)]
    products=[];sample=[]
    for i,j,di,mi in [(0,4,0,0),(3,5,3,1)]:
        values=[]
        for item,(rot,shift,rms),atom,partner in [(items[0],transforms[0],di,j),(items[1],transforms[1],mi,i)]:
            position=item['cube']['xyz'][atom]*Bohr@rot+shift
            direction=ts[partner]-position;direction/=np.linalg.norm(direction)
            point=position+.45*direction
            values.append(interpolate(item['cube'],(point-shift)@rot.T))
        assert min(map(abs,values))>.005
        products.append(float(np.prod(values)));sample.append(dict(pair=[i+1,j+1],orbital_amplitudes=values))
    assert products[0]*products[1]>0
    record=dict(orbitals=audit,phase_samples=sample,phase_products=products,
                same_relative_phase_at_both_contacts=True,fragment_fit_RMSD_A=[v[2] for v in transforms],
                sample_offset_A=.45,interpretation='Both contact phase products have the same sign; one global sign choice makes both contacts constructive. These are isolated-fragment orbitals, not interacting-TS orbitals.',
                source='https://www.faccts.de/docs/orca/6.1/manual/contents/utilitiesvisualization/plots.html')
    save_json(DIR/'audit.json',record)
    print(json.dumps(dict(norms={k:v['norm_integral'] for k,v in audit.items()},phase_products=products,
                         fit_RMSD_A=record['fragment_fit_RMSD_A']),indent=2))


if __name__=='__main__':main()
