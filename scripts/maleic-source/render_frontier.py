#!/usr/bin/env python3
"""Render Gaussian cubegen orbitals; never evaluate a new wavefunction."""
import hashlib
import json
from pathlib import Path
import numpy as np
from PIL import ImageDraw
from audit_stage1 import ROOT, BOHR_A, fchk_fields
from postprocess_static import canvas, font, save, INK, MUTED, GRID, view_matrix
from prepare_stage1 import expected_edges

RAW = ROOT/'received/frontier-direct-20261015/Gau3-orbital-direct-20261015-run2'
ISO = .03
COLORS = np.array([[50,111,154],[197,120,62]])


def read_cube(path):
    with path.open() as f:
        comments=[next(f).strip(),next(f).strip()]
        row=next(f).split(); n=int(row[0]); origin=np.array(row[1:4],float)
        rows=[next(f).split() for _ in range(3)]
        shape=tuple(int(r[0]) for r in rows); steps=np.array([r[1:4] for r in rows],float)
        assert all(v>1 for v in shape)
        atoms=[next(f).split() for _ in range(abs(n))]
        numbers=[int(r[0]) for r in atoms]; xyz=np.array([r[2:5] for r in atoms],float)
        orbital=list(map(int,next(f).split())) if n<0 else []
        assert len(orbital)==2 and orbital[0]==1
        values=np.fromstring(f.read().replace('D','E'),sep=' ')
    assert values.size==np.prod(shape) and np.isfinite(values).all()
    return dict(comments=comments,origin=origin,steps=steps,xyz=xyz,numbers=numbers,
                mo=orbital[1],grid=values.reshape(shape))


def surface(cube, iso):
    grid=cube['grid']
    corners=np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]])
    tets=[[0,1,3,7],[0,3,2,7],[0,2,6,7],[0,6,4,7],[0,4,5,7],[0,5,1,7]]
    shape=np.array(grid.shape)-1
    idx=np.array(np.meshgrid(*(np.arange(s) for s in shape),indexing='ij')).reshape(3,-1).T
    vals=np.array([grid[c[0]:c[0]+shape[0],c[1]:c[1]+shape[1],c[2]:c[2]+shape[2]].ravel() for c in corners]).T
    keep=(vals.min(axis=1)<iso)&(vals.max(axis=1)>iso)
    idx,vals=idx[keep],vals[keep]; triangles=[]
    for tet in tets:
        tv=vals[:,tet]
        loc=(idx[:,None,:]+corners[tet][None,:,:])@cube['steps']+cube['origin']
        masks=np.sum((tv>iso)*[1,2,4,8],axis=1)
        for mask in range(1,15):
            rows=masks==mask
            if not rows.any(): continue
            inside=[k for k in range(4) if mask&(1<<k)]
            outside=[k for k in range(4) if not mask&(1<<k)]
            def edge(a,b):
                fraction=(iso-tv[rows,a])/(tv[rows,b]-tv[rows,a])
                return loc[rows,a]+fraction[:,None]*(loc[rows,b]-loc[rows,a])
            if len(inside)==1:
                triangles.append(np.stack([edge(inside[0],j) for j in outside],axis=1))
            elif len(outside)==1:
                triangles.append(np.stack([edge(outside[0],j) for j in inside],axis=1))
            else:
                a,b=inside; c,d=outside
                ac,ad,bc,bd=edge(a,c),edge(a,d),edge(b,c),edge(b,d)
                triangles.extend([np.stack([ac,ad,bc],axis=1),np.stack([ad,bd,bc],axis=1)])
    return np.concatenate(triangles)*BOHR_A


def interpolate(cube, point_A):
    p=(np.asarray(point_A)/BOHR_A-cube['origin'])@np.linalg.inv(cube['steps'])
    i=np.floor(p).astype(int); fraction=p-i
    assert np.all(i>=0) and np.all(i+1<cube['grid'].shape)
    value=0.
    for x in [0,1]:
        for y in [0,1]:
            for z in [0,1]:
                corner=np.array([x,y,z])
                weight=np.prod(np.where(corner,fraction,1-fraction))
                value+=weight*cube['grid'][tuple(i+corner)]
    return float(value)


def fit(source,target):
    a=source.mean(axis=0); b=target.mean(axis=0)
    u,_,vt=np.linalg.svd((source-a).T@(target-b))
    fix=np.diag([1,1,np.linalg.det(u@vt)])
    rotation=u@fix@vt; shift=b-a@rotation
    rms=float(np.sqrt(np.mean(np.sum((source@rotation+shift-target)**2,axis=1))))
    return rotation,shift,rms


def render(im,items,center,scale,view=None):
    view=view_matrix(30) if view is None else view
    d=ImageDraw.Draw(im); primitives=[]
    def project(x): return x[...,:2]*[scale,-scale]+np.array(center)
    for item in items:
        rotation,shift=item.get('transform',(np.eye(3),np.zeros(3)))
        atoms=(item['cube']['xyz']*BOHR_A@rotation+shift)@view.T
        for phase,triangles in enumerate(item['meshes']):
            tri=(triangles@rotation+shift)@view.T
            for t in tri:
                normal=np.cross(t[1]-t[0],t[2]-t[0]); normal/=max(np.linalg.norm(normal),1e-20)
                light=.60+.40*abs(normal@np.array([.2,-.3,.9327379]))
                color=tuple(np.clip(COLORS[phase]*light+255*(1-light)*.15,0,255).astype(int))
                primitives.append((t[:,2].mean(),'tri',project(t),color))
        for a,b in item['edges']:
            a,b=atoms[a-1],atoms[b-1]
            # Clip sticks inside the atom radii before depth ordering.
            r=.11/np.linalg.norm(b-a)
            for t in np.linspace(r,1-r,18,endpoint=False):
                p=a+(b-a)*t; q=a+(b-a)*min(t+(1-2*r)/18,1-r)
                primitives.append(((p[2]+q[2])/2,'bond',project(np.array([p,q])),None))
        for i,p in enumerate(atoms):
            primitives.append((p[2],'atom',project(p),item['cube']['numbers'][i]))
    for _,kind,p,color in sorted(primitives,key=lambda v:v[0]):
        if kind=='tri': d.polygon([tuple(x) for x in p],fill=color)
        elif kind=='bond': d.line([tuple(x) for x in p],fill='#737c82',width=9)
        else:
            radius=18 if color==6 else 12
            d.ellipse((p[0]-radius,p[1]-radius,p[0]+radius,p[1]+radius),
                      fill='#596570' if color==6 else '#eff1f3',outline='#414d57',width=2)
    # A thin skeletal projection keeps atom identities legible behind the
    # otherwise opaque phase surfaces; surface geometry itself is unchanged.
    for item in items:
        rotation,shift=item.get('transform',(np.eye(3),np.zeros(3)))
        atoms=(item['cube']['xyz']*BOHR_A@rotation+shift)@view.T
        points=project(atoms)
        for a,b in item['edges']:
            if item['cube']['numbers'][a-1]==item['cube']['numbers'][b-1]==6:
                d.line([tuple(points[a-1]),tuple(points[b-1])],fill='white',width=9)
                d.line([tuple(points[a-1]),tuple(points[b-1])],fill='#596570',width=4)
        for index,p in enumerate(points):
            if item['cube']['numbers'][index]!=6: continue
            label=index+1 if item['name']=='butadiene-opt' else index+5
            r=20
            d.ellipse((p[0]-r,p[1]-r,p[0]+r,p[1]+r),fill='white',outline=INK,width=2)
            d.text(tuple(p),str(label),font=font(23,True),fill=INK,anchor='mm')


def main():
    state=json.loads((ROOT/'cluster-progress.json').read_text())
    records={}; items=[]
    for name,stem,mo in [('butadiene-opt','butadiene',15),('ethylene-opt','ethylene',9)]:
        role='homo' if mo==15 else 'lumo'
        original=ROOT/state['jobs'][name]['local_directory']/'calc.fchk'
        copy=RAW/f'{stem}.fchk'
        assert copy.read_bytes()==original.read_bytes()
        f=fchk_fields(copy); cube=read_cube(RAW/f'{stem}-{role}.cube')
        assert cube['mo']==mo and cube['numbers']==f['Atomic numbers']
        xyz=np.array(f['Current cartesian coordinates']).reshape(-1,3)
        coord_error=float(np.max(np.abs(cube['xyz']-xyz)))
        assert coord_error<1e-6
        norm=float(np.sum(cube['grid']**2)*abs(np.linalg.det(cube['steps'])))
        assert abs(norm-1)<.005, norm
        boundary=max(float(np.max(np.abs(np.take(cube['grid'],[0,-1],axis=axis)))) for axis in range(3))
        assert boundary<ISO/10
        meshes=[surface(cube,iso) for iso in [ISO,-ISO]]
        items.append(dict(name=name,cube=cube,meshes=meshes,edges=expected_edges(name)))
        records[name]=dict(source=str(original.relative_to(ROOT)),mo=mo,orbital_energy_eV=f['Alpha Orbital Energies'][mo-1]*27.211386245988,
                           grid_shape=list(cube['grid'].shape),grid_steps_bohr=cube['steps'].tolist(),norm_integral=norm,
                           max_coordinate_error_bohr=coord_error,boundary_amplitude_max=boundary,
                           surface_triangle_counts=[len(m) for m in meshes],cube_sha256=hashlib.sha256((RAW/f'{stem}-{role}.cube').read_bytes()).hexdigest())
    im=canvas('Frontier-orbital phase surfaces','B3LYP-D3(BJ)/def2-SVP  |  Isolated optimized reactants  |  Isovalue: ±0.030 bohr⁻³ᐟ²',(2400,1500))
    d=ImageDraw.Draw(im)
    for index,item in enumerate(items):
        x=90+index*1180; center=item['cube']['xyz'].mean(axis=0)*BOHR_A
        item['transform']=(np.eye(3),-center)
        title='(a) 1,3-Butadiene HOMO' if index==0 else '(b) Ethylene LUMO'
        r=records[item['name']]
        d.text((x,240),title,font=font(37,True),fill=INK)
        d.text((x,305),f'MO {r["mo"]}  |  {r["orbital_energy_eV"]:+.3f} eV',font=font(31),fill=MUTED)
        render(im,[item],(x+500,770),scale=162)
        d.text((x,1190),'C1-C2-C3-C4 torsion: 26.90°' if index==0 else 'Antibonding π* orbital',font=font(30),fill=INK)
    d.line((1190,230,1190,1270),fill=GRID,width=2)
    d.text((90,1360),'Blue: positive phase; orange: negative phase. The overall sign of each orbital is arbitrary.',font=font(29),fill=MUTED)
    d.text((90,1415),'Same scale and isovalue; numbered carbon skeletons are overlaid for orientation.',font=font(29),fill=MUTED)
    save(im,'frontier-orbital-surfaces.png')

    tsf=fchk_fields(ROOT/state['jobs']['ts-opt']['local_directory']/'calc.fchk')
    ts=np.array(tsf['Current cartesian coordinates']).reshape(-1,3)*BOHR_A
    for item,indices in zip(items,[[0,1,2,3,6,7,8,9,10,11],[4,5,12,13,14,15]]):
        source=item['cube']['xyz']*BOHR_A
        rotation,shift,rms=fit(source,ts[indices])
        item['transform']=(rotation,shift)
        records[item['name']]['placement_fit_RMSD_A']=rms
        records[item['name']]['placement_rotation']=rotation.tolist()
        records[item['name']]['placement_translation_A']=shift.tolist()
    amplitudes=[]
    for a,b,di,ei in [(0,4,0,0),(3,5,3,1)]:
        direction=(ts[b]-ts[a])/np.linalg.norm(ts[b]-ts[a])
        values=[]
        for item,idx,sign in [(items[0],di,1),(items[1],ei,-1)]:
            rot,shift=item['transform']; atom=item['cube']['xyz'][idx]*BOHR_A@rot+shift
            point=atom+sign*.65*direction
            values.append(interpolate(item['cube'],(point-shift)@rot.T))
        amplitudes.append(values)
    products=np.prod(amplitudes,axis=1)
    assert products[0]*products[1]>0
    sign=1 if products[0]>0 else -1
    if sign<0: items[1]['meshes']=list(reversed(items[1]['meshes']))
    mid=ts.mean(axis=0)
    for item in items:
        r,t=item['transform'];item['transform']=(r,t-mid)
    im=canvas('Orbital phases along the two bond-forming contacts',
              'Isolated-reactant orbitals rigidly placed in the TS approach geometry',(2200,1700))
    render(im,items,(1070,820),scale=170,view=view_matrix(42))
    d=ImageDraw.Draw(im)
    d.text((95,235),'Ethylene LUMO (C5, C6)',font=font(34,True),fill=INK)
    d.text((95,1295),'1,3-Butadiene HOMO (C1-C4)',font=font(34,True),fill=INK)
    d.text((95,1395),'C1-C5 and C4-C6: the facing lobes can have matching phase at both contacts.',font=font(30),fill=INK)
    d.text((95,1460),'Isovalue: ±0.030 bohr⁻³ᐟ². Molecular geometries are preserved during rigid placement.',font=font(28),fill=MUTED)
    d.text((95,1520),'The ethylene orbital is globally sign-aligned; this does not change any physical result.',font=font(28),fill=MUTED)
    d.text((95,1580),'This is a fragment-orbital comparison, not a molecular orbital of the transition state.',font=font(28),fill=MUTED)
    save(im,'frontier-phase-matching.png')
    result=dict(isovalue_au=ISO,orbitals=records,renderer='marching tetrahedra; depth-sorted orthographic projection',
                contact_samples_distance_A=.65,contact_amplitudes_original=amplitudes,ethylene_global_phase=sign,
                contacts_same_phase_after_global_alignment=True,
                cube_title_note='Inherited Gaussian job title contains initial structure; cube coordinates independently matched to final optimized fchk.')
    (ROOT/'analysis/postprocessing/orbital-validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
