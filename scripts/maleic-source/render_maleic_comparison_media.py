"""Source-backed endo/exo modes and isolated-fragment orbital comparisons.

Only display transforms are applied; no quantum calculations or path inference.
"""
import json, hashlib
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw
from ase.io import read
from ase.units import Bohr
from render_frontier import read_cube, surface, interpolate, fit, COLORS
from render_maleic_results import molecule
from render_maleic_multiview import reference_frame,rx,ry
from maleic_common import graph
from postprocess_static import font,INK

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiment-report/figures'
AUDIT=ROOT/'analysis/maleic-endo-exo-20261022'
def load(path):return json.loads((ROOT/path).read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def alignment(x,target):
    r,t,_=fit(x[:4],target[:4]);assert np.linalg.det(r)>.999999
    return r,t

def mode_media():
    audits={('Gaussian','endo'):load('analysis/maleic-ts-review-20261021/audit.json'),
      ('Gaussian','exo'):load('analysis/maleic-endo-exo-20261022/exo-gaussian-independent-audit.json'),
      ('ORCA','endo'):load('analysis/maleic-stage1-review-20261021/audit.json')['orca']['ts'],
      ('ORCA','exo'):load('analysis/maleic-endo-exo-20261022/exo-orca-independent-audit.json')}
    target=np.array(audits['Gaussian','endo']['positions_A']);basis=reference_frame(target)
    frames={};metadata={}
    for (program,approach),a in audits.items():
        if program=='Gaussian':x=np.array(a['positions_A']);numbers=a['atomic_numbers']
        else:
            run='Gau3-local-orca-maleic-da-20261021' if approach=='endo' else 'Gau3-local-orca-maleic-exo-20261022'
            atoms=read(ROOT/run/'ts-opt/calc.xyz');x=atoms.positions;numbers=atoms.numbers
        r,t=alignment(x,target);v=np.array(a['imaginary_mode_cartesian'])@r@basis
        x=(x@r+t)@basis
        # Normalize all modes to the same maximum atomic amplitude; choose bond shortening phase.
        derivative=sum(np.dot((x[j]-x[i])/np.linalg.norm(x[j]-x[i]),v[j]-v[i]) for i,j in [(0,4),(3,5)])
        if derivative>0:v=-v
        v=v/np.max(np.linalg.norm(v,axis=1))*.22
        freq=abs(a['frequency_cm_1'][0]);key=f'{program.lower()}-{approach}'
        fs=[]
        for phase in np.linspace(0,2*np.pi,32,endpoint=False):
            im=Image.new('RGB',(900,1000),'white');d=ImageDraw.Draw(im)
            d.text((450,38),f'{program} {approach}  {freq:.1f}i cm⁻¹',font=font(29,True),fill=INK,anchor='mm')
            for cy,camera in zip([280,740],[rx(35)@ry(12),rx(12)@ry(85)]):
                molecule(im,x+np.sin(phase)*v,numbers,(450,cy),scale=78,camera=camera,labels=False)
            fs.append(im)
        frames[key]=fs
        fs[0].save(OUT/f'maleic-mode-{key}.gif',save_all=True,append_images=fs[1:],duration=90,loop=0)
        metadata[key]={'frequency_cm_1':freq,'max_display_displacement_A':.22,'rotation_det':float(np.linalg.det(r))}
    groups={
      'gaussian-endo-exo':['gaussian-endo','gaussian-exo'],
      'orca-endo-exo':['orca-endo','orca-exo'],
      'endo-programs':['gaussian-endo','orca-endo'],
      'exo-programs':['gaussian-exo','orca-exo']}
    for name,keys in groups.items():
        joined=[]
        for i in range(32):
            im=Image.new('RGB',(1800,1000),'white')
            for j,key in enumerate(keys):im.paste(frames[key][i],(900*j,0))
            joined.append(im)
        joined[0].save(OUT/f'maleic-mode-{name}.gif',save_all=True,append_images=joined[1:],duration=90,loop=0)
    return metadata

def draw_orbitals(im,items,center,camera,scale):
    d=ImageDraw.Draw(im);primitives=[]
    def xy(p):return p[...,:2]*[scale,-scale]+center
    for item in items:
        atoms=item['positions']@camera.T
        for phase,mesh in enumerate(item['meshes']):
            triangles=mesh@camera.T
            normals=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])
            normals/=np.maximum(np.linalg.norm(normals,axis=1)[:,None],1e-20)
            light=.6+.4*np.abs(normals@np.array([.2,-.3,.9327379]))
            colors=np.clip(COLORS[phase]*light[:,None]+255*(1-light[:,None])*.15,0,255).astype(int)
            for z,p,c in zip(triangles[:,:,2].mean(1),xy(triangles),colors):primitives.append((z,'tri',p,tuple(c)))
        for i,j in item['edges']:
            for t in np.linspace(.05,.95,15):
                p=atoms[i]*(1-t)+atoms[j]*t;q=atoms[i]*(1-t-.04)+atoms[j]*(t+.04)
                primitives.append(((p[2]+q[2])/2,'bond',xy(np.array([p,q])),None))
        for p,n in zip(atoms,item['numbers']):primitives.append((p[2],'atom',xy(p),int(n)))
    for _,kind,p,c in sorted(primitives,key=lambda a:a[0]):
        if kind=='tri':d.polygon([tuple(q) for q in p],fill=c)
        elif kind=='bond':d.line([tuple(q) for q in p],fill='#647480',width=5)
        else:
            radius=8 if c==1 else 12;color='#eeeef0' if c==1 else '#b4544b' if c==8 else '#526775'
            d.ellipse((p[0]-radius,p[1]-radius,p[0]+radius,p[1]+radius),fill=color,outline=INK,width=1)

def orbital_media():
    exported=load('analysis/maleic-orbitals-20261021/export.json');sources=[];meta={}
    for name in ['butadiene','maleic-anhydride']:
        entry=exported[name];p=ROOT/entry['cube'];assert sha(p)==entry['cube_sha256']
        cube=read_cube(p);assert cube['mo']==entry['orbital_index_zero_based']
        atoms=read(ROOT/entry['source_directory']/'calc.xyz')
        assert np.max(abs(cube['xyz']*Bohr-atoms.positions))<6e-7
        sources.append(dict(cube=cube,meshes=[surface(cube,v) for v in [.03,-.03]],edges=graph(atoms),numbers=atoms.numbers))
        meta[name]={'cube_sha256':sha(p),'energy_eV':entry['energy_eV']}
    target=read(ROOT/'Gau3-local-orca-maleic-da-20261021/ts-opt/calc.xyz').positions
    basis=reference_frame(target);maps=[[0,1,2,3,11,12,13,14,15,16],[4,5,6,7,8,9,10,17,18]]
    oriented={}
    for approach,run in [('endo','Gau3-local-orca-maleic-da-20261021'),('exo','Gau3-local-orca-maleic-exo-20261022')]:
        ts=read(ROOT/run/'ts-opt/calc.xyz').positions;r,t=alignment(ts,target);ts=(ts@r+t)@basis;ts-=ts.mean(0)
        items=[];products=[];transforms=[]
        for source,mapping in zip(sources,maps):
            rot,shift,rms=fit(source['cube']['xyz']*Bohr,ts[mapping]);assert np.linalg.det(rot)>.999999
            transforms.append((rot,shift,rms))
            items.append(dict(positions=source['cube']['xyz']*Bohr@rot+shift,
                meshes=[mesh@rot+shift for mesh in source['meshes']],edges=source['edges'],numbers=source['numbers']))
        for di,mi,i,j in [(0,0,0,4),(3,1,3,5)]:
            values=[]
            for k,idx,partner in [(0,di,j),(1,mi,i)]:
                rot,shift,_=transforms[k];p=items[k]['positions'][idx];direction=ts[partner]-p;direction/=np.linalg.norm(direction)
                values.append(interpolate(sources[k]['cube'],(p+.45*direction-shift)@rot.T))
            assert min(map(abs,values))>.005
            products.append(float(np.prod(values)))
        assert products[0]*products[1]>0
        # Apply one global LUMO sign reversal, when needed, consistently to both contacts.
        if products[0]<0:items[1]['meshes']=items[1]['meshes'][::-1]
        meta[approach]={'phase_products_before_sign_choice':products,'fragment_RMSD_A':[v[2] for v in transforms],
                       'LUMO_global_sign_reversed':products[0]<0}
        oriented[approach]=items
    still=Image.new('RGB',(1800,1400),'white');d=ImageDraw.Draw(still)
    for cx,approach in [(470,'endo'),(1330,'exo')]:
        d.text((cx,55),approach,font=font(37,True),fill=INK,anchor='mm')
        for cy,camera in zip([400,1020],[rx(35)@ry(12),rx(12)@ry(85)]):draw_orbitals(still,oriented[approach],np.array([cx,cy]),camera,90)
    still.save(OUT/'maleic-orbital-endo-exo.png',dpi=(300,300))
    all_frames={}
    for approach in ['endo','exo']:
        still.crop((0 if approach=='endo' else 900,0,900 if approach=='endo' else 1800,1400)).save(OUT/f'maleic-orbital-{approach}.png')
        frames=[]
        for angle in np.linspace(0,360,24,endpoint=False):
            im=Image.new('RGB',(900,800),'white');d=ImageDraw.Draw(im)
            d.text((450,40),f'{approach}  HOMO / LUMO',font=font(29,True),fill=INK,anchor='mm')
            draw_orbitals(im,oriented[approach],np.array([450,420]),rx(25)@ry(angle),82)
            frames.append(im)
        frames[0].save(OUT/f'maleic-orbital-{approach}.gif',save_all=True,append_images=frames[1:],duration=130,loop=0)
        all_frames[approach]=frames
    joined=[]
    for i in range(24):
        im=Image.new('RGB',(1800,800),'white');im.paste(all_frames['endo'][i],(0,0));im.paste(all_frames['exo'][i],(900,0));joined.append(im)
    joined[0].save(OUT/'maleic-orbital-endo-exo.gif',save_all=True,append_images=joined[1:],duration=130,loop=0)
    return meta

if __name__=='__main__':
    modes=mode_media();print('Mode animations complete',flush=True)
    orbitals=orbital_media();print('Orbital figures and rotations complete',flush=True)
    record={'modes':modes,'orbitals':orbitals,'mode_note':'Synchronized phase and normalized illustrative amplitude; not real-time dynamics or IRC.',
            'orbital_note':'Rigidly fitted isolated ORCA orbitals at ±0.03 bohr^-3/2; rotation only, not a reaction trajectory.',
            'files':{p.name:sha(p) for pattern in ['maleic-mode-*.gif','maleic-orbital-*.png','maleic-orbital-*.gif'] for p in OUT.glob(pattern)}}
    (AUDIT/'comparison-media-audit.json').write_text(json.dumps(record,indent=2))
