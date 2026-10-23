"""Render the new reaction's measured structures, normal mode and energy data."""
import json
import math
from pathlib import Path
import numpy as np
from ase import Atoms
from PIL import Image, ImageDraw
from maleic_common import graph, FORMING
from postprocess_static import canvas, font, dashed, arrow, INK, MUTED, GRID, BLUE, ORANGE

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'analysis/maleic-postprocessing-20261021/data.json'
DEST = ROOT/'experiment-report/figures'
RED = '#ba554d'


def fit(x, target):
    x = np.asarray(x); target = np.asarray(target)
    u, _, vt = np.linalg.svd((x-x.mean(0)).T@(target-target.mean(0)))
    r = u@np.diag([1, 1, np.linalg.det(u@vt)])@vt
    assert np.linalg.det(r) > .999999
    return (x-x.mean(0))@r+target.mean(0)


def molecule(im, x, numbers, center, scale=120, mode=None, camera=None, labels=True):
    x = np.asarray(x); d = ImageDraw.Draw(im)
    # A common proper camera rotation is used for every aligned geometry.
    a = math.radians(35)
    r = np.array([[1,0,0],[0,math.cos(a),math.sin(a)],[0,-math.sin(a),math.cos(a)]])
    if camera is not None:
        r=np.asarray(camera)
        assert np.allclose(r@r.T,np.eye(3)) and np.linalg.det(r)>.999999
    p = (x-x.mean(0))@r.T
    xy = p[:,:2]*[scale,-scale]+center
    edges = graph(Atoms(numbers=numbers,positions=x))
    for i,j in sorted(edges,key=lambda e:p[list(e),2].mean()):
        d.line([tuple(xy[i]),tuple(xy[j])],fill='#84949e',width=9)
    for i,j in sorted(FORMING):
        if (i,j) not in edges:
            dashed(d,xy[i],xy[j],ORANGE,4,10,8)
    for i in np.argsort(p[:,2]):
        radius = 20 if numbers[i] == 1 else 28
        color = '#eef1f3' if numbers[i] == 1 else RED if numbers[i] == 8 else BLUE
        q=xy[i];d.ellipse((q[0]-radius,q[1]-radius,q[0]+radius,q[1]+radius),fill=color,outline='#536a79',width=2)
        if numbers[i] != 1 and labels:
            label = ('O' if numbers[i] == 8 else 'C')+str(i+1)
            d.text(tuple(q),label,font=font(16,True),fill='white',anchor='mm')
    if mode is not None:
        vectors=np.asarray(mode)@r.T
        for i in range(len(x)):
            if np.linalg.norm(vectors[i]) > .035:
                arrow(d,xy[i],xy[i]+vectors[i,:2]*[230,-230],fill=ORANGE,width=5,head=13)


def axis(d, box, xr, yr, xt, yt, xlabel, ylabel):
    x0,y0,x1,y1=box
    def xy(x,y):
        return (x0+(x-xr[0])/(xr[1]-xr[0])*(x1-x0),y1-(y-yr[0])/(yr[1]-yr[0])*(y1-y0))
    for val in yt:
        yy=xy(0,val)[1];d.line((x0,yy,x1,yy),fill=GRID,width=2)
        d.text((x0-20,yy),f'{val:g}',font=font(26),fill=MUTED,anchor='rm')
    for val in xt:
        xx=xy(val,0)[0];d.line((xx,y1,xx,y1+10),fill=INK,width=2)
        d.text((xx,y1+23),f'{val:g}',font=font(26),fill=MUTED,anchor='mt')
    d.line((x0,y0,x0,y1,x1,y1),fill=INK,width=3)
    d.text(((x0+x1)/2,y1+72),xlabel,font=font(28),fill=INK,anchor='mt')
    d.text((x0,y0-55),ylabel,font=font(28,True),fill=INK)
    return xy


def irc_figure(data):
    points=data['Gaussian_IRC']; x=[p['coordinate_sqrt_amu_bohr'] for p in points]
    im=canvas('Butadiene + maleic anhydride: Gaussian IRC',
              'B3LYP-D3(BJ)/def2-SVP | Electronic energies relative to separated reactants',size=(2400,1500))
    d=ImageDraw.Draw(im)
    xr=(min(x)-.2,max(x)+.2); ticks=list(range(math.ceil(xr[0]),math.floor(xr[1])+1,2))
    xy=axis(d,(170,280,2270,835),xr,(-55,8),ticks,[-50,-40,-30,-20,-10,0],
            'Mass-weighted IRC coordinate / (amu^1/2 bohr)','Electronic energy / (kcal mol−1)')
    line=[xy(s,p['relative_energy_kcal_mol']) for s,p in zip(x,points)]
    d.line(line,fill=BLUE,width=5)
    for px,py in line:d.ellipse((px-4,py-4,px+4,py+4),fill=BLUE)
    tx,ty=xy(0,data['relative_to_separated_reactants']['Gaussian']['TS']['E_kcal_mol'])
    d.text((tx+25,ty-15),'TS  +2.552',font=font(29,True),fill=ORANGE,anchor='lb')
    d.text((230,735),'Reactant side: 22 points',font=font(27),fill=MUTED)
    d.text((1630,350),'Product side: 30 points',font=font(27),fill=MUTED)
    xy2=axis(d,(170,1070,2270,1330),xr,(1.45,3.3),ticks,[1.5,2,2.5,3],
             'Mass-weighted IRC coordinate / (amu^1/2 bohr)','Forming-bond distance / Å')
    for k,color in [(0,BLUE),(1,ORANGE)]:
        line=[xy2(s,p['forming_distances_A'][k]) for s,p in zip(x,points)]
        if k==0:d.line(line,fill=color,width=6)
        else:
            for p,q in zip(line,line[1:]):dashed(d,p,q,color,3,8,6)
    d.text((1750,1010),'C1–C5     C4–C6',font=font(27),fill=MUTED)
    im.save(DEST/'maleic-gaussian-irc.png',dpi=(300,300))


def ladder(data):
    im=canvas('Butadiene + maleic anhydride: energy levels',
              'Each program uses its own separated-reactant reference | 298.15 K, 1 atm, symmetry-corrected RRHO',size=(2400,1550))
    d=ImageDraw.Draw(im);vals=data['relative_to_separated_reactants']
    for panel,key,title,yr,yt in [(0,'E_kcal_mol','Electronic energy',(-56,8),[-50,-40,-30,-20,-10,0]),
                                  (1,'G_symmetry_corrected_kcal_mol','Gibbs free energy',(-38,22),[-30,-20,-10,0,10,20])]:
        left=170+panel*1180; right=left+1010
        xy=axis(d,(left,320,right,1250),(0,3),yr,[],yt,'Reaction states (schematic spacing)',title+' / (kcal mol−1)')
        labels=['Diene + MA','Complex','TS','Cis product']; names=[None,'RC','TS','product']
        xx=[.26,1.08,1.90,2.72]; centers=[]
        for x,name,label in zip(xx,names,labels):
            val=vals['Gaussian'][name][key] if name else 0
            a,b=xy(x-.20,val),xy(x+.20,val);centers.append((a,b))
            d.line([a,b],fill=BLUE,width=9)
            d.text((xy(x,val)[0],xy(x,val)[1]-22),f'{val:+.3f}',font=font(27,True),fill=BLUE,anchor='mb')
            d.text((xy(x,yr[0])[0],1290),label,font=font(23),fill=INK,anchor='mt')
        # Association is not an IRC segment; the remaining links indicate states on the verified Gaussian route.
        for k in range(3):dashed(d,centers[k][1],centers[k+1][0],MUTED,2,8,8)
        for x,name in [(1.90,'TS'),(2.72,'product')]:
            val=vals['ORCA'][name][key];px,py=xy(x,val)
            d.ellipse((px-7,py-7,px+7,py+7),fill=ORANGE)
            d.text((px,py+23),f'ORCA {val:+.3f}',font=font(24),fill=ORANGE,anchor='mt')
    d.text((170,1450),'Blue: Gaussian stationary points    Orange: independent ORCA TS / product    Dashed links: level guide, not a continuous G path',font=font(27),fill=MUTED)
    im.save(DEST/'maleic-energy-levels.png',dpi=(300,300))


def structures_and_animations(data):
    g=data['species']['Gaussian'];target=np.array(g['TS']['positions_A']);numbers=g['TS']['atomic_numbers']
    ts=json.loads((ROOT/'analysis/maleic-ts-review-20261021/audit.json').read_text())
    mode=np.array(ts['imaginary_mode_cartesian'])
    im=canvas('Butadiene + maleic anhydride: structures and TS mode',
              'Gaussian optimized geometries | Orange arrows: unique cooperative forming-bond mode, 442.956i cm−1',size=(2400,1200))
    d=ImageDraw.Draw(im)
    for center,name,label in [((400,610),'RC','Reactant complex'),((1200,610),'TS','Cycloaddition TS'),((2000,610),'product','Cis adduct')]:
        x=fit(g[name]['positions_A'],target)
        molecule(im,x,numbers,center,scale=95,mode=mode if name=='TS' else None)
        d.text((center[0],270),label,font=font(35,True),fill=INK,anchor='mm')
        distances=[np.linalg.norm(x[j]-x[i]) for i,j in sorted(FORMING)]
        d.text((center[0],1000),f'{distances[0]:.4f} / {distances[1]:.4f} Å',font=font(32),fill=INK,anchor='mm')
    d.text((100,1110),'C: blue   O: red   H: pale gray | Both forming distances follow C1–C5 / C4–C6',font=font(28),fill=MUTED)
    im.save(DEST/'maleic-structures-mode.png',dpi=(300,300))
    frames=[]
    for phase in np.linspace(0,2*np.pi,32,endpoint=False):
        frame=canvas('Cycloaddition TS: imaginary mode','442.956i cm−1 | Exaggerated displacement, not a time trajectory',size=(1200,1000))
        molecule(frame,target+.65*np.sin(phase)*mode,numbers,(600,570),scale=110)
        frames.append(frame)
    frames[0].save(DEST/'maleic-ts-imaginary-mode.gif',save_all=True,append_images=frames[1:],duration=90,loop=0)
    frames=[]
    for point in data['Gaussian_IRC']:
        frame=canvas('Gaussian IRC: reactants to cis product',
                     f's = {point["coordinate_sqrt_amu_bohr"]:+.3f} amu^1/2 bohr | ΔE = {point["relative_energy_kcal_mol"]:+.3f} kcal/mol',size=(1200,1000))
        molecule(frame,fit(point['coordinates_A'],target),numbers,(600,575),scale=105)
        frames.append(frame)
    frames[0].save(DEST/'maleic-gaussian-irc.gif',save_all=True,append_images=frames[1:],duration=130,loop=0)


if __name__=='__main__':
    data=json.loads(DATA.read_text());irc_figure(data);ladder(data);structures_and_animations(data)
    print('Rendered 3 static figures and 2 source-backed animations.')
