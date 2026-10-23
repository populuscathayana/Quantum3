#!/usr/bin/env python3
"""Render source-backed static figures; no wavefunction calculation is run.

All coordinates, frequencies and eigenvalues come from this experiment's fchk.
G values are retained verbatim and are not silently symmetry-corrected.
"""
import json
import math
from pathlib import Path
import re

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from audit_stage1 import ROOT, BOHR_A, fchk_fields
from audit_ts import analyze as audit_ts
from prepare_stage1 import expected_edges

DEST = ROOT / 'experiment-report/figures'
RESULT = ROOT / 'analysis/postprocessing'
EV = 27.211386245988
KCAL = 627.5094740631
INK, MUTED, GRID = '#243747', '#566772', '#dce3e8'
BLUE, ORANGE = '#326f9a', '#c5783e'


def font(size, bold=False):
    name = '<local-user>/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/libreoffice-headless/libreoffice/LibreOfficeDev.app/Contents/Resources/fonts/truetype/DejaVuSans' + ('-Bold' if bold else '') + '.ttf'
    return ImageFont.truetype(name, size)


def canvas(title, subtitle, size=(2400, 1500)):
    im = Image.new('RGB', size, 'white')
    d = ImageDraw.Draw(im)
    d.text((85, 58), title, font=font(52, True), fill=INK)
    d.text((85, 133), subtitle, font=font(30), fill=MUTED)
    return im


def save(im, name):
    im.save(DEST / name, dpi=(300, 300))


def dashed(d, a, b, fill, width=4, dash=16, gap=11):
    a, b = np.asarray(a, float), np.asarray(b, float)
    vec = b-a
    length = np.linalg.norm(vec)
    if length < 1e-9:
        return
    for start in np.arange(0, length, dash+gap):
        d.line([tuple(a+vec*start/length), tuple(a+vec*min(start+dash, length)/length)], fill=fill, width=width)


def arrow(d, a, b, fill=ORANGE, width=6, head=17):
    a, b = np.asarray(a, float), np.asarray(b, float)
    v = b-a
    length = np.linalg.norm(v)
    if length < 2:
        return
    u = v/length
    side = np.array([-u[1], u[0]])
    d.line([tuple(a), tuple(b)], fill=fill, width=width)
    d.polygon([tuple(b), tuple(b-head*u+.45*head*side), tuple(b-head*u-.45*head*side)], fill=fill)


def view_matrix(angle=43):
    theta = math.radians(angle)
    return np.array([[1,0,0],[0,math.cos(theta),math.sin(theta)],
                     [0,-math.sin(theta),math.cos(theta)]])


def molecule(im, xyz, numbers, edges, center, scale=145, forming=(), mode=None,
             carbon_labels=None, view=None, fixed_center=None):
    d = ImageDraw.Draw(im)
    view = view_matrix() if view is None else view
    xyz = np.asarray(xyz)
    x3 = xyz @ view.T
    origin = (x3.max(axis=0)+x3.min(axis=0))/2 if fixed_center is None else fixed_center
    x3 = x3-origin
    xy = x3[:,:2]*[scale,-scale] + np.asarray(center)
    for a,b in forming:
        dashed(d, xy[a-1], xy[b-1], BLUE, 5)
    # Draw sticks first, then atoms and labels; bond endpoints cannot obscure
    # the common carbon numbering used throughout the report.
    for a,b in edges:
        a,b = a-1,b-1
        d.line([tuple(xy[a]),tuple(xy[b])],fill='#89939b',width=11)
    labels = carbon_labels or {i+1:str(i+1) for i,z in enumerate(numbers) if z==6}
    for i in np.argsort(x3[:,2]):
        p=xy[i]; carbon=numbers[i]==6; r=24 if carbon else 15
        d.ellipse((p[0]-r,p[1]-r,p[0]+r,p[1]+r),fill=INK if carbon else '#f2f4f5',outline='#89939b',width=2)
        if carbon:
            d.text(tuple(p),labels.get(i+1,str(i+1)),font=font(21,True),fill='white',anchor='mm')
    for i,z in enumerate(numbers):
        if z==6:
            d.text(tuple(xy[i]),labels.get(i+1,str(i+1)),font=font(21,True),fill='white',anchor='mm',stroke_width=1,stroke_fill=INK)
    if mode is not None:
        projected = np.asarray(mode) @ view.T
        for i,z in enumerate(numbers):
            if z==6:
                delta=projected[i,:2]*[scale,-scale]*1.9
                if np.linalg.norm(delta)>8:
                    start=xy[i]+delta/np.linalg.norm(delta)*29
                    arrow(d,start,start+delta,ORANGE,7,19)
    return xy


def load_sources():
    progress=json.loads((ROOT/'cluster-progress.json').read_text())
    sources={}
    for name in ['ethylene-opt','butadiene-opt','cyclohexene-opt','ts-opt']:
        directory=ROOT/progress['jobs'][name]['local_directory']
        f=fchk_fields(directory/'calc.fchk')
        out=(directory/'calc.out').read_text()
        assert 'Normal termination of Gaussian' in out
        assert 'Standard basis: def2SVP (5D, 7F)' in out
        assert 'Grimme-D3(BJ) Dispersion energy=' in out
        thermo={}
        patterns={'zpe':'Zero-point correction=',
                  'E_zpe':'Sum of electronic and zero-point Energies=',
                  'H':'Sum of electronic and thermal Enthalpies=',
                  'G':'Sum of electronic and thermal Free Energies='}
        for key,patt in patterns.items():
            matches=re.findall(re.escape(patt)+r'\s*([-\d.]+)',out)
            assert len(matches)==1
            thermo[key]=float(matches[0])
        thermo['E']=f['Total Energy']
        thermo['rotational_symmetry_number']=int(re.findall(r'Rotational symmetry number\s+(\d+)\.',out)[-1])
        assert re.search(r'Temperature\s+298\.150 Kelvin\.\s+Pressure\s+1\.00000 Atm',out)
        modes=f['Number of Normal Modes']
        freqs=f['Vib-E2'][:modes]
        sources[name]=dict(path=str(directory.relative_to(ROOT)),fields=f,
                           xyz=np.array(f['Current cartesian coordinates']).reshape(-1,3)*BOHR_A,
                           numbers=f['Atomic numbers'],thermo=thermo,frequencies=freqs)
    return sources


def structures(sources, ts):
    im=canvas('Optimized structures', 'B3LYP-D3(BJ)/def2-SVP  |  Gas phase  |  Common molecular scale')
    d=ImageDraw.Draw(im)
    names=['ethylene-opt','butadiene-opt','ts-opt','cyclohexene-opt']
    titles=['(a) Ethylene','(b) 1,3-Butadiene','(c) Transition state','(d) Cyclohexene']
    for k,(name,title) in enumerate(zip(names,titles)):
        x=90+(k%2)*1170; y=235+(k//2)*590
        d.text((x,y),title,font=font(38,True),fill=INK)
        s=sources[name]
        edges=expected_edges('ts-constrained' if name=='ts-opt' else name)
        molecule(im,s['xyz'],s['numbers'],edges,(x+520,y+230),scale=80,
                 forming=[(1,5),(4,6)] if name=='ts-opt' else [],
                 carbon_labels={1:'5',2:'6'} if name=='ethylene-opt' else None)
        modes=len(s['frequencies']); negative=sum(v<0 for v in s['frequencies'])
        if name=='ts-opt':
            text=f'42 modes; 1 imaginary mode: {abs(ts["imaginary_frequency_cm_1"]):.1f}i cm⁻¹'
            extra='C1-C5 = 2.2712 Å; C4-C6 = 2.2712 Å'
        elif name=='butadiene-opt':
            text='24 modes; no imaginary frequencies'
            extra='C1-C2-C3-C4 torsion = 26.90° (near s-cis)'
        else:
            text=f'{modes} modes; no imaginary frequencies'
            extra='C5=C6' if name=='ethylene-opt' else 'C2=C3 retained; six-membered ring formed'
        d.text((x,y+433),text,font=font(29),fill=INK)
        d.text((x,y+482),extra,font=font(27),fill=MUTED)
    d.line((1180,230,1180,1370),fill=GRID,width=2)
    d.line((85,795,2315,795),fill=GRID,width=2)
    d.text((85,1430),'Carbon labels follow C1-C4 (diene), C5-C6 (ethylene). Dashed lines denote the two forming bonds.',font=font(27),fill=MUTED)
    save(im,'optimized-structures.png')


def mode_figures(sources, record, mode):
    s=sources['ts-opt']; xyz=s['xyz']; numbers=s['numbers']; edges=expected_edges('ts-constrained')
    if record['mode_bond_derivatives']['1-5']>0:
        mode=-mode
    im=canvas('Transition-state imaginary mode',
              f'{abs(record["imaginary_frequency_cm_1"]):.1f}i cm⁻¹  |  B3LYP-D3(BJ)/def2-SVP  |  One imaginary mode among 42',
              (2200,1450))
    d=ImageDraw.Draw(im)
    molecule(im,xyz,numbers,edges,(1080,695),scale=185,forming=[(1,5),(4,6)],mode=mode)
    d.text((90,1180),'The two forming C-C bonds shorten together in the displayed displacement direction.',font=font(34),fill=INK)
    d.text((90,1245),'Arrow lengths are enlarged; the opposite sign describes the reverse displacement.',font=font(29),fill=MUTED)
    d.text((90,1335),'Normal-mode displacement, not an IRC trajectory.',font=font(29),fill=MUTED)
    save(im,'ts-imaginary-mode.png')
    view=view_matrix(); x3=xyz@view.T; center=(x3.max(axis=0)+x3.min(axis=0))/2
    frames=[]
    for i in range(40):
        amp=.6*math.sin(2*math.pi*i/40)
        frame=canvas('TS imaginary mode',f'{abs(record["imaginary_frequency_cm_1"]):.1f}i cm⁻¹  |  Linear normal-mode displacement',(1400,1000))
        coords=xyz+amp*mode
        molecule(frame,coords,numbers,edges,(680,525),scale=128,forming=[(1,5),(4,6)],fixed_center=center)
        dr=ImageDraw.Draw(frame)
        for j,(a,b) in enumerate([(1,5),(4,6)]):
            length=np.linalg.norm(coords[a-1]-coords[b-1])
            dr.text((85+600*j,875),f'C{a}-C{b}: {length:.3f} Å',font=font(28),fill=INK)
        dr.text((85,935),'Motion is amplified for visualization; this is not the IRC path.',font=font(25),fill=MUTED)
        frames.append(frame)
    frames[0].save(DEST/'ts-imaginary-mode.gif',save_all=True,append_images=frames[1:],duration=100,loop=0,optimize=False)
    preview=Image.new('RGB',(2800,1000),'white')
    preview.paste(frames[10],(0,0));preview.paste(frames[30],(1400,0))
    preview.save(RESULT/'ts-mode-extrema-qa.png')


def energy_profile(sources):
    ref=sources['ethylene-opt']['thermo']['E']+sources['butadiene-opt']['thermo']['E']
    values=[0.,(sources['ts-opt']['thermo']['E']-ref)*KCAL,(sources['cyclohexene-opt']['thermo']['E']-ref)*KCAL]
    im=canvas('Electronic energies of the stationary structures',
              'B3LYP-D3(BJ)/def2-SVP  |  Reference: isolated near-s-cis butadiene + ethylene',(2200,1400))
    d=ImageDraw.Draw(im); left,right,top,bottom=260,2040,275,1120
    low=10*math.floor((min(values)-10)/10); high=10*math.ceil((max(values)+10)/10)
    py=lambda v:bottom-(v-low)/(high-low)*(bottom-top)
    for value in np.arange(low,high+1,10):
        y=py(value)
        d.line((left,y,right,y),fill=GRID,width=2)
        d.text((left-25,y),f'{value:.0f}',font=font(29),fill=MUTED,anchor='rm')
    d.line((left,top,left,bottom),fill=INK,width=3)
    d.text((85,218),'Relative E / kcal mol⁻¹',font=font(31),fill=INK)
    xs=[490,1130,1770]; colors=[INK,ORANGE,BLUE]
    for x,value,color in zip(xs,values,colors):
        y=py(value);d.line((x-160,y,x+160,y),fill=color,width=9)
        d.text((x,y-54),f'{value:+.2f}' if value else '0.00',font=font(38,True),fill=color,anchor='mm')
    for x,text in zip(xs,['Separated reactants','Transition state','Cyclohexene']):
        d.text((x,1190),text,font=font(31),fill=INK,anchor='mm')
    d.text((85,1295),'Electronic energy only: no zero-point or thermal correction. The three levels are not an IRC curve.',font=font(28),fill=MUTED)
    save(im,'electronic-energy-profile.png')
    return dict(reference_E_hartree=ref,relative_E_kcal_mol=dict(zip(['reactants','TS','product'],values)))


def frontier_levels(sources):
    records={}
    for name in ['butadiene-opt','ethylene-opt']:
        f=sources[name]['fields']; n=f['Number of alpha electrons']; e=f['Alpha Orbital Energies']
        records[name]={'HOMO':{'mo':n,'Eh':e[n-1],'eV':e[n-1]*EV},'LUMO':{'mo':n+1,'Eh':e[n],'eV':e[n]*EV}}
    im=canvas('Frontier-orbital energy levels',
              'Isolated optimized reactants  |  B3LYP-D3(BJ)/def2-SVP  |  Gas phase',(2200,1400))
    d=ImageDraw.Draw(im); top,bottom=285,1120; left,right=245,2050
    py=lambda value:bottom-(value+9)/11*(bottom-top)
    for v in range(-9,3):
        y=py(v); d.line((left,y,right,y),fill=GRID,width=2)
        d.text((left-25,y),str(v),font=font(28),fill=MUTED,anchor='rm')
    d.line((left,top,left,bottom),fill=INK,width=3)
    d.text((85,218),'Orbital energy / eV',font=font(31),fill=INK)
    for k,name in enumerate(['butadiene-opt','ethylene-opt']):
        x=690+k*850
        for role in ['HOMO','LUMO']:
            r=records[name][role]; y=py(r['eV'])
            focus=(name=='butadiene-opt' and role=='HOMO') or (name=='ethylene-opt' and role=='LUMO')
            color=(BLUE if role=='HOMO' else ORANGE) if focus else '#727e87'
            if role=='HOMO': d.line((x-170,y,x+170,y),fill=color,width=8)
            else: dashed(d,(x-170,y),(x+170,y),color,8,24,12)
            d.text((x,y-50),f'{role}  (MO {r["mo"]})',font=font(33,True),fill=color,anchor='mm')
            d.text((x,y+48),f'{r["eV"]:+.3f} eV',font=font(32),fill=INK,anchor='mm')
    for x,label in [(690,'Near-s-cis 1,3-butadiene'),(1540,'Ethylene')]:
        d.text((x,1185),label,font=font(34),fill=INK,anchor='mm')
    d.text((85,1290),'Solid: occupied; dashed: unoccupied. Orbital eigenvalue differences are not activation or excitation energies.',font=font(27),fill=MUTED)
    save(im,'frontier-orbital-levels.png')
    records['course_pair_gap_eV']=records['ethylene-opt']['LUMO']['eV']-records['butadiene-opt']['HOMO']['eV']
    return records


def main():
    DEST.mkdir(parents=True,exist_ok=True); RESULT.mkdir(parents=True,exist_ok=True)
    sources=load_sources(); ts,numbers,xyz,mode=audit_ts()
    assert np.max(np.abs(xyz-sources['ts-opt']['xyz']))<1e-10
    structures(sources,ts)
    mode_figures(sources,ts,mode)
    profile=energy_profile(sources)
    frontier=frontier_levels(sources)
    table={name:{'source':s['path'],'energy_hartree':s['thermo'],'frequencies_cm_1':s['frequencies'],
                 'mode_count':len(s['frequencies']),'imaginary_count':sum(v<0 for v in s['frequencies'])}
           for name,s in sources.items()}
    ref={key:sources['ethylene-opt']['thermo'][key]+sources['butadiene-opt']['thermo'][key] for key in ['E','E_zpe','H','G']}
    derived={name:{key:(sources[name]['thermo'][key]-ref[key])*KCAL for key in ref} for name in ['ts-opt','cyclohexene-opt']}
    result={'method':'B3LYP-D3(BJ)/def2-SVP','temperature_K':298.15,'pressure_atm':1,
            'sources':table,'relative_energies_kcal_mol':derived,'electronic_profile':profile,
            'frontier_levels':frontier,'G_status':'raw Gaussian values preserved; separate symmetry-corrected values and audit in thermochemistry-audit.json',
            'IRC_status':'not represented in these plots','figures':[p.name for p in sorted(DEST.glob('*'))]}
    (RESULT/'static-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'figures':result['figures'],'relative_energies_kcal_mol':derived,'frontier_levels':frontier},indent=2))


if __name__=='__main__':
    main()
