"""Physical-size typography for the report's existing scientific renderers."""
import json
import re
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
QA=ROOT/'analysis/figure-typography-20261022'
QA.mkdir(exist_ok=True)
WIDTH_PT=146.4/25.4*72.27
RECORDS={}
ORIGINAL_TEXT=ImageDraw.ImageDraw.text
ORIGINAL_WRITE=Path.write_text
SHORT={
 'Separated reactants':'R', 'Separated R':'R', 'Encounter complex':'RC',
 'Weak complex':'RC', 'Cycloaddition TS':'TS',
 'Cyclohexene A':'A', 'Cyclohexene B':'B',
 'Reactant-side IRC end':'Reactant-side end',
 'Product-side IRC end':'Product-side end',
 'Electronic energy relative to TS / kcal mol⁻¹':'ΔE / kcal mol⁻¹',
 'Forming C-C distance / Å':'C–C / Å',
 'Oriented IRC coordinate, s / (amu¹/² bohr)':'s / (amu¹/² bohr)',
 'Mass-weighted IRC coordinate / (amu^1/2 bohr)':'s / (amu¹/² bohr)',
 'Cumulative Cartesian path length / Å':'Path length / Å',
 'Peak above weak complex / kcal mol⁻¹':'ΔE‡ / kcal mol⁻¹',
 'NEB sampled maximum':'NEB', 'CI-NEB climbing image':'CI-NEB',
 'Gaussian optimized TS':'Gaussian TS',
 'ΔE / kcal mol⁻¹; own-program B = 0':'ΔE / kcal mol⁻¹',
 'ORCA NEB (step 35)':'ORCA NEB', 'ORCA CI-NEB (step 50)':'ORCA CI-NEB',
 'Near-s-cis 1,3-butadiene':'1,3-Butadiene',
 'Endo transition state':'Endo TS', 'Reactant complex':'RC',
 'Electronic energy / (kcal mol−1)':'ΔE / kcal mol⁻¹',
 'Gibbs free energy / (kcal mol−1)':'ΔG / kcal mol⁻¹',
 'Reaction states':'', 'Diene + MA':'R', 'Complex':'RC', 'Cis product':'P',
 'ORCA refined TS':'ORCA TS', 'Gaussian TS3346':'Gaussian TS',
 'Saved result':'Method', 'Conformer TS':'TS(A/B)',
 'Conformer comparison':'A',
}

def record(draw,xy,text,*args,**kwargs):
    title=draw._image.info.get('report_title','')
    if not title:return ORIGINAL_TEXT(draw,xy,text,*args,**kwargs)
    assert not args
    kw=dict(kwargs);f=kw.pop('font')
    RECORDS.setdefault(title,[]).append(dict(xy=list(xy),text=text,
        font_path=str(f.path),font_size=f.size,kwargs=kw))

def wrap(text,font,width):
    lines=[]
    for para in text.split('\n'):
        line=''
        for word in para.split():
            test=(line+' '+word).strip()
            if line and font.getlength(test)>width:lines.append(line);line=word
            else:line=test
        lines.append(line)
    return '\n'.join(lines)

def concise_numbers(text,name):
    """Display-only rounding; original rows and plotted coordinates stay intact."""
    def replace(match):
        raw=match.group();value=float(raw);a=abs(value)
        if 'Å' in text and 'eV' not in text: places=3
        elif 'cm' in text or match.end()<len(text) and text[match.end()]=='i': places=1
        else: places=2
        if 0<a<.1 and places==2:
            return format(value,'.2g')
        return f'{value:.{places}f}'
    return re.sub(r'(?<![\w.])-?\d+\.\d+(?![\d.])',replace,text)

def render(im,rows,name):
    # Preserve the full canvas width: the font-to-point conversion is exact.
    im=im.copy();d=ImageDraw.Draw(im);W,H=im.size
    factor={'ts-imaginary-mode.png':.91,'frontier-phase-matching.png':.84}.get(name,1)
    target=8.5;placed=[];audit=[]
    for row in sorted(rows,key=lambda r:(r['xy'][1],r['xy'][0])):
        text=SHORT.get(row['text'],row['text']);x,y=row['xy'];kw=dict(row['kwargs'])
        text=concise_numbers(text,name)
        if not text:continue
        if text.startswith('Oriented IRC coordinate'):text='s / (amu¹/² bohr)'
        if name=='irc-full-structures.png' and text.startswith('ΔE ='):
            text=text.replace(' kcal mol⁻¹','')
        if name=='frontier-orbital-levels.png':
            text=re.sub(r'  \(MO (\d+)\)',r' (\1)',text)
        if name=='neb-cineb-irc-comparison.png':
            if row['text']=='Last IRC point':text='IRC end';y=780;x=1720
            if row['xy'][1]==1050:y=1080
            if row['xy'][1] in (1110,1165,1220):
                y={1110:1170,1165:1250,1220:1330}[row['xy'][1]]
        if name=='maleic-energy-levels.png' and text.startswith('ORCA '):
            text=text.replace('ORCA ','')
        if name=='electronic-energy-ladder-complete.png' and row['xy'][1]==1030:
            text='' # Job identifiers are already attached to the source tables.
        if not text:continue
        atom=bool(re.fullmatch(r'(?:[COH])?\d{1,2}',text)) and kw.get('anchor')=='mm' and (kw.get('fill')=='white' or kw.get('stroke_fill'))
        atom_origin=(x,y)
        if atom and name=='maleic-structures-multiview.png':
            if y>900:
                offsets={'O9':(-75,0),'O10':(85,-35),'O11':(-95,-15),
                         'C7':(-80,30),'C8':(85,15),'C5':(-65,35),'C6':(65,-15),
                         'C1':(75,25),'C4':(75,10),'C2':(-70,5),'C3':(70,-30)}
            else:offsets={'C1':(-35,28),'C5':(-35,-15)}
            dx,dy=offsets.get(text,(0,0));x+=dx;y+=dy
        pt=8.0 if atom else target
        size=round(pt*W/(WIDTH_PT*factor))
        font=ImageFont.truetype(row['font_path'],size)
        anchor=kw.pop('anchor','la');kw.pop('stroke_width',None);kw.pop('stroke_fill',None)
        # Titles and categorical labels wrap within their own column, not shrink.
        peers=[r['xy'][0] for r in rows if abs(r['xy'][1]-y)<18 and r['xy'][0]>x+50]
        leftpeers=[r['xy'][0] for r in rows if abs(r['xy'][1]-y)<18 and r['xy'][0]<x-50]
        maxw=W-80
        if anchor[0]=='l':maxw=min(maxw,W-x-35,(min(peers)-x-25) if peers else W)
        elif anchor[0]=='m':
            maxw=min(maxw,2*min(x-30,W-x-30))
            if peers:maxw=min(maxw,min(peers)-x-25)
            if leftpeers:maxw=min(maxw,x-max(leftpeers)-25)
        else:maxw=min(maxw,x-30)
        if not atom:text=wrap(text,font,max(250,maxw))
        # Use glyph bounds to keep the original anchor and avoid edge clipping.
        box=d.multiline_textbbox((0,0),text,font=font,spacing=8,align='center' if anchor[0]=='m' else 'left')
        w,h=box[2]-box[0],box[3]-box[1]
        l=x-({'l':0,'m':w/2,'r':w}.get(anchor[0],0))
        t=y-({'m':h/2,'b':h,'s':h,'t':0,'a':0}.get(anchor[-1],0))
        l=max(25,min(l,W-w-25));t=max(15,t)
        if not atom:
            # Deterministic label-only collision avoidance; data marks never move.
            for _ in range(20):
                hits=[b for b in placed if l<b[2]+10 and l+w>b[0]-10 and t<b[3]+8 and t+h>b[1]-8]
                if not hits:break
                t=max(b[3] for b in hits)+10
            placed.append((l,t,l+w,t+h))
        audit.append(dict(text=text,point_size=pt,bbox=[l,t,l+w,t+h],atom=atom))
        if t+h>im.height-20:
            expanded=Image.new('RGB',(W,int(t+h+45)),'white');expanded.paste(im,(0,0));im=expanded;d=ImageDraw.Draw(im)
        kw.pop('align',None)
        if atom:
            if atom_origin!=(x,y):d.line((atom_origin,(x,y)),fill='#657783',width=2)
            d.rounded_rectangle((l-4,t-3,l+w+4,t+h+3),radius=5,fill='white')
            kw.update(fill='#243747')
        with patch.object(ImageDraw.ImageDraw,'text',ORIGINAL_TEXT):
            ORIGINAL_TEXT(d,(l-box[0],t-box[1]),text,font=font,spacing=8,align='center' if anchor[0]=='m' else 'left',**kw)
    # Trim vertical white margins only, never the calibrated width.
    from PIL import ImageChops
    box=ImageChops.difference(im,Image.new('RGB',im.size,'white')).getbbox()
    im=im.crop((0,max(0,box[1]-35),W,min(im.height,box[3]+35)))
    return im,audit

def save_calibrated(im,path,original_save):
    name=Path(path).name;rows=RECORDS[im.info['report_title']]
    original_save(im,QA/(name+'.geometry.png'))
    ORIGINAL_WRITE(QA/(name+'.texts.json'),json.dumps(rows,ensure_ascii=False,indent=2))
    out,audit=render(im,rows,name);original_save(out,path,dpi=(300,300))
    ORIGINAL_WRITE(QA/(name+'.audit.json'),json.dumps(audit,ensure_ascii=False,indent=2))

if __name__=='__main__':
    for p in QA.glob('*.texts.json'):
        name=p.name.removesuffix('.texts.json')
        im=Image.open(QA/(name+'.geometry.png')).convert('RGB')
        out,audit=render(im,json.loads(p.read_text()),name)
        out.save(ROOT/'experiment-report/figures'/name,dpi=(300,300))
        (QA/(name+'.audit.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2))
