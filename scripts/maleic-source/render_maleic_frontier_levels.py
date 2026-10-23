"""Four frontier levels from existing optimized ORCA outputs; no new SCF."""
import json,re,hashlib
from pathlib import Path
from PIL import Image,ImageDraw
from postprocess_static import font,INK,BLUE,ORANGE,dashed

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    exported=json.loads((ROOT/'analysis/maleic-orbitals-20261021/export.json').read_text())
    records={}
    for kind,electrons in [('butadiene',30),('maleic-anhydride',50)]:
        source=ROOT/exported[kind]['source_directory'];out=(source/'calc.out').read_text()
        assert 'ORCA TERMINATED NORMALLY' in out[-3000:]
        assert sha(source/'calc.gbw')==exported[kind]['source_gbw_sha256']
        table=out.rsplit('ORBITAL ENERGIES',1)[1]
        rows={int(i):dict(occupation=float(o),energy_hartree=float(e),energy_eV=float(ev))
            for i,o,e,ev in re.findall(r'^\s*(\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$',table,re.M)}
        occupied=[i for i,v in rows.items() if v['occupation']>0]
        assert sum(v['occupation'] for v in rows.values())==electrons
        homo=max(occupied);lumo=min(i for i,v in rows.items() if v['occupation']==0)
        assert lumo==homo+1
        records[kind]={'HOMO':dict(rows[homo],index_zero_based=homo),'LUMO':dict(rows[lumo],index_zero_based=lumo),
                       'source_output':str((source/'calc.out').relative_to(ROOT)),'source_sha256':sha(source/'calc.out')}
    forward=records['maleic-anhydride']['LUMO']['energy_eV']-records['butadiene']['HOMO']['energy_eV']
    reverse=records['butadiene']['LUMO']['energy_eV']-records['maleic-anhydride']['HOMO']['energy_eV']
    assert abs(forward-3.1732)<1e-10 and abs(reverse-7.3212)<1e-10
    im=Image.new('RGB',(1800,1200),'white');d=ImageDraw.Draw(im)
    def y(e):return 900-(e+9)*80
    for e in range(-9,1):
        yy=y(e);d.line((180,yy,1650,yy),fill='#e5e9ed',width=2)
        d.text((145,yy),str(e),font=font(37),fill=INK,anchor='rm')
    d.line((180,y(0),180,y(-9)),fill=INK,width=3)
    d.text((180,80),'ε / eV',font=font(37,True),fill=INK)
    for cx,kind,title in [(600,'butadiene','1,3-Butadiene'),(1320,'maleic-anhydride','Maleic anhydride')]:
        for orbital in ['HOMO','LUMO']:
            value=records[kind][orbital]['energy_eV'];yy=y(value)
            color=BLUE if kind=='butadiene' and orbital=='HOMO' else ORANGE if kind=='maleic-anhydride' and orbital=='LUMO' else '#78858d'
            if orbital=='HOMO':d.line((cx-150,yy,cx+150,yy),fill=color,width=6)
            else:dashed(d,(cx-150,yy),(cx+150,yy),color,5,18,10)
            d.text((cx,yy-37),orbital,font=font(37,True),fill=color,anchor='mm')
            d.text((cx,yy+39),f'{value:.2f}',font=font(37),fill=INK,anchor='mm')
        d.text((cx,960),title,font=font(37,True),fill=INK,anchor='mm')
    d.text((900,1060),f'HOMO(diene) → LUMO(anhydride)     Δε = {forward:.2f} eV',font=font(37),fill=INK,anchor='mm')
    d.text((900,1130),f'HOMO(anhydride) → LUMO(diene)     Δε = {reverse:.2f} eV',font=font(37),fill=INK,anchor='mm')
    dest=ROOT/'experiment-report/figures/maleic-frontier-levels.png';im.save(dest,dpi=(300,300))
    data={'orbitals':records,'forward_gap_eV':forward,'reverse_gap_eV':reverse,
          'shared_by_approaches':['endo','exo'],'display_decimals':2,'figure':str(dest.relative_to(ROOT)),
          'figure_sha256':sha(dest),'note':'Isolated optimized ORCA fragment orbital energies, not activation barriers.'}
    (ROOT/'analysis/maleic-orbitals-20261021/levels.json').write_text(json.dumps(data,indent=2))
    print(json.dumps(data,indent=2))

if __name__=='__main__':main()
