import json
from PIL import ImageDraw
from figure_style import *
FIG = DEST / 'electronic-energy-ladder-complete.png'

def draw(data):
    im = canvas('Complete electronic-energy ladder',
        '1,3-Butadiene + ethylene → cyclohexene  |  B3LYP-D3(BJ)/def2-SVP  |  Electronic energy, not G',
        size=(2400,1900))
    d = ImageDraw.Draw(im)
    d.text((90,220),'A  Gaussian stationary levels: separated reactants = 0',font=font(36,True),fill=INK)
    l,r,t,b = 210,2310,330,950
    py = lambda e:b-(e+55)/75*(b-t)
    for e in range(-50,21,10):
        d.line((l,py(e),r,py(e)),fill=GRID,width=2)
        d.text((l-25,py(e)),str(e),font=font(30),fill=INK,anchor='rm')
    d.text((90,285),'ΔE / kcal mol⁻¹',font=font(30),fill=MUTED)
    d.line((l,t,l,b,r,b),fill=INK,width=3)
    xs = [350,700,1050,1400,1750,2100]; levels = data['Gaussian_global_states']
    # IRC-supported connections only. No joined curve through the conformer levels.
    for ia,ib in [(1,2),(2,3)]:
        d.line((xs[ia]+95,py(levels[ia]['E_kcal_mol']),xs[ib]-95,py(levels[ib]['E_kcal_mol'])),fill=BLUE,width=4)
    for x, row in zip(xs,levels):
        color = ORANGE if row['key']=='TS_conf_candidate' else BLUE if row['key'] in ('RC','TS_DA','B') else INK
        y=py(row['E_kcal_mol']);d.line((x-95,y,x+95,y),fill=color,width=7)
        d.text((x,y-36),f'{row["E_kcal_mol"]:+.4f}',font=font(33,True),fill=color,anchor='mb')
        d.text((x,990),row['label'],font=font(29,True),fill=color,anchor='mm')
        d.text((x,1030),row['id'],font=font(27),fill=MUTED,anchor='mm')
    d.text((1480,430),'A/B: CI-NEB supports a candidate path;',font=font(30),fill=INK)
    d.text((1480,480),'two-sided IRC remains unconfirmed.',font=font(30),fill=ORANGE)
    d.text((1480,540),'The two high product-side levels are',font=font(28),fill=MUTED)
    d.text((1480,582),'nearly degenerate; see the zoom below.',font=font(28),fill=MUTED)
    d.text((90,1100),'Solid guides: IRC-supported RC–TS–B connection. Unmapped association and unconfirmed A/B links are not drawn.',font=font(28),fill=MUTED)
    d.line((90,1160,2310,1160),fill=GRID,width=2)
    d.text((90,1200),'B  Conformer peak heights / kcal mol⁻¹',font=font(35,True),fill=INK)
    d.text((90,1250),'Each program uses its own optimized A and B references.',font=font(27),fill=MUTED)
    d.text((90,1330),'Saved result',font=font(30,True),fill=INK)
    d.text((850,1330),'Above A',font=font(30,True),fill=INK,anchor='ra')
    d.text((1150,1330),'Above B',font=font(30,True),fill=INK,anchor='ra')
    d.line((90,1380,1180,1380),fill=INK,width=2)
    for i,row in enumerate(data['conformer_rows']):
        y=1430+i*65
        label=row['label']+(f' (step {row["stage_step"]})' if 'stage_step' in row else '')
        d.text((90,y),label,font=font(29,i==2),fill=ORANGE if i==2 else INK)
        for key,x in [('A',850),('B',1150)]:d.text((x,y),f'{row[key]:.4f}',font=font(31,i==2),fill=ORANGE if i==2 else INK,anchor='ra')
    d.text((90,1720),'ORCA refined TS: one ring-twist mode, 67.56i cm⁻¹.',font=font(29),fill=ORANGE)
    # A quantitative, explicitly separate zoom: no visual claim of error bars.
    d.text((1370,1200),'C  Shallow B-side height (zoom)',font=font(34,True),fill=INK)
    d.text((1370,1250),'ΔE / kcal mol⁻¹; own-program B = 0',font=font(27),fill=MUTED)
    zx0,zx1,zt,zb=1490,2310,1320,1650
    zy=lambda e:zb-e/.04*(zb-zt)
    for e in [0,.01,.02,.03,.04]:
        d.line((zx0,zy(e),zx1,zy(e)),fill=GRID,width=2)
        d.text((zx0-22,zy(e)),f'{e:.2f}',font=font(27),fill=INK,anchor='rm')
    d.line((zx0,zt,zx0,zb,zx1,zb),fill=INK,width=3)
    for x,row,label,col in zip([1620,1890,2170],data['conformer_rows'][1:],
                              ['CI-NEB','Refined TS','Gaussian TS'],[BLUE,ORANGE,INK]):
        y=zy(row['B']);d.line((x-80,y,x+80,y),fill=col,width=7)
        d.text((x,y-25),f'{row["B"]:.4f}',font=font(30,True),fill=col,anchor='mb')
        d.text((x,1690),label,font=font(26),fill=col,anchor='mm')
    d.text((1370,1750),'Differences at this scale are not resolved rate barriers.',font=font(26),fill=MUTED)
    d.text((90,1810),'Discrete energy levels, not a continuous reaction coordinate. No ZPE/thermal correction; no mixed-program absolute energies.',font=font(28),fill=MUTED)
    im.save(FIG,dpi=(300,300))
    return FIG

if __name__ == '__main__':
    draw(json.loads((ROOT/'analysis/energy-ladder.json').read_text()))
