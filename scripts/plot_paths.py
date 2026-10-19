#!/usr/bin/env python3
"""Render only measured path points; preserve ordering and separate energy zeros."""
import json
import sys
from PIL import Image, ImageDraw
from figure_style import ROOT, DEST, canvas, font, dashed, INK, MUTED, GRID, BLUE, ORANGE

SOURCE=ROOT/'analysis/path-comparison.json'
FIG=DEST/'neb-cineb-irc-comparison.png'
GREY='#4c5365'


def marker(d,x,y,kind,color,r=7):
    if kind=='circle': d.ellipse((x-r,y-r,x+r,y+r),fill='white',outline=color,width=3)
    elif kind=='square': d.rectangle((x-r,y-r,x+r,y+r),fill=color,outline='white',width=1)
    else: d.polygon([(x,y-r),(x+r,y),(x,y+r),(x-r,y)],fill='white',outline=color,width=3)


def main():
    assert '--replace-render' in sys.argv or not FIG.exists(),'Preserve existing figure unless explicitly re-rendering this generated asset'
    a=json.loads(SOURCE.read_text());assert a['status']=='passed_with_stated_limits'
    im=canvas('NEB, CI-NEB and the original Gaussian IRC',
        'B3LYP-D3(BJ)/def2-SVP  |  Electronic energy relative to each program’s weak reactant complex',size=(2400,1700))
    d=ImageDraw.Draw(im)
    series=[('Gaussian IRC + TS',a['Gaussian_IRC']['rows'],GREY,'circle'),
            ('ORCA NEB (52)',a['NEB']['rows'],BLUE,'diamond'),
            ('ORCA CI-NEB (400)',a['CI_NEB']['rows'],ORANGE,'square')]
    for x,(name,_,color,kind) in zip([120,830,1540],series):
        if kind=='diamond':dashed(d,(x,235),(x+100,235),color,5,14,9)
        else:d.line((x,235,x+100,235),fill=color,width=5)
        marker(d,x+50,235,kind,color,9)
        d.text((x+120,235),name,font=font(30,True),fill=INK,anchor='lm')

    def panel(box,xlimits,ylimits,xticks,yticks):
        l,t,r,b=box;w,h=r-l,b-t;hi,lo=xlimits;ylo,yhi=ylimits
        xp=lambda q:(hi-q)/(hi-lo)*w
        yp=lambda e:h-(e-ylo)/(yhi-ylo)*h
        layer=Image.new('RGB',(w+1,h+1),'white');ld=ImageDraw.Draw(layer)
        for v in yticks:ld.line((0,yp(v),w,yp(v)),fill=GRID,width=2)
        for _,data,color,kind in series:
            pts=[(xp(p['q_mean_forming_A']),yp(p['relative_complex_E_kcal_mol'])) for p in data]
            for p,q in zip(pts,pts[1:]):
                if kind=='diamond':dashed(ld,p,q,color,4,13,9)
                else:ld.line((p,q),fill=color,width=4)
            for x,y in pts:marker(ld,x,y,kind,color)
        # Independent Gaussian optimized endpoints: crosses, never appended to IRC.
        for p,e in [(a['NEB']['rows'][0],0),
                    (a['NEB']['rows'][-1],(a['Gaussian_TS']['product_energy_hartree']-a['Gaussian_TS']['complex_energy_hartree'])*a['definitions']['kcal_per_hartree'])]:
            x,y=xp(p['q_mean_forming_A']),yp(e)
            ld.line((x-13,y-13,x+13,y+13),fill=GREY,width=4);ld.line((x-13,y+13,x+13,y-13),fill=GREY,width=4)
        im.paste(layer,(l,t));d=ImageDraw.Draw(im)
        d.line((l,t,l,b,r,b),fill=INK,width=3)
        for v in yticks:d.text((l-22,t+yp(v)),f'{v:g}',font=font(27),fill=INK,anchor='rm')
        for v in xticks:
            x=l+xp(v);d.line((x,b,x,b+12),fill=INK,width=2)
            d.text((x,b+22),f'{v:g}',font=font(27),fill=INK,anchor='mt')
        return lambda q:l+xp(q),lambda e:t+yp(e)

    d.text((110,305),'A  Full electronic-energy profile',font=font(34,True),fill=INK)
    d.text((1450,305),'B  Reactant-side enlargement',font=font(34,True),fill=INK)
    d.text((220,370),'ΔE / kcal mol⁻¹',font=font(29),fill=MUTED)
    d.text((1570,370),'ΔE / kcal mol⁻¹',font=font(29),fill=MUTED)
    px,py=panel((220,430,1320,1170),(3.85,1.45),(-50,20),[3.5,3,2.5,2,1.5],[-50,-40,-30,-20,-10,0,10,20])
    rx,ry=panel((1570,430,2260,890),(3.85,2.85),(-.25,4),[3.8,3.5,3.2,2.9],[0,1,2,3,4])
    d=ImageDraw.Draw(im)
    d.text((px(2.27),py(15.74)-35),'TS region',font=font(28,True),fill=INK,anchor='mb')
    d.text((px(1.80)-25,py(-43)),'Cyclohexene B',font=font(27),fill=INK,anchor='rm')
    d.text((px(3.64)+12,py(0)+24),'Weak complex',font=font(27),fill=INK)
    d.text((rx(3.3485)+16,ry(.50256)+35),'Last IRC point',font=font(26),fill=GREY)
    d.text((770,1250),'q = mean of the two forming C–C distances / Å',font=font(29),fill=INK,anchor='mt')
    d.text((1915,970),'q / Å  (reactant-side projection)',font=font(28),fill=INK,anchor='mt')
    d.text((1450,1050),'Peak above weak complex / kcal mol⁻¹',font=font(29,True),fill=INK)
    for y,key,label,color in [(1110,'NEB','NEB sampled maximum',BLUE),(1165,'CI_NEB','CI-NEB climbing image',ORANGE),(1220,'Gaussian_TS','Gaussian optimized TS',GREY)]:
        d.text((1450,y),label,font=font(28),fill=color)
        d.text((2250,y),f'{a[key]["forward_complex_barrier_kcal_mol"]:.4f}',font=font(29,True),fill=color,anchor='rt')
    d.text((1450,1280),'Extra digits compare outputs, not physical accuracy.',font=font(25),fill=MUTED)
    notes=[
        'q is a geometric projection, not an IRC coordinate or NEB path length. Original point order is retained.',
        'Lines connect calculated points only. Crosses mark separately optimized Gaussian endpoints, not extra IRC points.',
        'The IRC reuses the original TS; it is not a fourth independent saddle search. No Gibbs energies are plotted.',
        'Near-equal peak energies do not certify NEB stationarity or pointwise agreement of the complete paths.']
    for y,line in zip([1410,1465,1520,1575],notes):d.text((110,y),line,font=font(27),fill=MUTED)
    im.save(FIG,dpi=(300,300));print(FIG)


if __name__=='__main__':main()
