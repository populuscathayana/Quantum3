"""Measured Gaussian endo/exo TS geometries, common proper front/side views."""
import json
import numpy as np
from PIL import Image, ImageDraw
from render_maleic_results import ROOT, DEST, fit, molecule
from render_maleic_multiview import reference_frame, rx, ry
from postprocess_static import font, INK

def main():
    endo=json.loads((ROOT/'analysis/maleic-ts-review-20261021/audit.json').read_text())
    exo=json.loads((ROOT/'analysis/maleic-endo-exo-20261022/exo-gaussian-independent-audit.json').read_text())
    target=np.array(endo['positions_A']); basis=reference_frame(target)
    image=Image.new('RGB',(1800,1400),'white');d=ImageDraw.Draw(image)
    # 37 pixels at 1800 px / 416.6 TeX pt gives 8.56 pt in the report.
    for cx,record,label in [(470,endo,'endo'),(1320,exo,'exo')]:
        original=np.array(record['positions_A'])
        u,_,vt=np.linalg.svd((original[:4]-original[:4].mean(0)).T@(target[:4]-target[:4].mean(0)))
        rotation=u@np.diag([1,1,np.linalg.det(u@vt)])@vt
        assert np.linalg.det(rotation)>.999999
        x=((original-original[:4].mean(0))@rotation+target[:4].mean(0))@basis
        assert np.max(abs(np.linalg.norm(x[:,None]-x[None,:],axis=-1)-np.linalg.norm(original[:,None]-original[None,:],axis=-1)))<1e-10
        d.text((cx,65),label,font=font(37,True),fill=INK,anchor='mm')
        for cy,camera in zip([395,1000],[rx(35)@ry(12),rx(12)@ry(85)]):
            molecule(image,x,record['atomic_numbers'],(cx,cy),scale=105,camera=camera,labels=False)
        distances=[np.linalg.norm(original[j]-original[i]) for i,j in [(0,4),(3,5)]]
        d.text((cx,1335),f'{distances[0]:.3f} / {distances[1]:.3f} Å',font=font(37),fill=INK,anchor='mm')
    image.save(DEST/'maleic-endo-exo-ts.png',dpi=(300,300))

if __name__=='__main__':main()
