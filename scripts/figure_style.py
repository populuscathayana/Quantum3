"""Shared plotting style ported from the experiment's original figures."""
import math
import os
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont
ROOT = Path(__file__).resolve().parents[1]
DEST = Path(os.environ.get('QUANTUM3_OUTPUT', ROOT / 'build/figures'))
DEST.mkdir(parents=True, exist_ok=True)
INK, MUTED, GRID = '#243747', '#566772', '#dce3e8'
BLUE, ORANGE = '#326f9a', '#c5783e'

def font(size, bold=False):
    stem = 'DejaVuSans' + ('-Bold' if bold else '') + '.ttf'
    if os.environ.get('QUANTUM3_FONT_DIR'):
        return ImageFont.truetype(str(Path(os.environ['QUANTUM3_FONT_DIR']) / stem), size)
    from matplotlib import get_data_path
    return ImageFont.truetype(str(Path(get_data_path()) / 'fonts/ttf' / stem), size)

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
