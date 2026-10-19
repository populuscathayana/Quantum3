import json
import numpy as np
from PIL import Image, ImageDraw
from figure_style import *
def expected_edges(name):
    if name == 'ethylene-opt':
        return {(1,2),(1,3),(1,4),(2,5),(2,6)}
    if name == 'butadiene-opt':
        return {(1,2),(2,3),(3,4),(1,5),(1,6),(2,7),(3,8),(4,9),(4,10)}
    bonds = {(1,2),(2,3),(3,4),(5,6),(1,7),(1,8),(2,9),(3,10),(4,11),(4,12),(5,13),(5,14),(6,15),(6,16)}
    if name == 'cyclohexene-opt': bonds |= {(1,5),(4,6)}
    return bonds
def render(record, coords):
    points = record['points']
    s = np.array([p['s'] for p in points])
    energy = np.array([p['relative_electronic_energy_kcal_mol'] for p in points])
    bonds = np.array([p['forming_bonds_A'] for p in points])
    view = view_matrix()
    projected = np.concatenate(coords)@view.T
    center = (projected.max(axis=0)+projected.min(axis=0))/2
    numbers, edges = [6]*6+[1]*10, expected_edges('ts-constrained')
    im = canvas('Cycloaddition IRC: energy and bond formation',
                'B3LYP-D3(BJ)/def2-SVP  |  54 converged path points and the transition state', (2400, 1700))
    d = ImageDraw.Draw(im)
    left, right = 260, 2250
    px = lambda a: left+(a+9)/15*(right-left)
    panels = [(300, 790, -65, 5, range(-60, 1, 10), 'Electronic energy relative to TS / kcal mol⁻¹'),
              (1040, 1430, 1.4, 3.5, [1.5, 2., 2.5, 3., 3.5], 'Forming C-C distance / Å')]
    for k, (top, bottom, low, high, ticks, title) in enumerate(panels):
        py = lambda v: bottom-(v-low)/(high-low)*(bottom-top)
        d.text((90, top-65), title, font=font(32, True), fill=INK)
        for v in ticks:
            y = py(v)
            d.line((left, y, right, y), fill=GRID, width=2)
            d.text((left-25, y), f'{v:g}', font=font(27), fill=MUTED, anchor='rm')
        dashed(d, (px(0), top), (px(0), bottom), GRID, 3, 10, 10)
        d.line((left, top, left, bottom, right, bottom), fill=INK, width=3)
        for v in [-9, -6, -3, 0, 3, 6]:
            z = px(v)
            d.line((z, bottom, z, bottom+10), fill=INK, width=2)
            d.text((z, bottom+20), str(v), font=font(27), fill=MUTED, anchor='mt')
        if k == 0:
            plot = [(px(a), py(b)) for a, b in zip(s, energy)]
            d.line(plot, fill=BLUE, width=6)
            for x, y in plot: d.ellipse((x-5, y-5, x+5, y+5), fill=BLUE)
            d.text((px(0)+20, py(0)+15), 'TS', font=font(32, True), fill=INK)
            d.text((px(s[0])+20, py(energy[0])+28), f'{energy[0]:.2f}', font=font(30, True), fill=BLUE)
            d.text((px(s[-1])-20, py(energy[-1])-30), f'{energy[-1]:.2f}', font=font(30, True), fill=BLUE, anchor='rm')
            d.text((left+10, 900), 'Reactant side', font=font(29, True), fill=MUTED)
            d.text((right, 900), 'Cyclohexene B side', font=font(29, True), fill=MUTED, anchor='ra')
        else:
            for col, color in [(0, BLUE), (1, ORANGE)]:
                plot = [(px(a), py(b)) for a, b in zip(s, bonds[:, col])]
                if col == 0: d.line(plot, fill=color, width=6)
                else:
                    for a, b in zip(plot, plot[1:]): dashed(d, a, b, color, 6, 14, 10)
            d.line((1530, 975, 1600, 975), fill=BLUE, width=6)
            d.text((1620, 975), 'C1-C5', font=font(28), fill=INK, anchor='lm')
            dashed(d, (1870, 975), (1940, 975), ORANGE, 6)
            d.text((1960, 975), 'C4-C6', font=font(28), fill=INK, anchor='lm')
    d.text((1250, 1540), 'Oriented IRC coordinate, s / (amu¹ᐟ² bohr)', font=font(32), fill=INK, anchor='mm')
    d.text((90, 1630), 'Negative: reactant side. Positive: product side. Endpoint optimizations verified separately.', font=font(27), fill=MUTED)
    save(im, 'irc-full-profile.png')
    frames = []
    for i, (p, xyz) in enumerate(zip(points, coords)):
        im = canvas('Cycloaddition IRC: reactant side to product side',
                    'B3LYP-D3(BJ)/def2-SVP  |  Actual computed structures', (1800, 1350))
        d = ImageDraw.Draw(im)
        label = 'Transition state' if p['branch'] == 'TS' else ('Reactant-side branch' if p['branch'] == 'reactant' else 'Product-side branch')
        d.text((90, 210), f'{label}    |    Frame {i+1:02d} / 55', font=font(31, True), fill=INK)
        d.text((90, 265), f's = {p["s"]:+.5f} amu¹ᐟ² bohr     ΔE = {energy[i]:+.2f} kcal mol⁻¹', font=font(29), fill=INK)
        molecule(im, xyz, numbers, edges, (880, 710), scale=150, forming=[(1, 5), (4, 6)], fixed_center=center)
        d.text((90, 1140), f'C1-C5: {bonds[i, 0]:.4f} Å       C4-C6: {bonds[i, 1]:.4f} Å', font=font(32, True), fill=INK)
        d.text((90, 1210), 'Rigid alignment only; no interpolated structures. Playback is not a time trajectory.', font=font(28), fill=MUTED)
        frames.append(im)
    durations = [200]*55
    for i in [0, 34, 54]: durations[i] = 1500
    frames[0].save(DEST/'irc-full-path.gif', save_all=True,
                   append_images=frames[1:], duration=durations, loop=0, optimize=False)
    im = canvas('Structures along the cycloaddition IRC',
                'Common scale and carbon numbering  |  Electronic energies relative to the TS', (2400, 1120))
    d = ImageDraw.Draw(im)
    for k, (i, label) in enumerate([(0, 'Reactant-side IRC end'), (34, 'Transition state'), (54, 'Product-side IRC end')]):
        left = 65+k*790
        d.text((left, 245), label, font=font(33, True), fill=INK)
        molecule(im, coords[i], numbers, edges, (left+365, 590), scale=66,
                 forming=[(1, 5), (4, 6)], fixed_center=center)
        d.text((left, 835), f's = {s[i]:+.3f}', font=font(28), fill=INK)
        d.text((left, 890), f'ΔE = {energy[i]:+.2f} kcal mol⁻¹', font=font(28), fill=INK)
        d.text((left, 945), f'C-C: {bonds[i, 0]:.3f} / {bonds[i, 1]:.3f} Å', font=font(26), fill=MUTED)
        if k: d.line((left-25, 230, left-25, 1010), fill=GRID, width=2)
    d.text((65, 1055), 's in amu¹ᐟ² bohr. Dashed links mark C1-C5 and C4-C6. End structures shown are IRC points.', font=font(27), fill=MUTED)
    save(im, 'irc-full-structures.png')
    qa = Image.new('RGB', (5400, 1350), 'white')
    for k, i in enumerate([0, 34, 54]): qa.paste(frames[i], (k*1800, 0))
    qa.save(DEST/'full-irc-animation-qa.png')
    (DEST/'full-irc-data.json').write_text(json.dumps(record, indent=2)+'\n')

if __name__ == '__main__':
    record = json.loads((ROOT/'data/irc.json').read_text())
    target = np.array(next(p['coordinates_A'] for p in record['points'] if p['branch']=='TS'))
    aligned = []
    for point in record['points']:
        x = np.array(point['coordinates_A'])
        u, _, vt = np.linalg.svd((x-x.mean(0)).T @ (target-target.mean(0)))
        rot = u @ np.diag([1, 1, np.linalg.det(u@vt)]) @ vt
        aligned.append((x-x.mean(0)) @ rot+target.mean(0))
    render(record, np.array(aligned))
