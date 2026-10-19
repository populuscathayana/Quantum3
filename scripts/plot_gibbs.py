#!/usr/bin/env python3
"""Discrete, symmetry-corrected RRHO free-energy levels, not IRC free energy."""
import json
from PIL import ImageDraw
from figure_style import ROOT, DEST
from figure_style import canvas, font, save, dashed, INK, MUTED, GRID, BLUE, ORANGE


def main():
    data = json.loads((ROOT/'analysis/thermochemistry.json').read_text())
    relative = data['relative_to_separated_reactants_kcal_mol']
    levels = [dict(label='Separated reactants', job=None, x=440, G=0., G_raw=0.),
              dict(label='Encounter complex', job='reactant-endpoint-opt', x=880),
              dict(label='Cycloaddition TS', job='ts-opt', x=1320),
              dict(label='Cyclohexene B', job='product-endpoint-opt', x=1760),
              dict(label='Cyclohexene A', job='cyclohexene-opt', x=2240)]
    for level in levels[1:]:
        level['G'] = relative[level['job']]['G_symmetry_corrected']
        level['G_raw'] = relative[level['job']]['G']
    im = canvas('Gibbs free-energy levels of the cycloaddition',
                'B3LYP-D3(BJ)/def2-SVP  |  Gas phase, 298.15 K, 1 atm  |  RRHO with rotational symmetry correction', (2600, 1550))
    d = ImageDraw.Draw(im)
    top, bottom, left, right = 320, 1170, 240, 2460
    py = lambda g: bottom-(g+40)/75*(bottom-top)
    d.text((90, 235), 'Relative G / kcal mol⁻¹', font=font(34, True), fill=INK)
    for v in range(-40, 31, 10):
        y = py(v)
        d.line((left, y, right, y), fill=GRID, width=2)
        d.text((left-25, y), str(v), font=font(29), fill=MUTED, anchor='rm')
    d.line((left, top, left, bottom, right, bottom), fill=INK, width=3)
    dashed(d, (2000, top-15), (2000, bottom+130), GRID, 3)
    d.text((2240, 285), 'Conformer comparison', font=font(28, True), fill=MUTED, anchor='mm')
    for a, b in zip(levels[:3], levels[1:4]):
        dashed(d, (a['x']+105, py(a['G'])), (b['x']-105, py(b['G'])), BLUE, 4, 14, 11)
    for i, level in enumerate(levels):
        x, y = level['x'], py(level['G'])
        color = ORANGE if i == 4 else BLUE
        d.line((x-105, y, x+105, y), fill=color, width=9)
        d.text((x, y-35), f'{level["G"]:+.2f}', font=font(38, True), fill=color, anchor='mb')
        d.text((x, 1220), level['label'], font=font(28, True), fill=INK, anchor='mm')
    d.text((1100, 810), f'ΔG‡ = {levels[2]["G"]:.2f} kcal mol⁻¹', font=font(35, True), fill=INK, anchor='mm')
    d.text((2240, 1275), 'No A/B path shown', font=font(25), fill=MUTED, anchor='mm')
    d.text((90, 1350), 'Reference: isolated butadiene (C1-C2-C3-C4 torsion 26.9°) + ethylene; each species represented by one conformer.', font=font(27), fill=MUTED)
    d.text((90, 1410), 'Discrete RRHO free energies; connectors guide the eye. These are not the electronic energies along the IRC.', font=font(27), fill=MUTED)
    d.text((90, 1470), 'B is the IRC-connected product. The A/B interconversion pathway is under verification.', font=font(27), fill=MUTED)
    save(im, 'gibbs-free-energy-profile.png')
    metadata = dict(source='analysis/thermochemistry.json',
                    temperature_K=298.15, pressure_atm=1.,
                    levels=levels, figure='experiment-report/figures/gibbs-free-energy-profile.png',
                    quantity='Discrete RRHO Gibbs free energies with external rotational symmetry correction',
                    reference='One isolated twisted s-cis-region butadiene conformer and one isolated ethylene',
                    verified_connection='Encounter complex to cycloaddition TS to cyclohexene B',
                    A_B_connection_displayed=False, IRC_energies_used=False)
    (DEST/'free-energy-figure-data.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__': main()
