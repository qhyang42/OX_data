#!/usr/bin/env python3
"""Create a 16:9 committee-slide montage of the saved binary-mask overlap."""
import argparse
import shutil
from datetime import datetime
from pathlib import Path
import json
import nibabel as nib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--plane', choices=('axial', 'coronal'), default='axial')
parser.add_argument('--context', choices=('ODOR', 'PERSON', 'FOOD', 'LOCATION', 'COMBINED'), default='ODOR')
args = parser.parse_args()
plane = args.plane
context = args.context
axis = 2 if plane == 'axial' else 1
coordinate = 'z' if plane == 'axial' else 'y'
root = Path(__file__).resolve().parents[1]
folder = root / 'results/contrast_map'
key = 'ODOR_gt_REST' if context == 'ODOR' else context + '_gt_OTHER_CONTEXTS'
source = folder / f'{key}_n_significant_MNI.nii.gz'
if context == 'COMBINED':
    source = folder / 'PERSON_gt_OTHER_CONTEXTS_n_significant_MNI.nii.gz'
reference = root / 'ROIs/CZ_ROIs_MNI_1mm/MNI152_T1_1mm_brain.nii'
target = folder / f'{key}_n_significant_MNI_{plane}_committee.png'
if target.exists():
    archive = folder / 'old' / ('display_before_update_' + datetime.now().strftime('%Y%m%dT%H%M%S_%f'))
    archive.mkdir(parents=True)
    shutil.copy2(target, archive / target.name)
img = nib.as_closest_canonical(nib.load(source))
ref = nib.as_closest_canonical(nib.load(reference))
assert img.shape == ref.shape and np.allclose(img.affine, ref.affine)
counts = img.get_fdata(); brain = ref.get_fdata()
assert np.isin(counts, np.arange(6)).all()
# 8-mm sampling spans the observed -64 to +22 mm extent in 12 panels.
slices = list(range(-64, 25, 8)) if plane == 'axial' else list(range(-100, 71, 10))
interval = 8 if plane == 'axial' else 10
colors = ['#ffff66', '#ffcc33', '#ff881a', '#ee3b16', '#c90016']
context_colors = {
    'PERSON': ['#ff9999','#ff6666','#ff3333','#e60000','#b30000'],
    'FOOD': ['#99ccff','#66b3ff','#3399ff','#0073e6','#0059b3'],
    'LOCATION': ['#99ff99','#66ee66','#33cc33','#00aa00','#008000']}
minimum = 2 if context == 'ODOR' else 1
if context in context_colors:
    colors = context_colors[context]
cmap = ListedColormap(colors[minimum-1:])
norm = BoundaryNorm(np.arange(minimum-.5, 6, 1), 6-minimum)
combined = {}
if context == 'COMBINED':
    for name in context_colors:
        im = nib.as_closest_canonical(nib.load(folder / f'{name}_gt_OTHER_CONTEXTS_n_significant_MNI.nii.gz'))
        assert im.shape == ref.shape and np.allclose(im.affine, ref.affine)
        combined[name] = im.get_fdata()
        assert np.isin(combined[name], np.arange(6)).all()
fig, axes = plt.subplots(2 if plane == 'axial' else 3, 6, figsize=(16, 9), facecolor='#090909')
fig.subplots_adjust(left=.025, right=.975, top=.85, bottom=.19, wspace=.025, hspace=.14)
vmax = np.percentile(brain[brain > 0], 99)
for ax, z in zip(axes.flat, slices):
    k = int(round((z-img.affine[axis, 3])/img.affine[axis, axis]))
    ax.set_facecolor('#090909')
    ax.imshow(np.take(brain, k, axis=axis).T, origin='lower', cmap='gray', vmin=0, vmax=vmax, interpolation='bilinear')
    if context != 'COMBINED':
        ax.imshow(np.ma.masked_less(np.take(counts, k, axis=axis).T, minimum), origin='lower', cmap=cmap, norm=norm, interpolation='nearest')
    else:
        from matplotlib.colors import to_rgb
        shape = np.take(brain, k, axis=axis).T.shape
        rgb = np.zeros((*shape, 3)); number = np.zeros(shape)
        for name, data in combined.items():
            values = np.take(data, k, axis=axis).T.astype(int)
            palette = np.array([[0,0,0]] + [to_rgb(c) for c in context_colors[name]])
            rgb += palette[values]; number += values > 0
        rgb /= np.maximum(number, 1)[...,None]
        rgba = np.concatenate([rgb, (number>0)[...,None]], axis=-1)
        ax.imshow(rgba, origin='lower', interpolation='nearest')
    ax.set_axis_off()
    ax.set_title(f'{coordinate} = {z:+d} mm', fontsize=15, color='white', pad=7)
    ax.text(.02, .5, 'L', transform=ax.transAxes, color='#cccccc', fontsize=11, va='center')
    ax.text(.98, .5, 'R', transform=ax.transAxes, color='#cccccc', fontsize=11, ha='right', va='center')
fig.text(.5, .955, ('Odor > Rest' if context == 'ODOR' else 'Semantic context contrasts' if context == 'COMBINED' else context + ' > Other Contexts'), ha='center', color='white', fontsize=29, weight='bold')
fig.text(.5, .903, f'Overlap of participant-level FWE p < .001 maps  |  5 participants  |  Display: {minimum}–5', ha='center', color='#dddddd', fontsize=17)
legends = context_colors if context == 'COMBINED' else {context: colors}
for j, (name, palette) in enumerate(legends.items()):
    x, width = (.12 + j*.28, .20) if context == 'COMBINED' else (.32, .36)
    cm = ListedColormap(palette[minimum-1:])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cm), cax=fig.add_axes([x, .115, width, .027]), orientation='horizontal', ticks=range(minimum,6))
    cb.ax.tick_params(colors='white', labelsize=14, length=0, pad=6)
    cb.outline.set_visible(False)
    cb.set_label(name + ' | significant participants' if context == 'COMBINED' else 'Number of significant participants', color='white', fontsize=12, labelpad=8)
fig.text(.5, .025, f'MNI152  •  {plane.title()} slices every {interval} mm  •  Neurological orientation (L = left)  •  {'Overlapping contexts use blended colors' if context == 'COMBINED' else 'Descriptive overlap'}', ha='center', color='#bbbbbb', fontsize=12)
fig.savefig(target, dpi=240, facecolor=fig.get_facecolor())
plt.close(fig)
(folder / f'{key}_{plane}_committee_plot.json').write_text(json.dumps(dict(source=str(source),underlay=str(reference),output=str(target),plane=plane,slice_coordinate=coordinate,slices_mni_mm=slices,colors_1_to_5=colors,display_minimum=minimum,below_minimum='transparent',context=context,context_palettes=context_colors if context == 'COMBINED' else None,combined_rule='Mean RGB of active context count colors' if context == 'COMBINED' else None,dimensions_pixels=[3840,2160],orientation='neurological'),indent=2))
print(target)
