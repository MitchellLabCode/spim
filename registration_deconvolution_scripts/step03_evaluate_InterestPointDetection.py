#!/usr/bin/env python3
# Evaluate # of interest points per view
# MATLAB → Python port
# Requires: numpy, matplotlib, scipy

import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from datetime import datetime

# ----------------------------
# User parameters
# ----------------------------
datdir = 'F:\\PROJECTS\\LightMicroscopyBootcamp2025\\bynGAL4klar_UASmChCAAXUASGFPnls\\20251011181532_bynGAL4klar_UASmChCAAX_UASGFPnls_combined\\'
rootdir = os.path.join(datdir,'interestpoints\\')
# rootdir = './interestpoints/'           # directory containing *.beads.ip.txt files

# imfn = 'c*_t*_a*.ome.tif' 
imfn = 'tp_*_c*_angle_*.ome.tif' 
tca = 'tca'   # the order of time, channel, and angle in the image filenames

previewIPs = True
pauseTime = 0.05
previewEveryN = 30
densityThres = 0.02   # threshold for local density color scale
dz = 1.0              # um, axial sampling
dx = 0.2925  # 0.195            # um, in-plane sampling

# ----------------------------
# Helpers
# ----------------------------
def read_ip_file(fn):
    """Read interest point file with header 'id\\tx\\ty\\tz' -> (ids, xyz)."""
    data = np.loadtxt(fn, delimiter='\t', skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    ids = data[:, 0].astype(int)
    # assert size(ips,1) == max(id)+1
    if data.shape[0] != (np.max(ids) + 1):
        raise AssertionError(
            f"Row count {data.shape[0]} != max(id)+1 ({np.max(ids)+1}) in {fn}"
        )
    xyz = data[:, 1:4]
    return ids, xyz

import re

def build_regex_from_template(imfn: str):
    """
    Build a regex with named groups (?P<t>), (?P<c>), (?P<a>) from a template.
    Supports either token mode:  'tp_{t}_c{c}_angle_{a}.ome.tif'
    or heuristic wildcard mode: 'tp_*_c*_angle_*.ome.tif' (common MVB/Luxendo style).
    """
    # Token mode (preferred)
    if any(tok in imfn for tok in ('{t}', '{c}', '{a}')):
        pat = re.escape(imfn)
        pat = pat.replace(r'\{t\}', r'(?P<t>\d+)')
        pat = pat.replace(r'\{c\}', r'(?P<c>\d+)')
        pat = pat.replace(r'\{a\}', r'(?P<a>\d+)')
        # For any other '*' you may have in the template, treat as non-greedy wildcard
        pat = pat.replace(r'\*', r'.*?')
        return re.compile(r'^' + pat + r'$')

    # Heuristic mode for the specific style: 'tp_*_c*_angle_*.ome.tif'
    # We interpret the '*' after 'tp_' as {t}, after 'c' as {c}, after 'angle_' as {a}.
    # Any other '*' are treated as wildcards.
    s = re.escape(imfn)

    # Replace escaped anchors for the three key numeric fields
    s = s.replace(r'tp_\*', r'tp_(?P<t>\d+)')
    s = s.replace(r'c\*',   r'c(?P<c>\d+)')
    s = s.replace(r'angle_\*', r'angle_(?P<a>\d+)')

    # Remaining '*' (if any) become non-greedy wildcards
    s = s.replace(r'\*', r'.*?')
    return re.compile(r'^' + s + r'$')


def parse_filename(name: str, imfn: str, tca: str = 'tca'):
    """
    Parse filename using `imfn` template and return (c, t, a) following your `tca` ordering.
    - `tca` is a permutation of 'cta' describing the order of numeric fields in the *string*.
    """
    tca = tca.lower()
    if sorted(tca) != ['a', 'c', 't']:
        raise ValueError("tca must be a permutation of 'cta'.")

    pat = build_regex_from_template(imfn)
    m = pat.search(name)
    if not m:
        raise ValueError(f"Filename does not match template: {name}  vs  {imfn}")

    # We captured by *names*, so get canonical values:
    try:
        vals = {k:int(v) for k,v in m.groupdict().items()}
    except Exception as e:
        raise ValueError(f"Template must expose numeric fields for t,c,a. Error: {e}")

    # If you want to *enforce/reflect* the declared tca order, we can check consistency:
    # Build the order we *saw* in the filename by scanning for the named groups.
    # (This is optional; named groups already give correct values regardless of order.)
    # Here we just return the canonical triplet:
    return vals['c'], vals['t'], vals['a']


def local_density_knn(xyz, k=10, z_scale=1.0):
    """Return local density = 1 / mean(distance to kNN) after scaling z."""
    pts = xyz.copy()
    pts[:, 2] *= z_scale
    tree = cKDTree(pts)
    dists, _ = tree.query(pts, k=k)        # includes self (distance 0)
    mean_d = np.mean(dists, axis=1)
    with np.errstate(divide='ignore'):
        rho = 1.0 / mean_d
    # Handle inf if any zero means
    if np.any(~np.isfinite(rho)):
        finite = rho[np.isfinite(rho)]
        fill = np.max(finite) if finite.size else 0.0
        rho[~np.isfinite(rho)] = fill
    return rho

# ----------------------------
# Auto-detect timepoints and view-tiles from file names
# ----------------------------
# Expecting names like: c<chan>_t<tp>_a<angle>.ome.tif or similar. The order can be swapped based on tca.
tp_files = glob.glob(os.path.join(datdir, imfn))

if not tp_files:
    raise FileNotFoundError("No files matching 'c*_t*_a*.ome.tif' were found.")

tp_vals = []
vtile_pairs = []  # (channel, angle)
for fn in tp_files:
    base = os.path.basename(fn)
    try:
        c, t, a = parse_filename(base, imfn=imfn, tca=tca)
    except ValueError:
        # filename didn't match the template; skip
        continue
    tp_vals.append(t)
    vtile_pairs.append((c, a))

tps = sorted(set(tp_vals))

# Unique (channel, angle) tuples; then define view IDs as 0..N-1
unique_pairs = sorted(set(vtile_pairs))
vtiles = list(range(len(unique_pairs)))  # MATLAB’s vtiles = 0:size(unique_rows)-1

# ----------------------------
# Count interest points
# ----------------------------
os.makedirs(rootdir, exist_ok=True)

colors = plt.cm.tab20(np.linspace(0, 1, max(20, len(vtiles))))
nips = [np.zeros(len(tps), dtype=int) for _ in vtiles]
missing_files = []

z_scale = dz / dx

for tidx, tp in enumerate(tps, start=1):  # 1-based to mirror MATLAB modulo logic
    if tp % 10 == 0:
        print(f"counting interest points for tp={tp}/{max(tps)}")

    for vidx, vtile in enumerate(vtiles):
        fn = os.path.join(rootdir, f"tpId_{tp}_viewSetupId_{vtile}.beads.ip.txt")
        if os.path.isfile(fn):
            try:
                ids, xyz = read_ip_file(fn)
            except Exception as e:
                print(f"[WARN] Could not read '{fn}': {e}")
                missing_files.append(f"tp={tp} vtile={vtile}")
                continue

            nips[vidx][tidx - 1] = xyz.shape[0]

            # Preview every previewEveryN timepoints (when tidx==1, 1+N, ...)
            if previewIPs and ((tidx - 1) % previewEveryN == 0) and xyz.shape[0] > 0:
                localDensities = local_density_knn(xyz, k=10, z_scale=z_scale)

                plt.clf()
                ax = plt.axes(projection='3d')
                sc = ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2],
                                s=10, c=localDensities, cmap='bwr')
                ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
                ax.set_box_aspect([1, 1, 1])
                sc.set_clim(0, densityThres)
                plt.colorbar(sc)
                ax.view_init(elev=90, azim=-90)  # view(2) top-down
                plt.title(f"tp={tp}: view {vtile}")
                plt.pause(pauseTime)
        else:
            missing_files.append(f"tp={tp} vtile={vtile}")

# ----------------------------
# Report missing
# ----------------------------
if missing_files:
    print("Missing files:")
    for msg in missing_files:
        print(msg)
else:
    print("All expected files are present.")

# ----------------------------
# Plot summary over time
# ----------------------------
plt.close('all')
plt.figure(figsize=(8, 5))
leg = []
for vidx, vtile in enumerate(vtiles):
    plt.plot(tps, nips[vidx], '.-', color=colors[vidx % len(colors)])
    leg.append(f"View {vtile}")

plt.xlabel('timepoint')
plt.ylabel('#interest points')
plt.title('Interest points summary')
plt.legend(leg, loc='upper left', bbox_to_anchor=(1.02, 1.0))
plt.gcf().set_facecolor('white')
plt.tight_layout()

stamp = datetime.now().strftime('%Y%m%d%H%M')
outpng = f"interestpoint_counts_{stamp}.png"
plt.savefig(os.path.join(datdir,outpng), dpi=200)
print(f"Saved: {outpng}")
