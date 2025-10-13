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
rootdir = './interestpoints/'
previewIPs = True
pauseTime = 0.1
previewEveryN = 10
densityThres = 0.02   # threshold for local density color scale
dz = 1.0              # um, axial sampling
dx = 0.195            # um, in-plane sampling

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
# Expecting names like: c<chan>_t<tp>_a<angle>.ome.tif
tp_files = glob.glob('./c*_t*_a*.ome.tif')

if not tp_files:
    raise FileNotFoundError("No files matching './c*_t*_a*.ome.tif' were found.")

tp_vals = []
vtile_pairs = []  # (channel, angle)
pat = re.compile(r'c(\d+)_t(\d+)_a(\d+)\.ome\.tif$')

for fn in tp_files:
    m = pat.search(os.path.basename(fn))
    if not m:
        continue
    c, t, a = map(int, m.groups())
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
plt.savefig(outpng, dpi=200)
print(f"Saved: {outpng}")
