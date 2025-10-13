#!/usr/bin/env python3
# Remove interest points based on local density (via KDTree)
# Port of NPMitchell 2023 MATLAB script → Python
# Requires: numpy, matplotlib, scipy

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from datetime import datetime

# ----------------------------
# User parameters (edit here)
# ----------------------------
rootdir = './interestpoints/'           # directory containing *.beads.ip.txt files
tps = list(range(0, 119))               # timepoints to prune (0:118)
vtiles = list(range(0, 18))             # view tiles to prune (0:17)

# Experimental data spacing
dz = 1.0      # um, step between frames in z
dx = 0.195    # um, pixel size in xy

# Density-based pruning
useDensity = True
densityThres = 0.02
preview = True
pausetime = 0.1  # seconds

# ROI-based pruning (placeholder to mirror MATLAB behavior)
useROI = False
ROIs = [ [] ]        # if useROI is True, define as list of ROI specs
ROIstyle = 'cylinder'  # 'box' or 'cylinder' (not implemented below)

# Include or exclude interest points matching condition
includeExclude = 'exclude'  # 'exclude' or 'include'

# ----------------------------
# Helper functions
# ----------------------------

def read_ip_file(fn):
    """
    Read interest point file with header 'id\tx\ty\tz' and
    return (ids, xyz) with xyz shape (N,3).
    """
    try:
        data = np.loadtxt(fn, delimiter='\t', skiprows=1)
    except Exception as e:
        raise RuntimeError(f"Failed to read {fn}: {e}")

    if data.ndim == 1:
        data = data.reshape(1, -1)

    # Validate that #rows == max(id)+1 as in MATLAB assert
    ids = data[:, 0].astype(int)
    if data.shape[0] != int(np.max(ids)) + 1:
        raise AssertionError(
            f"Row count {data.shape[0]} != max(id)+1 ({int(np.max(ids))+1}) in {fn}"
        )
    xyz = data[:, 1:4]
    return ids, xyz


def write_ip_file(fn, xyz):
    """
    Write interest points back to file with header and zero-based ids.
    """
    header = "id\tx\ty\tz"
    n = xyz.shape[0]
    ids = np.arange(n, dtype=int)

    # Write text with specified formatting
    with open(fn, 'w') as f:
        f.write(header + '\n')
        for i in range(n):
            f.write(f"{ids[i]}\t{xyz[i,0]:0.14f}\t{xyz[i,1]:0.14f}\t{xyz[i,2]:0.14f}\n")


def local_density_knn(xyz, k=10, z_scale=1.0):
    """
    Compute a simple local density estimate = 1 / mean(distance to k-NN),
    after scaling z by z_scale to account for anisotropy (dz/dx).

    Returns:
        local_densities: (N,) array
        distances: (N, k) distances to k nearest neighbors
    """
    pts = xyz.copy()
    pts[:, 2] *= z_scale
    tree = cKDTree(pts)
    dists, idxs = tree.query(pts, k=k)  # includes self at distance 0
    # mean across neighbors (including self, matching MATLAB behavior)
    mean_d = np.mean(dists, axis=1)
    with np.errstate(divide='ignore'):
        rho = 1.0 / mean_d
    # replace inf (for mean_d==0) with large number
    rho[np.isinf(rho)] = np.nanmax(rho[np.isfinite(rho)]) if np.any(np.isfinite(rho)) else 0.0
    return rho, dists


def do_preview(tp, view_id, xyz, values, reject_mask, densityThres, pausetime):
    """
    Show side-by-side scatter plots like the MATLAB preview.
    """
    plt.clf()
    plt.colormaps()  # no-op, but keeps parity with MATLAB colormap usage
    cmap = plt.cm.bwr

    ax1 = plt.subplot(1, 2, 1, projection='3d')
    p1 = ax1.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], s=10, c=values, cmap=cmap)
    ax1.set_xlabel('x'); ax1.set_ylabel('y'); ax1.set_zlabel('z'); ax1.set_box_aspect([1,1,1])
    ax1.set_title('Local density')
    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_clim(0, densityThres)
    plt.colorbar(p1, ax=ax1, orientation='horizontal', fraction=0.046, pad=0.04)

    ax2 = plt.subplot(1, 2, 2, projection='3d')
    p2 = ax2.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], s=10, c=reject_mask.astype(float), cmap=cmap)
    ax2.set_xlabel('x'); ax2.set_ylabel('y'); ax2.set_zlabel('z'); ax2.set_box_aspect([1,1,1])
    ax2.set_title('Rejected (1=yes)')
    plt.suptitle(f"tp={tp}: view {view_id}")
    plt.gcf().set_facecolor('white')
    plt.pause(pausetime)


# ----------------------------
# Main
# ----------------------------

def main():
    os.makedirs(rootdir, exist_ok=True)

    # Color map for plotting counts per view
    colors = plt.cm.tab20(np.linspace(0, 1, max(20, len(vtiles))))  # enough distinct colors

    # Prealloc: counts removed/kept per view x time
    nipsRm = [np.zeros(len(tps), dtype=int) for _ in vtiles]
    nipsKeep = [np.zeros(len(tps), dtype=int) for _ in vtiles]

    for tidx, tp in enumerate(tps):
        if tp % 10 == 0:
            print(f"pruning interest points for tp={tp}/{max(tps)}")

        for vidx, vtile in enumerate(vtiles):
            fn = os.path.join(rootdir, f"tpId_{tp}_viewSetupId_{vtile}.beads.ip.txt")

            if not os.path.isfile(fn):
                raise FileNotFoundError(
                    f"No interest point file found for tp={tp}, view={vtile}: {fn}"
                )

            # Read
            _, ips_xyz = read_ip_file(fn)  # shape (N,3)

            # Default: keep all
            reject = np.zeros(ips_xyz.shape[0], dtype=bool)

            # Density-based pruning
            if useDensity and ips_xyz.shape[0] > 0:
                z_scale = dz / dx
                localDensities, _ = local_density_knn(ips_xyz, k=10, z_scale=z_scale)

                if includeExclude.lower() == 'exclude':
                    reject = localDensities > densityThres
                elif includeExclude.lower() == 'include':
                    reject = localDensities < densityThres
                else:
                    raise ValueError("includeExclude must be 'include' or 'exclude'")

                if preview:
                    do_preview(tp, vtile, ips_xyz, localDensities, reject, densityThres, pausetime)

            # ROI-based pruning (not implemented, matching MATLAB error)
            if useROI:
                raise NotImplementedError(
                    "ROI pruning step not yet implemented; independent of density pruning."
                )

            # Apply mask, re-id from 0..N-1, and write back
            kept_xyz = ips_xyz[~reject]
            write_ip_file(fn, kept_xyz)

            # Record counts
            nipsRm[vidx][tidx] = int(np.sum(reject))
            nipsKeep[vidx][tidx] = kept_xyz.shape[0]

    # Plot counts after pruning
    plt.close('all')
    plt.figure(figsize=(8, 5))
    legEntries = []
    for vidx, vtile in enumerate(vtiles):
        plt.plot(tps, nipsKeep[vidx], '.-', color=colors[vidx % len(colors)])
        legEntries.append(f"View {vtile}")
    plt.xlabel('timepoint')
    plt.ylabel('# interest points')
    plt.title('Counts after interest point pruning')
    plt.legend(legEntries, loc='center left', bbox_to_anchor=(1.02, 0.5))
    plt.tight_layout()

    stamp = datetime.now().strftime('%Y%m%d%H%M')
    outpng = f"interestpoint_counts_Pruned_{stamp}.png"
    plt.savefig(outpng, dpi=200)
    print(f"Saved: {outpng}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
