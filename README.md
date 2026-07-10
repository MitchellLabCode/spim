# spimCode

Code for processing SPIM (Selective Plane Illumination Microscopy) / lightsheet imaging data: unpacking raw acquisitions, registration and multi-view fusion, deconvolution, and generating preview movies/MIPs.

## Repository layout

| Directory | Purpose |
|---|---|
| [`registration_deconvolution_scripts/`](registration_deconvolution_scripts/) | Core pipeline: extract raw stacks, detect interest points, register/fuse multi-view data, deconvolve, and generate MIPs/movies at each stage. |
| [`luxendo/`](luxendo/) | Scripts for unpacking and processing data from Luxendo lightsheet microscopes (unpack to Fiji-readable stacks, max projections, masking/background subtraction, bundling hyperstacks). |
| [`zebrafishImaging/`](zebrafishImaging/) | Scripts specific to zebrafish heart imaging experiments (defining stacks/rotation angles, deleting bad timepoints, tracking acquisition progress). |
| [`SLURM/`](SLURM/) | Template for running batch deconvolution jobs on an HPC cluster (RCC) via `sbatch`. See [SLURM/README](SLURM/README) for the full protocol. |
| [`hash_md5_parsing/`](hash_md5_parsing/) | Utilities for verifying that a dataset transferred between machines is identical, by comparing MD5 checksum lists across different file paths. |

## Registration & deconvolution pipeline

Scripts in `registration_deconvolution_scripts/` are numbered by pipeline stage:

1. **`step00_extractStacksAndCropLive2.m`** — Disentangle raw Micro-Manager `.ome.tif` output (interleaved cameras/timepoints/z-slices) into per-timepoint, per-angle, per-channel stacks.
2. **`step01_makeMIPS_of_unpackedData.m`** — Generate maximum-intensity projections (MIPs) of the unpacked data for a quick preview.
3. **`step02_makeMoviesFromRawDataMIPs_2color.ijm`** — Build preview movies from the raw-data MIPs (2-color).
4. **`step03*`** — Interest point detection, registration, and fusion:
   - `detectIP_Register_Fuse.ijm`, `register_5x*.ijm` — detect interest points and register views (single timepoint, all-to-all, or per-timepoint; 1- or 2-color variants).
   - `step03_combineInterestPointXML.m` — merge interest points detected separately per channel/timepoint into one XML file.
   - `step03b_pruneInterestPoints.m` / `.py` — prune interest points by local density (KD-tree).
   - `step03_evaluate_InterestPointDetection.m` / `.py`, `step03_checkInterestPoints.py`, `step03_check_fused_images_for_boundingbox.ijm` — QC for interest point detection and fusion bounding box.
5. **`step04*`** — Deconvolution and final MIPs:
   - `deconvolve_timepoints.ijm` — run multi-view deconvolution (BigStitcher/Fiji) over a range of timepoints.
   - `step04_makeMIPS.m`, `step04_makeMIPs_fusedOrDeconvolved_2color.ijm`, `step04_makeMIPs_resliceFusedOrDeconvolved_2color.ijm` — generate MIPs (including resliced views) from fused or deconvolved output.
   - `step04_evaluateFusedorDeconvolvedIJ.ijm` — visual QC of fused/deconvolved output.

Supporting helpers: `combine_xml_files.m`, `bwr.m` (colormap), `extractStacks_Andor_Fiji_4view_series_2color_npm.m`, `equalize_histogram_for_each_timepoint.ijm`.

## Cluster deconvolution (SLURM)

`SLURM/` holds a template for running deconvolution as a batch job on RCC:

- `decon_main.sh` — top-level `sbatch` script; sets CPU/RAM allocation and the timepoint range to process, then calls `deconvolve_timepoints.sh`.
- `deconvolve_timepoints.sh` — launches Fiji headless and runs the ImageJ macro `batch_macro.bsh`.
- `batch_macro.bsh` — the deconvolution macro; edit this with your dataset's file path, channel, magnification, and bounding box.

**Usage:**

1. Copy this directory elsewhere (don't edit the template in place) and follow the full protocol on [LabResources](https://docs.google.com/document/d/1POgY4802YNJ5gKTr1Buk5kvwCBillGMHnD1IpRqOSRE/edit?usp=sharing).
2. Edit `batch_macro.bsh` with your parameters: XML dataset directory, processing channel, timepoint(s), bounding box limits, iterations, output directory, and PSF paths (specified twice — change both).
3. Edit `decon_main.sh` for memory allocation and timepoint selection.
4. Submit the job: `sbatch decon_main.sh` (or the equivalent `*.sh` script in your copy), then monitor/manage it with:
   - `squeue -u <cnetID>` — check job status
   - `scancel <jobID>` — cancel a job

Notes:
- Known-working shared-partition configs: `partition=caslake, mem=0G, cpus-per-task=40` or `partition=amd, mem=0G, cpus-per-task=64`. Also set `account=pi-npmitchell` and a job name. See the [RCC partitions guide](https://rcc-uchicago.github.io/user-guide/partitions/#configurations) for more options.
- The default Fiji install path is `/project/npmitchell/fiji-linux64/Fiji.app/`; only change it if you need a different Fiji installation.

## Checksum comparison (`hash_md5_parsing/`)

Verify that a dataset copied between two computers is actually identical, using MD5 checksums (e.g. from `md5deep`) generated separately on each machine.

Since the two machines mount the data at different paths, the checksum lists won't match directly by path. Fix that first, then compare:

```bash
cd hash_md5_parsing/

# 1. Rewrite paths in the second machine's hash list to match the first's
python replace_substring_in_hashfile.py examples/hashlist2.txt examples/hashlist2_fixed.txt "E:\Runt" "D:\Atlas_Data\Runt"

# 2. Compare the two hash lists
python compare_hashes.py <hashlist_file> <master_hashlist_file> <output_filename>
# e.g.
python compare_hashes.py examples/hashlist1.txt examples/hashlist_tweak.txt examples/hash_comparison.txt
```

See `hash_md5_parsing/examples/` for sample inputs/outputs.
