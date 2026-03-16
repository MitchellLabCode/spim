#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Simple script to mask a TIFF file based on thresholded probabilities from an H5 file.
"""

import os
import h5py
import numpy as np
import tifffile as tiff
from pathlib import Path
from scipy.ndimage import zoom as ndi_zoom

# ==== USER CONFIG ============================================================
datdir = Path("/Users/npmitchell/Dropbox/Soft_Matter_Biophysics/PAPER/gut_chirality/figures/byn_domain")  # Base directory for data
H5_PATH = os.path.join(datdir, "APDV_Time_000005_8x_Probabilities.h5")        # Path to HDF5 file with probabilities
H5_DATASET = "/exported_data"                # Dataset path within H5 file
PROB_CHANNEL = 0                             # Which channel/class to use (0 for first)
DOWNSAMPLE_FACTOR = 4                        # Downsampling factor of prob data (e.g., 4 = 4x downsampled)
THRESHOLD = 0.5                              # Probability threshold for mask (0.0-1.0)
TIFF_INPUT = os.path.join(datdir, "APDV_Time_000005_2x.tif")             # Input TIFF file
TIFF_OUTPUT = os.path.join(datdir, "APDV_Time_000005_2x_masked.tif")    # Output masked TIFF file
TIFF_CHANNEL_IDX = 0                         # Which channel to extract if 4D (0 for first)
TIFF_CHANNEL_AXIS = 1                       # Which axis holds channels if 4D (-1=last, 0=first, 1=second, etc)
PRESERVE_DTYPE = True                        # Cast output back to input dtype
# ============================================================================


def load_probabilities_from_h5(h5_path, dataset_path, channel_idx):
    """
    Load probability data from HDF5 file.
    
    Assumes 4D data with shape (Z, C, Y, X) where:
    - Z: depth/slices
    - C: channels/classes
    - Y, X: spatial dimensions
    
    Returns: 3D array (Z, Y, X) for specified channel
    """
    with h5py.File(h5_path, "r") as f:
        data = f[dataset_path][()]  # Load to memory
    
    if data.ndim != 4:
        raise ValueError(f"Expected 4D data (Z,C,Y,X), got shape {data.shape}")
    
    Z, C, Y, X = data.shape
    if channel_idx < 0 or channel_idx >= C:
        raise IndexError(f"Channel {channel_idx} out of range [0, {C-1}]")
    
    # Extract channel and ensure float32
    prob = data[:, channel_idx, :, :].astype(np.float32)
    return prob


def upsample_to_shape(arr, target_shape, order=1):
    """
    Upsample array (Z,Y,X) to target_shape using SciPy interpolation.
    Uses scipy.ndimage.zoom for smooth interpolation.
    """
    sz, sy, sx = arr.shape
    tz, ty, tx = target_shape
    
    # Calculate zoom factors
    zoom_factors = (tz / sz, ty / sy, tx / sx)
    
    # Upsample
    upsampled = ndi_zoom(arr, zoom=zoom_factors, order=order, prefilter=(order > 1))
    
    # Ensure exact shape (zoom may be off by 1 due to rounding)
    if upsampled.shape != target_shape:
        # Pad or crop to exact target shape
        upsampled = _fit_to_shape(upsampled, target_shape)
    
    return upsampled


def _fit_to_shape(arr, target_shape, pad_value=0):
    """
    Center-pad or center-crop array to exact target shape.
    """
    z, y, x = arr.shape
    tz, ty, tx = target_shape
    
    # Handle Z dimension
    if z < tz:
        pad_before = (tz - z) // 2
        pad_after = tz - z - pad_before
        arr = np.pad(arr, ((pad_before, pad_after), (0, 0), (0, 0)),
                     mode="constant", constant_values=pad_value)
    elif z > tz:
        start = (z - tz) // 2
        arr = arr[start:start + tz, :, :]
    
    z, y, x = arr.shape
    
    # Handle Y dimension
    if y < ty:
        pad_before = (ty - y) // 2
        pad_after = ty - y - pad_before
        arr = np.pad(arr, ((0, 0), (pad_before, pad_after), (0, 0)),
                     mode="constant", constant_values=pad_value)
    elif y > ty:
        start = (y - ty) // 2
        arr = arr[:, start:start + ty, :]
    
    z, y, x = arr.shape
    
    # Handle X dimension
    if x < tx:
        pad_before = (tx - x) // 2
        pad_after = tx - x - pad_before
        arr = np.pad(arr, ((0, 0), (0, 0), (pad_before, pad_after)),
                     mode="constant", constant_values=pad_value)
    elif x > tx:
        start = (x - tx) // 2
        arr = arr[:, :, start:start + tx]
    
    return arr


def apply_mask(volume, mask):
    """
    Apply binary mask to volume by multiplication.
    
    Args:
        volume: Input image array (any dtype)
        mask: Binary mask (0 or 1 values)
    
    Returns:
        Masked volume (float32)
    """
    return volume.astype(np.float32) * mask.astype(np.float32)


def main():
    print("[INFO] Loading probabilities from H5...")
    prob = load_probabilities_from_h5(H5_PATH, H5_DATASET, PROB_CHANNEL)
    print(f"[INFO] Probability shape: {prob.shape} (Z, Y, X)")
    print(f"[INFO] Probability downsampling factor: {DOWNSAMPLE_FACTOR}x")
    
    print("[INFO] Loading TIFF image...")
    tiff_data = tiff.imread(str(TIFF_INPUT))
    print(f"[INFO] TIFF shape: {tiff_data.shape}")
    orig_dtype = tiff_data.dtype
    
    # Determine if we have multi-channel data and reshape to (Z, Y, X, C) internally if needed
    if tiff_data.ndim == 3:
        # Single channel (Z, Y, X)
        spatial_shape = tiff_data.shape
        num_channels = 1
        tiff_data_reshaped = tiff_data[..., np.newaxis]  # Add channel dim: (Z, Y, X, 1)
    elif tiff_data.ndim == 4:
        # Multi-channel - determine which axis is the channel axis
        if TIFF_CHANNEL_AXIS == 1:  # (Z, C, Y, X) - need to reorder to (Z, Y, X, C)
            tiff_data = np.transpose(tiff_data, (0, 2, 3, 1))
        elif TIFF_CHANNEL_AXIS == 0:  # (C, Z, Y, X) - need to reorder to (Z, Y, X, C)
            tiff_data = np.transpose(tiff_data, (1, 2, 3, 0))
        # else axis -1 (3) is already last: (Z, Y, X, C)
        
        spatial_shape = tiff_data.shape[:3]  # (Z, Y, X)
        num_channels = tiff_data.shape[3]
        tiff_data_reshaped = tiff_data
    else:
        raise ValueError(f"Expected 3D or 4D TIFF, got {tiff_data.ndim}D shape {tiff_data.shape}")
    
    print(f"[INFO] Detected {num_channels} channel(s), spatial shape: {spatial_shape}")
    
    # Upsample probabilities to match TIFF spatial shape if needed
    if prob.shape != spatial_shape:
        print(f"[INFO] Upsampling probabilities from {prob.shape} to {spatial_shape}...")
        print(f"[INFO]   (accounting for {DOWNSAMPLE_FACTOR}x downsampling factor)")
        prob = upsample_to_shape(prob, spatial_shape, order=1)
    
    # Threshold probabilities to create binary mask
    print(f"[INFO] Thresholding probabilities at {THRESHOLD}...")
    mask = (prob >= THRESHOLD).astype(np.float32)
    print(f"[INFO] Mask coverage: {np.sum(mask > 0) / mask.size * 100:.1f}%")
    
    # Apply mask to all channels
    print(f"[INFO] Applying mask to {num_channels} channel(s)...")
    masked_all_channels = np.zeros_like(tiff_data_reshaped, dtype=np.float32)
    for ch in range(num_channels):
        masked_all_channels[:, :, :, ch] = apply_mask(tiff_data_reshaped[:, :, :, ch], mask)
    
    # Reshape back to original structure
    if num_channels == 1:
        masked = masked_all_channels[:, :, :, 0]  # Remove added channel dimension
    elif TIFF_CHANNEL_AXIS == 1:  # Need to return to (Z, C, Y, X)
        masked = np.transpose(masked_all_channels, (0, 3, 1, 2))
    elif TIFF_CHANNEL_AXIS == 0:  # Need to return to (C, Z, Y, X)
        masked = np.transpose(masked_all_channels, (3, 0, 1, 2))
    else:  # Keep as (Z, Y, X, C)
        masked = masked_all_channels
    
    # Cast back to original dtype if requested
    if PRESERVE_DTYPE:
        if np.issubdtype(orig_dtype, np.integer):
            info = np.iinfo(orig_dtype)
            masked = np.clip(np.rint(masked), info.min, info.max).astype(orig_dtype)
        else:
            masked = masked.astype(orig_dtype)
    
    # Create output directory if needed
    output_path = Path(TIFF_OUTPUT)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Write output - always reshape to (Z, Y, X, C) for proper tifffile/Fiji compatibility
    print(f"[INFO] Writing masked TIFF to {TIFF_OUTPUT}...")
    print(f"[INFO] Output shape before write: {masked.shape}")
    
    # Ensure output is in (Z, Y, X, C) format for tifffile
    if masked.ndim == 3:
        # Single channel, add channel dim
        masked_write = masked[..., np.newaxis]
    elif masked.ndim == 4:
        # Multi-channel - check if we need to reshape
        if TIFF_CHANNEL_AXIS == 1 or TIFF_CHANNEL_AXIS == 0:
            # Currently in (Z, C, Y, X) or (C, Z, Y, X), need to reshape to (Z, Y, X, C)
            if masked.shape[1] <= 4:  # Likely channels in axis 1
                masked_write = np.transpose(masked, (0, 2, 3, 1))
            else:  # Likely channels in axis 0 or it's already (Z, Y, X, C)
                masked_write = masked
        else:
            # Already in (Z, Y, X, C) format
            masked_write = masked
    else:
        masked_write = masked
    
    print(f"[INFO] Output shape to write: {masked_write.shape}")
    tiff.imwrite(str(output_path), masked_write, compression=None)
    
    print("[OK] Done!")


if __name__ == "__main__":
    main()
