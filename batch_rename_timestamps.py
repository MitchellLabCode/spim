from pathlib import Path

# --- Config ---
datdir = Path('F:\\PROJECTS\\LightMicroscopyBootcamp2025\\bynGAL4klar_UASmChCAAXUASGFPnls\\20251011181532_unpacked\\')

channels = [0, 1]
angles = list(range(21, 361, 60))   # MATLAB 0:60:360 (inclusive)
times = range(2, 3)              # MATLAB 0:199 (inclusive)
addt = 0

for ang in angles:
    for ch in channels:
        for tt in times[::-1]:
            # new timestamp
            tout = tt + addt   # floor(tt * 0.5) + 2

            # expected input filename
            filename = f"stack_*_channel_{ch}_tp_{tt:05d}_angle_{ang}.ome.tif"
            matches = list(datdir.glob(filename))

            if len(matches) > 1:
                raise RuntimeError(f"More than one match for filename: {filename}")
            elif len(matches) == 1:
                src = matches[0]
                dst = src.with_name(f"tp_{tout:06d}_c{ch}_angle_{ang}.ome.tif")

                print(f"Seeking: {filename}")
                print(f"renaming {src.name} > {dst.name}")

                if dst.exists():
                    raise FileExistsError(
                        f"Output filename already exists! Would overwrite: {dst}"
                    )

                try:
                    src.rename(dst)
                except Exception as e:
                    raise RuntimeError(f"Could not move file: {e}")
            else:
                print(f"Could not find: {filename}")
