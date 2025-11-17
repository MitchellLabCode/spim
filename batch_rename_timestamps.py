from pathlib import Path

# --- Config ---
datdir = Path('F:\\PROJECTS\\LightMicroscopyBootcamp2025\\bynGAL4klar_UASmChCAAXUASGFPnls\\20251011181532_unpacked\\')
datdir = Path('E:\\avistrok\\bapGAL4_UAShidUASstingerHiRFP\\20250806120735_bapGAL4_UAShidUASStingerHiRFP_5mpf_22x\\20250806121230_bapGAL4_UAShidUASStingerHiRFP_5mpf_22x\\unpacked\\')
datdir = Path('E:\\avistrok\\bapGAL4_UAShidUASstingerHiRFP\\20250806120735_bapGAL4_UAShidUASStingerHiRFP_5mpf_22x\\2025-08-06_120735_tp0\\unpacked\\')

channels = [0, 1]
angles = list(range(6, 361, 60))   # MATLAB 0:60:360 (inclusive)
times = range(0, 65)              # MATLAB 0:199 (inclusive)
addt = -1

for ang in angles:
    for ch in channels:
        for tt in times[::-1]:
            # new timestamp
            tout = tt + addt   # floor(tt * 0.5) + 2

            # expected input filename
            filename = f"stack_*_channel_{ch}_*_tp_{tt:05d}_angle_{ang}.ome.tif"
            matches = list(datdir.glob(filename))

            if len(matches) > 1:
                raise RuntimeError(f"More than one match for filename: {filename}")
            elif len(matches) == 1:
                src = matches[0]
                dst = src.with_name(f"t{tout:06d}_c{ch}_a{ang}.ome.tif")

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
