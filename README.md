# prym-direction-check

Author: Heywood Geblomi. Public snapshot built 24 November 2026.

## What this is

A preregistered check of a counting question: does the recorded count of "20 exact 2:1 closed geodesics" on the Prym / X(1,1) work double-count the same unoriented closed geodesic (rotation or reversal of the shear word)?

The criterion is in `CRITERION.md`, hashed in `CRITERION.sha256` before the orbit data was scored. The public repo [HeywoodGeblomi/prym-omega-ed6-h1plus](https://github.com/HeywoodGeblomi/prym-omega-ed6-h1plus) at `c3f09f506f51` was read only through the GitHub API.

## Recorded result

The per-orbit list of the 20 was not in that repo, and not in any of its 11 commits. `data/x11_affine_maps.json` is cited in a note there but was never committed. Only 3 of the 20 are identifiable from what is present (traces 52, 44, and 148). Those three traces are distinct, so they are three pairwise distinct unoriented closed geodesics. Counts up to rotation, or up to rotation and reversal, stay undetermined and are only known to be at least 3. The check used a trace-only proxy.

The example `A = [[31,10],[65,21]]` has trace 52 and ρ = (2+√3)^3. That example's `M_plus` spectrum is consistent with an exact 2:1 ratio. A naive shear product of the example word does not reconstruct that A. See `RESULTS.md`.

Later rebuild and square-tiled notes live in `rebuild/` and `square_tiled/`, each with their own prereg and hash files. Read those RESULTS files for what those follow-up checks did and did not show.

## How to check it

```
sha256sum -c CRITERION.sha256
sha256sum -c rebuild/HASHES.sha256
sha256sum -c rebuild/PREREG.sha256
sha256sum -c square_tiled/HASHES.sha256
sha256sum -c square_tiled/PREREG.sha256
python3 rebuild/checks.py
python3 square_tiled/validate.py
```

## Added data

`darryl_share.txt` and the literature text files `rebuild/lit/ln_g4.txt` and `rebuild/lit/ln_g34.txt` are in this repo. The two literature PDFs could not be uploaded: the GitHub connector used here accepts text only, and `$file:` rejection for non-text files blocked `rebuild/lit/ln_g4.pdf` and `rebuild/lit/ln_g34.pdf`. I did not use `gh` or a raw API token as a workaround.

## Restoring large and binary files

The literature PDFs could not be stored as binary blobs through the connector used to publish this repo. They are here as base64 text: `rebuild/lit/ln_g4.pdf.b64` and `rebuild/lit/ln_g34.pdf.b64`.

`BINARIES.sha256` lists the sha256 of each original file. To restore those originals, run:

```
python3 restore_binaries.py
```

That script uses only the Python standard library. It decodes the base64 files, writes the original paths, and checks each sha256.
