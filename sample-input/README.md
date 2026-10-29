# Sample input

No FCS files are committed. The app reads what `fcs-prep` writes:

```
INPUT_DIR/
  fcs/<sample_id>.fcs          FCS 3.1, float32, marker channels arcsinh-transformed
  metadata/panel.csv           channel, marker, type (marker|scatter|time|other), cofactor, ...
  metadata/samples.csv         (not used by this app)
```

## Demo data

The workflow's demo data are the CytoBatchFlagR spectral flow files (Koenig Z, Andorf S,
Tilburgs T, Eswar S. *CytoBatchFlagR: Spectral flow cytometry dataset.* Zenodo 2025,
doi:10.5281/zenodo.15388817, CC BY 4.0): 24 FCS files, Cytek Aurora, PBMC pre-gated in
FlowJo to CD45+CD14− live single cells, 28 markers, 25k–600k events per file.

Download them (1.46 GB) into a folder, e.g. `../workflow-immune-cytometry/data/zenodo-15388817/`,
and put a sample sheet `samples.csv` (`file, sample_id, batch, condition`, any extra columns) next to them; without one the sample id is the file name and the batch its `$DATE`:

```sh
mkdir -p ../workflow-immune-cytometry/data/zenodo-15388817 && cd $_
curl -sL https://zenodo.org/api/records/15388817 | python3 -c \
  'import json,sys; [print(f["links"]["self"], f["key"]) for f in json.load(sys.stdin)["files"]]' |
  while read url name; do curl -sL -o "$name" "$url"; done
```

Then either run `fcs-prep` over them, or use the stand-in that ships with this repo:

```sh
pip install -r requirements-dev.txt
make fixture              # 20,000 events per file -> data/fixture/{fcs,metadata}
make fixture EVENTS=0     # every event (1.4 GB)
make run                  # -> data/output/{fcs,qc}
```

`tests/make_fixture.py` does what `fcs-prep` does by contract: keeps the channels present
in every file (drops `BV785-A`/CD28, missing from batches 05–06), arcsinh-transforms the
channels that have a `$PnS` with cofactor 6000, leaves scatter, `AF-A` and `Time` as
acquired, and writes `panel.csv` and `samples.csv`. Subsampled events stay in acquisition
order. For Cleanet the usual input is peacoqc's output (`fcs/` from `peacoqc`, `metadata/` from `fcs-prep`); the stand-in works too.
