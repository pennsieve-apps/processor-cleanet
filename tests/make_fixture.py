#!/usr/bin/env python3
"""
Stand-in for the fcs-prep step: turn raw FCS files (the Zenodo 15388817 demo data) into
what fcs-prep writes, so this app can be run before fcs-prep exists.

    <out>/fcs/<sample_id>.fcs          FCS 3.1, float32; marker channels arcsinh(x / cofactor),
                                       scatter / Time / other channels as acquired
    <out>/metadata/samples.csv         sample_id, file, batch, condition, n_events (+ extra sheet columns)
    <out>/metadata/panel.csv           channel, marker, type, cofactor, clustering, state

It follows docs/conventions.md of workflow-immune-cytometry: channels missing from some
files are dropped from all (the demo's BV785-A / CD28), marker = $PnS, cofactor 6000.
--events N keeps N events per file, chosen at random (seeded) but kept in acquisition
order, so PeacoQC still sees the time course; --events 0 keeps every event.

    python3 tests/make_fixture.py <demo dir> <out dir> [--events 20000] [--files 4]

Needs flowio and numpy (pip install -r requirements-dev.txt).
"""
import argparse
import csv
import os
import re
import sys

import flowio
import numpy as np

SCATTER = re.compile(r"^(FSC|SSC)", re.I)
TIME = re.compile(r"^time$", re.I)
# Functional / activation markers of the demo panel (differential-state candidates);
# everything else with a $PnS is treated as a lineage marker.
STATE = {"CD25", "CD69", "PD1", "HLADR", "CD95", "TIGIT", "CD39", "GZMB", "PFN", "GNLY", "CD27", "CD28"}
KEEP_KEYWORDS = ("date", "btim", "etim", "cyt", "cytsn", "timestep", "op", "inst", "src", "fil")


def channel_type(pnn, pns):
    if TIME.match(pnn):
        return "time"
    if SCATTER.match(pnn):
        return "scatter"
    return "marker" if pns else "other"


def read(path):
    fd = flowio.FlowData(path)
    return fd, np.reshape(np.asarray(fd.events, dtype=np.float64), (-1, fd.channel_count))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source")
    ap.add_argument("out")
    ap.add_argument("--events", type=int, default=20000, help="events per file (0 = all)")
    ap.add_argument("--files", type=int, default=0, help="only the first N files of the sample sheet (0 = all)")
    ap.add_argument("--cofactor", type=float, default=6000.0)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args(argv)

    sheet_path = os.path.join(a.source, "samples.csv")
    if os.path.exists(sheet_path):
        with open(sheet_path, newline="") as fh:
            sheet = list(csv.DictReader(fh))
    else:  # no sample sheet: as fcs-prep does, sample_id from the file name, batch = $DATE
        sheet = []
        for f in sorted(os.listdir(a.source)):
            if f.lower().endswith(".fcs") and not f.startswith("."):
                fd = flowio.FlowData(os.path.join(a.source, f), only_text=True)
                sheet.append({"file": f, "sample_id": re.sub(r"[^A-Za-z0-9_-]+", "_", f[:-4]).strip("_"),
                              "batch": fd.text.get("date", "NA"), "condition": "NA"})
    sheet.sort(key=lambda r: r["sample_id"])
    if a.files:
        sheet = sheet[: a.files]

    # Panel = channels present in every file.
    headers = {}
    for row in sheet:
        fd = flowio.FlowData(os.path.join(a.source, row["file"]), only_text=True)
        headers[row["file"]] = list(zip(fd.pnn_labels, [s or "" for s in fd.pns_labels]))
    common = set.intersection(*(set(n for n, _ in h) for h in headers.values()))
    first = headers[sheet[0]["file"]]
    panel = [(n, s) for n, s in first if n in common]
    dropped = sorted(set(n for h in headers.values() for n, _ in h) - common)

    os.makedirs(os.path.join(a.out, "fcs"), exist_ok=True)
    os.makedirs(os.path.join(a.out, "metadata"), exist_ok=True)
    rng = np.random.default_rng(a.seed)
    samples = []
    for row in sheet:
        fd, ev = read(os.path.join(a.source, row["file"]))
        names = fd.pnn_labels
        idx = [names.index(n) for n, _ in panel]
        ev = ev[:, idx]
        if a.events and ev.shape[0] > a.events:
            keep = np.sort(rng.choice(ev.shape[0], a.events, replace=False))
            ev = ev[keep]
        types = [channel_type(n, s) for n, s in panel]
        for j, t in enumerate(types):
            if t == "marker":
                ev[:, j] = np.arcsinh(ev[:, j] / a.cofactor)
        meta = {k: v for k, v in fd.text.items() if k in KEEP_KEYWORDS}
        meta["fil"] = row["sample_id"] + ".fcs"
        meta["PENNSIEVE_PREP"] = "transform=arcsinh;cofactor=%g;source=make_fixture" % a.cofactor
        for j, (n, _) in enumerate(panel, start=1):
            src = names.index(n) + 1
            meta["p%dr" % j] = "8" if types[j - 1] == "marker" else fd.text.get("p%dr" % src, "262144")
        out = os.path.join(a.out, "fcs", row["sample_id"] + ".fcs")
        with open(out, "wb") as fh:
            flowio.create_fcs(fh, ev.astype(np.float32).ravel(),
                              [n for n, _ in panel], [s for _, s in panel], metadata_dict=meta)
        s = dict(row)
        s["n_events"] = str(ev.shape[0])
        samples.append(s)
        print("%s: %d events" % (row["sample_id"], ev.shape[0]), flush=True)

    cols = ["sample_id", "file", "batch", "condition", "n_events"]
    cols += [c for c in samples[0] if c not in cols]
    with open(os.path.join(a.out, "metadata", "samples.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(samples)
    with open(os.path.join(a.out, "metadata", "panel.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["channel", "marker", "type", "cofactor", "clustering", "state"])
        for n, s in panel:
            t = channel_type(n, s)
            marker = s or n
            w.writerow([n, marker, t, "%g" % a.cofactor if t == "marker" else "",
                        str(t == "marker" and marker not in STATE).lower(), str(t == "marker" and marker in STATE).lower()])
    print("panel: %d channels; dropped (not in every file): %s" % (len(panel), ", ".join(dropped) or "none"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
