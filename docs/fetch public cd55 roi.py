#!/usr/bin/env python3
"""
fetch_public_cd55_roi.py
========================
Slice public CpG methylation datasets down to the CD55 region of interest and
emit one small JSON file that docs/public-methylation.html can load directly.

Why it works without downloading whole genomes
----------------------------------------------
bigWig and bigBed are indexed, range-request-friendly formats. The UCSC command
line tools read them straight over HTTPS and pull only the bytes covering the
requested interval. A 2 GB ENCODE methylation file costs a few hundred KB of
traffic when you ask for 6.5 kb of chr1.

Requirements
------------
    conda create -n cd55roi -c bioconda -c conda-forge \
        ucsc-bigwigtobedgraph ucsc-bigbedtobed ucsc-liftover python=3.11
    conda activate cd55roi

Usage
-----
    python fetch_public_cd55_roi.py --out docs/data/public_cd55_roi.json
    python fetch_public_cd55_roi.py --only encode --min-cov 5 -v
    python fetch_public_cd55_roi.py --list          # show configured sources

Run it on Skynet, not on a laptop — nothing here is heavy, but the network path
to ENCODE/UCSC is far more reliable from the lab.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
from datetime import datetime, timezone

# --------------------------------------------------------------------------
# Region of interest. hg38 is the working assembly: almost every public
# methylome is aligned to it, and it is the assembly the site page quotes.
# Coordinates below are 1-based inclusive (browser style); BED/bigWig slicing
# converts to 0-based half-open internally.
# --------------------------------------------------------------------------
ROI = {
    "assembly": "hg38",
    "chrom": "chr1",
    "start": 207318058,   # 1-based inclusive
    "end":   207324558,
    "expected_cpgs": 137,
    "t2t_equivalent": "chr1:206583354-206589854",  # T2T-CHM13v2.0, for reference only
    # NOTE ON THE EDGE BASE: written 1-based inclusive, this span is 6,501 bp,
    # while the ROI is described everywhere as 6,500 bp. That means the numbers
    # were most likely recorded BED-style (0-based start). One base at the left
    # edge is at stake. Set BED_STYLE_START=True if the coordinates came
    # straight out of a BED file; the script reports how many CpGs it actually
    # found either way, so compare against expected_cpgs to settle it.
}

BED_STYLE_START = False

USER_AGENT = "epicme-cd55-roi/1.0 (research use; contact epicme.bio)"

# --------------------------------------------------------------------------
# Sources. Add rows here as new datasets are found — the JSON schema and the
# viewer do not need to change.
#
# kind:
#   "encode"  -> resolve files through the ENCODE REST API by experiment accession
#   "bigwig"  -> a direct bigWig URL carrying a methylation fraction (0-1 or 0-100)
#   "hub"     -> a UCSC track hub; every bigWig in trackDb is sliced
# assembly:
#   the assembly the FILE is in. Anything that is not hg38 gets lifted over.
# --------------------------------------------------------------------------
SOURCES = [
    {
        "id": "encode_k562_wgbs",
        "kind": "encode",
        "accession": "ENCSR765JPC",
        "label": "K562 · WGBS (ENCODE / HAIB)",
        "cell": "K562",
        "cell_group": "k562",
        "assay": "WGBS",
        "assembly": "hg38",
    },
    {
        "id": "encode_tcell_wgbs",
        "kind": "encode",
        "accession": "ENCSR663MXB",
        "label": "Primary CD3+ T cells · WGBS (ENCODE / Roadmap)",
        "cell": "Primary T cells (CD3+)",
        "cell_group": "tcell",
        "assay": "WGBS",
        "assembly": "hg38",
    },
    {
        "id": "encode_na12878_wgbs",
        "kind": "encode",
        "accession": "ENCSR890UQO",
        "label": "NA12878 · WGBS (ENCODE)",
        "cell": "NA12878 lymphoblastoid",
        "cell_group": "other",
        "assay": "WGBS",
        "assembly": "hg38",
    },
    {
        "id": "atlas_hub",
        "kind": "hub",
        "url": "https://files.cs.huji.ac.il/tommy/Meth_Atlas_hg38/hub.txt",
        "label": "Human Methylation Atlas (Loyfer 2023)",
        "cell_group": "tcell",
        "assay": "WGBS",
        "assembly": "hg38",
        # Keep the payload sane: only pull tracks whose name looks immune.
        # Drop `track_filter` entirely to take all 39 cell types.
        "track_filter": r"(?i)(t.?cell|cd4|cd8|naive|treg|b.?cell|nk|monocyte|granulo|macroph)",
    },
    # ---- hg19 sources: lifted over automatically, see --chain ----
    # Roadmap ships fraction and coverage as separate bigWigs per epigenome.
    # Uncomment and set the EIDs you want after checking which ones have WGBS at
    # https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/
    # {
    #     "id": "roadmap_E039",
    #     "kind": "bigwig",
    #     "url": "https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/"
    #            "FractionalMethylation_bigwig/E039_WGBS_FractionalMethylation.bigwig",
    #     "coverage_url": "https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/"
    #                     "ReadCoverage_bigwig/E039_WGBS_ReadCoverage.bigwig",
    #     "label": "E039 primary T helper naive · WGBS (Roadmap)",
    #     "cell": "Primary T helper naive cells",
    #     "cell_group": "tcell",
    #     "assay": "WGBS",
    #     "assembly": "hg19",
    # },
]

ENCODE_BASE = "https://www.encodeproject.org"


# ==========================================================================
# helpers
# ==========================================================================
def log(msg: str, verbose: bool = True) -> None:
    if verbose:
        print(f"  {msg}", file=sys.stderr)


def need(tool: str) -> str:
    path = shutil.which(tool)
    if not path:
        sys.exit(
            f"error: `{tool}` not found on PATH.\n"
            f"       install with: conda install -c bioconda ucsc-{tool.lower()}"
        )
    return path


def get_json(url: str) -> dict:
    req = urllib.request.Request(
        url, headers={"Accept": "application/json", "User-Agent": USER_AGENT}
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.load(resp)


def get_text(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read().decode("utf-8", errors="replace")


def run(cmd: list[str]) -> str:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed: {proc.stderr.strip()[:400]}")
    return proc.stdout


# ==========================================================================
# coordinate handling
# ==========================================================================
def roi_zero_based() -> tuple[str, int, int]:
    """ROI as 0-based half-open, which is what BED/bigWig tools want."""
    start = ROI["start"] if BED_STYLE_START else ROI["start"] - 1
    return ROI["chrom"], start, ROI["end"]


def lift_roi(chain: str | None, verbose: bool) -> tuple[str, int, int]:
    """Convert the hg38 ROI to hg19 with UCSC liftOver. No hard-coded numbers."""
    if not chain:
        sys.exit(
            "error: an hg19 source is enabled but --chain was not given.\n"
            "       get it with:\n"
            "         wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/"
            "liftOver/hg38ToHg19.over.chain.gz"
        )
    need("liftOver")
    chrom, start, end = roi_zero_based()
    with tempfile.TemporaryDirectory() as tmp:
        src = os.path.join(tmp, "roi.hg38.bed")
        dst = os.path.join(tmp, "roi.hg19.bed")
        unmapped = os.path.join(tmp, "unmapped.bed")
        with open(src, "w") as fh:
            fh.write(f"{chrom}\t{start}\t{end}\tCD55_ROI\n")
        run(["liftOver", src, chain, dst, unmapped])
        with open(dst) as fh:
            line = fh.readline().strip()
        if not line:
            sys.exit("error: the ROI did not lift over to hg19 — check the chain file.")
        f = line.split("\t")
        lifted = (f[0], int(f[1]), int(f[2]))
    log(f"ROI lifted to hg19: {lifted[0]}:{lifted[1] + 1}-{lifted[2]}", verbose)
    return lifted


# ==========================================================================
# slicing
# ==========================================================================
def slice_bigwig(url: str, chrom: str, start: int, end: int) -> dict[int, float]:
    """Remote bigWig -> {0-based position: value}. Only the ROI bytes are fetched."""
    need("bigWigToBedGraph")
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "roi.bedGraph")
        run([
            "bigWigToBedGraph",
            f"-chrom={chrom}", f"-start={start}", f"-end={end}",
            url, out,
        ])
        values: dict[int, float] = {}
        with open(out) as fh:
            for line in fh:
                f = line.split()
                if len(f) < 4:
                    continue
                a, b, v = int(f[1]), int(f[2]), float(f[3])
                for pos in range(a, b):
                    values[pos] = v
    return values


def slice_bigbed(url: str, chrom: str, start: int, end: int) -> list[list[str]]:
    """Remote bigBed (bedMethyl) -> list of fields for rows inside the ROI."""
    need("bigBedToBed")
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "roi.bed")
        run([
            "bigBedToBed",
            f"-chrom={chrom}", f"-start={start}", f"-end={end}",
            url, out,
        ])
        rows = []
        with open(out) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 11:
                    rows.append(f)
    return rows


def parse_bedmethyl(rows: list[list[str]], min_cov: int) -> list[dict]:
    """
    bedMethyl: chrom start end name score strand thickStart thickEnd rgb cov pct
    Strands are collapsed to CpG level: the plus-strand C and the minus-strand G
    of the same dinucleotide are pooled by read count.
    """
    by_cpg: dict[int, dict] = {}
    for f in rows:
        start, strand = int(f[1]), f[5]
        try:
            cov, pct = int(f[9]), float(f[10])
        except ValueError:
            continue
        # anchor both strands on the C of the plus strand
        key = start if strand != "-" else start - 1
        rec = by_cpg.setdefault(key, {"meth_reads": 0.0, "cov": 0})
        rec["meth_reads"] += cov * pct / 100.0
        rec["cov"] += cov

    out = []
    for pos in sorted(by_cpg):
        rec = by_cpg[pos]
        if rec["cov"] < min_cov:
            continue
        out.append({
            "pos": pos + 1,                                    # back to 1-based
            "meth": round(rec["meth_reads"] / rec["cov"], 4),
            "cov": rec["cov"],
        })
    return out


# ==========================================================================
# source handlers
# ==========================================================================
def encode_pick_files(accession: str, verbose: bool) -> list[dict]:
    """
    Ask the ENCODE API which released GRCh38 files carry per-CpG methylation.
    Preference: bigBed 'methylation state at CpG' (has coverage), then bigWig.
    """
    meta = get_json(f"{ENCODE_BASE}/experiments/{accession}/?format=json")
    picked = []
    for f in meta.get("files", []):
        if f.get("status") != "released":
            continue
        if f.get("assembly") != "GRCh38":
            continue
        otype = (f.get("output_type") or "").lower()
        ftype = (f.get("file_type") or "").lower()
        if "methylation state at cpg" in otype and "bigbed" in ftype:
            picked.append({
                "accession": f.get("accession"),
                "url": ENCODE_BASE + f["href"],
                "flavour": "bigbed",
                "biological_replicates": f.get("biological_replicates"),
            })
    if not picked:
        for f in meta.get("files", []):
            if f.get("status") != "released" or f.get("assembly") != "GRCh38":
                continue
            otype = (f.get("output_type") or "").lower()
            if "methylation state at cpg" in otype and "bigwig" in (f.get("file_type") or ""):
                picked.append({
                    "accession": f.get("accession"),
                    "url": ENCODE_BASE + f["href"],
                    "flavour": "bigwig",
                    "biological_replicates": f.get("biological_replicates"),
                })
    log(f"{accession}: {len(picked)} candidate file(s)", verbose)
    return picked


def handle_encode(src: dict, min_cov: int, verbose: bool) -> list[dict]:
    chrom, start, end = roi_zero_based()
    out = []
    for f in encode_pick_files(src["accession"], verbose):
        log(f"  slicing {f['accession']} ({f['flavour']})", verbose)
        if f["flavour"] == "bigbed":
            cpgs = parse_bedmethyl(slice_bigbed(f["url"], chrom, start, end), min_cov)
        else:
            vals = slice_bigwig(f["url"], chrom, start, end)
            cpgs = [{"pos": p + 1, "meth": round(v / 100.0 if v > 1 else v, 4), "cov": None}
                    for p, v in sorted(vals.items())]
        reps = f.get("biological_replicates") or []
        suffix = f" rep{','.join(map(str, reps))}" if reps else ""
        out.append(make_dataset(
            src, cpgs,
            dataset_id=f"{src['id']}_{f['accession']}",
            label=src["label"] + suffix,
            provenance=f"ENCODE {src['accession']} / {f['accession']}",
            file_url=f["url"],
        ))
    return out


def handle_bigwig(src: dict, min_cov: int, chain: str | None, verbose: bool) -> list[dict]:
    if src.get("assembly", "hg38") == "hg38":
        chrom, start, end = roi_zero_based()
    else:
        chrom, start, end = lift_roi(chain, verbose)

    frac = slice_bigwig(src["url"], chrom, start, end)
    cov = slice_bigwig(src["coverage_url"], chrom, start, end) if src.get("coverage_url") else {}

    cpgs = []
    for pos in sorted(frac):
        c = int(cov.get(pos, 0)) if cov else None
        if cov and c < min_cov:
            continue
        v = frac[pos]
        cpgs.append({"pos": pos + 1, "meth": round(v / 100.0 if v > 1 else v, 4), "cov": c})

    return [make_dataset(
        src, cpgs,
        dataset_id=src["id"],
        label=src["label"],
        provenance=src.get("provenance", src["url"]),
        file_url=src["url"],
        note=None if src.get("assembly") == "hg38"
             else "positions lifted from hg19 to the hg38 ROI window",
    )]


def handle_hub(src: dict, min_cov: int, verbose: bool) -> list[dict]:
    """Parse a UCSC track hub and slice every bigWig it declares."""
    hub_url = src["url"]
    base = hub_url.rsplit("/", 1)[0] + "/"
    hub = get_text(hub_url)

    genomes_rel = re.search(r"^genomesFile\s+(\S+)", hub, re.M)
    genomes = get_text(base + genomes_rel.group(1)) if genomes_rel else ""
    trackdb_rel = re.search(r"^trackDb\s+(\S+)", genomes, re.M)
    if not trackdb_rel:
        log("hub: no trackDb found, skipping", verbose)
        return []
    trackdb_path = trackdb_rel.group(1)
    trackdb = get_text(base + trackdb_path)
    tdb_base = (base + trackdb_path).rsplit("/", 1)[0] + "/"

    pattern = src.get("track_filter")
    rx = re.compile(pattern) if pattern else None

    chrom, start, end = roi_zero_based()
    results = []
    blocks = re.split(r"\n\s*\n", trackdb)
    for block in blocks:
        name = re.search(r"^\s*track\s+(\S+)", block, re.M)
        data = re.search(r"^\s*bigDataUrl\s+(\S+)", block, re.M)
        if not name or not data:
            continue
        if not data.group(1).endswith((".bw", ".bigWig", ".bigwig")):
            continue
        short = re.search(r"^\s*shortLabel\s+(.+)$", block, re.M)
        label = short.group(1).strip() if short else name.group(1)
        if rx and not rx.search(label + " " + name.group(1)):
            continue
        url = data.group(1)
        if not url.startswith("http"):
            url = tdb_base + url
        log(f"  hub track: {label}", verbose)
        try:
            vals = slice_bigwig(url, chrom, start, end)
        except RuntimeError as exc:
            log(f"    skipped ({exc})", verbose)
            continue
        cpgs = [{"pos": p + 1, "meth": round(v / 100.0 if v > 1 else v, 4), "cov": None}
                for p, v in sorted(vals.items())]
        results.append(make_dataset(
            src, cpgs,
            dataset_id=f"{src['id']}_{name.group(1)}",
            label=label,
            provenance=f"Methylation Atlas hub / {name.group(1)}",
            file_url=url,
            cell=label,
        ))
    return results


# ==========================================================================
# assembly
# ==========================================================================
def make_dataset(src, cpgs, dataset_id, label, provenance, file_url,
                 cell=None, note=None) -> dict:
    meths = [c["meth"] for c in cpgs if c["meth"] is not None]
    covs = [c["cov"] for c in cpgs if c.get("cov")]
    return {
        "id": dataset_id,
        "label": label,
        "cell": cell or src.get("cell", ""),
        "cell_group": src.get("cell_group", "other"),
        "assay": src.get("assay", "WGBS"),
        "assembly": src.get("assembly", "hg38"),
        "provenance": provenance,
        "file_url": file_url,
        "note": note,
        "summary": {
            "n_cpg": len(cpgs),
            "mean_meth": round(statistics.fmean(meths), 4) if meths else None,
            "median_cov": int(statistics.median(covs)) if covs else None,
        },
        "cpgs": cpgs,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="docs/data/public_cd55_roi.json")
    ap.add_argument("--min-cov", type=int, default=5,
                    help="drop CpGs below this read depth (default 5); "
                         "sources without a coverage track are kept as-is and flagged")
    ap.add_argument("--only", action="append", default=None,
                    help="restrict to source ids or kinds; repeatable")
    ap.add_argument("--chain", default=None, help="hg38ToHg19.over.chain.gz, for hg19 sources")
    ap.add_argument("--list", action="store_true", help="print configured sources and exit")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    if args.list:
        for s in SOURCES:
            print(f"{s['id']:28s} {s['kind']:8s} {s.get('assembly','hg38'):5s} {s['label']}")
        return 0

    selected = SOURCES
    if args.only:
        wanted = set(args.only)
        selected = [s for s in SOURCES if s["id"] in wanted or s["kind"] in wanted]
        if not selected:
            sys.exit(f"error: nothing matched --only {args.only}")

    datasets, failures = [], []
    for src in selected:
        print(f"[{src['id']}] {src['label']}", file=sys.stderr)
        try:
            if src["kind"] == "encode":
                datasets += handle_encode(src, args.min_cov, args.verbose)
            elif src["kind"] == "bigwig":
                datasets += handle_bigwig(src, args.min_cov, args.chain, args.verbose)
            elif src["kind"] == "hub":
                datasets += handle_hub(src, args.min_cov, args.verbose)
            else:
                failures.append((src["id"], f"unknown kind {src['kind']}"))
        except (urllib.error.URLError, RuntimeError, OSError) as exc:
            failures.append((src["id"], str(exc)[:300]))
            print(f"  FAILED: {exc}", file=sys.stderr)

    payload = {
        "roi": ROI,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "min_cov": args.min_cov,
        "datasets": datasets,
        "failures": [{"source": s, "error": e} for s, e in failures],
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(payload, fh, separators=(",", ":"))

    size_kb = os.path.getsize(args.out) / 1024
    print(f"\nwrote {args.out} — {len(datasets)} dataset(s), {size_kb:.1f} KB",
          file=sys.stderr)
    expected = ROI["expected_cpgs"]
    for d in datasets:
        s = d["summary"]
        mean = f"{s['mean_meth']:.3f}" if s["mean_meth"] is not None else "n/a"
        flag = "" if s["n_cpg"] == expected else f"  <- expected {expected}"
        print(f"  {d['id']:44s} n={s['n_cpg']:4d}  mean={mean}{flag}", file=sys.stderr)

    off = [d["id"] for d in datasets if d["summary"]["n_cpg"] not in (0, expected)]
    if off:
        print(f"\nCpG count differs from {expected} in {len(off)} dataset(s). Usual causes,\n"
              "in order of likelihood: CpGs dropped by --min-cov; the edge base (see\n"
              "BED_STYLE_START in this file); or a genuine reference difference between\n"
              "hg38 and the T2T coordinates the 137 sites were called on.", file=sys.stderr)
    if failures:
        print(f"\n{len(failures)} source(s) failed — see 'failures' in the JSON.",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
