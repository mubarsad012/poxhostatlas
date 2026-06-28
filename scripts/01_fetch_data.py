#!/usr/bin/env python3
"""Fetch GSE278320 count files and assemble raw count/metadata tables."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import GEOparse
import pandas as pd
import requests


ACCESSION = "GSE278320"
REPO_ROOT = Path(__file__).resolve().parents[1]
COUNT_DIR = REPO_ROOT / "data" / "raw" / "count_files"
METADATA_DIR = REPO_ROOT / "data" / "raw" / "geo_metadata"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
SOFT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278320/soft/"
    "GSE278320_family.soft.gz"
)


@dataclass(frozen=True)
class PrimarySample:
    sample_id: str
    infection: str
    replicate: int
    count_file: str

    @property
    def url(self) -> str:
        return (
            "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8544nnn/"
            f"{self.sample_id}/suppl/{self.count_file}"
        )


PRIMARY_SAMPLES = [
    PrimarySample("GSM8544817", "mock", 2, "GSM8544817_Par_Mock_total_2_gene_counts-ORIGINAL.txt.gz"),
    PrimarySample("GSM8544823", "mock", 3, "GSM8544823_Par_Mock_total_3_gene_counts-ORIGINAL.txt.gz"),
    PrimarySample("GSM8544828", "mock", 4, "GSM8544828_Par_Mock_total_4_gene_counts-ORIGINAL.txt.gz"),
    PrimarySample("GSM8544820", "VacV", 2, "GSM8544820_Par_WR_total_2_gene_counts-ORIGINAL.txt.gz"),
    PrimarySample("GSM8544826", "VacV", 3, "GSM8544826_Par_WR_total_3_gene_counts-ORIGINAL.txt.gz"),
    PrimarySample("GSM8544830", "VacV", 4, "GSM8544830_Par_WR_total_7_gene_counts-ORIGINAL.txt.gz"),
    PrimarySample("GSM8544836", "VacV", 5, "GSM8544836_Par_WR_total_4_gene_counts-ORIGINAL.txt.gz"),
    PrimarySample("GSM8544842", "VacV", 6, "GSM8544842_Par_WR_total_5_gene_counts-ORIGINAL.txt.gz"),
    PrimarySample("GSM8544848", "VacV", 7, "GSM8544848_Par_WR_total_6_gene_counts-ORIGINAL.txt.gz"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="Re-download existing raw count files.")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in (COUNT_DIR, METADATA_DIR, PROCESSED_DIR):
        path.mkdir(parents=True, exist_ok=True)


def download_url(url: str, destination: Path, force: bool = False) -> None:
    if destination.exists() and destination.stat().st_size > 0 and not force:
        return

    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = destination.with_suffix(destination.suffix + ".tmp")
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with tmp_path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    tmp_path.replace(destination)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_geo_metadata(force: bool = False) -> tuple[pd.DataFrame, dict[str, str]]:
    """Download GEO metadata and return selected sample metadata plus title lookup."""
    soft_path = METADATA_DIR / f"{ACCESSION}_family.soft.gz"
    try:
        # GEOparse provides an independently parsed copy of the GEO family record.
        GEOparse.get_GEO(geo=ACCESSION, destdir=str(METADATA_DIR), silent=True)
    except Exception as exc:  # noqa: BLE001
        print(f"Warning: GEOparse metadata fetch failed; using direct SOFT download. {exc}", file=sys.stderr)

    download_url(SOFT_URL, soft_path, force=force)

    title_by_sample: dict[str, str] = {}
    infection_by_sample: dict[str, str] = {}
    cell_line_by_sample: dict[str, str] = {}
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as handle:
        current_sample: str | None = None
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                current_sample = line.split("=", 1)[1].strip()
            elif current_sample and line.startswith("!Sample_title = "):
                title_by_sample[current_sample] = line.split("=", 1)[1].strip()
            elif current_sample and line.startswith("!Sample_characteristics_ch1 = "):
                characteristic = line.split("=", 1)[1].strip()
                if characteristic.startswith("infection: "):
                    infection_by_sample[current_sample] = characteristic.split(": ", 1)[1]
                elif characteristic.startswith("cell line: "):
                    cell_line_by_sample[current_sample] = characteristic.split(": ", 1)[1]

    metadata_rows = []
    for sample in PRIMARY_SAMPLES:
        metadata_rows.append(
            {
                "sample_id": sample.sample_id,
                "title": title_by_sample.get(sample.sample_id, ""),
                "infection": sample.infection,
                "geo_infection": infection_by_sample.get(sample.sample_id, ""),
                "cell_line": cell_line_by_sample.get(sample.sample_id, ""),
                "fraction": "total",
                "fraction_source": "supplementary_count_filename",
                "replicate": sample.replicate,
            }
        )

    metadata = pd.DataFrame(metadata_rows).set_index("sample_id")
    return metadata, title_by_sample


def infer_fraction_from_title(title: str) -> str:
    lowered = title.lower()
    if "total" in lowered:
        return "total"
    if "polysome" in lowered:
        return "polysome"
    if "80s" in lowered or "ribosome" in lowered:
        return "80s"
    return "unknown"


def infer_fraction_from_filename(filename: str) -> str:
    lowered = filename.lower()
    if "_total" in lowered:
        return "total"
    if "dipol" in lowered or "polysome" in lowered:
        return "polysome"
    if "_80s" in lowered:
        return "80s"
    return "unknown"


def read_count_file(path: Path, sample_id: str) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["gene_id", "gene_symbol", sample_id],
        dtype={"gene_id": "string", "gene_symbol": "string", sample_id: "int64"},
    )
    if frame["gene_id"].duplicated().any():
        duplicated = frame.loc[frame["gene_id"].duplicated(), "gene_id"].head().tolist()
        raise ValueError(f"Duplicate gene IDs in {path.name}: {duplicated}")
    return frame


def build_outputs(force: bool = False) -> None:
    ensure_dirs()
    metadata, title_by_sample = fetch_geo_metadata(force=force)

    manifest_rows = []
    count_frames: list[pd.DataFrame] = []
    for sample in PRIMARY_SAMPLES:
        destination = COUNT_DIR / sample.count_file
        print(f"Downloading/verifying {sample.sample_id}: {sample.count_file}")
        download_url(sample.url, destination, force=force)

        title = title_by_sample.get(sample.sample_id, "")
        filename_fraction = infer_fraction_from_filename(sample.count_file)
        title_fraction = infer_fraction_from_title(title)
        metadata_warning = ""
        if title_fraction != "unknown" and title_fraction != filename_fraction:
            metadata_warning = f"title_fraction={title_fraction};filename_fraction={filename_fraction}"

        manifest_rows.append(
            {
                "sample_id": sample.sample_id,
                "infection": sample.infection,
                "replicate": sample.replicate,
                "title": title,
                "count_filename": sample.count_file,
                "count_url": sample.url,
                "file_size_bytes": destination.stat().st_size,
                "sha256": sha256_file(destination),
                "fraction_from_title": title_fraction,
                "fraction_from_filename": filename_fraction,
                "metadata_warning": metadata_warning,
            }
        )
        count_frames.append(read_count_file(destination, sample.sample_id))

    counts = count_frames[0]
    for frame in count_frames[1:]:
        counts = counts.merge(frame, on=["gene_id", "gene_symbol"], how="outer", validate="one_to_one")

    sample_columns = [sample.sample_id for sample in PRIMARY_SAMPLES]
    counts[sample_columns] = counts[sample_columns].fillna(0).astype("int64")
    counts = counts[["gene_id", "gene_symbol", *sample_columns]].sort_values("gene_id")

    manifest = pd.DataFrame(manifest_rows)
    metadata = metadata.loc[sample_columns]

    counts.to_csv(PROCESSED_DIR / "counts.csv", index=False)
    metadata.to_csv(PROCESSED_DIR / "metadata.csv")
    manifest.to_csv(PROCESSED_DIR / "sample_manifest.csv", index=False)
    (METADATA_DIR / "selected_samples.json").write_text(
        json.dumps([row.__dict__ for row in PRIMARY_SAMPLES], indent=2) + "\n",
        encoding="utf-8",
    )

    expected_counts = {"mock": 3, "VacV": 6}
    actual_counts = metadata["infection"].value_counts().to_dict()
    if actual_counts != expected_counts:
        raise ValueError(f"Unexpected infection counts: {actual_counts}; expected {expected_counts}")

    print(f"Wrote {PROCESSED_DIR / 'counts.csv'} with {counts.shape[0]} genes and {len(sample_columns)} samples.")
    warnings = manifest["metadata_warning"].replace("", pd.NA).dropna()
    if not warnings.empty:
        print("Metadata warnings recorded in sample_manifest.csv:")
        for sample_id, warning in zip(manifest.loc[warnings.index, "sample_id"], warnings, strict=True):
            print(f"  {sample_id}: {warning}")


def main() -> None:
    args = parse_args()
    try:
        build_outputs(force=args.force)
    except (requests.RequestException, OSError, ValueError) as exc:
        raise SystemExit(f"Data acquisition failed: {exc}") from exc


if __name__ == "__main__":
    main()
