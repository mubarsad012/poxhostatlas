#!/usr/bin/env python3
"""Pack the full project into a single zip with at most 10 top-level files."""

from __future__ import annotations

import hashlib
import shutil
import subprocess
import tarfile
import zipfile
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
STAGING = REPO_ROOT / "build" / "project_10files_staging"
OUTPUT_ZIP = REPO_ROOT / "poxvirus_transcriptomics_project_10files.zip"
DESKTOP_ZIP = REPO_ROOT.parent / "poxvirus_transcriptomics_project_10files.zip"

EXCLUDE_DIR_NAMES = {".git", "venv", "__pycache__", ".mplcache", "build"}
EXCLUDE_FILE_NAMES = {".DS_Store"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def should_skip(path: Path) -> bool:
    return any(part in EXCLUDE_DIR_NAMES for part in path.parts) or path.name in EXCLUDE_FILE_NAMES


def make_tar(source_dir: Path, archive_path: Path, arc_prefix: str) -> None:
    with tarfile.open(archive_path, "w:gz") as tar:
        if not source_dir.exists():
            return
        for item in sorted(source_dir.rglob("*")):
            if should_skip(item) or not item.is_file():
                continue
            tar.add(item, arcname=str(Path(arc_prefix) / item.relative_to(source_dir)))


def make_root_config_tar(archive_path: Path) -> None:
    root_files = ["README.md", "run_all.sh", "requirements.txt", "requirements-lock.txt", ".gitignore"]
    with tarfile.open(archive_path, "w:gz") as tar:
        for name in root_files:
            path = REPO_ROOT / name
            if path.exists():
                tar.add(path, arcname=f"project_root/{name}")


def write_readme(path: Path) -> None:
    path.write_text(
        f"""Poxvirus Transcriptomics — Full Project Bundle (10 files)
Date: {date.today().isoformat()}

This zip contains the complete project in 10 top-level files (excluding venv).
Upload or extract anywhere, then unpack each .tar.gz into a single folder.

CONTENTS
  1. README.txt                      — this file
  2. manuscript_overleaf.zip         — Overleaf-ready LaTeX + figures
  3. scripts.tar.gz                  — all pipeline scripts (scripts/)
  4. docs.tar.gz                     — documentation and manuscript sources
  5. data.tar.gz                     — raw, processed, and external data
  6. results.tar.gz                  — figures, tables, meta-analysis outputs
  7. release.tar.gz                  — manuscript-ready release package
  8. project_root.tar.gz             — README.md, run_all.sh, requirements
  9. MANIFEST.txt                    — file listing with sizes
 10. SHA256SUMS.txt                   — checksums for integrity verification

QUICK START (after extract)
  mkdir poxvirus-transcriptomics && cd poxvirus-transcriptomics
  tar -xzf project_root.tar.gz --strip-components=1
  tar -xzf scripts.tar.gz
  tar -xzf docs.tar.gz
  tar -xzf data.tar.gz
  tar -xzf results.tar.gz
  tar -xzf release.tar.gz
  python3 -m venv venv && source venv/bin/activate
  pip install -r requirements.txt
  bash run_all.sh

OVERLEAF
  Upload manuscript_overleaf.zip directly to Overleaf (separate from this bundle).

NOTE
  The Python virtual environment (venv/) is excluded; recreate locally with pip.
""",
        encoding="utf-8",
    )


def main() -> None:
    if STAGING.exists():
        shutil.rmtree(STAGING)
    STAGING.mkdir(parents=True)

    bundles = {
        "scripts.tar.gz": lambda p: make_tar(REPO_ROOT / "scripts", p, "scripts"),
        "docs.tar.gz": lambda p: make_tar(REPO_ROOT / "docs", p, "docs"),
        "data.tar.gz": lambda p: make_tar(REPO_ROOT / "data", p, "data"),
        "results.tar.gz": lambda p: make_tar(REPO_ROOT / "results", p, "results"),
        "release.tar.gz": lambda p: make_tar(REPO_ROOT / "release", p, "release"),
        "project_root.tar.gz": make_root_config_tar,
    }

    write_readme(STAGING / "README.txt")

    overleaf_src = REPO_ROOT / "docs" / "manuscript" / "poxvirus_manuscript_overleaf.zip"
    if overleaf_src.exists():
        shutil.copy2(overleaf_src, STAGING / "manuscript_overleaf.zip")

    for name, builder in bundles.items():
        print(f"Creating {name}...")
        builder(STAGING / name)

    top_level = sorted(f for f in STAGING.iterdir() if f.is_file())
    manifest_lines = ["Top-level bundle files", "=" * 40, ""]
    for path in top_level:
        manifest_lines.append(f"{path.name}\t{path.stat().st_size:,} bytes")
    (STAGING / "MANIFEST.txt").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    hash_lines = []
    for path in sorted(STAGING.iterdir()):
        if path.is_file():
            hash_lines.append(f"{sha256_file(path)}  {path.name}")
    (STAGING / "SHA256SUMS.txt").write_text("\n".join(hash_lines) + "\n", encoding="utf-8")

    final_files = sorted(f.name for f in STAGING.iterdir() if f.is_file())
    if len(final_files) != 10:
        raise SystemExit(f"Expected exactly 10 top-level files, got {len(final_files)}: {final_files}")

    if OUTPUT_ZIP.exists():
        OUTPUT_ZIP.unlink()
    with zipfile.ZipFile(OUTPUT_ZIP, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name in final_files:
            zf.write(STAGING / name, arcname=name)

    shutil.copy2(OUTPUT_ZIP, DESKTOP_ZIP)
    print(f"\nWrote: {OUTPUT_ZIP}")
    print(f"Copied: {DESKTOP_ZIP}")
    print(f"Top-level files ({len(final_files)}):")
    for name in final_files:
        print(f"  - {name}")
    print(f"Total zip size: {OUTPUT_ZIP.stat().st_size / (1024 * 1024):.1f} MB")


if __name__ == "__main__":
    main()
