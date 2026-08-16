#!/usr/bin/env python3
"""Audit the bundled 05_fibroblast notebooks and source data without executing them."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="List notebooks, kernels, code-cell counts, and other bundled files."
    )
    parser.add_argument("source", type=Path, help="Path to the upstream 05_fibroblast folder")
    parser.add_argument("--json", action="store_true", help="Emit JSON instead of a table")
    args = parser.parse_args()

    if not args.source.is_dir():
        parser.error(f"not a directory: {args.source}")

    records = []
    failures = []
    for path in sorted(args.source.iterdir()):
        record = {"name": path.name, "bytes": path.stat().st_size, "kind": "data"}
        if path.suffix == ".ipynb":
            record["kind"] = "notebook"
            try:
                notebook = json.loads(path.read_text(encoding="utf-8"))
                metadata = notebook.get("metadata", {})
                record["kernel"] = (
                    metadata.get("kernelspec", {}).get("language")
                    or metadata.get("language_info", {}).get("name")
                    or "unknown"
                )
                record["code_cells"] = sum(
                    cell.get("cell_type") == "code" for cell in notebook.get("cells", [])
                )
            except (OSError, UnicodeError, json.JSONDecodeError) as exc:
                record["error"] = str(exc)
                failures.append(path.name)
        records.append(record)

    if args.json:
        print(json.dumps(records, indent=2))
    else:
        print("name\tkind\tkernel\tcode_cells\tbytes")
        for item in records:
            print(
                f"{item['name']}\t{item['kind']}\t{item.get('kernel', '-')}\t"
                f"{item.get('code_cells', '-')}\t{item['bytes']}"
            )

    notebooks = sum(item["kind"] == "notebook" for item in records)
    data_files = sum(item["kind"] == "data" for item in records)
    if not args.json:
        print(f"summary\tnotebooks={notebooks}\tdata_files={data_files}\tfailures={len(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
