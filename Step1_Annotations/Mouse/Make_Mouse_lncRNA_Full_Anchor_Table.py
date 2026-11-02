#!/usr/bin/env python3
"""
Merge mouse lncRNA genomic coordinates with their upstream/downstream PCG anchors.

Inputs
------
1. outputmousegencodestep1.txt
2. Mouse_lncRNA_PCG_anchors.txt

Output
------
Mouse_lncRNA_Full_Anchor_Table.tsv
"""

from __future__ import annotations

import ast
import csv
import os
import sys
from collections import defaultdict
from typing import Dict, List


# ---------------------------------------------------------------------------
# EDIT THESE FILENAMES IF YOUR FILES HAVE DIFFERENT NAMES
# ---------------------------------------------------------------------------

STEP1_FILE = "outputmousegencodestep1.txt"
ANCHOR_FILE = "Mouse_lncRNA_PCG_anchors.txt"
OUTPUT_FILE = "Mouse_lncRNA_Full_Anchor_Table.tsv"


def empty_anchor() -> Dict[str, object]:
    return {
        "pcg_name": "",
        "distance": "",
        "chromosome": "",
        "start": "",
        "end": "",
        "strand": "",
    }


def parse_anchor(text: str) -> Dict[str, object]:
    text = text.strip()

    if not text or text == "None":
        return empty_anchor()

    try:
        value = ast.literal_eval(text)
    except (SyntaxError, ValueError):
        return empty_anchor()

    if not isinstance(value, tuple) or len(value) < 6:
        return empty_anchor()

    return {
        "pcg_name": str(value[0]),
        "distance": value[1],
        "chromosome": str(value[2]),
        "start": value[3],
        "end": value[4],
        "strand": str(value[5]),
    }


def read_mouse_lncrna_coordinates(path: str) -> Dict[str, List[Dict[str, object]]]:
    """
    Reads lncRNA gene-block coordinate lines from the Step 1 output.

    Expected coordinate-line structure:
        lncRNA_name chromosome start end strand isoform_count

    The parser begins after 'lncRNA dictionary gene blocks' and stops when
    the PCG section begins.
    """

    coordinates: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    inside_lncrna_section = False

    with open(path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            lower = line.lower()

            if lower.startswith("lncrna dictionary gene blocks"):
                inside_lncrna_section = True
                continue

            if inside_lncrna_section and (
                lower.startswith("pcg dictionary gene blocks")
                or lower.startswith("protein coding")
                or lower.startswith("protein-coding")
            ):
                break

            if not inside_lncrna_section:
                continue

            if lower.startswith("total lncrna"):
                continue

            fields = line.split()

            # Ignore name/count summary rows such as:
            # Gm12345 4
            if len(fields) < 5:
                continue

            name = fields[0]
            chromosome = fields[1]

            if not chromosome.startswith("chr"):
                continue

            try:
                start = int(fields[2])
                end = int(fields[3])
            except ValueError:
                continue

            strand = fields[4]

            if strand not in {"+", "-", "."}:
                continue

            isoform_count = ""
            if len(fields) >= 6:
                try:
                    isoform_count = int(fields[5])
                except ValueError:
                    isoform_count = fields[5]

            coordinates[name].append(
                {
                    "lncrna_name": name,
                    "lncrna_chr": chromosome,
                    "lncrna_start": start,
                    "lncrna_end": end,
                    "lncrna_strand": strand,
                    "lncrna_isoform_count": isoform_count,
                }
            )

    return dict(coordinates)


def read_mouse_anchor_records(path: str) -> List[Dict[str, object]]:
    records: List[Dict[str, object]] = []

    current_name = ""
    current_upstream = empty_anchor()

    with open(path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith("lncRNA:"):
                current_name = line.split(":", 1)[1].strip()
                current_upstream = empty_anchor()

            elif line.startswith("upstream PCG:"):
                current_upstream = parse_anchor(line.split(":", 1)[1].strip())

            elif line.startswith("downstream PCG:"):
                current_downstream = parse_anchor(line.split(":", 1)[1].strip())

                if current_name:
                    records.append(
                        {
                            "lncrna_name": current_name,
                            "upstream": current_upstream,
                            "downstream": current_downstream,
                        }
                    )

                current_name = ""
                current_upstream = empty_anchor()

    return records


def choose_coordinate_match(
    anchor_record: Dict[str, object],
    candidates: List[Dict[str, object]],
) -> tuple[List[Dict[str, object]], str]:
    """
    Resolve duplicate lncRNA names using the chromosome of their PCG anchors.

    Most names will have one coordinate record. When a name occurs at multiple
    loci, retain candidates whose chromosome agrees with either anchor.
    """

    if not candidates:
        return [], "NO_COORDINATE_MATCH"

    if len(candidates) == 1:
        return candidates, "UNIQUE_MATCH"

    upstream = anchor_record["upstream"]
    downstream = anchor_record["downstream"]

    anchor_chromosomes = {
        str(upstream.get("chromosome", "")),
        str(downstream.get("chromosome", "")),
    }
    anchor_chromosomes.discard("")

    chromosome_matches = [
        record
        for record in candidates
        if str(record["lncrna_chr"]) in anchor_chromosomes
    ]

    if len(chromosome_matches) == 1:
        return chromosome_matches, "CHROMOSOME_RESOLVED_MATCH"

    if chromosome_matches:
        return chromosome_matches, "MULTIPLE_COORDINATE_MATCHES"

    return candidates, "MULTIPLE_NAME_MATCHES_UNRESOLVED"


def main() -> None:
    missing = [
        path for path in (STEP1_FILE, ANCHOR_FILE)
        if not os.path.isfile(path)
    ]

    if missing:
        print("ERROR: These input files were not found:", file=sys.stderr)
        for path in missing:
            print(f"  - {path}", file=sys.stderr)
        print(
            "\nPut this script in the same folder as both files, or edit the "
            "filenames at the top of the script.",
            file=sys.stderr,
        )
        raise SystemExit(1)

    print("Reading mouse lncRNA coordinates...")
    coordinate_index = read_mouse_lncrna_coordinates(STEP1_FILE)

    print("Reading mouse lncRNA PCG anchors...")
    anchor_records = read_mouse_anchor_records(ANCHOR_FILE)

    output_columns = [
        "mouse_lncrna_name",
        "mouse_lncrna_chr",
        "mouse_lncrna_start",
        "mouse_lncrna_end",
        "mouse_lncrna_strand",
        "mouse_lncrna_isoform_count",
        "mouse_upstream_pcg",
        "mouse_upstream_distance",
        "mouse_upstream_chr",
        "mouse_upstream_start",
        "mouse_upstream_end",
        "mouse_upstream_strand",
        "mouse_downstream_pcg",
        "mouse_downstream_distance",
        "mouse_downstream_chr",
        "mouse_downstream_start",
        "mouse_downstream_end",
        "mouse_downstream_strand",
        "coordinate_match_status",
    ]

    rows_written = 0
    no_coordinate_match = 0
    ambiguous_match = 0

    with open(OUTPUT_FILE, "w", newline="") as outfile:
        writer = csv.DictWriter(
            outfile,
            fieldnames=output_columns,
            delimiter="\t",
        )
        writer.writeheader()

        for anchor_record in anchor_records:
            name = str(anchor_record["lncrna_name"])
            candidates = coordinate_index.get(name, [])

            selected, status = choose_coordinate_match(
                anchor_record,
                candidates,
            )

            if not selected:
                selected = [
                    {
                        "lncrna_name": name,
                        "lncrna_chr": "",
                        "lncrna_start": "",
                        "lncrna_end": "",
                        "lncrna_strand": "",
                        "lncrna_isoform_count": "",
                    }
                ]
                no_coordinate_match += 1

            if status.startswith("MULTIPLE"):
                ambiguous_match += 1

            upstream = anchor_record["upstream"]
            downstream = anchor_record["downstream"]

            for coordinate in selected:
                writer.writerow(
                    {
                        "mouse_lncrna_name": coordinate["lncrna_name"],
                        "mouse_lncrna_chr": coordinate["lncrna_chr"],
                        "mouse_lncrna_start": coordinate["lncrna_start"],
                        "mouse_lncrna_end": coordinate["lncrna_end"],
                        "mouse_lncrna_strand": coordinate["lncrna_strand"],
                        "mouse_lncrna_isoform_count":
                            coordinate["lncrna_isoform_count"],

                        "mouse_upstream_pcg": upstream["pcg_name"],
                        "mouse_upstream_distance": upstream["distance"],
                        "mouse_upstream_chr": upstream["chromosome"],
                        "mouse_upstream_start": upstream["start"],
                        "mouse_upstream_end": upstream["end"],
                        "mouse_upstream_strand": upstream["strand"],

                        "mouse_downstream_pcg": downstream["pcg_name"],
                        "mouse_downstream_distance": downstream["distance"],
                        "mouse_downstream_chr": downstream["chromosome"],
                        "mouse_downstream_start": downstream["start"],
                        "mouse_downstream_end": downstream["end"],
                        "mouse_downstream_strand": downstream["strand"],

                        "coordinate_match_status": status,
                    }
                )
                rows_written += 1

    coordinate_count = sum(
        len(records) for records in coordinate_index.values()
    )

    print("\nFinished.")
    print(f"Mouse coordinate records read: {coordinate_count}")
    print(f"Mouse anchor records read: {len(anchor_records)}")
    print(f"Rows written: {rows_written}")
    print(f"Anchor records without coordinates: {no_coordinate_match}")
    print(f"Ambiguous duplicate-name matches: {ambiguous_match}")
    print(f"Output file: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
