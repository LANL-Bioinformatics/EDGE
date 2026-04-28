#!/usr/bin/env python3
"""
## FASTQ Header Examples Detected by Platform

The FASTQ inspector checks read header lines to infer the sequencing platform. The examples below show header formats that are recognized.

#### Illumina

Modern Illumina / CASAVA-style headers:

```text
@M00176:17:000000000-A3JHG:1:1101:16713:1403 1:N:0:1
@A016U:114:H7NYMDSX7:1:1101:16853:1322 2:N:0:ATCACG
@:114:A016U:1:1:16853:1322 1:N:0:1
```

Older Illumina-style headers with /1 or /2 read markers:
```text
@HWUSI-EAS100R:6:73:941:1973#0/1
@HWI-ST1234:88:1101:1234:5678/2
```

#### Oxford Nanopore
Nanopore headers commonly include UUID-style read IDs and metadata such as `runid`, `read`, `ch`, or `start_time`:
```text
@b3a6f4f2-89aa-4a4f-8b84-f6fd2c1f74bd runid=1a2b3c4d5e6f read=42 ch=128 start_time=2024-03-21T12:34:56Z
@e0c10aa1-9b7a-4d2a-9a6f-4e6b3f86d58f runid=abcdef1234567890
@9f2a1c3b-1111-4222-8333-abcdefabcdef
```

#### PacBio
PacBio headers commonly use movie/ZMW/read-coordinate formats:
```text
@m64011_190830_220126/14/0_12345
@m64011_190830_220126/14/ccs
@m54238_180901_011437/4194312/100_2500
```

These examples are not exhaustive, but they represent the header shapes currently matched by the FASTQ inspector.
"""

import argparse
import gzip
import re
from collections import Counter


ILLUMINA_PATTERNS = [
    # Modern Illumina CASAVA 1.8+
    # @MACHINE:RUN:FLOWCELL:LANE:TILE:X:Y READ:FILTER:CONTROL:INDEX
    re.compile(
        r"^@[^:\s]*:\d+:[^:\s]+:\d+:\d+:\d+:\d+\s+[12]:[YN]:\d+:\S+",
        re.IGNORECASE,
    ),

    # Older Illumina format
    # @MACHINE:LANE:TILE:X:Y#INDEX/READ
    re.compile(
        r"^@[^:\s]*:\d+:\d+:\d+:\d+#?[A-Z0-9]*\/[12]",
        re.IGNORECASE,
    ),

    # Illumina-like 7 colon-separated fields
    re.compile(
        r"^@[^:\s]*:\d+:[^:\s]+:\d+:\d+:\d+:\d+",
        re.IGNORECASE,
    ),
]


NANOPORE_PATTERNS = [
    # Common Guppy / Dorado / MinKNOW style
    # @read_id runid=... read=... ch=... start_time=...
    re.compile(r"^@[a-f0-9-]{20,}\s+.*\brunid=", re.IGNORECASE),
    re.compile(r"\brunid=[a-f0-9]+", re.IGNORECASE),
    re.compile(r"\bread=\d+", re.IGNORECASE),
    re.compile(r"\bch=\d+", re.IGNORECASE),
    re.compile(r"\bstart_time=", re.IGNORECASE),

    # ONT UUID-like read IDs
    re.compile(
        r"^@[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}",
        re.IGNORECASE,
    ),
]


PACBIO_PATTERNS = [
    # PacBio subread / CCS / HiFi style
    # @movie/zmw/start_end
    # @movie/zmw/ccs
    re.compile(r"^@[^/\s]+/\d+/\d+_\d+", re.IGNORECASE),
    re.compile(r"^@[^/\s]+/\d+/ccs", re.IGNORECASE),
    re.compile(r"/ccs\b", re.IGNORECASE),
    re.compile(r"^@[^/\s]+/\d+/", re.IGNORECASE),
]


def open_sequence_file(path):
    """
    Open sequence file, supporting plain text and gzip-compressed files.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def normalize_header(header):
    """
    Convert FASTA headers beginning with '>' to FASTQ-like headers beginning
    with '@' so the same regex patterns can be reused.
    """
    header = header.strip()

    if header.startswith(">"):
        return "@" + header[1:]

    return header


def score_header(header):
    """
    Return platform scores for a single FASTQ or FASTA header.
    """
    header = normalize_header(header)

    scores = {
        "Illumina": 0,
        "Nanopore": 0,
        "PacBio": 0,
    }

    for pattern in ILLUMINA_PATTERNS:
        if pattern.search(header):
            scores["Illumina"] += 1

    for pattern in NANOPORE_PATTERNS:
        if pattern.search(header):
            scores["Nanopore"] += 1

    for pattern in PACBIO_PATTERNS:
        if pattern.search(header):
            scores["PacBio"] += 1

    return scores


def detect_file_format(path):
    """
    Guess whether the input is FASTQ or FASTA based on the first non-empty line.
    """
    with open_sequence_file(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith("@"):
                return "fastq"
            if line.startswith(">"):
                return "fasta"

            return "unknown"

    return "unknown"


def iter_fastq_headers(handle, max_reads):
    """
    Yield headers from a FASTQ file.
    """
    reads_seen = 0

    while reads_seen < max_reads:
        header = handle.readline()
        if not header:
            break

        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()

        if not quality:
            break

        header = header.strip()

        if header.startswith("@"):
            yield header
            reads_seen += 1


def iter_fasta_headers(handle, max_reads):
    """
    Yield headers from a FASTA file.
    """
    reads_seen = 0

    for line in handle:
        if reads_seen >= max_reads:
            break

        line = line.strip()

        if line.startswith(">"):
            yield line
            reads_seen += 1


def detect_platform(path, max_reads=1000, file_format="auto"):
    """
    Inspect up to max_reads FASTQ or FASTA records and infer sequencing platform.
    """
    if file_format == "auto":
        file_format = detect_file_format(path)

    total_scores = Counter()
    headers_checked = 0
    example_headers = {}

    with open_sequence_file(path) as handle:
        if file_format == "fastq":
            header_iter = iter_fastq_headers(handle, max_reads)
        elif file_format == "fasta":
            header_iter = iter_fasta_headers(handle, max_reads)
        else:
            return {
                "platform": "Unknown",
                "confidence": "none",
                "file_format": file_format,
                "headers_checked": 0,
                "scores": dict(total_scores),
                "examples": example_headers,
            }

        for header in header_iter:
            scores = score_header(header)

            for platform, score in scores.items():
                total_scores[platform] += score
                if score > 0 and platform not in example_headers:
                    example_headers[platform] = header

            headers_checked += 1

    if headers_checked == 0:
        return {
            "platform": "Unknown",
            "confidence": "none",
            "file_format": file_format,
            "headers_checked": 0,
            "scores": dict(total_scores),
            "examples": example_headers,
        }

    # Ensure all platforms are represented, even if their score is zero.
    for platform in ["Illumina", "Nanopore", "PacBio"]:
        total_scores.setdefault(platform, 0)

    best_platform, best_score = total_scores.most_common(1)[0]

    if best_score == 0:
        platform = "Unknown"
        confidence = "none"
    else:
        sorted_scores = total_scores.most_common()
        second_score = sorted_scores[1][1] if len(sorted_scores) > 1 else 0

        if best_score >= 5 and best_score >= 2 * max(second_score, 1):
            confidence = "high"
        elif best_score > second_score:
            confidence = "medium"
        else:
            confidence = "low"

        platform = best_platform

    return {
        "platform": platform,
        "confidence": confidence,
        "file_format": file_format,
        "headers_checked": headers_checked,
        "scores": dict(total_scores),
        "examples": example_headers,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Detect sequencing platform from FASTQ or FASTA headers."
    )
    parser.add_argument(
        "input",
        help="Input FASTQ or FASTA file, optionally gzipped"
    )
    parser.add_argument(
        "-n",
        "--max-reads",
        type=int,
        default=1000,
        help="Maximum number of records to inspect. Default: 1000"
    )
    parser.add_argument(
        "--format",
        choices=["auto", "fastq", "fasta"],
        default="auto",
        help="Input format. Default: auto"
    )
    parser.add_argument(
        "--show-example",
        action="store_true",
        help="Show example matching headers"
    )

    args = parser.parse_args()

    result = detect_platform(
        args.input,
        max_reads=args.max_reads,
        file_format=args.format,
    )

    print(f"File: {args.input}")
    print(f"Detected input format: {result['file_format']}")
    print(f"Detected platform: {result['platform']}")
    print(f"Confidence: {result['confidence']}")
    print(f"Headers checked: {result['headers_checked']}")
    print("Scores:")

    for platform in ["Illumina", "Nanopore", "PacBio"]:
        print(f"  {platform}: {result['scores'].get(platform, 0)}")

    if args.show_example and result["examples"]:
        print("\nExample matching headers:")
        for platform, header in result["examples"].items():
            print(f"  {platform}: {header}")


if __name__ == "__main__":
    main()
