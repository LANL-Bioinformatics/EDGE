#!/usr/bin/env python3

import argparse
import gzip
import re
from collections import Counter


ILLUMINA_PATTERNS = [
    # Modern Illumina CASAVA 1.8+
    # @MACHINE:RUN:FLOWCELL:LANE:TILE:X:Y READ:FILTER:CONTROL:INDEX
    re.compile(
        r"^@[^:\s]+:\d+:[^:\s]+:\d+:\d+:\d+:\d+\s+[12]:[YN]:\d+:\S+",
        re.IGNORECASE,
    ),

    # Older Illumina format
    # @MACHINE:LANE:TILE:X:Y#INDEX/READ
    re.compile(
        r"^@[^:\s]+:\d+:\d+:\d+:\d+#?[A-Z0-9]*\/[12]",
        re.IGNORECASE,
    ),

    # Illumina-like 7 colon-separated fields
    re.compile(
        r"^@[^:\s]+:\d+:[^:\s]+:\d+:\d+:\d+:\d+",
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
    # PacBio subread / CCS style
    # @movie/zmw/start_end
    # @movie/zmw/ccs
    re.compile(r"^@[^/\s]+/\d+/\d+_\d+", re.IGNORECASE),
    re.compile(r"^@[^/\s]+/\d+/ccs", re.IGNORECASE),

    # PacBio HiFi often has /ccs
    re.compile(r"/ccs\b", re.IGNORECASE),

    # PacBio ZMW names
    re.compile(r"^@[^/\s]+/\d+/", re.IGNORECASE),
]


def open_fastq(path):
    """
    Open FASTQ file, supporting both plain text and gzip-compressed files.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def score_header(header):
    """
    Return platform scores for a single FASTQ header.
    """
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


def detect_platform(path, max_reads=1000):
    """
    Inspect up to max_reads FASTQ records and infer sequencing platform.
    """
    total_scores = Counter()
    headers_checked = 0
    example_headers = {}

    with open_fastq(path) as handle:
        while headers_checked < max_reads:
            header = handle.readline()
            if not header:
                break

            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()

            if not quality:
                break

            header = header.strip()

            if not header.startswith("@"):
                continue

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
            "headers_checked": 0,
            "scores": dict(total_scores),
            "examples": example_headers,
        }

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
        "headers_checked": headers_checked,
        "scores": dict(total_scores),
        "examples": example_headers,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Detect sequencing platform from FASTQ headers."
    )
    parser.add_argument(
        "fastq",
        help="Input FASTQ file, optionally gzipped"
    )
    parser.add_argument(
        "-n",
        "--max-reads",
        type=int,
        default=1000,
        help="Maximum number of FASTQ records to inspect. Default: 1000"
    )
    parser.add_argument(
        "--show-example",
        action="store_true",
        help="Show example matching headers"
    )

    args = parser.parse_args()

    result = detect_platform(args.fastq, max_reads=args.max_reads)

    print(f"File: {args.fastq}")
    print(f"Detected platform: {result['platform']}")
    print(f"Confidence: {result['confidence']}")
    print(f"Headers checked: {result['headers_checked']}")
    print("Scores:")

    for platform, score in result["scores"].items():
        print(f"  {platform}: {score}")

    if args.show_example and result["examples"]:
        print("\nExample matching headers:")
        for platform, header in result["examples"].items():
            print(f"  {platform}: {header}")


if __name__ == "__main__":
    main()
