"""Adapter trimming for a successfully paired FASTQ task.

The published ``output/R1.fastq`` / ``output/R2.fastq`` files of a succeeded
task are opened read-only.  Every pair is validated again while streaming
(paired read ids, four-line FASTQ structure, sequence/quality length), adapter
contamination is removed from the 3' end of both mates, and the surviving
pairs are written to a fresh ZIP archive together with a JSON summary.

Design points:

* One request never touches another request's files.  Every call creates an
  isolated ``trims/<uuid>/`` directory below the task directory and, on any
  failure, removes only that directory.
* Validation and output happen in a single streaming pass over the inputs, so
  the total input size may exceed the process memory cap.
* A partial ZIP is never exposed: the archive is built under a temporary name
  inside the request directory and renamed only after it is complete and
  fsynced.  Any error removes the request directory entirely.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import uuid
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence

from . import config
from .fastq import FastqFormatError, _iter_raw_records, extract_read_id

_ADAPTER_RE = re.compile(
    r"^[ACGTacgt]{%d,%d}$"
    % (config.TRIM_ADAPTER_MIN_LENGTH, config.TRIM_ADAPTER_MAX_LENGTH)
)

ZIP_NAME = "trimmed_fastq.zip"
_META_NAME = "trim.json"


class TrimSourceError(Exception):
    """Source files are missing, inconsistent or malformed; no ZIP is made."""


class TrimLimitError(TrimSourceError):
    """The source has more pairs than one trim request may process."""


class TrimReadLengthError(TrimSourceError):
    """A read is outside the supported 40..1000 base range."""


@dataclass(frozen=True)
class TrimStats:
    total_pairs: int
    retained_pairs: int
    dropped_pairs: int


@dataclass(frozen=True)
class TrimResult:
    trim_id: str
    path: Path
    stats: TrimStats
    adapters: tuple[str, ...]


def normalize_adapters(adapters: Sequence[str]) -> list[str]:
    """Validate 1..4 adapters of 16..32 A/C/G/T bases; return upper case."""
    if not 1 <= len(adapters) <= config.TRIM_MAX_ADAPTERS:
        raise ValueError(
            f"expected 1..{config.TRIM_MAX_ADAPTERS} adapters, got {len(adapters)}"
        )
    normalized: list[str] = []
    for index, adapter in enumerate(adapters, 1):
        if not isinstance(adapter, str) or _ADAPTER_RE.fullmatch(adapter) is None:
            raise ValueError(
                f"adapter {index} must be "
                f"{config.TRIM_ADAPTER_MIN_LENGTH}.."
                f"{config.TRIM_ADAPTER_MAX_LENGTH} characters of A/C/G/T"
            )
        normalized.append(adapter.upper())
    return normalized


class AdapterMatcher:
    """Longest 3'-overlap matcher with at most one mismatch.

    Adapters are 16..32 bases long, so every hit involves the read's last
    32 bases.  A prefilter indexes adapter 8-mers at each overlap length; a
    read's final 8 bases are compared against
    all 24 single-substitution variants plus the exact 8-mer.  Any alignment
    with at most one mismatch must have its last 8 aligned bases in that
    index, so the prefilter cannot hide a hit.  Candidates are verified on the
    full overlap, longest first.
    """

    def __init__(self, adapters: Sequence[str]) -> None:
        self.adapters = [adapter.encode("ascii") for adapter in adapters]
        minimum = config.TRIM_MIN_OVERLAP
        self._index: dict[bytes, list[tuple[int, int]]] = {}
        for adapter_index, adapter in enumerate(self.adapters):
            for overlap in range(minimum, len(adapter) + 1):
                key = adapter[overlap - minimum:overlap]
                self._index.setdefault(key, []).append((adapter_index, overlap))
        for hits in self._index.values():
            hits.sort(key=lambda item: item[1], reverse=True)

    @staticmethod
    def _variants(mer: bytes) -> Iterator[bytes]:
        yield mer
        for position in range(len(mer)):
            for base in b"ACGT":
                if mer[position] != base:
                    yield mer[:position] + bytes((base,)) + mer[position + 1:]

    def cut_position(self, seq: bytes) -> int:
        """Index where adapter contamination starts (``len(seq)`` = no hit)."""
        length = len(seq)
        minimum = config.TRIM_MIN_OVERLAP
        if length < minimum:
            return length

        last = seq[-minimum:].upper()
        index = self._index
        candidates: set[tuple[int, int]] = set()
        for variant in self._variants(last):
            for hit in index.get(variant, ()):
                candidates.add(hit)

        for adapter_index, overlap in sorted(
            candidates, key=lambda item: item[1], reverse=True
        ):
            if overlap > length:
                continue
            left = seq[length - overlap:].upper()
            right = self.adapters[adapter_index][:overlap]
            if left == right:
                return length - overlap
            mismatches = 0
            for a, b in zip(left, right):
                if a != b:
                    mismatches += 1
                    if mismatches > config.TRIM_MAX_MISMATCHES:
                        break
            else:
                return length - overlap
        return length


def _split_eol(line: bytes) -> tuple[bytes, bytes]:
    """Split a raw FASTQ line into content and its original line terminator."""
    if line.endswith(b"\r\n"):
        return line[:-2], b"\r\n"
    if line.endswith(b"\n"):
        return line[:-1], b"\n"
    return line, b""


def _paired_records(
    r1_path: Path, r2_path: Path
) -> Iterator[
    tuple[
        tuple[bytes, bytes, bytes, bytes],
        tuple[bytes, bytes, bytes, bytes],
    ]
]:
    """Yield (record_r1, record_r2) in lockstep, validating both files.

    ``_iter_raw_records`` enforces the four-line FASTQ structure, non-empty
    sequence/quality, equal sequence/quality lengths and detects a truncated
    tail.  Here we additionally require both files to end at the same pair.
    """
    records_r1 = _iter_raw_records(r1_path, "R1")
    records_r2 = _iter_raw_records(r2_path, "R2")
    while True:
        try:
            record_r1 = next(records_r1)
        except StopIteration:
            record_r1 = None
        try:
            record_r2 = next(records_r2)
        except StopIteration:
            record_r2 = None
        if record_r1 is None and record_r2 is None:
            return
        if record_r1 is None:
            raise TrimSourceError("R2 contains more records than R1")
        if record_r2 is None:
            raise TrimSourceError("R1 contains more records than R2")
        yield record_r1, record_r2


def _trim_to_temp(
    r1_path: Path,
    r2_path: Path,
    matcher: AdapterMatcher,
    out_r1: Path,
    out_r2: Path,
    expected_pairs: int | None,
) -> TrimStats:
    """Stream both source files, trimming contamination into temp files."""
    total = 0
    dropped = 0
    pending_r1: list[bytes] = []
    pending_r2: list[bytes] = []
    pending_bytes = 0

    try:
        with out_r1.open("wb") as f1, out_r2.open("wb") as f2:
            for record_r1, record_r2 in _paired_records(r1_path, r2_path):
                (header1, seq1_line, plus1, qual1_line) = record_r1
                (header2, seq2_line, plus2, qual2_line) = record_r2
                total += 1
                if total > config.TRIM_MAX_PAIRS:
                    raise TrimLimitError(
                        f"source has more than {config.TRIM_MAX_PAIRS} pairs; "
                        "split the task before trimming"
                    )
                if expected_pairs is not None and total > expected_pairs:
                    raise TrimSourceError(
                        "published result files contain more pairs than the "
                        f"{expected_pairs} recorded during pairing"
                    )

                seq1, seq1_eol = _split_eol(seq1_line)
                qual1, qual1_eol = _split_eol(qual1_line)
                seq2, seq2_eol = _split_eol(seq2_line)
                qual2, qual2_eol = _split_eol(qual2_line)

                read_id1 = extract_read_id(header1.rstrip(b"\r\n"))
                read_id2 = extract_read_id(header2.rstrip(b"\r\n"))
                if read_id1 != read_id2:
                    raise TrimSourceError(
                        f"pair {total} has mismatched read ids: "
                        f"R1 '{read_id1}' vs R2 '{read_id2}'"
                    )
                for side, seq, qual in (("R1", seq1, qual1), ("R2", seq2, qual2)):
                    if len(seq) != len(qual):
                        raise TrimSourceError(
                            f"pair {total} ({side}) sequence length {len(seq)} "
                            f"!= quality length {len(qual)}"
                        )
                    if not (
                        config.TRIM_MIN_READ_LENGTH
                        <= len(seq)
                        <= config.TRIM_MAX_READ_LENGTH
                    ):
                        raise TrimReadLengthError(
                            f"pair {total} ({side}) read '{read_id1}' has "
                            f"{len(seq)} bases; trim supports "
                            f"{config.TRIM_MIN_READ_LENGTH}.."
                            f"{config.TRIM_MAX_READ_LENGTH} bases per end"
                        )

                cut1 = matcher.cut_position(seq1)
                cut2 = matcher.cut_position(seq2)
                if (
                    cut1 < config.TRIM_MIN_RETAINED_LENGTH
                    or cut2 < config.TRIM_MIN_RETAINED_LENGTH
                ):
                    dropped += 1
                    continue

                blob1 = b"".join((
                    header1, seq1[:cut1], seq1_eol,
                    plus1, qual1[:cut1], qual1_eol,
                ))
                blob2 = b"".join((
                    header2, seq2[:cut2], seq2_eol,
                    plus2, qual2[:cut2], qual2_eol,
                ))
                pending_r1.append(blob1)
                pending_r2.append(blob2)
                pending_bytes += len(blob1) + len(blob2)
                if pending_bytes >= config.TRIM_WRITE_BATCH_BYTES:
                    f1.write(b"".join(pending_r1))
                    f2.write(b"".join(pending_r2))
                    pending_r1.clear()
                    pending_r2.clear()
                    pending_bytes = 0
            if pending_r1:
                f1.write(b"".join(pending_r1))
                f2.write(b"".join(pending_r2))
            f1.flush()
            f2.flush()
            os.fsync(f1.fileno())
            os.fsync(f2.fileno())
    except FastqFormatError as exc:
        raise TrimSourceError(f"{exc}") from exc
    except OSError as exc:
        raise TrimSourceError(f"cannot read source files: {exc}") from exc

    if expected_pairs is not None and total != expected_pairs:
        raise TrimSourceError(
            f"published result files contain {total} pairs but pairing "
            f"recorded {expected_pairs}"
        )
    return TrimStats(
        total_pairs=total,
        retained_pairs=total - dropped,
        dropped_pairs=dropped,
    )


def _build_zip(
    zip_path: Path,
    trimmed_r1: Path,
    trimmed_r2: Path,
    summary: dict,
) -> None:
    """Write the ZIP by streaming the two trimmed files in fixed chunks."""
    with zipfile.ZipFile(
        zip_path,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=1,
        allowZip64=True,
    ) as archive:
        for arcname, source in (
            ("R1.fastq", trimmed_r1),
            ("R2.fastq", trimmed_r2),
        ):
            with source.open("rb") as src, archive.open(arcname, "w") as dest:
                while True:
                    chunk = src.read(config.TRIM_WRITE_BATCH_BYTES)
                    if not chunk:
                        break
                    dest.write(chunk)
        archive.writestr(
            "summary.json",
            (json.dumps(summary, ensure_ascii=False, indent=2) + "\n").encode(),
        )


def discard_trim_dir(task_id: str, trim_id: str,
                     tasks_dir: Path | None = None) -> None:
    """Remove one trim request directory and nothing else."""
    root = tasks_dir if tasks_dir is not None else config.TASKS_DIR
    shutil.rmtree(root / task_id / "trims" / trim_id, ignore_errors=True)


def create_trim_zip(
    task_id: str,
    r1_path: Path,
    r2_path: Path,
    adapters: Sequence[str],
    *,
    expected_pairs: int | None = None,
    tasks_dir: Path | None = None,
) -> TrimResult:
    """Trim one succeeded task and return the new ZIP plus its statistics.

    The caller is responsible for checking that the task exists and succeeded.
    Every request is isolated in its own UUID directory; only that directory
    is removed when the request fails.
    """
    normalized = normalize_adapters(adapters)
    matcher = AdapterMatcher(normalized)

    root = tasks_dir if tasks_dir is not None else config.TASKS_DIR
    trim_id = uuid.uuid4().hex
    work_dir = root / task_id / "trims" / trim_id
    work_dir.mkdir(parents=True, exist_ok=False)

    trimmed_r1 = work_dir / ".R1.fastq.part"
    trimmed_r2 = work_dir / ".R2.fastq.part"
    zip_tmp = work_dir / (".%s.part" % ZIP_NAME)
    final_zip = work_dir / ZIP_NAME

    try:
        stats = _trim_to_temp(
            r1_path,
            r2_path,
            matcher,
            trimmed_r1,
            trimmed_r2,
            expected_pairs,
        )
        summary = {
            "source_task_id": task_id,
            "trim_id": trim_id,
            "adapters": normalized,
            "min_overlap": config.TRIM_MIN_OVERLAP,
            "max_mismatches": config.TRIM_MAX_MISMATCHES,
            "min_retained_length": config.TRIM_MIN_RETAINED_LENGTH,
            "total_pairs": stats.total_pairs,
            "retained_pairs": stats.retained_pairs,
            "dropped_pairs": stats.dropped_pairs,
        }
        _build_zip(zip_tmp, trimmed_r1, trimmed_r2, summary)
        with zip_tmp.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(zip_tmp, final_zip)
        (work_dir / _META_NAME).write_text(
            json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        for part in (trimmed_r1, trimmed_r2):
            try:
                part.unlink(missing_ok=True)
            except OSError:
                pass
        return TrimResult(
            trim_id=trim_id,
            path=final_zip,
            stats=stats,
            adapters=tuple(normalized),
        )
    except BaseException:
        # Remove only this request's directory.  Concurrent requests for the
        # same task have different UUID directories and are never touched.
        discard_trim_dir(task_id, trim_id, root)
        raise
