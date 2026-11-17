"""Memory test: trimming a source larger than the 256 MiB process cap.

Builds a synthetic *succeeded* task with 200,000 pairs of 400-base reads
(~330 MiB across both published files), runs ``trim.create_trim_zip`` in a
child process and measures the child's peak RSS (VmHWM).  The run must:

  * produce a complete ZIP (R1/R2 plus summary.json),
  * drop pairs whose either end falls below 20 retained bases and keep the
    rest in R1 order,
  * stay below the memory limit.

Usage:
    .venv/bin/python -I tests/test_trim_memory.py [pairs]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import textwrap
import time
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAIRS = int(sys.argv[1]) if len(sys.argv) > 1 else 200_000
RSS_LIMIT_KIB = 230 * 1024
SEQ_LEN = 400
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTC"
DROPPED_EVERY = 1000  # every 1000th pair is trimmed below 20 bases

CHILD = textwrap.dedent("""
    import sys
    sys.path.insert(0, %r)
    from pathlib import Path
    from app import config, trim

    data_dir = Path(sys.argv[1])
    config.DATA_DIR = data_dir
    config.DB_PATH = data_dir / "tasks.db"
    config.TASKS_DIR = data_dir / "tasks"
    task_id = "memorystress"
    base = data_dir / "tasks" / task_id
    result = trim.create_trim_zip(
        task_id,
        base / "output" / "R1.fastq",
        base / "output" / "R2.fastq",
        [%r],
    )
    print("RESULT=" + result.path.as_posix())
    print("TOTAL=%%d" %% result.stats.total_pairs)
    print("RETAINED=%%d" %% result.stats.retained_pairs)
    print("DROPPED=%%d" %% result.stats.dropped_pairs)
""")


def write_inputs(base: Path, n: int) -> tuple[Path, Path]:
    out = base / "output"
    out.mkdir(parents=True)
    r1 = out / "R1.fastq"
    r2 = out / "R2.fastq"
    # Every DROPPED_EVERY-th pair carries the adapter right after 10 bases, so
    # trimming leaves fewer than 20 bases and the pair must be dropped whole.
    truncated = b"A" * 10 + ADAPTER.encode()
    clean = b"A" * SEQ_LEN
    qual_truncated = b"~" * len(truncated)
    qual = b"~" * SEQ_LEN
    with r1.open("wb") as f1, r2.open("wb") as f2:
        buf1 = bytearray()
        buf2 = bytearray()
        for i in range(n):
            dropped = i % DROPPED_EVERY == 0
            record = (b"@read-%09d/1\n%s\n+\n%s\n"
                      % (i, truncated if dropped else clean,
                         qual_truncated if dropped else qual))
            buf1 += record
            buf2 += (b"@read-%09d/2\n%s\n+\n%s\n" % (i, clean, qual))
            if len(buf1) >= 4 * 1024 * 1024:
                f1.write(buf1)
                f2.write(buf2)
                buf1.clear()
                buf2.clear()
        f1.write(buf1)
        f2.write(buf2)
    return r1, r2


def read_vmhwm(pid: int) -> int:
    for line in Path(f"/proc/{pid}/status").read_text().splitlines():
        if line.startswith("VmHWM:"):
            return int(line.split()[1])
    raise RuntimeError("VmHWM not available")


def verify(zip_path: Path, n: int) -> tuple[bool, str]:
    """Stream both ZIP members and verify order, trimming and counts."""
    clean_line = b"A" * SEQ_LEN + b"\n"
    dropped = sum(1 for i in range(n) if i % DROPPED_EVERY == 0)
    expected_retained = n - dropped
    with zipfile.ZipFile(zip_path) as archive:
        names = set(archive.namelist())
        if names != {"R1.fastq", "R2.fastq", "summary.json"}:
            return False, f"unexpected ZIP members {sorted(names)}"
        summary = json.loads(archive.read("summary.json"))
        if (
            summary["total_pairs"] != n
            or summary["retained_pairs"] != expected_retained
            or summary["dropped_pairs"] != dropped
        ):
            return False, f"summary mismatch: {summary}"
        with archive.open("R1.fastq") as f1, archive.open("R2.fastq") as f2:
            kept = 0
            expected_id = 0
            while True:
                h1 = f1.readline()
                h2 = f2.readline()
                if not h1 and not h2:
                    break
                s1 = f1.readline()
                f1.readline()  # plus line
                q1 = f1.readline()
                s2 = f2.readline()
                f2.readline()  # plus line
                q2 = f2.readline()
                # Skip source ids that were expected to be dropped.
                while expected_id % DROPPED_EVERY == 0:
                    expected_id += 1
                if expected_id >= n:
                    return False, f"unexpected extra pair after source id {n - 1}"
                want1 = b"@read-%09d/1\n" % expected_id
                want2 = b"@read-%09d/2\n" % expected_id
                if h1 != want1 or h2 != want2:
                    return False, (
                        f"order mismatch at source id {expected_id}: "
                        f"{h1!r} {h2!r}"
                    )
                if s1 != clean_line or s2 != clean_line:
                    return False, f"sequence mismatch at pair {expected_id}"
                if len(q1) != len(s1) or len(q2) != len(s2):
                    return False, f"quality length mismatch at pair {expected_id}"
                kept += 1
                expected_id += 1
            if kept != expected_retained:
                return False, f"ZIP has {kept} pairs, expected {expected_retained}"
            if f1.readline() or f2.readline():
                return False, "trailing data after expected pairs"
    return True, f"verified {kept} retained pairs, {dropped} dropped"


def main() -> int:
    scratch = Path(tempfile.mkdtemp(prefix="fastq-trim-mem-",
                                    dir=str(ROOT / "tests" / "output")))
    data_dir = scratch / "data"
    base = data_dir / "tasks" / "memorystress"
    print(f"Generating {PAIRS:,} pairs (~"
          f"{2 * PAIRS * (SEQ_LEN * 2 + 30) / (1024**2):.0f} MiB total)...")
    t0 = time.time()
    r1, r2 = write_inputs(base, PAIRS)
    print(f"  inputs written in {time.time() - t0:.1f}s: "
          f"R1={r1.stat().st_size / (1024**2):.1f} MiB, "
          f"R2={r2.stat().st_size / (1024**2):.1f} MiB")

    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    proc = subprocess.Popen(
        [sys.executable, "-I", "-c", CHILD % (str(ROOT), ADAPTER), str(data_dir)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
        text=True,
    )
    t0 = time.time()
    peak = 0
    while proc.poll() is None:
        try:
            peak = max(peak, read_vmhwm(proc.pid))
        except (FileNotFoundError, ProcessLookupError):
            break
        time.sleep(0.1)
    out, err = proc.communicate()
    runtime = time.time() - t0
    print(out.strip())
    if proc.returncode != 0:
        print(err[-4000:])
        print("RESULT: FAIL (child exited %d)" % proc.returncode)
        return 1

    zip_path = next(
        (Path(line.split("=", 1)[1]) for line in out.splitlines()
         if line.startswith("RESULT=")), None,
    )
    if zip_path is None or not zip_path.is_file():
        print("RESULT: FAIL (no ZIP produced)")
        return 1

    ok, detail = verify(zip_path, PAIRS)
    print(detail)
    print(f"peak RSS: {peak / 1024:.1f} MiB "
          f"(limit {RSS_LIMIT_KIB / 1024:.0f} MiB), runtime {runtime:.1f}s, "
          f"ZIP {zip_path.stat().st_size / (1024**2):.1f} MiB")
    if peak > RSS_LIMIT_KIB:
        print("RESULT: FAIL (memory cap exceeded)")
        return 1
    if not ok:
        print("RESULT: FAIL (ZIP verification)")
        return 1
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
