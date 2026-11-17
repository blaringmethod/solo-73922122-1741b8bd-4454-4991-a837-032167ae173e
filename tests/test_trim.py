"""Self-tests for adapter trimming.

Run with the project virtualenv:
    .venv/bin/python -I tests/test_trim.py

The tests cover the trimming algorithm, the HTTP API, validation of the
source FASTQ stream, isolation between concurrent requests for one task, and
the documented limits (1..4 adapters of 16..32 A/C/G/T, 200k pairs, 40..1000
bases per end).
"""

from __future__ import annotations

import asyncio
import hashlib
import io
import json
import os
import shutil
import sys
import tempfile
import time
import uuid
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import httpx  # noqa: E402

from app import config, db, trim, worker  # noqa: E402
from app.main import app  # noqa: E402

RESULTS: list[tuple[str, bool, str]] = []

ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTC"       # 32 bases
ADAPTER_B = "CTGTCTCTTATACACATCT"                  # 19 bases


def record(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))


def make_record(rid: str, seq: bytes, quality: bytes | None = None,
                suffix: str = "", desc: bytes = b"") -> bytes:
    if quality is None:
        quality = b"~" * len(seq)
    header = b"@" + rid.encode() + suffix.encode()
    if desc:
        header += b" " + desc
    return b"\n".join((header, seq, b"+", quality)) + b"\n"


def fastq(records: list[bytes]) -> io.BytesIO:
    return io.BytesIO(b"".join(records))


async def upload(client: httpx.AsyncClient, r1: io.BytesIO, r2: io.BytesIO):
    r1.seek(0)
    r2.seek(0)
    return await client.post(
        "/tasks",
        files={"r1": ("r1.fastq", r1, "application/octet-stream"),
               "r2": ("r2.fastq", r2, "application/octet-stream")},
    )


async def wait_terminal(client: httpx.AsyncClient, task_id: str,
                        timeout: float = 60.0) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        body = (await client.get(f"/tasks/{task_id}")).json()
        if body["status"] in ("succeeded", "failed", "cancelled"):
            return body
        await asyncio.sleep(0.05)
    raise AssertionError(f"task {task_id} did not finish in {timeout}s")


def parse_fastq_bytes(data: bytes) -> list[tuple[str, bytes, bytes, bytes, bytes]]:
    lines = data.splitlines(keepends=True)
    assert len(lines) % 4 == 0, f"not four-line FASTQ: {len(lines)} lines"
    out = []
    for i in range(0, len(lines), 4):
        header = lines[i].rstrip(b"\r\n")
        rid = header[1:].split(None, 1)[0]
        if rid.endswith(b"/1") or rid.endswith(b"/2"):
            rid = rid[:-2]
        seq = lines[i + 1].rstrip(b"\r\n")
        plus = lines[i + 2].rstrip(b"\r\n")
        qual = lines[i + 3].rstrip(b"\r\n")
        out.append((rid.decode("latin-1"), lines[i], seq, plus, qual))
    return out


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_zip(zip_path: Path) -> tuple[dict, bytes, bytes]:
    with zipfile.ZipFile(zip_path) as archive:
        names = set(archive.namelist())
        assert names == {"R1.fastq", "R2.fastq", "summary.json"}, names
        summary = json.loads(archive.read("summary.json"))
        return summary, archive.read("R1.fastq"), archive.read("R2.fastq")


def dirs_under(path: Path) -> list[str]:
    if not path.is_dir():
        return []
    return sorted(p.name for p in path.iterdir())


def matcher_unit_tests() -> None:
    matcher = trim.AdapterMatcher([ADAPTER])
    record("cut: full adapter at 3' end",
           matcher.cut_position(b"A" * 18 + ADAPTER.encode()) == 18)
    record("cut: 8-base overlap",
           matcher.cut_position(b"C" * 42 + b"AGATCGGA") == 42)
    record("cut: 8-base overlap with one mismatch",
           matcher.cut_position(b"C" * 42 + b"AGATCGGT") == 42)
    record("no cut: 8-base overlap with two mismatches",
           matcher.cut_position(b"C" * 42 + b"AGATCGGT".replace(b"T", b"A")
                                .replace(b"G", b"C", 1)) == 50)
    record("no cut: only 7 matching bases",
           matcher.cut_position(b"C" * 43 + b"AGATCGG") == 50)
    record("no cut: adapter only internally",
           matcher.cut_position(b"AGATCGGA" + b"C" * 42) == 50)
    record("no cut: clean read ends with N",
           matcher.cut_position(b"C" * 49 + b"N") == 50)

    multi = trim.AdapterMatcher([ADAPTER, ADAPTER_B])
    overlap_b = ADAPTER_B[:11].encode()
    record("multi adapter: second adapter found",
           multi.cut_position(b"T" * 39 + overlap_b) == 39)
    record("multi adapter: longest hit wins",
           multi.cut_position(b"G" * 40 + ADAPTER.encode()[:20]) == 40)


def direct_trim_tests(tmp: Path) -> None:
    """Validation tests driven directly against crafted source files."""
    task_id = "direct"
    out = tmp / task_id / "output"
    out.mkdir(parents=True)
    r1 = out / "R1.fastq"
    r2 = out / "R2.fastq"

    def put(r1_bytes: bytes, r2_bytes: bytes) -> None:
        r1.write_bytes(r1_bytes)
        r2.write_bytes(r2_bytes)
        trims = tmp / task_id / "trims"
        shutil.rmtree(trims, ignore_errors=True)

    clean1 = make_record("a", b"C" * 60)
    clean2 = make_record("a", b"G" * 60)

    put(make_record("a", b"C" * 60), make_record("b", b"G" * 60))
    try:
        trim.create_trim_zip(task_id, r1, r2, [ADAPTER], tasks_dir=tmp)
        record("read check: mismatched ids rejected", False, "no error raised")
    except trim.TrimSourceError as exc:
        record("read check: mismatched ids rejected", "mismatched read ids" in str(exc),
               str(exc))
    record("read check: failed request leaves no trim dir",
           dirs_under(tmp / task_id / "trims") == [])

    put(clean1, make_record("a", b"G" * 59, b"~" * 58))
    try:
        trim.create_trim_zip(task_id, r1, r2, [ADAPTER], tasks_dir=tmp)
        record("read check: seq/quality length mismatch rejected", False)
    except trim.TrimSourceError as exc:
        record("read check: seq/quality length mismatch rejected",
               "length" in str(exc), str(exc))

    put(clean1, b"@a\n" + b"G" * 60 + b"\n+\n")
    try:
        trim.create_trim_zip(task_id, r1, r2, [ADAPTER], tasks_dir=tmp)
        record("read check: truncated tail rejected", False)
    except trim.TrimSourceError as exc:
        record("read check: truncated tail rejected",
               "four lines" in str(exc) or "missing lines" in str(exc), str(exc))
    record("read check: truncated request leaves no partial ZIP",
           dirs_under(tmp / task_id / "trims") == [])

    put(clean1, make_record("a", b"G" * 60))
    missing = tmp / task_id / "output" / "missing.fastq"
    try:
        trim.create_trim_zip(task_id, r1, missing, [ADAPTER], tasks_dir=tmp)
        record("read check: missing source file rejected", False)
    except trim.TrimSourceError as exc:
        record("read check: missing source file rejected",
               "cannot read source" in str(exc) or "No such file" in str(exc),
               str(exc))
    record("read check: missing file request leaves no trim dir",
           dirs_under(tmp / task_id / "trims") == [])

    # one extra record on one side
    put(clean1 + make_record("b", b"C" * 60), clean2)
    try:
        trim.create_trim_zip(task_id, r1, r2, [ADAPTER], tasks_dir=tmp)
        record("read check: unpaired tail rejected", False)
    except trim.TrimSourceError as exc:
        record("read check: unpaired tail rejected",
               "more records" in str(exc), str(exc))

    # read-length limits at the trim layer
    put(make_record("a", b"C" * 39), make_record("a", b"G" * 39))
    try:
        trim.create_trim_zip(task_id, r1, r2, [ADAPTER], tasks_dir=tmp)
        record("limit: 39-base read rejected", False)
    except trim.TrimReadLengthError as exc:
        record("limit: 39-base read rejected", "39 bases" in str(exc), str(exc))

    put(make_record("a", b"C" * 1001), make_record("a", b"G" * 1001))
    try:
        trim.create_trim_zip(task_id, r1, r2, [ADAPTER], tasks_dir=tmp)
        record("limit: 1001-base read rejected", False)
    except trim.TrimReadLengthError as exc:
        record("limit: 1001-base read rejected", "1001 bases" in str(exc), str(exc))

    # exactly the limits are accepted: 40 and 1000 bases
    put(make_record("a", b"C" * 40), make_record("a", b"G" * 40))
    result = trim.create_trim_zip(task_id, r1, r2, [ADAPTER], tasks_dir=tmp)
    record("limit: 40-base read accepted", result.stats.retained_pairs == 1)
    record("limit: 40-base read unchanged",
           read_zip(result.path)[1] == make_record("a", b"C" * 40))

    put(make_record("a", b"C" * 1000), make_record("a", b"G" * 1000))
    result = trim.create_trim_zip(task_id, r1, r2, [ADAPTER], tasks_dir=tmp)
    record("limit: 1000-base read accepted", result.stats.retained_pairs == 1)

    # 200k pairs are accepted, 200_001 are rejected at the pair-count limit.
    n = config.TRIM_MAX_PAIRS
    put(b"", b"")
    with r1.open("wb") as f1, r2.open("wb") as f2:
        for i in range(n + 1):
            f1.write(b"@r%d\n%s\n+\n%s\n" % (i, b"A" * 40, b"~" * 40))
            f2.write(b"@r%d\n%s\n+\n%s\n" % (i, b"T" * 40, b"~" * 40))
    try:
        trim.create_trim_zip(task_id, r1, r2, [ADAPTER], tasks_dir=tmp)
        record("limit: 200001 pairs rejected", False)
    except trim.TrimLimitError as exc:
        record("limit: 200001 pairs rejected", "more than" in str(exc), str(exc))
    record("limit: over-limit request leaves no trim dir",
           dirs_under(tmp / task_id / "trims") == [])

    # 200k pairs exactly: accepted with correct counts.
    put(b"", b"")
    with r1.open("wb") as f1, r2.open("wb") as f2:
        for i in range(n):
            f1.write(b"@r%d\n%s\n+\n%s\n" % (i, b"A" * 40, b"~" * 40))
            f2.write(b"@r%d\n%s\n+\n%s\n" % (i, b"T" * 40, b"~" * 40))
    result = trim.create_trim_zip(task_id, r1, r2, [ADAPTER], tasks_dir=tmp)
    record("limit: exactly 200000 pairs accepted",
           result.stats.total_pairs == n
           and result.stats.retained_pairs == n
           and result.stats.dropped_pairs == 0,
           f"{result.stats}")
    summary, _, _ = read_zip(result.path)
    record("limit: 200000-pair summary counts match",
           summary["total_pairs"] == n and summary["retained_pairs"] == n)


async def api_tests(client: httpx.AsyncClient, tmp: Path) -> None:
    # ---- happy path with trimming, dropping, order and byte fidelity
    p1_r1 = b"C" * 50 + ADAPTER[:12].encode()
    p1_r2 = b"G" * 45 + b"AGATCGGT"          # 8-base overlap, 1 mismatch
    p2_r1 = b"A" * 60
    p2_r2 = b"T" * 60
    p3_r1 = b"T" * 18 + ADAPTER.encode()      # cut to 18 < 20 -> drop the pair
    p3_r2 = b"G" * 60
    p4_r1 = b"C" * 43 + b"AGATCGG"            # 7 matching bases -> no cut
    p4_r2 = b"G" * 40 + ADAPTER[:9].encode()     # 9 aligned bases -> cut
    recs1 = [
        make_record("p1", p1_r1, suffix="/1", desc=b"sample=x"),
        make_record("p2", p2_r1, suffix="/1"),
        make_record("p3", p3_r1, suffix="/1"),
        make_record("p4", p4_r1, suffix="/1"),
    ]
    recs2 = [
        make_record("p1", p1_r2, suffix="/2", desc=b"sample=x"),
        make_record("p2", p2_r2, suffix="/2"),
        make_record("p3", p3_r2, suffix="/2"),
        make_record("p4", p4_r2, suffix="/2"),
    ]
    r = await upload(client, fastq(recs1), fastq(recs2))
    task = await wait_terminal(client, r.json()["id"])
    record("trim source task succeeded", task["status"] == "succeeded",
           task.get("message", ""))
    record("task exposes pair_count", task["pair_count"] == 4,
           str(task.get("pair_count")))
    tid = task["id"]
    src1 = config.TASKS_DIR / tid / "output" / "R1.fastq"
    src2 = config.TASKS_DIR / tid / "output" / "R2.fastq"
    before = (sha(src1), sha(src2))

    resp = await client.post(f"/tasks/{tid}/trims",
                             json={"task_id": tid, "adapters": [ADAPTER]})
    record("trim POST returns 201", resp.status_code == 201, str(resp.status_code))
    meta = resp.json()
    record("trim metadata counts",
           (meta["total_pairs"], meta["retained_pairs"], meta["dropped_pairs"])
           == (4, 3, 1), str(meta))
    record("trim metadata download url",
           meta["download_url"] == f"/tasks/{tid}/trims/{meta['trim_id']}/download")

    zr = await client.get(meta["download_url"])
    record("trim ZIP download 200", zr.status_code == 200, str(zr.status_code))
    record("trim ZIP content type",
           zr.headers.get("content-type", "").startswith("application/zip"),
           zr.headers.get("content-type", ""))
    zip_path = config.TASKS_DIR / tid / "trims" / meta["trim_id"] / trim.ZIP_NAME
    summary, z1, z2 = read_zip(zip_path)
    record("ZIP summary retained/dropped",
           (summary["total_pairs"], summary["retained_pairs"],
            summary["dropped_pairs"]) == (4, 3, 1), json.dumps(summary))
    r1_out = parse_fastq_bytes(z1)
    r2_out = parse_fastq_bytes(z2)
    record("trimmed output keeps R1 order", [x[0] for x in r1_out] == ["p1", "p2", "p4"],
           str([x[0] for x in r1_out]))
    record("trimmed R2 aligned with R1", [x[0] for x in r2_out] == ["p1", "p2", "p4"])
    record("sequence trimmed at adapter start",
           [len(x[2]) for x in r1_out] == [50, 60, 50],
           str([len(x[2]) for x in r1_out]))
    record("R2 trimmed with one mismatch",
           [len(x[2]) for x in r2_out] == [45, 60, 40],
           str([len(x[2]) for x in r2_out]))
    record("quality trimmed in lockstep",
           all(len(x[2]) == len(x[4]) for x in r1_out + r2_out))
    record("header and plus lines preserved",
           r1_out[0][1] == recs1[0].split(b"\n")[0] + b"\n"
           and r1_out[0][3] == b"+")
    record("source R1 unchanged", sha(src1) == before[0])
    record("source R2 unchanged", sha(src2) == before[1])
    record("get trim metadata",
           (await client.get(f"/tasks/{tid}/trims/{meta['trim_id']}")).json()
           ["retained_pairs"] == 3)
    record("trim listing under other task -> 404",
           (await client.get(f"/tasks/{uuid.uuid4().hex}/trims/"
                             f"{meta['trim_id']}")).status_code == 404)

    # lowercase adapters are normalized
    resp = await client.post(f"/tasks/{tid}/trims",
                             json={"task_id": tid,
                                   "adapters": [ADAPTER.lower()]})
    record("lowercase adapter accepted and normalized",
           resp.status_code == 201
           and resp.json()["adapters"] == [ADAPTER], str(resp.status_code))

    # ---- request validation
    bad_requests = {
        "no adapters": [],
        "five adapters": [ADAPTER] * 5,
        "15-base adapter": ["A" * 15],
        "33-base adapter": ["A" * 33],
        "non-ACGT adapter": ["A" * 16 + "N"],
        "adapter with space": ["A" * 16 + " "],
    }
    for name, adapters in bad_requests.items():
        resp = await client.post(f"/tasks/{tid}/trims",
                                 json={"task_id": tid, "adapters": adapters})
        record(f"adapter validation: {name} -> 422", resp.status_code == 422,
               str(resp.status_code))
    resp = await client.post(f"/tasks/{tid}/trims",
                             json={"task_id": "different", "adapters": [ADAPTER]})
    record("body/path task id mismatch -> 422", resp.status_code == 422,
           str(resp.status_code))
    unknown_id = uuid.uuid4().hex
    resp = await client.post(f"/tasks/{unknown_id}/trims",
                             json={"task_id": unknown_id,
                                   "adapters": [ADAPTER]})
    record("trim unknown task -> 404", resp.status_code == 404,
           str(resp.status_code))

    # ---- read-length limits and non-succeeded tasks through the API
    short1 = fastq([make_record("s", b"A" * 39)])
    short2 = fastq([make_record("s", b"T" * 39)])
    r = await upload(client, short1, short2)
    short_task = await wait_terminal(client, r.json()["id"])
    resp = await client.post(
        f"/tasks/{short_task['id']}/trims",
        json={"task_id": short_task["id"], "adapters": [ADAPTER]})
    record("API: 39-base read -> 422", resp.status_code == 422,
           str(resp.status_code))

    long1 = fastq([make_record("l", b"A" * 1001)])
    long2 = fastq([make_record("l", b"T" * 1001)])
    r = await upload(client, long1, long2)
    long_task = await wait_terminal(client, r.json()["id"])
    resp = await client.post(
        f"/tasks/{long_task['id']}/trims",
        json={"task_id": long_task["id"], "adapters": [ADAPTER]})
    record("API: 1001-base read -> 422", resp.status_code == 422,
           str(resp.status_code))

    # failed task cannot be trimmed
    bad = io.BytesIO(b"@x\nAAAA\n+\n")
    r = await upload(client, bad, fastq([make_record("x", b"aaaa")]))
    failed_task = await wait_terminal(client, r.json()["id"])
    resp = await client.post(
        f"/tasks/{failed_task['id']}/trims",
        json={"task_id": failed_task["id"], "adapters": [ADAPTER]})
    record("API: failed task -> 409", resp.status_code == 409,
           str(resp.status_code))

    # ---- source damage: no partial ZIP, no source modification
    damaged_tid = tid
    src2.write_bytes(src2.read_bytes()[:-10])  # corrupt the published tail
    damaged_hash = sha(src2)
    before_dirs = dirs_under(config.TASKS_DIR / damaged_tid / "trims")
    resp = await client.post(f"/tasks/{damaged_tid}/trims",
                             json={"task_id": damaged_tid,
                                   "adapters": [ADAPTER]})
    record("API: corrupt source tail -> 409", resp.status_code == 409,
           str(resp.status_code))
    record("API: corrupt source leaves no partial ZIP",
           dirs_under(config.TASKS_DIR / damaged_tid / "trims") == before_dirs)
    record("API: corrupt source not modified by trim", sha(src2) == damaged_hash)

    # missing result file
    missing_tid = long_task["id"]
    (config.TASKS_DIR / missing_tid / "output" / "R2.fastq").unlink()
    resp = await client.post(f"/tasks/{missing_tid}/trims",
                             json={"task_id": missing_tid,
                                   "adapters": [ADAPTER]})
    record("API: missing source file -> 409", resp.status_code == 409,
           str(resp.status_code))
    record("API: missing source leaves no trim dir",
           dirs_under(config.TASKS_DIR / missing_tid / "trims") == [])

    # ---- concurrent requests for the same task: one fails, one succeeds
    conc_recs1 = [make_record(f"c{i}", b"C" * 60) for i in range(20)]
    conc_recs2 = [make_record(f"c{i}", b"G" * 60) for i in range(20)]
    r = await upload(client, fastq(conc_recs1), fastq(conc_recs2))
    conc_task = await wait_terminal(client, r.json()["id"])
    conc_tid = conc_task["id"]
    conc_src1 = config.TASKS_DIR / conc_tid / "output" / "R1.fastq"
    conc_src2 = config.TASKS_DIR / conc_tid / "output" / "R2.fastq"
    conc_before = (sha(conc_src1), sha(conc_src2))

    real_build_zip = trim._build_zip

    def flaky_build_zip(zip_path, trimmed_r1, trimmed_r2, summary):
        if "TTTTTTTTTTTTTTTT" in summary["adapters"]:
            raise OSError("injected interruption while writing archive")
        return real_build_zip(zip_path, trimmed_r1, trimmed_r2, summary)

    trim._build_zip = flaky_build_zip
    try:
        good_req = client.post(f"/tasks/{conc_tid}/trims",
                               json={"task_id": conc_tid,
                                     "adapters": [ADAPTER]})
        bad_req = client.post(f"/tasks/{conc_tid}/trims",
                              json={"task_id": conc_tid,
                                    "adapters": ["T" * 16]})
        good_resp, bad_resp = await asyncio.gather(good_req, bad_req)
    finally:
        trim._build_zip = real_build_zip

    record("concurrent: good request succeeds", good_resp.status_code == 201,
           f"{good_resp.status_code} {good_resp.text[:120]}")
    record("concurrent: interrupted request fails",
           bad_resp.status_code == 500, str(bad_resp.status_code))
    good_meta = good_resp.json()
    good_zip = (config.TASKS_DIR / conc_tid / "trims" / good_meta["trim_id"]
                / trim.ZIP_NAME)
    summary, z1, z2 = read_zip(good_zip)
    record("concurrent: surviving ZIP intact",
           summary["retained_pairs"] == 20
           and len(parse_fastq_bytes(z1)) == 20
           and len(parse_fastq_bytes(z2)) == 20)
    record("concurrent: failed request left no directory",
           good_meta["trim_id"] in dirs_under(
               config.TASKS_DIR / conc_tid / "trims")
           and len(dirs_under(config.TASKS_DIR / conc_tid / "trims")) == 1,
           str(dirs_under(config.TASKS_DIR / conc_tid / "trims")))
    record("concurrent: source files untouched",
           (sha(conc_src1), sha(conc_src2)) == conc_before)
    dl = await client.get(good_meta["download_url"])
    record("concurrent: surviving ZIP still downloadable",
           dl.status_code == 200 and len(dl.content) > 0)
    record("concurrent: failed trim id not registered",
           (await client.get(
               f"/tasks/{conc_tid}/trims/{uuid.uuid4().hex}")).status_code == 404)

    # ---- restart cleanup of interrupted trim directories
    stale = config.TASKS_DIR / conc_tid / "trims" / "staleinterrupted"
    stale.mkdir(parents=True)
    (stale / "leftover.part").write_text("partial")
    removed = worker.recover_trim_dirs()
    record("restart removes unregistered trim dir", removed == 1
           and not stale.exists(), str(removed))
    record("restart keeps registered trim result", good_zip.is_file())
    record("restart keeps source results", conc_src1.is_file() and conc_src2.is_file())

    # ---- docs and OpenAPI
    docs = await client.get("/docs")
    record("/docs is served", docs.status_code == 200
           and b"swagger" in docs.content.lower(), str(docs.status_code))
    spec = (await client.get("/openapi.json")).json()
    record("OpenAPI lists trim endpoints",
           "/tasks/{task_id}/trims" in spec["paths"]
           and "/tasks/{task_id}/trims/{trim_id}/download" in spec["paths"])


async def main_async() -> int:
    tmp_root = Path(tempfile.mkdtemp(prefix="fastq-trim-",
                                     dir=str(ROOT / "tests" / "output")))
    data_dir = tmp_root / "data"
    os.environ["FASTQ_DATA_DIR"] = str(data_dir)
    config.DATA_DIR = data_dir
    config.DB_PATH = data_dir / "tasks.db"
    config.TASKS_DIR = data_dir / "tasks"
    config.PROGRESS_UPDATE_EVERY_BYTES = 64 * 1024

    db._initialized = False
    db.init_db()
    worker.recover_on_startup()

    try:
        matcher_unit_tests()
        direct_trim_tests(tmp_root / "direct-root")

        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(transport=transport,
                                     base_url="http://test",
                                     timeout=120) as client:
            await api_tests(client, tmp_root)
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)

    failed = [name for name, ok, _ in RESULTS if not ok]
    print("\n%d/%d checks passed" % (len(RESULTS) - len(failed), len(RESULTS)))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(asyncio.run(main_async()))
