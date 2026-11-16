"""End-to-end API and workflow tests."""
from __future__ import annotations

import io
import time
from pathlib import Path

from tests.conftest import DEMO_DIR


def upload_files(client, names):
    files = []
    handles = []
    try:
        for name in names:
            handle = open(DEMO_DIR / name, "rb")
            handles.append(handle)
            files.append(("files", (name, handle, "text/plain")))
        return client.post("/api/batches", files=files)
    finally:
        for handle in handles:
            handle.close()


def wait_for_terminal(client, batch_id, timeout=10.0):
    deadline = time.time() + timeout
    while time.time() < deadline:
        data = client.get(f"/api/batches/{batch_id}").json()
        gen = data["latest_generation"]
        if gen and gen["status"] in ("done", "failed", "superseded", "aborted"):
            return data
        time.sleep(0.05)
    raise AssertionError("generation did not reach a terminal state")


def test_health(client):
    assert client.get("/api/health").json() == {"status": "ok"}


def test_upload_count_must_be_2_to_3(client):
    handle = open(DEMO_DIR / "batch_a.vcf", "rb")
    resp = client.post(
        "/api/batches",
        files=[("files", ("batch_a.vcf", handle, "text/plain"))],
    )
    handle.close()
    assert resp.status_code == 400
    assert "2~3" in resp.json()["detail"]


def test_full_demo_workflow(client, app_env):
    resp = upload_files(client, ["batch_a.vcf", "batch_b.vcf", "batch_c.vcf"])
    assert resp.status_code == 200, resp.text
    data = resp.json()
    batch_id = data["id"]
    assert data["batch_valid"] is True
    assert data["reference_version"] == "GRCh38.demo.1.0"

    # Every sample has a source group; S2 (A/B) and S3 (B/C) are duplicated.
    groups = {g["sample_name"]: g for g in data["dup_groups"]}
    assert set(groups) == {"S1", "S2", "S3", "S4"}
    dup = {name: g for name, g in groups.items() if g["is_duplicate"]}
    assert set(dup) == {"S2", "S3"}
    # Unique-source samples resolve automatically; duplicated ones block.
    assert not data["ready_to_generate"]
    assert set(data["unresolved_samples"]) == {"S2", "S3"}

    # Generation is blocked while duplicate samples are unresolved.
    blocked = client.post(f"/api/batches/{batch_id}/generate")
    assert blocked.status_code == 409
    assert "尚未选择" in blocked.json()["detail"]

    def choose(sample, file_name_part):
        group = dup[sample]
        chosen = next(
            c["file_id"] for c in group["choices"] if file_name_part in c["file_name"]
        )
        form = {"file_id": chosen}
        r = client.put(
            f"/api/batches/{batch_id}/selection/{group['id']}", data=form
        )
        assert r.status_code == 200, r.text
        return r.json()

    data = choose("S2", "batch_a")
    data = choose("S3", "batch_c")
    assert data["ready_to_generate"] is True
    assert data["selection_version"] == 3  # initial 1 + two changes

    gen_resp = client.post(f"/api/batches/{batch_id}/generate")
    assert gen_resp.status_code == 200
    final = wait_for_terminal(client, batch_id)
    gen = final["latest_generation"]
    assert gen["status"] == "done", gen["error"]
    assert gen["downloadable"] is True

    download = client.get(f"/api/batches/{batch_id}/download")
    assert download.status_code == 200
    content = download.content.decode()
    # Cohort samples must include S1-S4 exactly once each.
    header_line = next(line for line in content.splitlines() if line.startswith("#CHROM"))
    header_samples = header_line.split("\t")[9:]
    assert sorted(header_samples) == ["S1", "S2", "S3", "S4"]
    assert len(header_samples) == len(set(header_samples))
    # Union site chr1:100 has 2 ALT alleles G,T.
    site100 = next(line for line in content.splitlines() if line.startswith("chr1\t100"))
    assert site100.split("\t")[4] == "G,T"
    # S2 came from batch_a where chr2:300 is missing -> ./. not 0/0.
    site300 = next(line for line in content.splitlines() if line.startswith("chr2\t300"))
    columns = site300.split("\t")
    s2_col = columns[9 + header_samples.index("S2")]
    assert s2_col.startswith("./.")


def test_ref_conflict_fails_whole_batch(client):
    resp = upload_files(client, ["batch_a.vcf", "bad_ref_conflict.vcf"])
    assert resp.status_code == 200
    data = resp.json()
    assert data["batch_valid"] is False
    assert any("REF 冲突" in error for error in data["batch_errors"])
    assert data["ready_to_generate"] is False
    blocked = client.post(f"/api/batches/{data['id']}/generate")
    assert blocked.status_code == 409


def test_invalid_single_file_reports_reason(client):
    # Upload the structurally valid batch A plus a truncated file.
    bad = io.BytesIO(b"##fileformat=VCFv4.2\ngarbage content without header columns\n")
    with open(DEMO_DIR / "batch_a.vcf", "rb") as good:
        resp = client.post(
            "/api/batches",
            files=[
                ("files", ("batch_a.vcf", good, "text/plain")),
                ("files", ("broken.vcf", bad, "text/plain")),
            ],
        )
    assert resp.status_code == 200
    data = resp.json()
    broken = next(f for f in data["files"] if f["original_name"] == "broken.vcf")
    assert broken["valid"] is False
    assert broken["error"]


def test_selection_change_invalidates_completed_result(client):
    data = upload_files(client, ["batch_a.vcf", "batch_b.vcf"]).json()
    batch_id = data["id"]
    dup = {g["sample_name"]: g for g in data["dup_groups"]}
    group = dup["S2"]
    a_choice = next(c["file_id"] for c in group["choices"] if "batch_a" in c["file_name"])
    b_choice = next(c["file_id"] for c in group["choices"] if "batch_b" in c["file_name"])

    client.put(f"/api/batches/{batch_id}/selection/{group['id']}",
               data={"file_id": a_choice})
    client.post(f"/api/batches/{batch_id}/generate")
    final = wait_for_terminal(client, batch_id)
    assert final["latest_generation"]["status"] == "done"
    output_path = Path(final["latest_generation"]["output_path"])
    assert output_path.is_file()

    # Change selection: the done generation must be superseded immediately and
    # its output removed; download must be refused.
    changed = client.put(
        f"/api/batches/{batch_id}/selection/{group['id']}",
        data={"file_id": b_choice},
    ).json()
    assert changed["latest_generation"]["status"] == "superseded"
    assert not output_path.exists()
    refused = client.get(f"/api/batches/{batch_id}/download")
    assert refused.status_code == 409

    # Changing back must not restore the old generation.
    back = client.put(
        f"/api/batches/{batch_id}/selection/{group['id']}",
        data={"file_id": a_choice},
    ).json()
    assert back["latest_generation"]["status"] == "superseded"
    assert back["selection_version"] == 4
    refused2 = client.get(f"/api/batches/{batch_id}/download")
    assert refused2.status_code == 409

    # A fresh generation under the current selection works and is downloadable.
    client.post(f"/api/batches/{batch_id}/generate")
    final2 = wait_for_terminal(client, batch_id)
    assert final2["latest_generation"]["status"] == "done"
    assert client.get(f"/api/batches/{batch_id}/download").status_code == 200


def test_selection_change_aborts_inflight_generation(client):
    data = upload_files(client, ["batch_a.vcf", "batch_b.vcf", "batch_c.vcf"]).json()
    batch_id = data["id"]
    dup = {g["sample_name"]: g for g in data["dup_groups"]}
    # Resolve both groups first.
    for sample, part in (("S2", "batch_a"), ("S3", "batch_b")):
        group = dup[sample]
        choice = next(c["file_id"] for c in group["choices"] if part in c["file_name"])
        client.put(f"/api/batches/{batch_id}/selection/{group['id']}",
                   data={"file_id": choice})

    client.post(f"/api/batches/{batch_id}/generate")
    # Immediately change a selection while the worker runs.
    group = dup["S2"]
    b_choice = next(c["file_id"] for c in group["choices"] if "batch_b" in c["file_name"])
    client.put(f"/api/batches/{batch_id}/selection/{group['id']}",
               data={"file_id": b_choice})

    final = wait_for_terminal(client, batch_id)
    gen = final["latest_generation"]
    assert gen["status"] == "superseded"
    assert gen["output_path"] is None or not Path(gen["output_path"]).exists()
    # No partial file is left behind.
    from app import config
    partials = list(Path(config.TMP_DIR).glob("*.partial"))
    assert partials == []


def test_restart_marks_running_generation_failed_and_blocks_download(client, app_env):
    data = upload_files(client, ["batch_a.vcf", "batch_b.vcf"]).json()
    batch_id = data["id"]
    group = next(g for g in data["dup_groups"] if g["sample_name"] == "S2")
    choice = group["choices"][0]["file_id"]
    client.put(f"/api/batches/{batch_id}/selection/{group['id']}",
               data={"file_id": choice})

    # Insert a generation directly in running state (simulating a crash).
    db = app_env["db"]
    config = app_env["config"]
    with db.get_db() as conn:
        gen_id = db.create_generation(conn, batch_id, 1)
        db.mark_generation(conn, gen_id, config.STATUS_RUNNING, started=True)

    # Re-run startup recovery as a process restart would.
    db.init_db()
    refreshed = client.get(f"/api/batches/{batch_id}").json()
    gen = refreshed["latest_generation"]
    assert gen["status"] == "failed"
    assert "重启" in gen["error"]
    refused = client.get(f"/api/batches/{batch_id}/download")
    assert refused.status_code == 409


def test_sample_source_displayed_per_file(client):
    data = upload_files(client, ["batch_a.vcf", "batch_b.vcf"]).json()
    names = {f["original_name"] for f in data["files"]}
    assert names == {"batch_a.vcf", "batch_b.vcf"}
    for file_row in data["files"]:
        assert file_row["samples"]
        assert file_row["record_count"] > 0


def test_gzipped_vcf_upload(client, tmp_path):
    from tests.bgzf_util import bgzip

    gz_path = tmp_path / "batch_a.vcf.gz"
    bgzip(str(DEMO_DIR / "batch_a.vcf"), str(gz_path))
    with open(gz_path, "rb") as gz, open(DEMO_DIR / "batch_b.vcf", "rb") as plain:
        resp = client.post(
            "/api/batches",
            files=[
                ("files", ("batch_a.vcf.gz", gz, "application/gzip")),
                ("files", ("batch_b.vcf", plain, "text/plain")),
            ],
        )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    gz_file = next(f for f in data["files"] if f["original_name"] == "batch_a.vcf.gz")
    assert gz_file["valid"] is True
    assert gz_file["record_count"] == 4
    assert data["batch_valid"] is True


def test_bad_extension_rejected(client):
    resp = client.post(
        "/api/batches",
        files=[
            ("files", ("a.vcf", open(DEMO_DIR / "batch_a.vcf", "rb"), "text/plain")),
            ("files", ("notes.txt", b"not a vcf", "text/plain")),
        ],
    )
    assert resp.status_code == 400
    assert "扩展名" in resp.json()["detail"]
