import io
import time
import zipfile

import pytest

pytest.importorskip("fastapi", reason="the GUI tests need the [web] extra (fastapi, uvicorn, httpx)")
pytest.importorskip("httpx", reason="the GUI tests need httpx")

from fastapi.testclient import TestClient  # noqa: E402

from msldprep.gui.app import create_app  # noqa: E402

from .test_gui_service import ANALOG_TEXT, _reference_text  # noqa: E402

TOKEN = "test-token"


@pytest.fixture
def client():
    return TestClient(create_app(TOKEN), headers={"X-Token": TOKEN})


def _wait_for_job(client):
    for _ in range(100):
        state = client.get("/api/state").json()
        if not state["job"]["running"]:
            return state
        time.sleep(0.1)
    raise AssertionError("the alignment job did not finish")


def test_every_route_needs_the_token():
    bare = TestClient(create_app(TOKEN))
    assert bare.get("/").status_code == 403
    assert bare.get("/api/state").status_code == 403
    assert bare.get("/api/state", headers={"X-Token": "wrong"}).status_code == 403
    assert bare.get("/", params={"token": TOKEN}).status_code == 200
    assert bare.get("/api/state", params={"token": TOKEN}).status_code == 200


def test_the_page_and_its_scripts_are_served_without_outside_links(client):
    page = client.get("/").text
    assert "/static/3Dmol-min.js" in page and "/static/app.js" in page
    assert "http://" not in page and "https://" not in page
    assert client.get("/static/app.js").status_code == 200
    assert client.get("/static/3Dmol-min.js").status_code == 200


def test_whole_flow_through_the_api(client):
    assert client.put("/api/reference", json={"name": "ref", "text": _reference_text()}).status_code == 200
    assert client.put("/api/analogs", json={"text": ANALOG_TEXT + "\nCCO ethanol\n"}).status_code == 200

    assert client.post("/api/align").status_code == 202
    state = _wait_for_job(client)
    assert [a["status"] for a in state["analogs"]] == ["aligned"] * 3 + ["failed"]

    view = client.get("/api/analogs/10749").json()
    assert view["status"] == "aligned" and view["core_atoms"] and "molblock" in view

    again = client.post("/api/analogs/10749/realign", json={"seed": 5}).json()
    assert again["seed"] == 5

    decided = client.post("/api/analogs/10749/decision", json={"decision": "accepted"}).json()
    assert decided["analogs"][0]["decision"] == "accepted"

    sdf = client.get("/api/export", params={"which": "accepted"})
    assert sdf.status_code == 200 and sdf.text.count("$$$$") == 2
    assert 'filename="aligned_accepted.sdf"' in sdf.headers["content-disposition"]
    mol2 = client.get("/api/export", params={"format": "mol2", "which": "accepted", "reference": "false"})
    assert mol2.status_code == 200 and mol2.text.count("@<TRIPOS>MOLECULE") == 1
    archive = client.get("/api/export", params={"format": "mol2", "which": "accepted", "per_file": "true"})
    assert archive.headers["content-type"] == "application/zip"
    with zipfile.ZipFile(io.BytesIO(archive.content)) as zf:
        assert zf.namelist() == ["ref.mol2", "10749.mol2"]
    assert len(client.get("/api/overlay").json()) == 3


def test_errors_come_back_as_messages_not_crashes(client):
    assert client.post("/api/align").status_code == 422  # nothing loaded
    junk = client.put("/api/reference", json={"name": "x", "text": "junk"})
    assert junk.status_code == 422 and "Expected one molecule" in junk.json()["detail"]
    assert client.get("/api/analogs/nope").status_code == 404
    assert client.put("/api/params", json={"timeout": 0}).status_code == 422
    assert client.put("/api/analogs", json={"text": "# nothing here"}).status_code == 422
    assert client.get("/api/export").status_code == 422  # no reference
    assert client.get("/api/export", params={"format": "pdb"}).status_code == 422


def test_start_over(client):
    client.put("/api/reference", json={"name": "ref", "text": _reference_text()})
    state = client.post("/api/reset").json()
    assert state["reference"] is None


def test_the_matching_settings_round_trip_through_the_api(client):
    assert client.get("/api/state").json()["params"]["ignore_bond_orders"] is True
    state = client.put("/api/params", json={"ignore_bond_orders": False, "ring_matches_ring_only": True}).json()
    assert state["params"]["ignore_bond_orders"] is False and state["params"]["ring_matches_ring_only"] is True
