#!/usr/bin/env python3
"""
test_harness_server.py — smoke tests for the four-column generation concurrency
in harness_server.py. No network, no gpt-image-2 calls: exercises only the pure
job-state helpers (_sanitize_tag, _reserve, _start).

Run:  python3 test_harness_server.py   ->  prints "all tests passed" or asserts.

FAILURE MODE UNDER TEST — index collision. Four generate columns run
concurrently. Filenames are namespaced by tag (`gpt_<tag>_NNN.png`), and the
next index is handed out by _reserve(). If _reserve were replaced by a bare
_disk_next() read (as the single-worker version effectively was), two jobs with
the same tag would both compute the same start index and silently overwrite each
other's PNGs. test_reserve_non_overlapping guards exactly that.
"""
import os
import sys
import time
import json
import tempfile
import threading

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness_server as hs  # noqa: E402  (path insert must precede import)


def test_sanitize_tag():
    assert hs._sanitize_tag("room1") == "room1"
    assert hs._sanitize_tag("Room 1") == "room_1"
    assert hs._sanitize_tag("a--b__c") == "a_b_c"
    assert hs._sanitize_tag("") == "gen"
    assert hs._sanitize_tag(None) == "gen"
    assert hs._sanitize_tag("__weird__") == "weird"


def test_reserve_non_overlapping():
    """Consecutive reservations for one prefix must not overlap (the collision guard)."""
    prefix = "gpt_zzztest_"  # unlikely to exist on disk
    hs.RESERVED.pop(prefix, None)
    a = hs._reserve(prefix, 3)      # hands out a, a+1, a+2
    b = hs._reserve(prefix, 2)      # must start at or after a+3
    c = hs._reserve(prefix, 1)
    assert b >= a + 3, (a, b)
    assert c >= b + 2, (b, c)
    hs.RESERVED.pop(prefix, None)


def test_reserve_distinct_tags_independent():
    for p in ("gpt_taga_", "gpt_tagb_"):
        hs.RESERVED.pop(p, None)
    a = hs._reserve("gpt_taga_", 4)
    b = hs._reserve("gpt_tagb_", 4)
    # distinct prefixes each start fresh; no interference
    assert a >= 1 and b >= 1
    for p in ("gpt_taga_", "gpt_tagb_"):
        hs.RESERVED.pop(p, None)


def test_start_busy_reject_and_concurrent_slots():
    for s in ("t1", "t2"):
        hs.JOBS.pop(s, None)
    gate = threading.Event()

    def block():
        gate.wait(2)

    assert hs._start("t1", "test", block, 1) is True       # slot t1 now active
    assert hs._start("t1", "test", block, 1) is False       # same slot busy -> rejected
    assert hs._start("t2", "test", block, 1) is True        # different slot -> allowed concurrently
    assert hs.JOBS["t1"]["active"] is True
    assert hs.JOBS["t2"]["active"] is True
    gate.set()
    time.sleep(0.05)
    # bare block() target doesn't clear active (only _run_* do); tidy up for isolation
    for s in ("t1", "t2"):
        hs.JOBS.pop(s, None)


def test_status_idle_default_shape():
    # the shape /api/status?slot=<unknown> returns
    assert hs._IDLE["active"] is False
    assert set(hs._IDLE) >= {"active", "kind", "done", "total", "outputs", "error", "tag"}


# --- commit-to-room ("Send to room") ---------------------------------------
# FAILURE MODE UNDER TEST — a broken commit. "Send to room" must copy the closed
# base AND its `_open` partner under STABLE names (scene.png / scene_open.png) into
# a directory INSIDE the escape_rooms tree, carrying the re-keyed wrap + hotspots —
# and must reject a roomDir that escapes the tree. A half-committed or escaping room
# would silently break the door swap or write outside the site.

def _png(path):
    with open(path, "wb") as f:
        f.write(b"\x89PNG\r\n")


def test_commit_room_full_pair():
    with tempfile.TemporaryDirectory() as scene, tempfile.TemporaryDirectory() as root:
        hs.SCENE, hs.ESCAPE_ROOT = scene, root
        _png(os.path.join(scene, "gpt_r_3.png"))
        _png(os.path.join(scene, "gpt_r_3_open.png"))
        json.dump({"gpt_r_3.png": {"haov": 360, "vaov": 90, "hfov": 120, "vOffset": -5, "pitch": -6}},
                  open(os.path.join(scene, "wrap.json"), "w"))
        json.dump({"image": "gpt_r_3.png", "haov": 360, "vaov": 90, "vOffset": -5,
                   "hotspots": [{"id": "laptop", "box": [0.4, 0.5, 0.5, 0.6], "action": "puzzle"}]},
                  open(os.path.join(scene, "hotspots.json"), "w"))
        written, rd = hs._commit_room("gpt_r_3.png", "data_vis/case/room1")
        dest = os.path.join(root, "data_vis", "case", "room1")
        assert rd == "data_vis/case/room1"
        assert set(written) == {"scene.png", "scene_open.png", "wrap.json", "hotspots.json"}
        assert os.path.exists(os.path.join(dest, "scene.png"))
        assert os.path.exists(os.path.join(dest, "scene_open.png"))       # door pair kept together
        assert json.load(open(os.path.join(dest, "wrap.json")))["scene.png"]["vaov"] == 90  # re-keyed
        assert json.load(open(os.path.join(dest, "hotspots.json")))["image"] == "scene.png"  # re-keyed


def test_commit_room_no_open_partner():
    with tempfile.TemporaryDirectory() as scene, tempfile.TemporaryDirectory() as root:
        hs.SCENE, hs.ESCAPE_ROOT = scene, root
        _png(os.path.join(scene, "a.png"))
        written, _ = hs._commit_room("a.png", "d/r")
        assert "scene.png" in written and "scene_open.png" not in written


def test_commit_room_rejects_escape_and_missing():
    with tempfile.TemporaryDirectory() as scene, tempfile.TemporaryDirectory() as root:
        hs.SCENE, hs.ESCAPE_ROOT = scene, root
        _png(os.path.join(scene, "a.png"))
        for bad in ("../../etc", "..", ""):
            try:
                hs._commit_room("a.png", bad)
                raise AssertionError(f"should have rejected roomDir={bad!r}")
            except ValueError:
                pass
        try:
            hs._commit_room("nope.png", "d/r")
            raise AssertionError("should have rejected a missing image")
        except ValueError:
            pass


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"  ok  {t.__name__}")
    print(f"all tests passed ({len(tests)})")
