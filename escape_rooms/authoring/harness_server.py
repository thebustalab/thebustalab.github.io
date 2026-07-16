#!/usr/bin/env python3
"""
harness_server.py — in-repo authoring server for the gpt-image-2 escape-room scenes.

Replaces the retired DiT360 harness (was ~/dit360_bench/harness_server.py). We
committed to the gpt-image-2 pseudo-360 wrap (2026-07-15), so this server carries
NONE of the old GPU machinery: no run.sh, no GPU eviction, no conda-for-generation.
gpt-image-2 is a cloud API and Real-ESRGAN is light — the server just shells out to
the authoring scripts and reports progress.

Serves the alaska_pano/ folder on :8751 and adds:
  GET  /api/status[?slot=N]        -> job state for one slot, or all slots if omitted
  GET  /api/scenes                 -> list of scene/*.png (base candidates)
  POST /api/generate {prompt,n,quality,size,slot,tag}
                                   -> gpt-image-2 -> scene/gpt_<tag>_NNN.png (N candidates;
                                      size e.g. 1536x576 for a native wide panorama).
                                      Four columns author a whole room series in one go:
                                      each column is a `slot` with its own `tag`, and slots
                                      run concurrently (distinct tag = distinct filename
                                      prefix, so parallel jobs never collide on the index).
  POST /api/suggest  {prompt,n}    -> Claude Haiku prompt variants (AAPI)
  POST /api/save-wrap {image,haov,vaov,hfov,vOffset,pitch}
                                   -> writes scene/wrap.json (frozen viewer defaults)
  POST /api/save-hotspots {image,haov,vaov,vOffset,hotspots:[...]}
                                   -> writes scene/hotspots.json (viewer reads it)
  POST /api/dooropen {image,box,prompt}
                                   -> masked gpt-image-2 edit of the box region ->
                                      scene/<image>_open.png (door-swap target)
  POST /api/commit-room {image,roomDir}
                                   -> copy the chosen base + its _open partner into
                                      <escape_rooms>/<roomDir>/ under STABLE names
                                      (scene.png / scene_open.png), + that image's wrap
                                      + hotspots re-keyed. Keeps the door pair together
                                      and makes a self-contained, playable room dir.

Keys come from the environment (OPENAI_API_KEY for gpt, AAPI for Claude); launch
through a login shell so ~/.bashrc is sourced:
  bash -lic 'python3 <this>/harness_server.py'
"""
import os
import re
import json
import glob
import shutil
import threading
import subprocess
import http.server
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "alaska_pano"))  # served web root
SCENE = os.path.join(ROOT, "scene")
ESCAPE_ROOT = os.path.abspath(os.path.join(HERE, ".."))  # escape_rooms/ (commit targets live under here)
GEN = os.path.join(HERE, "generate_scene.py")
PORT = 8751

# In-process job state, one entry per slot (the four generate columns each own a
# slot; door-open uses its own "door" slot). Guarded by LOCK. Slots run
# concurrently, so filenames are namespaced by tag and indices are reserved
# atomically to keep parallel jobs from colliding.
JOBS = {}          # slot(str) -> {active,kind,done,total,outputs,error,tag}
RESERVED = {}      # filename prefix -> highest index handed out so far
LOCK = threading.Lock()

_IDLE = {"active": False, "kind": None, "done": 0, "total": 0,
         "outputs": [], "error": None, "tag": None}


def _sanitize_tag(tag):
    tag = re.sub(r"[^a-z0-9]+", "_", (tag or "").lower()).strip("_")
    return tag or "gen"


def _disk_next(prefix):
    n = 0
    for p in glob.glob(os.path.join(SCENE, prefix + "*.png")):
        m = re.search(re.escape(prefix) + r"(\d+)\.png$", os.path.basename(p))
        if m:
            n = max(n, int(m.group(1)))
    return n + 1


def _reserve(prefix, n):
    """Atomically hand out n consecutive indices for prefix (disk + in-flight)."""
    with LOCK:
        start = max(_disk_next(prefix), RESERVED.get(prefix, 0) + 1)
        RESERVED[prefix] = start + n - 1
        return start


def _run_generate(slot, tag, prompt, n, quality, size):
    os.makedirs(SCENE, exist_ok=True)
    prefix = f"gpt_{tag}_"
    ptmp = os.path.join(SCENE, f".prompt_{slot}.txt")
    with open(ptmp, "w", encoding="utf-8") as f:
        f.write(prompt)
    start = _reserve(prefix, n)
    for i in range(n):
        out = os.path.join(SCENE, f"{prefix}{start + i}.png")
        try:
            subprocess.run(["python3", GEN, "gen", "--prompt-file", ptmp,
                            "--out", out, "--quality", quality, "--size", size],
                           check=True, capture_output=True, text=True)
            with LOCK:
                JOBS[slot]["outputs"].append(os.path.basename(out))
                JOBS[slot]["done"] += 1
        except subprocess.CalledProcessError as e:
            with LOCK:
                JOBS[slot]["error"] = (e.stderr or e.stdout or str(e)).strip()[-500:]
            break
    with LOCK:
        JOBS[slot]["active"] = False


def _run_dooropen(slot, image, box, prompt):
    inp = os.path.join(SCENE, os.path.basename(image))
    stem = os.path.splitext(os.path.basename(image))[0]
    out = os.path.join(SCENE, stem + "_open.png")
    boxstr = ",".join(str(x) for x in box)
    try:
        subprocess.run(["python3", GEN, "dooropen", "--input", inp, "--box", boxstr,
                        "--prompt", prompt, "--out", out],
                       check=True, capture_output=True, text=True)
        with LOCK:
            JOBS[slot]["outputs"].append(os.path.basename(out))
            JOBS[slot]["done"] = 1
    except subprocess.CalledProcessError as e:
        with LOCK:
            JOBS[slot]["error"] = (e.stderr or e.stdout or str(e)).strip()[-500:]
    with LOCK:
        JOBS[slot]["active"] = False


def _start(slot, kind, target, total, tag=None):
    with LOCK:
        j = JOBS.get(slot)
        if j and j["active"]:
            return False
        JOBS[slot] = {"active": True, "kind": kind, "done": 0, "total": total,
                      "outputs": [], "error": None, "tag": tag}
    threading.Thread(target=target, daemon=True).start()
    return True


def suggest_prompts(base, n):
    key = os.environ.get("AAPI")
    if not key:
        raise RuntimeError("AAPI not set")
    meta = (
        f"Here is a base text-to-image prompt for gpt-image-2 describing an escape-room "
        f"scene. Write {n} varied alternatives that keep the same subject and mood but "
        f"explore genuinely different renders; use natural descriptive sentences (never "
        f"tag lists); front-load the most important elements; about 70-110 words each. "
        f"Return ONLY a JSON array of {n} strings and nothing else.\n\nBASE:\n{base}"
    )
    body = {"model": "claude-haiku-4-5-20251001", "max_tokens": 2500,
            "messages": [{"role": "user", "content": meta}]}
    req = urllib.request.Request(
        "https://api.anthropic.com/v1/messages",
        data=json.dumps(body).encode(),
        headers={"x-api-key": key, "anthropic-version": "2023-06-01",
                 "content-type": "application/json"})
    r = json.load(urllib.request.urlopen(req, timeout=60))
    text = r["content"][0]["text"].strip()
    if text.startswith("```"):
        text = text.strip("`")
        text = text[4:] if text.lower().startswith("json") else text
    return [str(p).strip() for p in json.loads(text)][:n]


def _commit_room(image, room_dir):
    """Copy a chosen closed base + its `_open` partner into a room directory under
    STABLE names (scene.png / scene_open.png), so the pair travels together and
    corresponds. Also carries that image's wrap + hotspots (re-keyed to scene.png)
    so the room dir is self-contained and playable. room_dir is relative to the
    escape_rooms/ tree. Returns (written_files, normalised_room_dir)."""
    room_dir = (room_dir or "").strip().strip("/").replace("\\", "/")
    if not room_dir:
        raise ValueError("roomDir is empty")
    dest = os.path.abspath(os.path.join(ESCAPE_ROOT, room_dir))
    if dest == ESCAPE_ROOT or not dest.startswith(ESCAPE_ROOT + os.sep):
        raise ValueError("roomDir must be a subdirectory inside the escape_rooms tree")
    src = os.path.join(SCENE, os.path.basename(image))
    if not os.path.exists(src):
        raise ValueError(f"image not found in scene/: {image}")
    os.makedirs(dest, exist_ok=True)
    written = []
    shutil.copyfile(src, os.path.join(dest, "scene.png"))
    written.append("scene.png")
    stem = os.path.splitext(os.path.basename(image))[0]
    openp = os.path.join(SCENE, stem + "_open.png")
    if os.path.exists(openp):
        shutil.copyfile(openp, os.path.join(dest, "scene_open.png"))
        written.append("scene_open.png")
    # wrap for this image, re-keyed to the stable name
    try:
        wrap = json.load(open(os.path.join(SCENE, "wrap.json")))
        wi = wrap.get(image) or (wrap if wrap.get("image") == image else None)
        if wi:
            keep = {k: wi[k] for k in ("haov", "vaov", "hfov", "vOffset", "pitch") if k in wi}
            with open(os.path.join(dest, "wrap.json"), "w") as f:
                json.dump({"scene.png": keep}, f, indent=2)
            written.append("wrap.json")
    except Exception:
        pass
    # hotspots, only if they were authored for this image (re-key image -> scene.png)
    try:
        hs = json.load(open(os.path.join(SCENE, "hotspots.json")))
        if hs.get("image") == image:
            hs["image"] = "scene.png"
            with open(os.path.join(dest, "hotspots.json"), "w") as f:
                json.dump(hs, f, indent=2)
            written.append("hotspots.json")
    except Exception:
        pass
    return written, room_dir


class H(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, **k):
        super().__init__(*a, directory=ROOT, **k)

    def log_message(self, *a):
        pass

    def _json(self, obj, code=200):
        b = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(b)))
        self.end_headers()
        self.wfile.write(b)

    def _body(self):
        n = int(self.headers.get("Content-Length", 0))
        return json.loads(self.rfile.read(n) or b"{}")

    def do_GET(self):
        route = self.path.split("?")[0]
        if route == "/api/status":
            qs = urllib.parse.parse_qs(self.path.split("?", 1)[1] if "?" in self.path else "")
            with LOCK:
                if "slot" in qs:
                    return self._json(dict(JOBS.get(qs["slot"][0], _IDLE)))
                return self._json({s: dict(j) for s, j in JOBS.items()})
        if route == "/api/scenes":
            files = sorted(os.path.basename(p) for p in glob.glob(os.path.join(SCENE, "*.png")))
            return self._json({"scenes": files})
        return super().do_GET()

    def do_POST(self):
        route = self.path.split("?")[0]
        try:
            if route == "/api/generate":
                req = self._body()
                prompt = (req.get("prompt") or "").strip()
                if not prompt:
                    return self._json({"ok": False, "error": "empty prompt"}, 400)
                n = max(1, min(int(req.get("n", 1)), 6))
                q = req.get("quality", "medium")
                size = req.get("size", "1536x1024")
                slot = str(req.get("slot", "0"))
                tag = _sanitize_tag(req.get("tag"))
                if not _start(slot, "generate",
                              lambda: _run_generate(slot, tag, prompt, n, q, size),
                              n, tag=tag):
                    return self._json({"ok": False, "error": "this column is already running"}, 409)
                return self._json({"ok": True, "total": n, "slot": slot, "tag": tag})
            if route == "/api/dooropen":
                req = self._body()
                img = req.get("image")
                box = req.get("box")
                prompt = (req.get("prompt") or "").strip()
                if not img or not isinstance(box, list) or len(box) != 4:
                    return self._json({"ok": False, "error": "need image + box[4]"}, 400)
                if not prompt:
                    return self._json({"ok": False, "error": "empty door prompt"}, 400)
                if not _start("door", "dooropen",
                              lambda: _run_dooropen("door", img, box, prompt), 1):
                    return self._json({"ok": False, "error": "a door-open job is already running"}, 409)
                return self._json({"ok": True})
            if route == "/api/commit-room":
                req = self._body()
                img = req.get("image")
                room_dir = req.get("roomDir")
                if not img or not room_dir:
                    return self._json({"ok": False, "error": "need image + roomDir"}, 400)
                try:
                    written, rd = _commit_room(img, room_dir)
                except ValueError as ve:
                    return self._json({"ok": False, "error": str(ve)}, 400)
                return self._json({"ok": True, "dest": rd, "written": written})
            if route == "/api/save-hotspots":
                req = self._body()
                if not req.get("image"):
                    return self._json({"ok": False, "error": "no image"}, 400)
                keep = {k: req[k] for k in ("image", "haov", "vaov", "vOffset",
                                            "hotspots") if k in req}
                keep.setdefault("hotspots", [])
                with open(os.path.join(SCENE, "hotspots.json"), "w") as f:
                    json.dump(keep, f, indent=2)
                return self._json({"ok": True, "count": len(keep["hotspots"])})
            if route == "/api/suggest":
                req = self._body()
                prompts = suggest_prompts(req.get("prompt", ""),
                                          max(1, min(int(req.get("n", 3)), 8)))
                return self._json({"ok": True, "prompts": prompts})
            if route == "/api/save-wrap":
                req = self._body()
                img = req.get("image")
                if not img:
                    return self._json({"ok": False, "error": "no image"}, 400)
                path = os.path.join(SCENE, "wrap.json")
                data = {}
                if os.path.exists(path):
                    try:
                        data = json.load(open(path))
                    except Exception:
                        data = {}
                # migrate a legacy flat {image, haov, ...} record into the per-image map
                if isinstance(data, dict) and "image" in data and "haov" in data:
                    old = data.pop("image")
                    data = {old: {k: v for k, v in data.items()}}
                if not isinstance(data, dict):
                    data = {}
                data[img] = {k: req[k] for k in ("haov", "vaov", "hfov",
                                                 "vOffset", "pitch") if k in req}
                with open(path, "w") as f:
                    json.dump(data, f, indent=2)
                return self._json({"ok": True, "image": img})
        except Exception as e:  # noqa: BLE001 — report to the UI
            return self._json({"ok": False, "error": str(e)}, 500)
        return self._json({"error": "not found"}, 404)


if __name__ == "__main__":
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", PORT), H)
    print(f"gpt harness -> http://127.0.0.1:{PORT}/harness_gpt.html", flush=True)
    httpd.serve_forever()
