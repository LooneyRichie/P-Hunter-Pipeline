"""
p_hunter_web.py — local web UI for p_hunter_deepen.

Runs a tiny HTTP server with no external dependencies, drives
p_hunter_deepen as a library, and streams the scan log to the browser.

    python p_hunter_web.py            # opens http://localhost:8000
    python p_hunter_web.py --port 9000 --no-browser
"""
import argparse
import io
import json
import os
import sys
import threading
import time
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse

import p_hunter_deepen as pipeline

# ============================================================
# Shared state
# ============================================================

_LOCK = threading.Lock()
_CANCEL = threading.Event()
_STATE = {
    "status":     "idle",       # idle | running | done | error | cancelled
    "config":     None,
    "started":    None,
    "finished":   None,
    "error":      None,
    "survivors":  [],
    "known":      [],
    "history":    [],
    "close_out":  None,
    "log":        "",
}


class _LogBuffer(io.TextIOBase):
    """Thread-safe rolling log.  Keeps the last N lines."""
    def __init__(self, max_lines=4000):
        self._lines = []
        self._max = max_lines
        self._lock = threading.Lock()

    def write(self, s):
        with self._lock:
            self._lines.extend(s.splitlines(keepends=True))
            if len(self._lines) > self._max:
                self._lines = self._lines[-self._max:]
        return len(s)

    def flush(self):
        pass

    def getvalue(self):
        with self._lock:
            return "".join(self._lines)


# ============================================================
# Scan runner (executes in a worker thread)
# ============================================================

def _run_scan(config):
    log = _LogBuffer()
    old_stdout = sys.stdout
    sys.stdout = log
    try:
        with _LOCK:
            _STATE.update({
                "status":    "running",
                "config":    config,
                "started":   time.time(),
                "finished":  None,
                "error":     None,
                "survivors": [],
                "known":     [],
                "history":   [],
                "close_out": None,
            })

        survivors = pipeline.initial_filter(config["start"], config["end"])
        survivors, history = pipeline.deepen(
            survivors,
            start_k=config["k"],
            multiplier=config["mult"],
            interactive=False,
            max_rounds=config.get("rounds"),
            workers=config.get("workers"),
        )
        known = sorted(p for p in survivors
                       if p in pipeline.KNOWN_MERSENNE_EXPONENTS)
        rest = sorted(p for p in survivors
                      if p not in pipeline.KNOWN_MERSENNE_EXPONENTS)
        with _LOCK:
            _STATE["known"] = known
            _STATE["survivors"] = rest
            _STATE["history"] = history
            _STATE["status"] = "done"
    except Exception as e:
        import traceback
        traceback.print_exc()
        with _LOCK:
            _STATE["error"] = f"{type(e).__name__}: {e}"
            _STATE["status"] = "error"
    finally:
        sys.stdout = old_stdout
        with _LOCK:
            _STATE["finished"] = time.time()
            _STATE["log"] = log.getvalue()


def _run_close_out():
    try:
        with _LOCK:
            survivors = list(_STATE["survivors"])
        import time as _t
        # Inline close-out so we don't touch pipeline.globals
        try:
            import gmpy2
            test = lambda n: gmpy2.is_prime(n) == 2
            backend = "gmpy2"
        except ImportError:
            try:
                import numthy as _nt
                test = _nt.is_prime
                backend = "numthy"
            except ImportError:
                from sympy import isprime as _sp
                test = _sp
                backend = "sympy"
        out = []
        for p in survivors:
            t0 = _t.time()
            try:
                v = test((1 << p) - 1)
            except Exception as e:
                v = f"error: {e}"
            out.append({
                "p": p,
                "verdict": ("prime" if v is True else
                            "composite" if v is False else str(v)),
                "ms": (_t.time() - t0) * 1000,
            })
        with _LOCK:
            _STATE["close_out"] = {"backend": backend, "results": out}
    except Exception as e:
        with _LOCK:
            _STATE["close_out"] = {"error": f"{type(e).__name__}: {e}"}


# ============================================================
# HTTP handler
# ============================================================

class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a, **k):
        pass  # quiet

    def _json(self, obj, status=200):
        body = json.dumps(obj).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _html(self, html):
        body = html.encode()
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        path = urlparse(self.path).path
        if path in ("/", "/index.html"):
            self._html(HTML)
        elif path == "/status":
            with _LOCK:
                snap = dict(_STATE)
            if snap["status"] == "running" and snap["started"]:
                snap["elapsed"] = time.time() - snap["started"]
            else:
                snap["elapsed"] = (
                    (snap["finished"] - snap["started"])
                    if snap["started"] and snap["finished"] else None
                )
            self._json(snap)
        else:
            self.send_error(404)

    def do_POST(self):
        path = urlparse(self.path).path
        length = int(self.headers.get("Content-Length", 0))
        raw = self.rfile.read(length) if length else b"{}"
        try:
            body = json.loads(raw)
        except Exception as e:
            return self._json({"error": f"bad json: {e}"}, 400)

        if path == "/start":
            with _LOCK:
                if _STATE["status"] == "running":
                    return self._json(
                        {"error": "a scan is already running"}, 409)
            threading.Thread(target=_run_scan, args=(body,),
                             daemon=True).start()
            self._json({"ok": True})
        elif path == "/close-out":
            with _LOCK:
                if _STATE["status"] == "running":
                    return self._json(
                        {"error": "wait for scan to finish"}, 409)
                if not _STATE["survivors"]:
                    return self._json(
                        {"error": "no survivors to test"}, 409)
            threading.Thread(target=_run_close_out, daemon=True).start()
            self._json({"ok": True})
        else:
            self.send_error(404)


# ============================================================
# HTML/CSS/JS — single page
# ============================================================

HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>p_hunter</title>
<style>
  :root {
    --bg:#0f1115; --panel:#171a21; --border:#2a2f3a;
    --fg:#d8dee9; --dim:#8a94a6; --accent:#7aa2f7;
    --ok:#9ece6a; --warn:#e0af68; --err:#f7768e;
  }
  * { box-sizing: border-box; }
  body {
    margin:0; padding:16px; background:var(--bg); color:var(--fg);
    font:14px/1.45 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
    max-width:760px; margin:0 auto;
  }
  h1 { font-size:18px; margin:0 0 12px; }
  h1 span { color: var(--dim); font-weight:400; }
  .card {
    background:var(--panel); border:1px solid var(--border);
    border-radius:8px; padding:14px; margin-bottom:12px;
  }
  label { display:block; color:var(--dim); font-size:12px; margin-bottom:4px; }
  input[type=number] {
    width:100%; padding:10px 12px; font:inherit;
    background:#0c0e12; color:var(--fg);
    border:1px solid var(--border); border-radius:6px;
  }
  .grid { display:grid; grid-template-columns:1fr 1fr; gap:10px; }
  .grid3 { display:grid; grid-template-columns:1fr 1fr 1fr; gap:10px; }
  button {
    padding:12px 16px; font:inherit; font-weight:600;
    background:var(--accent); color:#0b0d12; border:none;
    border-radius:6px; cursor:pointer; width:100%;
  }
  button:disabled { opacity:.4; cursor:not-allowed; }
  button.secondary { background:#2a2f3a; color:var(--fg); }
  button + button { margin-top:8px; }
  .row { display:flex; gap:10px; }
  .row > * { flex:1; }
  .status { display:flex; align-items:center; gap:10px; }
  .dot {
    width:10px; height:10px; border-radius:50%;
    background:var(--dim); flex-shrink:0;
  }
  .dot.running { background:var(--warn); animation:pulse 1s infinite; }
  .dot.done { background:var(--ok); }
  .dot.error { background:var(--err); }
  .dot.cancelled { background:var(--dim); }
  @keyframes pulse { 0%,100%{opacity:1} 50%{opacity:.3} }
  .meta { color:var(--dim); font-size:12px; }
  pre.log {
    background:#0c0e12; border:1px solid var(--border); border-radius:6px;
    padding:10px; height:280px; overflow:auto;
    white-space:pre; font-size:12px; margin:0;
  }
  .survivors { display:flex; flex-wrap:wrap; gap:6px; }
  .chip {
    padding:4px 10px; background:#0c0e12; border:1px solid var(--border);
    border-radius:99px; font-size:12px;
  }
  .chip.prime { border-color:var(--ok); color:var(--ok); }
  .chip.composite { border-color:var(--err); color:var(--err); }
  .hidden { display:none; }
  .tag-ok { color:var(--ok); }
  .tag-err { color:var(--err); }
</style>
</head>
<body>

<h1>p_hunter <span>· deepening sieve</span></h1>

<!-- ---------- config form ---------- -->
<div class="card" id="form-card">
  <div class="grid">
    <div>
      <label>start p</label>
      <input type="number" id="start" value="10000">
    </div>
    <div>
      <label>end p</label>
      <input type="number" id="end" value="10050">
    </div>
  </div>
  <div class="grid3" style="margin-top:10px;">
    <div>
      <label>initial k</label>
      <input type="number" id="k" value="1000">
    </div>
    <div>
      <label>multiplier</label>
      <input type="number" id="mult" value="10">
    </div>
    <div>
      <label>workers</label>
      <input type="number" id="workers" placeholder="auto">
    </div>
  </div>
  <div style="margin-top:12px;">
    <button id="start-btn">Start scan</button>
  </div>
</div>

<!-- ---------- status card ---------- -->
<div class="card" id="status-card">
  <div class="status">
    <span class="dot" id="dot"></span>
    <span id="status-text">idle</span>
  </div>
  <div class="meta" id="meta" style="margin-top:6px;"></div>
</div>

<!-- ---------- log card ---------- -->
<div class="card hidden" id="log-card">
  <pre class="log" id="log"></pre>
</div>

<!-- ---------- results card ---------- -->
<div class="card hidden" id="results-card">
  <div style="margin-bottom:10px;">
    <strong>Survivors</strong>
    <span class="meta" id="survivor-count"></span>
  </div>
  <div class="survivors" id="survivors"></div>
  <div style="margin-top:14px;">
    <button id="closeout-btn" class="secondary">Run close-out primality test</button>
    <button id="download-btn" class="secondary">Download results JSON</button>
    <button id="reset-btn" class="secondary">New scan</button>
  </div>
</div>

<div class="card hidden" id="closeout-card">
  <div style="margin-bottom:10px;">
    <strong>Close-out results</strong>
    <span class="meta" id="closeout-backend"></span>
  </div>
  <div class="survivors" id="closeout-list"></div>
</div>

<script>
const $ = id => document.getElementById(id);
let pollHandle = null;

function fmtElapsed(s) {
  if (s == null) return "";
  if (s < 60) return s.toFixed(1) + "s";
  const m = Math.floor(s / 60), r = Math.floor(s % 60);
  return m + "m " + r + "s";
}

async function post(path, body) {
  const r = await fetch(path, {
    method: "POST",
    headers: {"Content-Type": "application/json"},
    body: JSON.stringify(body || {}),
  });
  const j = await r.json().catch(() => ({}));
  if (!r.ok) throw new Error(j.error || ("HTTP " + r.status));
  return j;
}

async function startScan() {
  const cfg = {
    start:   parseInt($("start").value, 10) || 0,
    end:     parseInt($("end").value, 10) || 0,
    k:       parseInt($("k").value, 10) || 1000,
    mult:    parseInt($("mult").value, 10) || 10,
    workers: $("workers").value ? parseInt($("workers").value, 10) : null,
  };
  $("start-btn").disabled = true;
  $("log-card").classList.remove("hidden");
  $("results-card").classList.add("hidden");
  $("closeout-card").classList.add("hidden");
  $("log").textContent = "";
  try {
    await post("/start", cfg);
    startPolling();
  } catch (e) {
    alert("Error: " + e.message);
    $("start-btn").disabled = false;
  }
}

function renderStatus(s) {
  const dot = $("dot"), txt = $("status-text"), meta = $("meta");
  dot.className = "dot " + s.status;
  txt.textContent = s.status;
  const bits = [];
  if (s.config) {
    bits.push(`p ∈ [${s.config.start}, ${s.config.end}]`);
    bits.push(`k0=${s.config.k} ×${s.config.mult}`);
  }
  if (s.elapsed != null) bits.push("elapsed " + fmtElapsed(s.elapsed));
  if (s.status === "error") bits.push(s.error);
  meta.textContent = bits.join("  ·  ");
}

function renderLog(text) {
  const el = $("log");
  if (el.textContent !== text) {
    // Auto-scroll if user is near bottom (don't fight manual scroll)
    const near = (el.scrollHeight - el.scrollTop - el.clientHeight) < 60;
    el.textContent = text;
    if (near) el.scrollTop = el.scrollHeight;
  }
}

function renderSurvivors(list, known) {
  const el = $("survivors");
  el.innerHTML = "";
  for (const p of known) {
    const c = document.createElement("span");
    c.className = "chip prime";
    c.textContent = p + " ✓";
    el.appendChild(c);
  }
  for (const p of list) {
    const c = document.createElement("span");
    c.className = "chip";
    c.textContent = p;
    el.appendChild(c);
  }
  $("survivor-count").textContent =
      ` (${list.length} unverified, ${known.length} known)`;
}

function renderCloseOut(co) {
  if (!co) return;
  $("closeout-card").classList.remove("hidden");
  $("closeout-backend").textContent = co.backend ? " — " + co.backend : "";
  const el = $("closeout-list");
  el.innerHTML = "";
  if (co.error) {
    el.textContent = "error: " + co.error;
    return;
  }
  for (const r of co.results) {
    const c = document.createElement("span");
    c.className = "chip " + (r.verdict === "prime" ? "prime" :
                             r.verdict === "composite" ? "composite" : "");
    c.textContent = `${r.p} ${r.verdict} (${r.ms.toFixed(1)}ms)`;
    el.appendChild(c);
  }
}

async function poll() {
  try {
    const r = await fetch("/status");
    const s = await r.json();
    renderStatus(s);
    renderLog(s.log || "");
    if (s.status === "done" || s.status === "error" ||
        s.status === "cancelled") {
      $("start-btn").disabled = false;
      if (s.status === "done") {
        $("results-card").classList.remove("hidden");
        renderSurvivors(s.survivors, s.known);
      }
      if (s.close_out) renderCloseOut(s.close_out);
      if (s.status !== "done") stopPolling();
      if (s.status === "done" && !s.close_out) {
        // keep polling for close-out result
      } else if (s.close_out) {
        stopPolling();
      }
    }
  } catch (e) { /* server gone; ignore */ }
}

function startPolling() {
  if (pollHandle) return;
  pollHandle = setInterval(poll, 800);
  poll();
}
function stopPolling() {
  if (pollHandle) { clearInterval(pollHandle); pollHandle = null; }
}

async function runCloseOut() {
  $("closeout-btn").disabled = true;
  try {
    await post("/close-out");
    startPolling();
  } catch (e) {
    alert("Error: " + e.message);
  } finally {
    $("closeout-btn").disabled = false;
  }
}

async function downloadJSON() {
  const r = await fetch("/status");
  const s = await r.json();
  const blob = new Blob([JSON.stringify(s, null, 2)],
                        {type: "application/json"});
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `p_hunter_${s.config?.start}_${s.config?.end}.json`;
  a.click();
}

function reset() {
  stopPolling();
  $("form-card").classList.remove("hidden");
  $("results-card").classList.add("hidden");
  $("closeout-card").classList.add("hidden");
  $("log").textContent = "";
  $("survivors").innerHTML = "";
  $("closeout-list").innerHTML = "";
}

$("start-btn").addEventListener("click", startScan);
$("closeout-btn").addEventListener("click", runCloseOut);
$("download-btn").addEventListener("click", downloadJSON);
$("reset-btn").addEventListener("click", reset);

// On page load, adopt any run already in progress
startPolling();
</script>
</body>
</html>
"""


# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--no-browser", action="store_true")
    args = ap.parse_args()

    url = f"http://{args.host}:{args.port}"
    print(f"p_hunter web UI: {url}")
    if args.host == "127.0.0.1":
        print("  (localhost only — safe to leave running)")
    print("Ctrl+C to stop.")

    if not args.no_browser:
        try:
            webbrowser.open(url)
        except Exception:
            pass

    server = ThreadingHTTPServer((args.host, args.port), Handler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nshutting down.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()