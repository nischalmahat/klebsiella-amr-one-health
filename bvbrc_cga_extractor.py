
#!/usr/bin/env python3
"""
BV-BRC Comprehensive Genome Analysis (CGA) extractor — balanced 210 K. pneumoniae
==================================================================================
Sample set : Kp_sampling_design_v2_1.xlsx, sheet "Balanced_210"
             3 Isolation_Groups x 5 continents x 14 genomes = 210 (all QC KEEP).
             Sample_ID, Source_Group, Continent, Collection_year and the QC
             fingerprint (genome size, contigs, GC) all come from that workbook.
Account / workspace : <your-username>@bvbrc  ->  /home/KP_CGA/CGA Results/
                      ├── Animal CGA/         (KP001–KP100)
                      ├── Clinical CGA/       (KP101–KP200)
                      └── Environmental CGA/  (KP201–KP300)
Each CGA job is named by Sample ID (KP001, KP002, ...).

PIPELINE
  [0] Authenticate and check that the token belongs to <your-username>@bvbrc
      (you have two BV-BRC accounts, so a token for the wrong one is caught here).
  [1] Resolve Sample ID -> BV-BRC Genome ID with two independent routes:
        Route A (workspace): CGA job folder .KPxxx -> annotated.genome (GTO) -> "id"
                             (the method validated in your 90-genome script)
        Route B (data API):  /api/genome owned by your account, Sample ID parsed
                             from genome_name
      The two routes are cross-checked against each other.
  [2] Check the design (unique IDs, all KEEP/INCLUDE_BALANCED, 14 per cell,
      Balanced_210 consistent with Sampling_v2_all300), then keep only those 210.
      Every other workspace job is listed with its reason (QC fail, duplicate,
      removed for balancing).
  [3] Download per genome: AMR specialty genes, virulence-factor specialty genes,
      AMR phenotype predictions (genome_amr) and subsystems. Results are cached,
      so an interrupted run resumes where it stopped.
  [4] Validate: genome length / contig count / GC from BV-BRC against your own QC
      table (a fingerprint that catches swapped or mislabelled genomes), row
      counts against the server's reported totals, one genome ID per sample,
      known anchor (KP001 = 573.104016), live re-query of random samples.
  [5] Write ONE Excel workbook with all three isolation groups. Every row carries
      Isolation_Group, Continent, Collection_year and Stratum (Group|Continent);
      presence/absence matrices; prevalence per group and per group x continent;
      a Design_Balance sheet proving all 15 cells still hold 14 genomes.

RUN
  pip install requests pandas openpyxl
  # put Kp_sampling_design_v2_1.xlsx in the same folder, then
  python bvbrc_cga_extractor.py --only KP001 KP002     # quick test first
  python bvbrc_cga_extractor.py                        # all 210

  Credentials: set BVBRC_USERNAME / BVBRC_PASSWORD (or BVBRC_TOKEN) as environment
  variables, or the script asks for them. Never written to disk or the outputs.

Exit code 0 = all checks PASS/WARN; 2 = at least one FAIL (outputs are still
written, but the README sheet is marked NOT READY).

v2.1.6: the one-off "swap KP256 -> KP216 in the sampling design" cell is now part of
this script (apply_design_swap, run at the [design] step before the design is read).
It is idempotent: if KP216 is already in Balanced_210 and KP256 is not, it changes
nothing. Everything else is identical to v2.1.5.
"""

from __future__ import annotations
import argparse
import getpass
import gzip
import json
import logging
import os
import random
import re
import shutil
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

__version__ = "2.1.6"

# ════════════════════════════════════════════════════════════════════════
# CREDENTIALS  — never type them into this file
#   Leave both empty. get_token() then uses, in this order:
#     1) the BVBRC_TOKEN environment variable
#     2) a saved CLI token (~/.patric_token or ~/.bvbrc_token) that belongs to WS_OWNER
#     3) the BVBRC_USERNAME / BVBRC_PASSWORD environment variables
#     4) a prompt (the password is typed hidden with getpass, never echoed or saved)
# ════════════════════════════════════════════════════════════════════════
BVBRC_USERNAME = ""
BVBRC_PASSWORD = ""

# ════════════════════════════════════════════════════════════════════════
# RUN SETTINGS  (edit here when running in Jupyter / Spyder / IDLE)
# ════════════════════════════════════════════════════════════════════════
DESIGN_FILE   = "Kp_sampling_design_v2_1.xlsx"   # path to the sampling-design workbook
DESIGN_SHEET  = "Balanced_210"                   # samples to extract
ALL_SHEET     = "Sampling_v2_all300"             # all 300 genomes (fingerprints, reasons)
ONLY_SAMPLES  = []          # e.g. ["KP001", "KP002", "KP003"] for a quick test; [] = all 210
OUTPUT_FILE   = "KP_CGA_Balanced210_metadata.xlsx"
CACHE_DIR     = "bvbrc_cache"
SKIP_GTO      = False       # True = genome-API ids only (faster, less strict)
REFRESH_CACHE = False       # True = re-download everything
LIVE_RECHECK  = 5           # random samples re-queried live for validation
APPLY_SWAP    = True        # True = make sure the design has SWAP_ADD in and SWAP_REMOVE out (safe to leave on)
SWAP_REMOVE   = "KP256"     # assembly-quality outlier, taken out of Balanced_210
SWAP_ADD      = "KP216"     # replacement (same Source_Group x Continent cell, QC KEEP)

# ════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ════════════════════════════════════════════════════════════════════════
WS_OWNER = "<your-username>@bvbrc"
WS_ROOT = f"/{WS_OWNER}/home/KP_CGA/CGA Results"
GROUP_FOLDERS = {                       # Isolation_Group -> workspace folder
    "Animal": "Animal CGA",
    "Clinical": "Clinical CGA",
    "Environmental": "Environmental CGA",
}
GROUP_ORDER = ["Animal", "Clinical", "Environmental"]
# Sample IDs whose Genome ID you have confirmed by hand on the website.
KNOWN_ANCHORS = {"KP001": "573.104016"}

KP_TAXON_ID = "573"                     # Klebsiella pneumoniae NCBI taxon

# Endpoints (env overrides exist only for offline testing)
AUTH_URLS = [u for u in [os.environ.get("BVBRC_AUTH_URL"),
                         "https://user.patricbrc.org/authenticate",
                         "https://user.bv-brc.org/authenticate"] if u]
WS_URLS = [u for u in [os.environ.get("BVBRC_WS_URL"),
                       "https://www.bv-brc.org/services/Workspace",
                       "https://p3.theseed.org/services/Workspace"] if u]
API_BASE = os.environ.get("BVBRC_API_BASE", "https://www.bv-brc.org/api")
if os.environ.get("BVBRC_WS_URL"):      # testing: do not fall back to real servers
    WS_URLS, AUTH_URLS = WS_URLS[:1], AUTH_URLS[:1]

PAGE_SIZE = 10000                       # rows per API page
SLEEP = float(os.environ.get("BVBRC_SLEEP", 0.15))   # politeness delay between calls (s)

# "KP001", ".KP001", "KP_001", "KP001_CGA" -> KP001 ; not "KP0011"
SAMPLE_RE = re.compile(r"(?<![A-Za-z0-9])KP[_\- ]?(\d{3})(?!\d)", re.I)
GENOME_ID_RE = re.compile(r"^\d+\.\d+$")

GENOME_FIELDS = ["genome_id", "genome_name", "owner", "public", "taxon_id",
                 "species", "genome_status", "genome_length", "contigs",
                 "gc_content", "contig_n50", "cds", "mlst", "genome_quality",
                 "checkm_completeness", "checkm_contamination",
                 "coarse_consistency", "fine_consistency", "date_inserted"]
SP_FIELDS = ["genome_id", "patric_id", "refseq_locus_tag", "gene", "product",
             "property", "source", "source_id", "evidence", "classification",
             "antibiotics_class", "antibiotics", "pmid", "query_coverage",
             "subject_coverage", "identity", "e_value"]
AMR_FIELDS = ["genome_id", "antibiotic", "resistant_phenotype", "evidence",
              "computational_method", "computational_method_performance",
              "measurement", "measurement_unit", "measurement_sign",
              "laboratory_typing_method", "laboratory_typing_platform",
              "testing_standard", "source"]
SUB_FIELDS = ["genome_id", "patric_id", "gene", "product", "role_name",
              "subsystem_id", "subsystem_name", "superclass", "class", "subclass"]

PASS, WARN, FAIL = "PASS", "WARN", "FAIL"
log = logging.getLogger("cga")


class ExtractionError(RuntimeError):
    pass


# ════════════════════════════════════════════════════════════════════════
# HTTP CLIENT (auth, workspace JSON-RPC, data API with verified paging)
# ════════════════════════════════════════════════════════════════════════
def make_session() -> requests.Session:
    s = requests.Session()
    retry = Retry(total=6, connect=6, read=6, backoff_factor=1.5,
                  status_forcelist=(429, 502, 503, 504),   # not 500: workspace uses 500 for "object not found"
                  allowed_methods=frozenset({"GET", "POST"}),
                  raise_on_status=False)
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.mount("http://", HTTPAdapter(max_retries=retry))
    return s
def token_user(token: str) -> str:
    m = re.search(r"(?:^|\|)un=([^|]+)", token or "")
    return m.group(1) if m else ""
class BVBRC:
    def __init__(self, token: str):
        self.token = token
        self.s = make_session()
        self.s.headers.update({"Authorization": token})
        self.ws_url = None
        self._rpc_id = 0
        self.n_calls = 0

    # ── authentication ──────────────────────────────────────────────────
    @staticmethod
    def login(username: str, password: str) -> str:
        s = make_session()
        last = ""
        for url in AUTH_URLS:
            try:
                r = s.post(url, data={"username": username, "password": password},
                           headers={"Content-Type": "application/x-www-form-urlencoded"},
                           timeout=30)
            except requests.RequestException as e:
                last = f"{url}: {e}"
                continue
            if r.status_code == 200 and "un=" in r.text:
                return r.text.strip()
            last = f"{url}: HTTP {r.status_code} {r.text[:150]}"
        raise ExtractionError(f"Authentication failed ({last})")

    # ── workspace ───────────────────────────────────────────────────────
    def rpc(self, method: str, params: list):
        urls = [self.ws_url] if self.ws_url else WS_URLS
        last = None
        for url in urls:
            self._rpc_id += 1
            try:
                r = self.s.post(url, json={"version": "1.1", "id": self._rpc_id,
                                           "method": f"Workspace.{method}",
                                           "params": params},
                                timeout=180)
                self.n_calls += 1
            except requests.RequestException as e:
                last = e
                continue
            if r.status_code != 200:
                last = f"HTTP {r.status_code}: {r.text[:300]}"
                # a JSON-RPC error from a working server is not a reason to switch URL
                try:
                    err = r.json().get("error")
                except ValueError:
                    err = None
                if err:
                    raise ExtractionError(f"Workspace.{method}: {err.get('message', err)}")
                continue
            data = r.json()
            if data.get("error"):
                raise ExtractionError(f"Workspace.{method}: {data['error'].get('message', data['error'])}")
            self.ws_url = url
            time.sleep(SLEEP)
            return data.get("result", [])
        raise ExtractionError(f"Workspace service unreachable ({last})")

    def ws_ls(self, path: str, recursive: bool = False) -> list[dict]:
        res = self.rpc("ls", [{"paths": [path], "excludeDirectories": 0,
                               "excludeObjects": 0, "recursive": 1 if recursive else 0}])
        raw = res[0].get(path, []) if res else []
        out = []
        for e in raw:
            # entry tuple: name, type, parent path, created, id, owner, size,
            #              user_meta, auto_meta, user_perm, global_perm, shock_url
            if not isinstance(e, list) or len(e) < 3:
                continue
            # e[2] is the PARENT folder (lesson from the 90-genome script),
            # so the full path is built as parent + "/" + name
            parent = str(e[2] or path).rstrip("/")
            out.append({"name": e[0], "type": e[1],
                        "path": f"{parent}/{e[0]}",
                        "created": e[3] if len(e) > 3 else "",
                        "size": e[6] if len(e) > 6 else None})
        return out

    @staticmethod
    def _unwrap_get(res):
        """Workspace.get returns [[ [meta, data] ]] (JSON-RPC result list -> list of
        objects -> (meta, data) tuple). Walk down until the (meta, data) pair is found."""
        x = res
        for _ in range(5):
            if (isinstance(x, list) and len(x) == 2 and isinstance(x[0], list)
                    and not isinstance(x[1], list)):
                return x[0], x[1]
            if isinstance(x, list) and x:
                x = x[0]
            else:
                break
        return None, None

    last_diag = ""                      # what the server sent for the last failed object

    @staticmethod
    def _decode(raw) -> dict:
        """bytes/str/dict -> dict. Handles gzip, BOM, double-encoded JSON; if the text is
        not valid JSON (e.g. cut off), recovers the genome id with a regex."""
        if isinstance(raw, dict):
            return raw
        b = raw.encode("utf-8") if isinstance(raw, str) else bytes(raw)
        if b[:2] == b"\x1f\x8b":                       # gzip
            b = gzip.decompress(b)
        t = b.decode("utf-8", "replace").lstrip("\ufeff").strip()
        try:
            obj = json.loads(t)
            if isinstance(obj, str):                   # JSON text stored as a JSON string
                obj = json.loads(obj)
        except ValueError:
            m = (re.search(r'"id"\s*:\s*"(\d+\.\d+)"', t[:500000])
                 or re.search(r"fig\|(\d+\.\d+)\.(?:peg|rna|CDS)", t[:2000000]))
            if m:
                return {"id": m.group(1), "_recovered_by_regex": True}
            raise ValueError(f"not JSON ({len(b)} bytes, starts {t[:120]!r})")
        if isinstance(obj, dict) and "id" not in obj and {"status", "data"} <= set(obj):
            raise ValueError(f"Shock returned node info, not the file: {str(obj)[:200]}")
        return obj

    def ws_get_json(self, path: str) -> dict:
        res = self.rpc("get", [{"objects": [path]}])
        meta, content = self._unwrap_get(res)
        if meta is None:
            BVBRC.last_diag = f"PATH {path}\nRAW RESPONSE {str(res)[:2000]}"
            raise ExtractionError(f"unexpected workspace response: {str(res)[:200]}")
        shock = meta[11] if len(meta) > 11 else ""
        diag = [f"PATH {path}", f"META {str(meta)[:1000]}",
                f"INLINE CONTENT type={type(content).__name__} len={len(content) if content else 0} "
                f"start={str(content)[:300]!r}"]
        attempts = []
        if content not in (None, "", {}):
            attempts.append(("inline", content))
        if shock:
            for url in (f"{shock}?download", f"{shock}?download_raw"):
                try:
                    r = self.s.get(url, headers={"Authorization": f"OAuth {self.token}"}, timeout=600)
                except requests.RequestException as e:
                    diag.append(f"SHOCK {url} -> {e}")
                    continue
                diag.append(f"SHOCK {url} -> HTTP {r.status_code} type={r.headers.get('Content-Type')} "
                            f"bytes={len(r.content)} final_url={r.url} start={r.content[:300]!r}")
                if r.status_code == 200 and r.content:
                    attempts.append((url, r.content))
                    break
        errors = []
        for where, raw in attempts:
            try:
                return self._decode(raw)
            except (ValueError, OSError) as e:
                errors.append(f"{where}: {e}")
        BVBRC.last_diag = "\n".join(diag)
        raise ExtractionError("; ".join(errors) or f"no content and no download link (shock='{shock}')")

    # ── data API ────────────────────────────────────────────────────────
    API_RETRY_WAITS = tuple(int(x) for x in os.environ.get("BVBRC_RETRY_WAITS", "5,15,30,60,120").split(","))

    def _get_with_retry(self, url: str, collection: str) -> requests.Response:
        """Data-API GET that survives BV-BRC's temporary database errors
        (HTTP 500 'socket hang up', timeouts, dropped connections)."""
        last = ""
        for attempt, wait in enumerate((0,) + self.API_RETRY_WAITS):
            if wait:
                log.warning(f"    BV-BRC busy ({last[:90]}) — retry {attempt}/{len(self.API_RETRY_WAITS)} in {wait}s")
                time.sleep(wait)
            try:
                r = self.s.get(url, headers={"Accept": "application/json"}, timeout=180)
                self.n_calls += 1
            except requests.RequestException as e:
                last = f"{type(e).__name__}: {e}"
                continue
            if r.status_code == 200:
                try:
                    r.json()
                    return r
                except ValueError:
                    last = f"HTTP 200 but not JSON: {r.text[:100]}"
                    continue
            last = f"HTTP {r.status_code}: {r.text[:150]}"
            if r.status_code in (400, 401, 403, 404):      # real errors: retrying will not help
                break
        raise ExtractionError(f"{collection} failed after retries — {last}. "
                              f"Re-run the same command: finished genomes are cached.")

    def query(self, collection: str, rql: str, fields: list[str]) -> list[dict]:
        """Fetch ALL rows of an RQL query; verify count against Content-Range."""
        rows, offset, total = [], 0, None
        select = ",".join(fields)
        while True:
            url = f"{API_BASE}/{collection}/?{rql}&select({select})&limit({PAGE_SIZE},{offset})"
            r = self._get_with_retry(url, collection)
            batch = r.json()
            if not isinstance(batch, list):
                raise ExtractionError(f"{collection}: unexpected response {str(batch)[:200]}")
            m = re.search(r"/(\d+)\s*$", r.headers.get("Content-Range", ""))
            if m:
                total = int(m.group(1))
            rows.extend(batch)
            offset += len(batch)
            time.sleep(SLEEP)
            if not batch or len(batch) < PAGE_SIZE or (total is not None and offset >= total):
                break
        if total is not None and len(rows) != total:
            raise ExtractionError(f"{collection} [{rql}]: got {len(rows)} rows, "
                                  f"server reports {total}")
        return rows


# ════════════════════════════════════════════════════════════════════════
# HELPERS
# ════════════════════════════════════════════════════════════════════════
def norm_sample(text: str) -> str | None:
    m = SAMPLE_RE.search(str(text))
    return f"KP{m.group(1)}" if m else None


def eqv(field: str, value: str) -> str:
    """RQL eq(). Values with spaces/commas/brackets MUST be double-quoted, otherwise
    Solr splits "Antibiotic Resistance" into two terms and matches nothing."""
    v = str(value)
    if re.search(r"[\s,()]", v):
        v = f'"{v}"'
    return f"eq({field},{quote(v, safe='')})"


AMR_PROPERTY, VF_PROPERTY = "Antibiotic Resistance", "Virulence Factor"
CACHE_SCHEMA = 3                        # bump -> old cached genome files are re-downloaded


def fetch_sp_genes(api: "BVBRC", gid: str, fields=None) -> tuple[list, list, str, dict]:
    """Specialty genes split into AMR and virulence.
    1st try: server-side filter on the quoted property value.
    Fallback (both empty): download all specialty genes of the genome and filter here
    (case-insensitive), so a naming/encoding difference can never silently give 0 rows."""
    fields = fields or SP_FIELDS
    g = eqv("genome_id", gid)
    amr = api.query("sp_gene", f"{g}&{eqv('property', AMR_PROPERTY)}", fields)
    vf = api.query("sp_gene", f"{g}&{eqv('property', VF_PROPERTY)}", fields)
    if amr or vf:
        return amr, vf, "server filter", {}
    allrows = api.query("sp_gene", g, list(dict.fromkeys(list(fields) + ["property"])))
    props = dict(sorted({p: sum(1 for x in allrows if x.get("property") == p)
                         for p in {str(x.get("property")) for x in allrows}}.items()))
    amr = [x for x in allrows if "resist" in str(x.get("property", "")).lower()]
    vf = [x for x in allrows if "virulen" in str(x.get("property", "")).lower()]
    return amr, vf, "local filter (all sp_gene rows)", props


def read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(p)
    if p.suffix.lower() in (".txt", ".list"):
        ids = [norm_sample(x) for x in p.read_text().split()]
        return pd.DataFrame({"Sample_ID": [i for i in ids if i]})
    return pd.read_csv(p)


def find_col(df: pd.DataFrame, *patterns: str) -> str | None:
    for pat in patterns:
        for c in df.columns:
            if re.search(pat, str(c), re.I):
                return c
    return None


class Checks:
    """Collects validation results -> 'Validation' sheet."""
    def __init__(self):
        self.rows = []

    def add(self, check, status, detail="", sample_id=""):
        self.rows.append({"Check": check, "Sample_ID": sample_id,
                          "Status": status, "Detail": detail})
        lvl = {PASS: logging.DEBUG, WARN: logging.WARNING, FAIL: logging.ERROR}[status]
        log.log(lvl, f"  [{status}] {check} {sample_id} {detail}")
    def df(self):
        return pd.DataFrame(self.rows, columns=["Check", "Sample_ID", "Status", "Detail"])

    def n(self, status):
        return sum(r["Status"] == status for r in self.rows)


# ════════════════════════════════════════════════════════════════════════
# STEP 1 — DISCOVER CGA JOBS AND RESOLVE GENOME IDs
# ════════════════════════════════════════════════════════════════════════
def discover_jobs(api: BVBRC, chk: Checks) -> pd.DataFrame:
    rows, unmatched = [], []
    for group in GROUP_ORDER:
        path = f"{WS_ROOT}/{GROUP_FOLDERS[group]}"
        entries = api.ws_ls(path)
        folders = [e for e in entries if e["type"] == "folder"]
        log.info(f"  {group:13s} {path}: {len(folders)} job folders, "
                 f"{sum(e['type'] == 'job_result' for e in entries)} job_result objects")
        for e in folders:
            sid = norm_sample(e["name"])
            if not sid:
                unmatched.append(e["path"])
                continue
            rows.append({"Sample_ID": sid, "Workspace_Group": group,
                         "Job_Folder": e["name"], "Job_Path": e["path"],
                         "Job_Created": e["created"]})
    for p in unmatched:
        chk.add("workspace_folder_named_by_sample_id", WARN, f"no KPxxx in name: {p}")
    df = pd.DataFrame(rows, columns=["Sample_ID", "Workspace_Group", "Job_Folder",
                                     "Job_Path", "Job_Created"])
    dup = df[df.duplicated("Sample_ID", keep=False)]
    for sid, g in dup.groupby("Sample_ID"):
        chk.add("one_cga_job_per_sample", FAIL,
                f"{len(g)} job folders: " + "; ".join(g["Job_Path"]), sid)
    return df


def gto_genome_id(api: BVBRC, job_path: str) -> tuple[str, str, dict]:
    """Route A: return (genome_id, gto_path, fingerprint) from the annotated GTO."""
    candidates = [f"{job_path}/annotated.genome"]
    try:
        gto = api.ws_get_json(candidates[0])
        gto_path = candidates[0]
    except ExtractionError:
        # fallback: look anywhere inside the job folder for a *.genome object
        found = [e["path"] for e in api.ws_ls(job_path, recursive=True)
                 if e["name"].endswith(".genome") or e["type"] == "genome"]
        if not found:
            raise ExtractionError(f"no annotated.genome under {job_path}")
        gto_path = sorted(found, key=lambda p: ("annotated" not in p, len(p)))[0]
        gto = api.ws_get_json(gto_path)
    gid = str(gto.get("id", ""))
    contigs = gto.get("contigs") or []
    fp = {"GTO_Name": gto.get("scientific_name", ""),
          "GTO_Contigs": len(contigs) if contigs else None,
          "GTO_Length": sum(len(c.get("dna", "")) for c in contigs) if contigs else None}
    return gid, gto_path, fp
def resolve_genome_ids(api: BVBRC, jobs: pd.DataFrame, cache: Path,
                       use_gto: bool, chk: Checks) -> pd.DataFrame:
    # ── Route B: every genome owned by the account ──────────────────────
    owned = api.query("genome", eqv("owner", WS_OWNER), GENOME_FIELDS)
    by_sample = defaultdict(list)
    for g in owned:
        sid = norm_sample(g.get("genome_name", ""))
        if sid:
            by_sample[sid].append(str(g["genome_id"]))
    log.info(f"  Route B: {len(owned)} genomes owned by {WS_OWNER}, "
             f"{len(by_sample)} with a KPxxx in genome_name")

    # ── Route A: GTO inside each CGA job (cached; GTOs are large) ───────
    map_file = cache / "gto_genome_ids.json"
    gto_cache = json.loads(map_file.read_text()) if map_file.exists() else {}
    debug_file = cache / "gto_debug.txt"
    fails_in_row, MAX_FAILS = 0, 5
    out = []
    for i, r in enumerate(jobs.itertuples(index=False), 1):
        sid = r.Sample_ID
        rec = {"Sample_ID": sid, "GTO_Genome_ID": "", "GTO_Path": "",
               "GTO_Name": "", "GTO_Contigs": None, "GTO_Length": None}
        if use_gto:
            key = r.Job_Path
            if not gto_cache.get(key):          # failures are retried on the next run
                try:
                    gid, p, fp = gto_genome_id(api, r.Job_Path)
                    gto_cache[key] = {"genome_id": gid, "gto_path": p, **fp}
                    map_file.write_text(json.dumps(gto_cache, indent=1))
                    fails_in_row = 0
                except (ExtractionError, requests.RequestException, ValueError) as e:
                    fails_in_row += 1
                    log.warning(f"  {sid}: GTO not readable ({str(e)[:250]})")
                    if fails_in_row == 1:
                        debug_file.write_text(BVBRC.last_diag or str(e), encoding="utf-8")
                        log.warning(f"  server reply saved to {debug_file} — send this file for diagnosis")
                    if fails_in_row >= MAX_FAILS:
                        use_gto = False
                        msg = (f"route A switched off after {MAX_FAILS} failures in a row; Genome IDs "
                               f"come from the genome API and are validated by the length/contig/GC "
                               f"fingerprint against your QC. Details: {debug_file}")
                        log.warning("  " + msg)
                        chk.add("route_a_annotated_genome", WARN, msg)
            c = gto_cache.get(key)
            if c:
                rec.update({"GTO_Genome_ID": c["genome_id"], "GTO_Path": c["gto_path"],
                            "GTO_Name": c.get("GTO_Name", ""),
                            "GTO_Contigs": c.get("GTO_Contigs"),
                            "GTO_Length": c.get("GTO_Length")})
            log.info(f"  [{i:3d}/{len(jobs)}] {sid} -> GTO id {rec['GTO_Genome_ID'] or '—'}")
        api_ids = sorted(set(by_sample.get(sid, [])))
        rec["API_Genome_IDs"] = ";".join(api_ids)
        gto_id = rec["GTO_Genome_ID"]

        if gto_id and not GENOME_ID_RE.match(gto_id):
            rec["Genome_ID"], rec["ID_Status"] = "", "BAD_GTO_ID"
        elif gto_id and gto_id in api_ids:
            rec["Genome_ID"] = gto_id
            rec["ID_Status"] = "CONFIRMED" if len(api_ids) == 1 else "CONFIRMED_MULTI_API"
        elif gto_id and not api_ids:
            rec["Genome_ID"], rec["ID_Status"] = gto_id, "GTO_ONLY"
        elif gto_id:
            rec["Genome_ID"], rec["ID_Status"] = gto_id, "CONFLICT"
        elif len(api_ids) == 1:
            rec["Genome_ID"], rec["ID_Status"] = api_ids[0], "API_ONLY"
        else:
            rec["Genome_ID"] = ""
            rec["ID_Status"] = "UNRESOLVED" if not api_ids else "AMBIGUOUS_API"
        out.append(rec)
    res = jobs.merge(pd.DataFrame(out), on="Sample_ID", how="left")
    return res, owned


# ════════════════════════════════════════════════════════════════════════
# STEP 2 — SAMPLING DESIGN (Kp_sampling_design_v2_1.xlsx -> Balanced_210)
# ════════════════════════════════════════════════════════════════════════
DESIGN_COLS = ["Sample_ID", "Assembly_Number", "Source_Group", "Continent", "Collection_year",
               "Final_Decision_v2", "Final_Decision", "Exclusion_Reason", "Duplicate_Of",
               "Completeness", "Contamination", "Contigs", "Genome_Size_bp", "GC_percent",
               "N50_own", "Kleborate_Species"]

SWAP_REASON = ("Replaced after BV-BRC CGA: assembly-quality outlier (likely uncorrected long-read "
               "assembly) - 7,487 CDS vs 4,866-6,220 in the other 209, BV-BRC fine consistency 80.5% "
               "(others >=94.8%), CheckM contamination 17.1%, CheckM2 completeness 91.2%. "
               "Replacement = Environmental x South America QC-KEEP genome with fewest contigs "
               "(rule fixed before looking at gene content).")
SWAP_DECISION_REMOVED = "EXCLUDE_ASSEMBLY_QUALITY"


def apply_design_swap(path: str, sheet: str, all_sheet: str, out_stem: str):
    """Swap SWAP_REMOVE -> SWAP_ADD in the sampling design workbook (was a separate one-off cell).

    Edits the workbook in place with openpyxl, so its formatting is kept. If SWAP_ADD is
    already in the balanced sheet and SWAP_REMOVE is not, nothing is changed and nothing is
    backed up. Otherwise the design and the current outputs are copied to a backup folder first.
    """
    from openpyxl import load_workbook

    design = Path(path)
    remove, add = SWAP_REMOVE, SWAP_ADD

    def need(cond, msg):
        if not cond:
            raise ExtractionError(f"design swap {remove} -> {add}: {msg}")

    wb = load_workbook(design)
    ALL, BAL, EXC = wb[all_sheet], wb[sheet], wb["Excluded_Balancing_63"]

    def header(ws):
        return [c.value for c in ws[1]]
    def find_row(ws, sid):
        for row in ws.iter_rows(min_row=2):
            if row[0].value == sid:
                return row
        return None

    H = header(ALL)
    need(header(BAL) == H and header(EXC) == H, "sheets do not share the same columns")
    col = {name: i for i, name in enumerate(H)}

    # ── safety checks (also makes this safe to run every time) ──────────────
    if find_row(BAL, add) and not find_row(BAL, remove):
        log.info(f"  Design swap: already done ({add} is in {sheet}, {remove} is not).")
    else:
        r_rem, r_add = find_row(ALL, remove), find_row(ALL, add)
        need(find_row(BAL, remove), f"{remove} not in {sheet}")
        need(find_row(EXC, add), f"{add} not in Excluded_Balancing_63")
        for k in ("Source_Group", "Continent"):
            need(r_rem[col[k]].value == r_add[col[k]].value, f"{k} differs - not the same design cell")
        need(r_add[col["Final_Decision"]].value == "KEEP", f"{add} did not pass QC")

        # backup of the design and of the current outputs (made before anything is changed)
        stamp = datetime.now().strftime("%Y%m%d_%H%M")
        backup = Path(f"backup_before_{add}_swap_{stamp}")
        backup.mkdir(exist_ok=True)
        shutil.copy2(design, backup / design.name)
        for f in Path(".").glob(f"{out_stem}*"):
            if f.is_file():
                shutil.copy2(f, backup / f.name)
        log.info(f"  Design swap {remove} -> {add}: backups -> {backup}/ ({len(list(backup.iterdir()))} files)")
        log.info(f"  Same cell: {r_add[col['Source_Group']].value} x {r_add[col['Continent']].value}; "
                 f"{add} QC: completeness {r_add[col['Completeness']].value}, "
                 f"contamination {r_add[col['Contamination']].value}, contigs {r_add[col['Contigs']].value}")

        # 1) master sheet: update the two decisions
        r_add[col["Final_Decision_v2"]].value = "INCLUDE_BALANCED"
        r_add[col["Exclusion_Reason"]].value = f"Replacement for {remove}"
        r_rem[col["Final_Decision_v2"]].value = SWAP_DECISION_REMOVED
        r_rem[col["Exclusion_Reason"]].value = SWAP_REASON

        # 2) balanced sheet: the SWAP_REMOVE row is overwritten with SWAP_ADD (same cell, same position)
        for cell, src in zip(find_row(BAL, remove), r_add):
            cell.value = src.value
        # 3) excluded sheet: the SWAP_ADD row is overwritten with SWAP_REMOVE
        for cell, src in zip(find_row(EXC, add), r_rem):
            cell.value = src.value

        # 4) README note
        rd = wb["README"]
        rd.append([])
        rd.append([f"Update {datetime.now():%Y-%m-%d}: {remove} replaced by {add} in {sheet}. {SWAP_REASON}"])
        rd.append([f"Excluded_Balancing_63 now holds 62 balancing exclusions + {remove} ({SWAP_DECISION_REMOVED})."])
        wb.save(design)

    # ── verify the saved workbook ───────────────────────────────────────────
    b = pd.read_excel(design, sheet)
    a = pd.read_excel(design, all_sheet).set_index("Sample_ID")
    need(len(b) == 210 and b.Sample_ID.is_unique, f"{sheet} is not 210 unique samples")
    need(add in set(b.Sample_ID) and remove not in set(b.Sample_ID), f"{add} in / {remove} out check failed")
    need((b.Final_Decision_v2 == "INCLUDE_BALANCED").all() and (b.Final_Decision == "KEEP").all(),
         f"not every {sheet} row is KEEP + INCLUDE_BALANCED")
    cells = pd.crosstab(b.Source_Group, b.Continent)
    need((cells.values == 14).all(), f"cells are not all 14:\n{cells}")
    need(a.at[remove, "Final_Decision_v2"] == SWAP_DECISION_REMOVED, f"{remove} decision not {SWAP_DECISION_REMOVED}")
    log.info(f"  ✓ Design checked: 210 genomes, all 15 cells = 14, {add} in, {remove} out")


def load_design(path: str, sheet: str, all_sheet: str, chk: Checks):
    """Read the balanced design + the all-300 sheet and check the design is intact."""
    bal = pd.read_excel(path, sheet_name=sheet)
    allx = pd.read_excel(path, sheet_name=all_sheet)
    for name, df in ((sheet, bal), (all_sheet, allx)):
        miss = [c for c in ("Sample_ID", "Source_Group", "Continent", "Genome_Size_bp",
                            "Contigs", "GC_percent") if c not in df]
        if miss:
            raise ExtractionError(f"sheet '{name}' lacks columns {miss}")
        df["Sample_ID"] = df["Sample_ID"].map(norm_sample)

    # 1. unique, valid IDs
    dup = bal[bal.Sample_ID.duplicated(keep=False)].Sample_ID.unique()
    chk.add("design_unique_sample_ids", FAIL if len(dup) or bal.Sample_ID.isna().any() else PASS,
            f"{len(bal)} rows" + (f"; duplicated {list(dup)}" if len(dup) else ""))
    # 2. every balanced sample passed QC and is marked INCLUDE_BALANCED
    bad = bal[(bal.get("Final_Decision", "KEEP").astype(str).str.upper() != "KEEP") |
              (bal.get("Final_Decision_v2", "INCLUDE_BALANCED").astype(str) != "INCLUDE_BALANCED")]
    chk.add("design_only_qc_kept_samples", FAIL if len(bad) else PASS,
            f"not KEEP/INCLUDE_BALANCED: {list(bad.Sample_ID)}" if len(bad) else "all KEEP + INCLUDE_BALANCED")
    # 3. balanced cells (Source_Group x Continent)
    cells = pd.crosstab(bal.Source_Group, bal.Continent)
    sizes = set(cells.values.ravel())
    chk.add("design_cells_balanced", PASS if len(sizes) == 1 else FAIL,
            f"{cells.shape[0]} groups x {cells.shape[1]} continents, cell size(s) {sorted(int(v) for v in sizes)}")
    unknown = set(bal.Source_Group) - set(GROUP_ORDER)
    if unknown:
        chk.add("design_known_groups", FAIL, f"unknown Source_Group {unknown}")
    # 4. balanced sheet agrees with the all-300 sheet
    m = bal.merge(allx, on="Sample_ID", how="left", suffixes=("", "_all"), indicator=True)
    n_bad = 0
    for c in ("Source_Group", "Continent", "Genome_Size_bp", "Contigs"):
        diff = m[m["_merge"].eq("both") & (m[c].astype(str) != m[f"{c}_all"].astype(str))]
        n_bad += len(diff)
        for r in diff.itertuples(index=False):
            chk.add("design_sheets_consistent", FAIL, f"{c}: {sheet} vs {all_sheet} differ", r.Sample_ID)
    not_in_all = m.loc[m["_merge"].eq("left_only"), "Sample_ID"].tolist()
    if not_in_all:
        chk.add("design_sheets_consistent", FAIL, f"missing from {all_sheet}: {not_in_all}")
    elif not n_bad:
        chk.add("design_sheets_consistent", PASS, f"{sheet} matches {all_sheet}")
    keep = [c for c in DESIGN_COLS if c in bal]
    return bal[keep].copy(), allx, cells


def build_sample_table(jobs: pd.DataFrame, design: pd.DataFrame, allx: pd.DataFrame,
                       only: list[str] | None, chk: Checks):
    final_ids = list(design.Sample_ID)
    if only:
        final_ids = [s for s in final_ids if s in set(only)]
        for s in sorted(set(only) - set(design.Sample_ID)):
            chk.add("only_ids_in_design", WARN, "requested with --only but not in the design", s)
    final_set = set(final_ids)

    # workspace jobs that are not part of the design -> reason from the all-300 sheet
    info = allx.set_index("Sample_ID")
    excluded = jobs[~jobs["Sample_ID"].isin(final_set)].copy()
    def why(sid):
        if only and sid in set(design.Sample_ID):
            return "in Balanced_210 but not requested in --only (test run)"
        if sid not in info.index:
            return "not in sampling design"
        r = info.loc[sid]
        dec = r.get("Final_Decision_v2", r.get("Final_Decision", ""))
        extra = r.get("Exclusion_Reason", "")
        if str(dec) == "EXCLUDE_BALANCING":
            extra = f"QC KEEP, removed to balance {r['Source_Group']} x {r['Continent']}"
        return f"{dec}: {extra}" if pd.notna(extra) and str(extra) else str(dec)
    excluded["Reason_Not_Extracted"] = excluded["Sample_ID"].map(why)

    missing = sorted(final_set - set(jobs["Sample_ID"]))
    for s in missing:
        chk.add("final_sample_has_cga_job", FAIL, "no CGA job folder found in workspace", s)

    samp = (pd.DataFrame({"Sample_ID": final_ids})
            .merge(jobs, on="Sample_ID", how="inner")
            .merge(design, on="Sample_ID", how="left"))
    samp.insert(1, "Isolation_Group", samp.pop("Source_Group"))
    samp.insert(2, "Continent", samp.pop("Continent"))
    samp.insert(3, "Collection_year", samp.pop("Collection_year") if "Collection_year" in samp else pd.NA)
    samp.insert(4, "Stratum", samp["Isolation_Group"] + "|" + samp["Continent"])
    for r in samp.itertuples(index=False):
        if r.Isolation_Group != r.Workspace_Group:
            chk.add("workspace_folder_matches_design_group", FAIL,
                    f"design group {r.Isolation_Group} but job is in '{GROUP_FOLDERS[r.Workspace_Group]}'",
                    r.Sample_ID)
    samp["Isolation_Group"] = pd.Categorical(samp["Isolation_Group"], GROUP_ORDER, ordered=True)
    samp = samp.sort_values(["Isolation_Group", "Continent", "Sample_ID"]).reset_index(drop=True)
    samp["Isolation_Group"] = samp["Isolation_Group"].astype(str)
    return samp, excluded, missing


# ════════════════════════════════════════════════════════════════════════
# STEP 3 — DOWNLOAD PER-GENOME DATA (cached / resumable)
# ════════════════════════════════════════════════════════════════════════
def fetch_genome_data(api: BVBRC, sid: str, gid: str, cache: Path, refresh: bool) -> dict:
    f = cache / "genomes" / f"{sid}_{gid}.json"
    if f.exists() and not refresh:
        d = json.loads(f.read_text())
        if d.get("schema") == CACHE_SCHEMA:     # files from older versions are re-fetched
            return d
    g = eqv("genome_id", gid)
    amr, vf, method, props = fetch_sp_genes(api, gid)
    data = {
        "schema": CACHE_SCHEMA, "sample_id": sid, "genome_id": gid,
        "fetched": datetime.now().isoformat(timespec="seconds"),
        "sp_method": method, "sp_properties_seen": props,
        "amr_genes": amr,
        "vf_genes": vf,
        "amr_phenotypes": api.query("genome_amr", g, AMR_FIELDS),
        "subsystems": api.query("subsystem", g, SUB_FIELDS),
    }
    f.parent.mkdir(parents=True, exist_ok=True)
    tmp = f.with_suffix(".tmp")
    tmp.write_text(json.dumps(data))
    tmp.replace(f)                      # atomic: a crash never leaves a half file
    return data


# ════════════════════════════════════════════════════════════════════════
# STEP 4 — VALIDATION
# ════════════════════════════════════════════════════════════════════════
def validate_identity(samp: pd.DataFrame, genomes: dict, qc_all: pd.DataFrame, chk: Checks):
    """Mapping checks + genome fingerprint (BV-BRC vs your own QC numbers)."""
    size_index = defaultdict(list)
    for r in qc_all.itertuples(index=False):
        size_index[int(r.Genome_Size_bp)].append(r.Sample_ID)

    status_level = {"CONFIRMED": PASS, "CONFIRMED_MULTI_API": WARN, "GTO_ONLY": WARN,
                    "API_ONLY": WARN, "CONFLICT": FAIL, "UNRESOLVED": FAIL,
                    "AMBIGUOUS_API": FAIL, "BAD_GTO_ID": FAIL}
    notes = {"CONFIRMED": "workspace GTO id == genome API id",
             "CONFIRMED_MULTI_API": "GTO id found, but account owns >1 genome with this KP id (re-run?)",
             "GTO_ONLY": "GTO id OK; genome_name has no KP id so route B could not confirm",
             "API_ONLY": "annotated.genome not readable; id from genome API only",
             "CONFLICT": "GTO id differs from genome API id",
             "UNRESOLVED": "no genome id found", "AMBIGUOUS_API": "several API ids, no GTO",
             "BAD_GTO_ID": "GTO id is not a genome id"}
    for r in samp.itertuples(index=False):
        chk.add("genome_id_resolution", status_level[r.ID_Status],
                f"{r.ID_Status}: {notes[r.ID_Status]} ({r.Genome_ID or r.API_Genome_IDs})",
                r.Sample_ID)

    for sid, want in KNOWN_ANCHORS.items():
        got = samp.loc[samp.Sample_ID == sid, "Genome_ID"]
        if len(got):
            chk.add("known_anchor_genome_id", PASS if got.iloc[0] == want else FAIL,
                    f"expected {want}, got {got.iloc[0]}", sid)

    dups = samp[samp.Genome_ID.ne("") & samp.duplicated("Genome_ID", keep=False)]
    for gid, g in dups.groupby("Genome_ID"):
        chk.add("genome_id_unique_across_samples", FAIL,
                f"{gid} assigned to {', '.join(g.Sample_ID)}")

    for r in samp.itertuples(index=False):
        if not r.Genome_ID:
            continue
        # GTO's own contigs (route A) vs your QC assembly size
        if pd.notna(r.GTO_Length) and r.GTO_Length and pd.notna(r.Genome_Size_bp):
            ok = int(r.GTO_Length) == int(r.Genome_Size_bp)
            chk.add("gto_length_vs_qc", PASS if ok else FAIL,
                    f"GTO {int(r.GTO_Length)} bp vs QC {int(r.Genome_Size_bp)} bp", r.Sample_ID)
        g = genomes.get(r.Genome_ID)
        if not g:
            chk.add("genome_visible_in_api", FAIL, f"{r.Genome_ID} not returned by /api/genome",
                    r.Sample_ID)
            continue
        if not str(g.get("taxon_id", "")).startswith(KP_TAXON_ID) and \
                "pneumoniae" not in str(g.get("species", "")).lower():
            chk.add("taxon_is_k_pneumoniae", WARN,
                    f"taxon {g.get('taxon_id')} / {g.get('species')}", r.Sample_ID)
        L, C, GC = g.get("genome_length"), g.get("contigs"), g.get("gc_content")
        if pd.isna(r.Genome_Size_bp):
            chk.add("fingerprint_vs_qc", WARN, "no QC numbers for this sample", r.Sample_ID)
            continue
        qL, qC, qGC = int(r.Genome_Size_bp), int(r.Contigs), float(r.GC_percent)
        if L is None:
            chk.add("fingerprint_vs_qc", WARN, "BV-BRC returned no genome_length", r.Sample_ID)
            continue
        # Two different K. pneumoniae genomes can differ by only ~0.2 % in size, so a
        # percentage tolerance alone cannot detect a swap. Rules:
        #   exact length (+contigs)               -> PASS
        #   length equals ANOTHER sample's length -> FAIL (swapped / mislabelled)
        #   BV-BRC slightly SMALLER, fewer contigs-> WARN (annotation dropped tiny contigs)
        #   anything else (bigger, or >0.2 % off) -> FAIL
        L = int(L)
        C = None if C is None else int(C)
        other = [s for s in size_index.get(L, []) if s != r.Sample_ID]
        if L == qL and (C is None or C == qC):
            gc_ok = GC is None or abs(float(GC) - qGC) <= 0.1
            chk.add("fingerprint_vs_qc", PASS if gc_ok else WARN,
                    f"length {L} bp, {C} contigs" + ("" if gc_ok else f", GC {GC} vs {qGC}"),
                    r.Sample_ID)
        elif other:
            chk.add("fingerprint_vs_qc", FAIL,
                    f"BV-BRC {L} bp/{C} contigs ≠ QC {qL} bp/{qC} contigs; this length is "
                    f"exactly QC sample {', '.join(other)} → genomes swapped or mislabelled",
                    r.Sample_ID)
        elif L < qL and (C is None or C < qC) and (qL - L) / qL <= 0.002:
            chk.add("fingerprint_vs_qc", WARN,
                    f"BV-BRC {L} bp/{C} contigs vs QC {qL} bp/{qC} contigs "
                    f"(−{qL - L} bp; small contigs probably filtered)", r.Sample_ID)
        else:
            chk.add("fingerprint_vs_qc", FAIL,
                    f"BV-BRC {L} bp/{C} contigs vs QC {qL} bp/{qC} contigs — wrong genome?",
                    r.Sample_ID)


def validate_content(sid: str, gid: str, d: dict, chk: Checks):
    for key in ("amr_genes", "vf_genes", "amr_phenotypes", "subsystems"):
        bad = {str(x.get("genome_id")) for x in d[key]} - {gid}
        if bad:
            chk.add("rows_belong_to_genome", FAIL, f"{key}: foreign genome ids {bad}", sid)
    # every K. pneumoniae chromosome carries intrinsic AMR genes (blaSHV/LEN/OKP,
    # fosA, oqxAB) and core virulence loci (e.g. ent, fim, mrk), so 0 rows is suspicious
    for key, label in (("amr_genes", "AMR specialty genes"),
                       ("vf_genes", "virulence specialty genes"),
                       ("subsystems", "subsystem rows")):
        if not d[key]:
            chk.add("non_empty_results", WARN, f"0 {label} — check CGA job finished", sid)


def live_recheck(api: BVBRC, samp: pd.DataFrame, data: dict, n: int, chk: Checks):
    """Re-query a random subset bypassing the cache; counts must be identical."""
    pool = [r for r in samp.itertuples(index=False) if r.Genome_ID]
    for r in random.Random(42).sample(pool, min(n, len(pool))):
        g = eqv("genome_id", r.Genome_ID)
        amr, vf, _, _ = fetch_sp_genes(api, r.Genome_ID, ["genome_id", "patric_id"])
        live = {"amr_genes": len(amr), "vf_genes": len(vf),
                "subsystems": len(api.query("subsystem", g, ["genome_id"]))}
        have = {k: len(data[r.Sample_ID][k]) for k in live}
        chk.add("live_requery_matches_output", PASS if live == have else FAIL,
                f"live {live} vs written {have}", r.Sample_ID)


# ════════════════════════════════════════════════════════════════════════
# STEP 5 — TABLES FOR ANALYSIS
# ════════════════════════════════════════════════════════════════════════
META = ["Sample_ID", "Isolation_Group", "Continent", "Collection_year", "Stratum"]


def long_table(samp, data, key, first_cols):
    """One row per BV-BRC record, each carrying the sample's design labels."""
    meta = samp.set_index("Sample_ID")[META[1:]]
    rows = []
    for sid, d in data.items():
        m = meta.loc[sid].to_dict()
        for x in d[key]:
            rows.append({"Sample_ID": sid, **m, **x})
    cols = META + first_cols
    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=cols)
    for c in cols:
        if c not in df:
            df[c] = ""
    df = df[cols + [c for c in df.columns if c not in cols]]
    for c in df.columns:                # lists (pmid, antibiotics) -> "a; b"
        if df[c].map(lambda v: isinstance(v, list)).any():
            df[c] = df[c].map(lambda v: "; ".join(map(str, v)) if isinstance(v, list) else v)
    df["Isolation_Group"] = pd.Categorical(df["Isolation_Group"], GROUP_ORDER, ordered=True)
    df = df.sort_values(["Isolation_Group", "Continent", "Sample_ID"]).reset_index(drop=True)
    df["Isolation_Group"] = df["Isolation_Group"].astype(str)
    return df
def gene_label(df):
    g = df["gene"].fillna("").astype(str).str.strip()
    return g.where(g != "", df["product"].fillna("").astype(str).str.strip())
def presence_tables(long_df: pd.DataFrame, samp: pd.DataFrame, extra_col: str):
    """Returns (binary sample x gene matrix, prevalence by group, prevalence by group x continent)."""
    base = samp[META + ["Genome_ID"]].copy()
    if long_df.empty:
        return base, pd.DataFrame(), pd.DataFrame()
    d = long_df.copy()
    d["Gene"] = gene_label(d)
    d["Feature"] = d["Gene"] + " [" + d["source"].astype(str) + "]"
    pres = d.groupby(["Sample_ID", "Feature"])["patric_id"].nunique().unstack(fill_value=0)
    pres = (pres > 0).astype(int)
    pres = pres[sorted(pres.columns, key=str.lower)]
    matrix = base.merge(pres, left_on="Sample_ID", right_index=True, how="left")
    feat_cols = list(pres.columns)
    matrix[feat_cols] = matrix[feat_cols].fillna(0).astype(int)

    info = d.groupby("Feature").agg(
        Gene=("Gene", "first"), Source=("source", "first"),
        **{extra_col: (extra_col, lambda s: "; ".join(sorted({str(v) for v in s if str(v) not in ("", "nan")})))},
        Product=("product", "first")).reset_index()

    # wide: one row per gene, n / N / % for each Isolation_Group
    N = samp.groupby("Isolation_Group")["Sample_ID"].nunique()
    prev = info.copy()
    for g in GROUP_ORDER:
        n = matrix.loc[matrix.Isolation_Group == g, feat_cols].sum()
        prev[f"{g}_n"] = prev["Feature"].map(n).fillna(0).astype(int)
        Ng = int(N.get(g, 0))
        prev[f"{g}_N"] = Ng
        prev[f"{g}_%"] = (100 * prev[f"{g}_n"] / Ng).round(1) if Ng else float("nan")
    prev["All_n"] = sum(prev[f"{g}_n"] for g in GROUP_ORDER)
    prev["All_N"] = int(N.sum())
    prev["All_%"] = (100 * prev["All_n"] / prev["All_N"]).round(1)
    prev = prev.sort_values(["Source", "All_n"], ascending=[True, False]).reset_index(drop=True)

    # long: one row per gene x group x continent (stratified tests, faceted heatmaps)
    cell = matrix.groupby(["Isolation_Group", "Continent"])[feat_cols].agg(["sum", "count"])
    rows = []
    for (g, c), r in cell.iterrows():
        for f in feat_cols:
            n, Nc = int(r[(f, "sum")]), int(r[(f, "count")])
            rows.append({"Feature": f, "Isolation_Group": g, "Continent": c, "n": n, "N": Nc,
                         "%": round(100 * n / Nc, 1) if Nc else float("nan")})
    strat = (pd.DataFrame(rows).merge(info[["Feature", "Gene", "Source", extra_col]], on="Feature")
             [["Feature", "Gene", "Source", extra_col, "Isolation_Group", "Continent", "n", "N", "%"]])
    strat["Isolation_Group"] = pd.Categorical(strat["Isolation_Group"], GROUP_ORDER, ordered=True)
    strat = strat.sort_values(["Source", "Feature", "Isolation_Group", "Continent"]).reset_index(drop=True)
    strat["Isolation_Group"] = strat["Isolation_Group"].astype(str)
    return matrix, prev, strat


def subsystem_summary(sub_long: pd.DataFrame) -> pd.DataFrame:
    if sub_long.empty:
        return sub_long
    keys = META + ["superclass", "class", "subclass", "subsystem_name"]
    s = sub_long.copy()
    for k in keys:
        s[k] = s[k].astype(str)
    return (s.groupby(keys, observed=True)["patric_id"].nunique()
            .rename("n_genes").reset_index())


# ════════════════════════════════════════════════════════════════════════
# EXCEL OUTPUT
# ════════════════════════════════════════════════════════════════════════
def write_excel(path: Path, sheets: dict[str, pd.DataFrame]):
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter
    hdr_fill = PatternFill("solid", start_color="1F4E79")
    hdr_font = Font(name="Arial", bold=True, color="FFFFFF", size=10)
    grp_fill = {"Animal": "FCE4D6", "Clinical": "DDEBF7", "Environmental": "E2EFDA"}
    stat_fill = {PASS: "C6EFCE", WARN: "FFEB9C", FAIL: "FFC7CE"}
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for name, df in sheets.items():
            df.to_excel(xw, sheet_name=name, index=False)
            ws = xw.sheets[name]
            for ci, col in enumerate(df.columns, 1):
                c = ws.cell(row=1, column=ci)
                c.fill, c.font = hdr_fill, hdr_font
                c.alignment = Alignment(wrap_text=True, vertical="center")
                sample = df[col].head(200).map(lambda v: len(str(v)))
                width = max([len(str(col))] + list(sample)) if len(df) else len(str(col))
                ws.column_dimensions[get_column_letter(ci)].width = min(max(width + 2, 8), 50)
            ws.freeze_panes = "C2" if "Isolation_Group" in df.columns[:2] else "A2"
            if len(df.columns):
                ws.auto_filter.ref = ws.dimensions
            # colour the group / status columns (cheap, only one column each)
            for key, fills in (("Isolation_Group", grp_fill), ("Status", stat_fill),
                               ("Overall_Status", stat_fill)):
                if key in df.columns and len(df) <= 200000:
                    ci = list(df.columns).index(key) + 1
                    for ri, v in enumerate(df[key].map(str), 2):
                        if v in fills:
                            ws.cell(row=ri, column=ci).fill = PatternFill("solid", start_color=fills[v])


# ════════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════════


def get_token() -> str:
    tok = os.environ.get("BVBRC_TOKEN", "").strip()
    typed = BVBRC_USERNAME not in {"", "your_username", "your@email.com"}
    if not tok and not typed:
        for f in (Path.home() / ".patric_token", Path.home() / ".bvbrc_token"):
            if f.exists():
                t = f.read_text().strip()
                if token_user(t) == WS_OWNER:
                    log.info(f"  using saved CLI token {f}")
                    tok = t
                    break
                log.info(f"  ignoring {f}: it belongs to '{token_user(t)}', not {WS_OWNER}")
    if not tok:
        placeholders = {"", "your_username", "yourpassword", "your@email.com"}
        user = BVBRC_USERNAME if BVBRC_USERNAME not in placeholders else os.environ.get("BVBRC_USERNAME", "")
        pw = BVBRC_PASSWORD if BVBRC_PASSWORD not in placeholders else os.environ.get("BVBRC_PASSWORD", "")
        if not user:
            user = input(f"BV-BRC username [{WS_OWNER.split('@')[0]}]: ").strip() or WS_OWNER.split("@")[0]
        if not pw:
            pw = getpass.getpass("BV-BRC password: ")
        log.info(f"  logging in as {user} ...")
        tok = BVBRC.login(user, pw)
    u = token_user(tok)
    if u != WS_OWNER:
        raise ExtractionError(f"Logged in as '{u}', but the CGA results belong to '{WS_OWNER}'. "
                              "Use the credentials of the new account.")
    return tok


def in_notebook() -> bool:
    try:
        from IPython import get_ipython
        return get_ipython() is not None and "IPKernelApp" in get_ipython().config
    except Exception:
        return False


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--design", default=DESIGN_FILE,
                    help="sampling-design workbook (default: %(default)s)")
    ap.add_argument("--design-sheet", default=DESIGN_SHEET,
                    help="sheet with the samples to extract (default: %(default)s)")
    ap.add_argument("--all-sheet", default=ALL_SHEET,
                    help="sheet with all 300 genomes, for fingerprints/exclusion reasons")
    ap.add_argument("--only", nargs="+", default=ONLY_SAMPLES or None,
                    help="test run on these Sample IDs, e.g. KP001 KP002")
    ap.add_argument("--out", default=OUTPUT_FILE)
    ap.add_argument("--cache", default=CACHE_DIR)
    ap.add_argument("--skip-gto", action="store_true", default=SKIP_GTO,
                    help="resolve ids via genome API only (fast, less strict)")
    ap.add_argument("--refresh", action="store_true", default=REFRESH_CACHE,
                    help="ignore cached genome data")
    ap.add_argument("--recheck", type=int, default=LIVE_RECHECK, help="random samples to re-query live")
    if argv is None and in_notebook():
        argv = []                       # Jupyter passes "-f kernel-xxxx.json": ignore it
    a, unknown = ap.parse_known_args(argv)
    if unknown:
        print(f"(ignoring unrecognised arguments: {' '.join(unknown)})")

    out = Path(a.out)
    if a.only and a.out == ap.get_default("out"):
        out = Path("KP_CGA_Balanced210_TEST.xlsx")
    cache = Path(a.cache)
    cache.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(message)s", force=True,
                        handlers=[logging.StreamHandler(sys.stdout),
                                  logging.FileHandler(out.with_suffix(".log"), mode="w", encoding="utf-8")])
    only = [norm_sample(x) for x in a.only] if a.only else None
    t0 = time.time()
    chk = Checks()
    log.info("=" * 70)
    log.info(f" BV-BRC CGA extractor v{__version__}   {datetime.now():%Y-%m-%d %H:%M}")
    log.info(f" Workspace: {WS_ROOT}")
    log.info("=" * 70)
    log.info(f"\n[design] {a.design} -> sheet '{a.design_sheet}'")
    if APPLY_SWAP:
        apply_design_swap(a.design, a.design_sheet, a.all_sheet, Path(a.out).stem)
    design, allx, cells = load_design(a.design, a.design_sheet, a.all_sheet, chk)
    log.info(f"  {len(design)} samples; cell size(s) {sorted(int(v) for v in set(cells.values.ravel()))}\n"
             + cells.to_string().replace("\n", "\n  ").join(["  ", ""]))
    if chk.n(FAIL):
        raise ExtractionError("sampling design failed its integrity checks: "
                              + "; ".join(f"{r['Check']} {r['Detail']}" for r in chk.rows if r["Status"] == FAIL))

    log.info("\n[0] Authenticating")
    api = BVBRC(get_token())
    log.info(f"  ✓ token user: {token_user(api.token)}")
    log.info("\n[1] Listing CGA jobs in the workspace")
    jobs = discover_jobs(api, chk)
    if jobs.empty:
        raise ExtractionError(f"No KPxxx job folders under {WS_ROOT}")
    log.info(f"  {len(jobs)} job folders named by Sample ID")

    log.info(f"\n[2] Matching workspace jobs to '{a.design_sheet}'")
    samp, excluded, missing = build_sample_table(jobs, design, allx, only, chk)
    log.info(f"  final set: {len(samp)} samples "
             f"({samp.Isolation_Group.value_counts().reindex(GROUP_ORDER).fillna(0).astype(int).to_dict()}), "
             f"{len(excluded)} workspace jobs left out, {len(missing)} missing")

    log.info("\n[3] Resolving Genome IDs (route A: annotated.genome, route B: genome API)")
    res, owned = resolve_genome_ids(api, samp[["Sample_ID", "Job_Path"]], cache,
                                    not a.skip_gto, chk)
    samp = samp.merge(res.drop(columns=["Job_Path"]), on="Sample_ID", how="left")
    ids = [g for g in samp.Genome_ID if g]
    genomes = {str(g["genome_id"]): g for g in owned}
    extra = [g for g in ids if g not in genomes]   # e.g. genomes not owned (shared)
    for i in range(0, len(extra), 100):
        chunk = ",".join(extra[i:i + 100])
        for g in api.query("genome", f"in(genome_id,({chunk}))", GENOME_FIELDS):
            genomes[str(g["genome_id"])] = g
    gdf = pd.DataFrame([genomes[g] for g in ids if g in genomes])
    if not gdf.empty:
        gdf = gdf.rename(columns={c: f"BVBRC_{c}" for c in gdf.columns if c != "genome_id"})
        samp = samp.merge(gdf, left_on="Genome_ID", right_on="genome_id", how="left").drop(columns="genome_id")
    validate_identity(samp, genomes, allx, chk)

    log.info("\n[4] Downloading CGA results per genome (cached in %s)" % cache)
    data = {}
    for i, r in enumerate(samp.itertuples(index=False), 1):
        if not r.Genome_ID:
            continue
        d = fetch_genome_data(api, r.Sample_ID, r.Genome_ID, cache, a.refresh)
        validate_content(r.Sample_ID, r.Genome_ID, d, chk)
        data[r.Sample_ID] = d
        log.info(f"  [{i:3d}/{len(samp)}] {r.Sample_ID} {r.Stratum:28s} {r.Genome_ID:14s} "
                 f"AMR {len(d['amr_genes']):3d} | VF {len(d['vf_genes']):3d} | "
                 f"pheno {len(d['amr_phenotypes']):3d} | subsys {len(d['subsystems']):5d}"
                 + ("" if d.get("sp_method") == "server filter" else f"  [{d.get('sp_method')}]"))
    if a.recheck:
        log.info(f"\n[5] Live re-query of {min(a.recheck, len(data))} random samples")
        live_recheck(api, samp[samp.Sample_ID.isin(data)], data, a.recheck, chk)

    log.info("\n[6] Building tables")
    samp_out = samp[samp.Sample_ID.isin(data) | samp.Genome_ID.eq("")].copy()
    amr = long_table(samp_out, data, "amr_genes",
                     ["genome_id", "gene", "product", "classification", "antibiotics_class",
                      "antibiotics", "source", "source_id", "evidence", "identity",
                      "query_coverage", "subject_coverage", "e_value", "patric_id"])
    vf = long_table(samp_out, data, "vf_genes",
                    ["genome_id", "gene", "product", "classification", "source", "source_id",
                     "evidence", "identity", "query_coverage", "subject_coverage",
                     "e_value", "patric_id"])
    pheno = long_table(samp_out, data, "amr_phenotypes",
                       ["genome_id", "antibiotic", "resistant_phenotype", "evidence",
                        "computational_method", "computational_method_performance"])
    sub = long_table(samp_out, data, "subsystems",
                     ["genome_id", "superclass", "class", "subclass", "subsystem_name",
                      "role_name", "gene", "product", "patric_id"])
    # denominators/matrix rows = genomes that actually have extracted data
    samp_ok = samp_out[samp_out.Sample_ID.isin(data)]
    amr_mat, amr_prev, amr_strat = presence_tables(amr, samp_ok, "antibiotics_class")
    vf_mat, vf_prev, vf_strat = presence_tables(vf, samp_ok, "classification")
    sub_sum = subsystem_summary(sub)

    # is the balanced design still intact after extraction?
    exp = cells.stack().rename("Design_n")
    got = samp_ok.groupby(["Isolation_Group", "Continent"]).size().rename("Extracted_n")
    if only:   # test run: compare only against the requested samples
        exp = design[design.Sample_ID.isin(only)].groupby(["Source_Group", "Continent"]).size().rename("Design_n")
    exp.index = exp.index.set_names(["Isolation_Group", "Continent"])
    balance = pd.concat([exp, got], axis=1).fillna(0).astype(int).reset_index()
    balance["Status"] = [PASS if d == e else FAIL for d, e in zip(balance.Design_n, balance.Extracted_n)]
    short = balance[balance.Status == FAIL]
    chk.add("balanced_design_complete", FAIL if len(short) else PASS,
            "all cells complete" if not len(short) else
            "; ".join(f"{r.Isolation_Group}|{r.Continent} {r.Extracted_n}/{r.Design_n}" for r in short.itertuples()))

    counts = pd.DataFrame({
        "Sample_ID": list(data),
        "AMR_gene_rows": [len(d["amr_genes"]) for d in data.values()],
        "VF_gene_rows": [len(d["vf_genes"]) for d in data.values()],
        "AMR_phenotype_rows": [len(d["amr_phenotypes"]) for d in data.values()],
        "Subsystem_rows": [len(d["subsystems"]) for d in data.values()],
        "SpecialtyGene_method": [d.get("sp_method", "") for d in data.values()],
        "SpecialtyGene_properties_seen": ["; ".join(f"{k}={v}" for k, v in d.get("sp_properties_seen", {}).items())
                                          for d in data.values()],
    })
    samp_out = samp_out.merge(counts, on="Sample_ID", how="left")
    vdf = chk.df()
    worst = (vdf.groupby("Sample_ID")["Status"]
             .agg(lambda s: FAIL if (s == FAIL).any() else WARN if (s == WARN).any() else PASS))
    samp_out["Overall_Status"] = samp_out["Sample_ID"].map(worst).fillna(PASS)
    front = META + ["Genome_ID", "ID_Status", "Overall_Status", "Assembly_Number",
             "BVBRC_genome_name", "BVBRC_mlst", "Genome_Size_bp", "BVBRC_genome_length",
             "Contigs", "BVBRC_contigs", "GC_percent", "BVBRC_gc_content"]
    front = [c for c in front if c in samp_out]
    samp_out = samp_out[front + [c for c in samp_out.columns if c not in front]]

    # excluded / missing
    ex = excluded[["Sample_ID", "Workspace_Group", "Job_Path", "Reason_Not_Extracted"]] \
        if len(excluded) else pd.DataFrame(columns=["Sample_ID", "Workspace_Group", "Job_Path", "Reason_Not_Extracted"])
    ex = pd.concat([ex, pd.DataFrame({"Sample_ID": missing, "Workspace_Group": "",
                                      "Job_Path": "", "Reason_Not_Extracted":
                                          "MISSING: in design but no CGA job in workspace"})],
                   ignore_index=True)
    n_fail, n_warn = chk.n(FAIL), chk.n(WARN)
    ready = n_fail == 0
    summary = [
        ("Status", "READY for analysis" if ready else f"NOT READY — {n_fail} FAIL checks (see Validation)"),
        ("Generated", datetime.now().strftime("%Y-%m-%d %H:%M")),
        ("Extractor version", __version__),
        ("BV-BRC account", WS_OWNER),
        ("Workspace root", WS_ROOT),
        ("Test run (--only)", ", ".join(only) if only else "no (full run)"),
        ("Sampling design", f"{Path(a.design).name} -> {a.design_sheet} ({len(design)} genomes, "
                            f"{cells.shape[0]} groups x {cells.shape[1]} continents x {cells.values.min()})"),
        ("Samples extracted", len(data)),
        *[(f"  {g}", int((samp_ok.Isolation_Group == g).sum())) for g in GROUP_ORDER],
        ("Workspace jobs not extracted", len(excluded)),
        ("Final samples missing in workspace", len(missing)),
        ("Validation PASS / WARN / FAIL", f"{chk.n(PASS)} / {n_warn} / {n_fail}"),
        ("AMR gene rows / VF gene rows", f"{len(amr)} / {len(vf)}"),
        ("AMR phenotype rows / subsystem rows", f"{len(pheno)} / {len(sub)}"),
        ("API calls this run", api.n_calls),
        ("Run time (min)", round((time.time() - t0) / 60, 1)),
        ("", ""),
        ("SHEETS", ""),
        ("Labels on every row", "Sample_ID, Isolation_Group, Continent, Collection_year, Stratum (= Group|Continent)"),
        ("Samples", "one row per genome: design labels, Genome ID, how it was resolved, "
                    "your QC metrics next to BV-BRC genome stats (length/contigs/GC/MLST), row counts"),
        ("Design_Balance", "genomes extracted per Isolation_Group x Continent vs the design (14 each)"),
        ("AMR_Genes", "BV-BRC specialty genes, property = Antibiotic Resistance (long format)"),
        ("Virulence_Genes", "BV-BRC specialty genes, property = Virulence Factor (long format)"),
        ("AMR_Phenotypes", "genome_amr records (computational predictions); may be empty for private genomes"),
        ("AMR_Gene_Matrix / VF_Gene_Matrix", "1/0 presence per sample; columns = 'gene [source]'"),
        ("AMR_Prevalence / VF_Prevalence", "per gene: n carriers, N genomes, % per Isolation_Group"),
        ("AMR_Prev_Group_Continent / VF_...", "long format: gene x Isolation_Group x Continent with n, N, % "
                                             "(for faceted heatmaps and continent-stratified tests)"),
        ("Subsystem_Summary", "genes per subsystem per sample (full gene list in *_subsystems.csv.gz)"),
        ("Validation", "every check, per sample, PASS/WARN/FAIL"),
        ("Not_Extracted", "workspace jobs left out (QC fail, duplicate, not selected) and missing jobs"),
        ("", ""),
        ("NOTES", ""),
        ("Gene sources", "Same gene can appear from CARD, NDARO and PATRIC under different names. "
                         "For prevalence heatmaps pick ONE source (filter Source) so genes are not double-counted."),
        ("Prevalence denominators", "N = number of extracted genomes in that Isolation_Group "
                                    "(70) or Group x Continent cell (14)."),
        ("Genome IDs are TEXT", "Read with pd.read_excel(f, sheet, dtype={'Genome_ID': str, 'genome_id': str}); "
                                "otherwise pandas turns 573.200100 into the float 573.2001."),
    ]
    readme = pd.DataFrame(summary, columns=["Item", "Value"])

    sheets = {"README": readme, "Samples": samp_out, "Design_Balance": balance,
              "AMR_Genes": amr, "Virulence_Genes": vf, "AMR_Phenotypes": pheno,
              "AMR_Gene_Matrix": amr_mat, "AMR_Prevalence": amr_prev, "AMR_Prev_Group_Continent": amr_strat,
              "VF_Gene_Matrix": vf_mat, "VF_Prevalence": vf_prev, "VF_Prev_Group_Continent": vf_strat,
              "Subsystem_Summary": sub_sum, "Validation": vdf, "Not_Extracted": ex}
    log.info(f"\n[7] Writing {out}")
    write_excel(out, sheets)
    sub_file = out.with_name(out.stem + "_subsystems.csv.gz")
    sub.to_csv(sub_file, index=False, compression="gzip")
    samp_out.to_csv(out.with_name(out.stem + "_samples.csv"), index=False)
    log.info("\n" + "=" * 70)
    for k, v in summary[:17]:
        log.info(f" {k:38s} {v}")
    if n_fail:
        log.info("\n FAIL checks:")
        for r in vdf[vdf.Status == FAIL].itertuples(index=False):
            log.info(f"   ✗ {r.Sample_ID:6s} {r.Check}: {r.Detail}")
    log.info(f"\n Outputs: {out}  |  {sub_file.name}  |  {out.with_suffix('.log').name}")
    log.info("=" * 70)
    return 0 if ready else 2


def run(argv=None) -> int:
    """Run and return an exit code (0 READY, 2 FAIL checks, 1 error) without killing Jupyter."""
    try:
        return main(argv)
    except ExtractionError as e:
        logging.getLogger("cga").error(f"\n✗ {e}")
        return 1
    except KeyboardInterrupt:
        print("\nInterrupted — run again to resume from the cache.")
        return 130


if __name__ == "__main__":
    code = run()
    if not in_notebook():               # in Jupyter, sys.exit() shows "SystemExit" noise
        sys.exit(code)
