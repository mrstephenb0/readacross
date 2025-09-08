#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# ra_core/java_tools.py
import os, subprocess, tempfile, shutil
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parents[1]

def _java_bin():
    # Prefer env override; else system java
    jb = os.environ.get("JAVA_BIN") or shutil.which("java")
    if not jb:
        raise RuntimeError("Java not available. Set JAVA_BIN or install OpenJDK.")
    return jb

def run_biotransformer(smiles: str, *, xmx="2G", timeout=120) -> list[str]:
    """Example: run your JAR via CLI and return raw lines (adjust flags/parse as needed)."""
    jar = APP_ROOT / "tools" / "biotransformer" / "BioTransformer.jar"
    if not jar.exists():
        raise RuntimeError(f"JAR not found: {jar}")

    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "bt_out.tsv"
        cmd = [
            _java_bin(), f"-Xmx{xmx}", "-jar", str(jar),
            "--smiles", smiles,
            "--task", "metabolism",
            "--output", str(out),
        ]
        res = subprocess.run(cmd, text=True, capture_output=True, timeout=timeout)
        if res.returncode != 0:
            # You can log res.stderr
            return []
        if not out.exists():
            return []
        return out.read_text(encoding="utf-8", errors="ignore").splitlines()

