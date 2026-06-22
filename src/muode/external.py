"""Helpers for shelling out to external bioinformatics tools.

The upstream phases (gene calling, CarveMe, gapseq, CheckM2, MetaPathPredict,
DLKcat, memote) are not Python libraries we import; they are programs invoked in
their own conda environments by the Snakemake workflow.  These helpers give the
wrappers a consistent way to locate a binary and run it with clear errors, so a
missing tool fails loudly with an actionable message instead of a cryptic
traceback.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Sequence


class MuodeToolError(RuntimeError):
    """Raised when an external tool is missing or exits non-zero."""


def require(binary: str, env_hint: str = "") -> str:
    """Return the resolved path to ``binary`` or raise with guidance."""
    path = shutil.which(binary)
    if path is None:
        hint = f" (install it, e.g. `{env_hint}`)" if env_hint else ""
        raise MuodeToolError(f"required tool '{binary}' not found on PATH{hint}")
    return path


def run(cmd: Sequence[str], log: str | Path | None = None, check: bool = True) -> subprocess.CompletedProcess:
    """Run a subprocess, optionally teeing combined output to ``log``."""
    proc = subprocess.run(
        [str(c) for c in cmd], capture_output=True, text=True
    )
    if log is not None:
        Path(log).write_text((proc.stdout or "") + (proc.stderr or ""))
    if check and proc.returncode != 0:
        raise MuodeToolError(
            f"command failed ({proc.returncode}): {' '.join(map(str, cmd))}\n{proc.stderr[-2000:]}"
        )
    return proc
