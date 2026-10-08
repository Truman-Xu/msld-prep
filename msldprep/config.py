"""System-level configuration for an MSLD build: block counts, replica/node
topology, and the handful of CHARMM scalar defaults (fnex, temp, box, ...).

Ported from values previously scattered across msld_crn.py's nsubs/nblocks/
nreps files, msld_wrt.py's name/nnodes/ncentral files, and literals baked
directly into msld_wrt.py's generated .inp script (msld_wrt.py:83,96-102).
Consolidated into one JSON-dumpable dataclass; `write_legacy_config_files`
still emits the individual flat files some out-of-scope consumers read
directly (Lg_Solvate.sh reads `outdir/name`; the external ALF driver
conventionally reads `outdir/nnodes`/`outdir/ncentral` as flat files).

Note: `nreps`/`nnodes`/`ncentral` are recorded here for the JSON dump and the
legacy flat files, but must NOT be baked as literal values into the
generated .inp script -- they stay as unresolved `@nreps`/`@ncentral` CHARMM
variable references there, since the external ALF driver sets them per-run
(see writers/alf_script.py).
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import asdict, dataclass, field

logger = logging.getLogger(__name__)


@dataclass
class SystemConfig:
    sysname: str
    outdir: str
    cgenff: bool
    nsubs: list[int] = field(default_factory=list)
    nreps: int = 1
    nnodes: int = 1
    ncentral: int = 0
    fnex: float = 5.5
    box: str = "????"
    temp: float = 298.15
    ligseg: str = "lig"
    resnum: int = 1
    builddir: str = "prep"

    @property
    def nsites(self) -> int:
        return len(self.nsubs)

    @property
    def nblocks(self) -> int:
        return sum(self.nsubs)


def write_config_json(config: SystemConfig, outdir: str, filename: str = "config.json") -> None:
    path = os.path.join(outdir, filename)
    with open(path, "w") as fp:
        json.dump(asdict(config), fp, indent=2)
    logger.debug("Wrote system config to %s", path)


def read_config_json(path: str) -> SystemConfig:
    with open(path) as fp:
        data = json.load(fp)
    return SystemConfig(**data)


def write_legacy_config_files(config: SystemConfig, outdir: str) -> None:
    """Backward-compat flat files: outdir/name, outdir/nnodes, outdir/ncentral.

    Ported from msld_wrt.py:48-57. `outdir/name` is read directly by the
    out-of-scope Lg_Solvate.sh; `nnodes`/`ncentral` are the conventional
    flat-file inputs to the external ALF driver.
    """
    with open(os.path.join(outdir, "name"), "w") as fp:
        fp.write(config.sysname)
    with open(os.path.join(outdir, "nnodes"), "w") as fp:
        fp.write(str(config.nnodes))
    with open(os.path.join(outdir, "ncentral"), "w") as fp:
        fp.write(str(config.ncentral))
