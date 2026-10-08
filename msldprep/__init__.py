"""msldprep -- RDKit-based, in-memory ligand preparation for MSLD.

Pipeline: align.py (load/type/align ligands, find shared core) ->
decompose.py (anchor/site/fragment decomposition) -> crn.py (charge
renormalization) -> writers/ (CHARMM/PDB output, optional).

No logging is configured here beyond a NullHandler, per standard library
convention -- attach handlers in the application that uses this package.
"""

import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())
