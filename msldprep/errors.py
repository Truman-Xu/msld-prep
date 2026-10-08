"""Typed exception hierarchy for msldprep."""


class MsldPrepError(Exception):
    """Base class for all msldprep errors."""


class AlignmentError(MsldPrepError):
    """Raised when loading, typing, or 3D-aligning a ligand fails."""


class DecompositionError(MsldPrepError):
    """Raised when core/anchor/site/fragment decomposition fails or is inconsistent."""


class CRNError(MsldPrepError):
    """Raised when charge renormalization produces an inconsistent result."""


class WriterError(MsldPrepError):
    """Raised when a CHARMM/PDB writer cannot produce valid output from its input records."""


class PrmError(MsldPrepError):
    """Raised when .prm parameter parsing/merging/derivation fails or is inconsistent."""


class ExternalToolError(MsldPrepError):
    """Raised when an external tool (cgenff) is missing or exits with an error."""
