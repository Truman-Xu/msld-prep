"""CHARMM control/bookkeeping writers: lpsites.inp, translation.txt, nsubs/nblocks/nreps.

Ported from msld_crn.py:797-798,841-847,896-903 (lpsites.inp, interleaved
with the rtf writers there; independent here per-writer-callable-alone),
:739-755 (translation.txt), :1001-1012 (nsubs/nblocks/nreps).
"""

from __future__ import annotations

import os
from typing import Sequence

from ..datamodel import CoreDecomposition, LigandRecord
from .records import (
    AtomNameTranslation,
    core_translation_for_ligand,
    compute_core_atom_block,
    compute_fragment_atom_block,
    resolve_translated_name,
)


def _lonepair_block(lp_translated: str, host1_translated: str, host2_translated: str, distance: float) -> str:
    # legacy always hardcodes SCAL 0.00 regardless of the source rtf's actual
    # scale value (msld_crn.py:847,903) -- preserved here deliberately.
    return (
        f"LONEPAIR COLI sele atom @ligseg @resnum {lp_translated} end -\n"
        f"              sele atom @ligseg @resnum {host1_translated} end -\n"
        f"              sele atom @ligseg @resnum {host2_translated} end -\n"
        f"         DIST {distance} SCAL 0.00\ncoor shake\n\n"
    )


def write_lpsites_inp(
    decomposition: CoreDecomposition,
    ligands: Sequence[LigandRecord],
    core_charges: dict[str, float],
    translation: AtomNameTranslation,
    outdir: str,
) -> None:
    """Always creates the file, even with zero lone pairs (msld_wrt.py streams it unconditionally)."""
    ligand_by_name = {lig.name: lig for lig in ligands}
    reference = ligand_by_name[decomposition.reference_name]

    lines = ["* Load LP Site Definitions (if applicable)\n*\n\n"]

    core_entries = compute_core_atom_block(decomposition, reference, core_charges, translation)
    for entry in core_entries:
        if not entry.is_lp:
            continue
        lp = entry.lp_record
        host1 = translation.core[lp.host_atoms[0]]
        host2 = translation.core[lp.host_atoms[1]]
        lines.append(_lonepair_block(entry.translated_name, host1, host2, lp.distance))

    for site in decomposition.sites:
        for ligand_name, fragment in site.fragments.items():
            ligand = ligand_by_name[ligand_name]
            # fragment_charges not needed for lpsites.inp content, only names/host refs
            dummy_charges = {name: 0.0 for name in [*fragment.atom_names, *fragment.lp_atom_names]}
            entries = compute_fragment_atom_block(site.index, ligand, fragment, dummy_charges, translation)
            local_core = core_translation_for_ligand(decomposition, translation, ligand_name)
            frag_translation = translation.fragments[(site.index, ligand_name)]
            for entry in entries:
                if not entry.is_lp:
                    continue
                lp = entry.lp_record
                host1 = resolve_translated_name(
                    lp.host_atoms[0], local_core_translation=local_core, fragment_translation=frag_translation
                )
                host2 = resolve_translated_name(
                    lp.host_atoms[1], local_core_translation=local_core, fragment_translation=frag_translation
                )
                lines.append(_lonepair_block(entry.translated_name, host1, host2, lp.distance))

    with open(os.path.join(outdir, "lpsites.inp"), "w") as fp:
        fp.writelines(lines)


def write_translation_table(
    decomposition: CoreDecomposition,
    translation: AtomNameTranslation,
    outdir: str,
) -> None:
    lines = ["Original Atom Name -> New Atom Name\n", "CORE\n"]
    for ref_name, new_name in translation.core.items():
        lines.append(f"{ref_name} {new_name}\n")
    lines.append("\n")
    for site in decomposition.sites:
        for frag_number, ligand_name in enumerate(site.fragments, start=1):
            lines.append(
                f"SITE {site.index + 1} site{site.index + 1}_sub{frag_number} (from {ligand_name})\n"
            )
            frag_translation = translation.fragments[(site.index, ligand_name)]
            for raw_name, new_name in frag_translation.items():
                lines.append(f"{raw_name} {new_name}\n")
            lines.append("\n")
        lines.append("\n")

    with open(os.path.join(outdir, "translation.txt"), "w") as fp:
        fp.writelines(lines)


def write_block_counts(decomposition: CoreDecomposition, outdir: str) -> None:
    nsubs = [len(site.fragments) for site in decomposition.sites]
    nblocks = sum(nsubs)
    with open(os.path.join(outdir, "nsubs"), "w") as fp:
        fp.write(" ".join(str(n) for n in nsubs) + " ")
    with open(os.path.join(outdir, "nblocks"), "w") as fp:
        fp.write(str(nblocks))
    with open(os.path.join(outdir, "nreps"), "w") as fp:
        fp.write("1")
