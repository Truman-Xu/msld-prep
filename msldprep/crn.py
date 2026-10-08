"""Charge renormalization (CRN): make each site's substituents -- and the
shared core -- carry a chemically consistent, mutually-compatible set of
partial charges, given each ligand's original CGenFF charges.

Ported from msld_crn.py:303-521, but `redistribute_charge_proportional`
replaces legacy's cyclic micro-charge-redistribution loop (add one +-1e-6
quantum at a time to successive non-lone-pair atoms, wrapping around the
atom list, until the target total is reached -- an O(qdiff/1e-6) loop that
also biased the correction toward whichever atoms happened to come first in
the list) with a magnitude-weighted proportional correction computed in
fixed-point Decimal arithmetic: atoms that already carry a larger partial
charge absorb proportionally more of the correction, and the result is
still guaranteed to sum to the target exactly after rounding. This one
function is shared between the per-site fragment CRN and the whole-core
CRN, which legacy duplicated.

Note (fixes msld_crn.py:482's inconsistency): the core-level lone-pair skip
in legacy checked for an atom named exactly "LP", while the fragment-level
skip checked for a name *prefix* of "LP" (matching real names like "LP1").
The core-level check was therefore inert. Both now use the same prefix
check.
"""

from __future__ import annotations

import logging
from decimal import Decimal, ROUND_HALF_UP
from typing import Sequence

import numpy as np

from .datamodel import CoreDecomposition, CRNResult, LigandRecord
from .errors import CRNError

logger = logging.getLogger(__name__)

DEFAULT_DECIMALS = 6


def _is_lone_pair(atom_name: str) -> bool:
    return atom_name.startswith("LP")


def _to_decimal(value: float, quantum: Decimal) -> Decimal:
    # str(value) round-trips through the same decimal text the charge was
    # originally read from (e.g. rtf "-0.170"), avoiding the long noisy
    # expansion Decimal(value) would give for a binary float.
    return Decimal(str(value)).quantize(quantum, rounding=ROUND_HALF_UP)


def redistribute_charge_proportional(
    atom_order: Sequence[str],
    charges: dict[str, float],
    target_total: float,
    *,
    decimals: int = DEFAULT_DECIMALS,
) -> dict[str, float]:
    """Adjust `charges` so their sum over `atom_order` equals `target_total`,
    distributing the correction across non-lone-pair atoms in proportion to
    each atom's own charge magnitude -- an atom that already carries a large
    partial charge absorbs proportionally more of the correction than one
    close to zero, keeping the *relative* distortion of every atom's charge
    small. Lone pairs (name prefix "LP") are never adjusted.

    Uses fixed-point Decimal arithmetic quantized to `decimals` places, so
    the *rounded* result is guaranteed to sum to `target_total` exactly --
    a plain float division would leave each atom's independently-rounded
    share summing to a precision unit off from the target.

    Returns a new dict (does not mutate `charges`). Raises CRNError if a
    nonzero correction is needed but `atom_order` has no non-lone-pair atoms.
    """
    quantum = Decimal(1).scaleb(-decimals)
    dec_charges = {name: _to_decimal(value, quantum) for name, value in charges.items()}
    current_total = sum((dec_charges[name] for name in atom_order), Decimal(0))
    target = _to_decimal(target_total, quantum)
    qdiff = target - current_total

    real_atoms = [name for name in atom_order if not _is_lone_pair(name)]
    if qdiff != 0 and not real_atoms:
        msg = "Cannot redistribute charge: atom_order has no non-lone-pair atoms."
        logger.error(msg)
        raise CRNError(msg)

    if real_atoms and qdiff != 0:
        weights = {name: abs(dec_charges[name]) for name in real_atoms}
        weight_total = sum(weights.values())
        if weight_total == 0:
            # every real atom is exactly zero -- split the correction evenly instead
            weights = {name: Decimal(1) for name in real_atoms}
            weight_total = Decimal(len(real_atoms))

        deltas = {
            name: (qdiff * weights[name] / weight_total).quantize(quantum, rounding=ROUND_HALF_UP)
            for name in real_atoms
        }

        # Rounding each proportional share to `decimals` can leave a residual
        # of a few quanta; sweep it onto the largest-magnitude atoms first
        # (bounded by len(real_atoms) steps, unlike legacy's O(qdiff/1e-6)
        # walk) so the total still lands exactly on target.
        residual = qdiff - sum(deltas.values())
        if residual != 0:
            ordered = sorted(real_atoms, key=lambda n: weights[n], reverse=True)
            step = quantum if residual > 0 else -quantum
            n_steps = int(abs(residual) / quantum)
            for i in range(n_steps):
                deltas[ordered[i % len(ordered)]] += step

        for name, delta in deltas.items():
            dec_charges[name] += delta

    final_total = sum((dec_charges[name] for name in atom_order), Decimal(0))
    if final_total != target:
        msg = (
            f"Charge redistribution failed: total ({float(final_total)}) doesn't match "
            f"target ({float(target)})."
        )
        logger.error(msg)
        raise CRNError(msg)

    return {name: float(value) for name, value in dec_charges.items()}


def _log_charge_diff(label: str, orig: float, new_total: float, q_cut_percent: float) -> None:
    diff = round(new_total - orig, DEFAULT_DECIMALS)
    pct_diff = 0.0 if orig == 0.0 else round(diff / orig * 100.0, 2)
    if abs(pct_diff) > q_cut_percent:
        logger.warning(
            "%s: CRN charge diff %.2f%% exceeds %.1f%% threshold (orig=%.6f, new=%.6f)",
            label, pct_diff, q_cut_percent, orig, new_total,
        )
    else:
        logger.info(
            "%s: Q(diff from avg)=%.6f Q(orig)=%.6f (%.2f%% diff)",
            label, diff, orig, pct_diff,
        )


def compute_crn(
    ligands: Sequence[LigandRecord],
    decomposition: CoreDecomposition,
    *,
    chk_q_change: bool = True,
    q_cut_percent: float = 5.0,
    decimals: int = DEFAULT_DECIMALS,
) -> CRNResult:
    """Compute charge-renormalized core + per-site fragment charges.

    Pure function: does not mutate `ligands`' Mol objects (see
    `apply_crn_charges` to write results back onto RDKit atoms).
    """
    ligand_by_name = {lig.name: lig for lig in ligands}
    reference = ligand_by_name[decomposition.reference_name]
    ref_charge = reference.net_charge

    charge_changed: dict[str, bool] = {}
    for ligand in ligands:
        changed = ligand.net_charge != ref_charge
        charge_changed[ligand.name] = changed
        if changed:
            logger.warning(
                "Charge change between reference (%s, QNET=%s) and %s (QNET=%s).",
                reference.name, ref_charge, ligand.name, ligand.net_charge,
            )
            if not chk_q_change:
                logger.warning("Recommended to enable chk_q_change when net charges differ.")
    any_charge_changed = any(charge_changed.values())

    # (1) core charge mean/stdev per reference-canonical core position (real atoms + LPs)
    core_charge_stats: dict[str, tuple[float, float]] = {}
    core_mean_charges: dict[str, float] = {}
    core_names_for_sum: list[str] = list(decomposition.core_atom_order)

    for pos, ref_name in enumerate(decomposition.core_atom_order):
        values = [
            ligand.atom_charges[decomposition.core_atom_names[ligand.name][pos]]
            for ligand in ligands
        ]
        mean = round(float(np.mean(values)), decimals)
        stdev = round(float(np.std(values)), decimals)
        core_charge_stats[ref_name] = (mean, stdev)
        core_mean_charges[ref_name] = mean

    for ref_host_name, lp_by_ligand in decomposition.core_lp_atom_names.items():
        lp_name = lp_by_ligand[reference.name]
        core_names_for_sum.append(lp_name)
        values = [
            ligand.atom_charges[lp_by_ligand[ligand.name]]
            for ligand in ligands
            if ligand.name in lp_by_ligand
        ]
        mean = round(float(np.mean(values)), decimals)
        stdev = round(float(np.std(values)), decimals)
        core_charge_stats[lp_name] = (mean, stdev)
        core_mean_charges[lp_name] = mean

    # (2) per-site fragment CRN
    fragment_charges: list[dict[str, dict[str, float]]] = []
    site_average: list[float] = []
    site_integer_charge: list[int] = []

    for site in decomposition.sites:
        orig_charges_by_ligand: dict[str, dict[str, float]] = {}
        orig_sums: dict[str, float] = {}
        for ligand_name, fragment in site.fragments.items():
            all_names = [*fragment.atom_names, *fragment.lp_atom_names]
            charges = {name: ligand_by_name[ligand_name].atom_charges[name] for name in all_names}
            orig_charges_by_ligand[ligand_name] = charges
            orig_sums[ligand_name] = round(sum(charges.values()), decimals)

        int_charges = {name: 0 for name in orig_sums}
        adjusted_sums = dict(orig_sums)
        site_int = 0
        if chk_q_change and any_charge_changed:
            int_charges = {name: int(round(value)) for name, value in orig_sums.items()}
            adjusted_sums = {
                name: round(orig_sums[name] - int_charges[name], decimals) for name in orig_sums
            }
            site_int = int(round(float(np.mean(list(int_charges.values())))))

        site_avg = round(float(np.mean(list(adjusted_sums.values()))), decimals)
        site_average.append(site_avg)
        site_integer_charge.append(site_int)

        site_fragment_charges: dict[str, dict[str, float]] = {}
        for ligand_name, fragment in site.fragments.items():
            all_names = [*fragment.atom_names, *fragment.lp_atom_names]
            target = round(site_avg + int_charges[ligand_name], decimals)
            new_charges = redistribute_charge_proportional(
                all_names, orig_charges_by_ligand[ligand_name], target,
                decimals=decimals,
            )
            site_fragment_charges[ligand_name] = new_charges
            _log_charge_diff(
                f"Site {site.index} fragment {ligand_name}",
                adjusted_sums[ligand_name], site_avg, q_cut_percent,
            )
        fragment_charges.append(site_fragment_charges)

    # (3) core CRN to neutralize the whole (reference) ligand
    site_sum = round(sum(site_average), decimals)
    int_sum = round(sum(site_integer_charge), decimals)
    core_sum = round(sum(core_mean_charges[name] for name in core_names_for_sum), decimals)
    target_core_sum = round(ref_charge - site_sum - int_sum, decimals)

    core_charges = redistribute_charge_proportional(
        core_names_for_sum, core_mean_charges, target_core_sum,
        decimals=decimals,
    )
    final_core_sum = round(sum(core_charges[name] for name in core_names_for_sum), decimals)
    _log_charge_diff("Core", core_sum, final_core_sum, q_cut_percent)

    return CRNResult(
        core_charges=core_charges,
        core_charge_stats=core_charge_stats,
        fragment_charges=fragment_charges,
        site_average=site_average,
        site_integer_charge=site_integer_charge,
        charge_changed=charge_changed,
    )


def apply_crn_charges(
    ligands: Sequence[LigandRecord], decomposition: CoreDecomposition, crn: CRNResult
) -> dict[str, dict[str, float]]:
    """Write CRN'd charges back onto RDKit atoms as `CRNPartialCharge`.

    Lone pairs have no RDKit Atom to attach a prop to; their CRN'd charges
    (unperturbed by redistribution, but present in `crn.core_charges` /
    `crn.fragment_charges`) are instead returned as {ligand_name: {lp_name: charge}}.
    """
    lp_charges: dict[str, dict[str, float]] = {lig.name: {} for lig in ligands}
    name_to_atom = {
        ligand.name: {
            atom.GetProp("_TriposAtomName"): atom
            for atom in ligand.mol.GetAtoms()
            if atom.HasProp("_TriposAtomName")
        }
        for ligand in ligands
    }

    for pos, ref_name in enumerate(decomposition.core_atom_order):
        charge = crn.core_charges[ref_name]
        for ligand in ligands:
            local_name = decomposition.core_atom_names[ligand.name][pos]
            name_to_atom[ligand.name][local_name].SetDoubleProp("CRNPartialCharge", charge)

    for ref_host_name, lp_by_ligand in decomposition.core_lp_atom_names.items():
        for ligand_name, lp_name in lp_by_ligand.items():
            lp_charges[ligand_name][lp_name] = crn.core_charges[lp_name]

    for site, site_charges in zip(decomposition.sites, crn.fragment_charges):
        for ligand_name, charges in site_charges.items():
            fragment = site.fragments[ligand_name]
            for name in fragment.atom_names:
                name_to_atom[ligand_name][name].SetDoubleProp("CRNPartialCharge", charges[name])
            for name in fragment.lp_atom_names:
                lp_charges[ligand_name][name] = charges[name]

    return lp_charges
