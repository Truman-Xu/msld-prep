import pytest

from msldprep.naming import NamingError, allocate_unique_atom_name, new_counters


def test_single_char_element_naming():
    counters = new_counters()
    assert allocate_unique_atom_name("C", counters) == "C001"
    assert allocate_unique_atom_name("C", counters) == "C002"
    # a different symbol draws from its own counter, not the shared
    # position-in-system count -- it starts back at 001.
    assert allocate_unique_atom_name("N", counters) == "N001"


def test_two_char_element_naming_uppercases():
    counters = new_counters()
    assert allocate_unique_atom_name("Cl", counters) == "CL0A"
    assert allocate_unique_atom_name("Br", counters) == "BR0A"
    assert allocate_unique_atom_name("cl", counters) == "CL0B"  # lowercase input also uppercased


def test_lp_naming():
    counters = new_counters()
    assert allocate_unique_atom_name("LP", counters) == "LP0A"
    assert allocate_unique_atom_name("lp", counters) == "LP0B"


def test_names_are_unique_across_calls():
    counters = new_counters()
    names = [allocate_unique_atom_name(sym, counters) for sym in ["C", "C", "Cl", "N", "Br", "C"]]
    assert len(names) == len(set(names))


def test_single_char_overflow_raises():
    counters = new_counters()
    counters["C"] = [1000]
    with pytest.raises(NamingError):
        allocate_unique_atom_name("C", counters)


def test_single_char_counters_are_independent_per_symbol():
    # regression: an earlier version shared one counter across every
    # single-character element, so a system with >999 atoms overall could
    # fail even with far fewer than 999 of any individual element.
    counters = new_counters()
    counters["C"] = [1000]
    with pytest.raises(NamingError):
        allocate_unique_atom_name("C", counters)
    assert allocate_unique_atom_name("N", counters) == "N001"


def test_three_char_symbol_raises():
    counters = new_counters()
    with pytest.raises(NamingError):
        allocate_unique_atom_name("Xyz", counters)
