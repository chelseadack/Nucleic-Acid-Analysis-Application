"""Unit tests for GUI_logic. Run with:  python3 -m pytest"""
import pytest

from GUI_logic import (StrandInput, complex_label, complex_sort_key, convert_sequence, default_strand_name,
    dot_bracket_pairs, find_input_problems, has_only_valid_bases)


# ------------------------------
# SEQUENCES
# ------------------------------
@pytest.mark.parametrize("sequence, material, expected", [
    ("ATGC", "DNA", True),
    ("AUGC", "RNA", True),
    ("AUGC", "DNA", False),  # U is not a DNA base
    ("ATGC", "RNA", False),  # T is not an RNA base
    ("ATGN", "DNA", False),
    ("", "DNA", True),  # empty rows are not flagged while typing
])
def test_has_only_valid_bases(sequence, material, expected):
    assert has_only_valid_bases(sequence, material) == expected


def test_convert_sequence_round_trip():
    assert convert_sequence("atgcTT", "RNA") == "AUGCUU"
    assert convert_sequence("AUGCUU", "DNA") == "ATGCTT"
    assert convert_sequence(convert_sequence("GATTACA", "RNA"), "DNA") == "GATTACA"


def test_default_strand_name_is_one_based():
    assert default_strand_name(0) == "strand1"


# ------------------------------
# INPUT VALIDATION
# ------------------------------
def test_valid_input_has_no_problems():
    strands = [StrandInput("a", "ATGC", 1.0), StrandInput("b", "GCAT", 0.5)]
    assert find_input_problems(strands, "DNA") == []


def test_no_strands():
    assert find_input_problems([], "DNA") == ["Add at least one strand."]


def test_empty_sequence_is_rejected():
    # an empty sequence crashes NUPACK, so it must be caught before the analysis
    problems = find_input_problems([StrandInput("a", "", 1.0)], "DNA")
    assert problems == ["'a' has no sequence."]


def test_invalid_bases_are_rejected():
    problems = find_input_problems([StrandInput("a", "AUGC", 1.0)], "DNA")
    assert len(problems) == 1 and "bases other than" in problems[0]


def test_duplicate_names_are_rejected():
    strands = [StrandInput("a", "ATGC", 1.0), StrandInput("a", "GCAT", 1.0)]
    assert find_input_problems(strands, "DNA") == ["The name 'a' is used more than once."]


@pytest.mark.parametrize("name", ["a+b", "a(1)", "<a>"])
def test_names_that_break_complex_labels_are_rejected(name):
    problems = find_input_problems([StrandInput(name, "ATGC", 1.0)], "DNA")
    assert len(problems) == 1 and "names cannot contain" in problems[0]


# ------------------------------
# COMPLEXES & STRUCTURES
# ------------------------------
def test_complex_label():
    assert complex_label("<Complex (a+b)>") == "a+b"
    assert complex_label("<Complex (strand1)>") == "strand1"


def test_complexes_sort_by_size_then_name():
    labels = ["b+c", "a", "a+b+c", "c", "a+b"]
    assert sorted(labels, key=complex_sort_key) == ["a", "c", "a+b", "b+c", "a+b+c"]


def test_dot_bracket_pairs_hairpin():
    assert dot_bracket_pairs("((..))") == {0: 5, 5: 0, 1: 4, 4: 1}


def test_dot_bracket_pairs_across_strands():
    # '+' is removed, so indices run continuously across both strands
    assert dot_bracket_pairs("((+))") == {0: 3, 3: 0, 1: 2, 2: 1}


def test_dot_bracket_pairs_unpaired():
    assert dot_bracket_pairs("....") == {}
