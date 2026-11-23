"""Core logic for the GUI.

These functions do not import Qt, NUPACK or NucDraw, so they can be unit tested
(see test_GUI_logic.py) without installing the full application.
"""
from __future__ import annotations

from dataclasses import dataclass

DNA_BASES = "ATGC"
RNA_BASES = "AUGC"

# characters that would make NUPACK complex names (e.g. "a+b") ambiguous
FORBIDDEN_NAME_CHARS = "+()<>"

UM_TO_M = 1e-6  # µM -> M (NUPACK works in molar)
M_TO_UM = 1e6  # M -> µM (for plotting)
MIN_CONCENTRATION_M = 1e-8  # complexes below this equilibrium concentration are not shown

@dataclass(frozen=True)
class StrandInput:
    """One strand entered by the user."""
    name: str
    sequence: str
    conc_um: float

def valid_bases(material: str) -> str:
    """Return the allowed bases for 'DNA' or 'RNA'."""
    return RNA_BASES if material == "RNA" else DNA_BASES

def has_only_valid_bases(sequence: str, material: str) -> bool:
    """True if every character is a valid base. An empty sequence passes, so a
    row is not flagged red before the user has typed anything."""
    return set(sequence) <= set(valid_bases(material))

def convert_sequence(sequence: str, material: str) -> str:
    """Convert a sequence to DNA (U -> T) or RNA (T -> U)."""
    sequence = sequence.upper()
    if material == "RNA":
        return sequence.replace("T", "U")
    return sequence.replace("U", "T")

def default_strand_name(index: int) -> str:
    """Name given to an unnamed strand in row `index` (0-based)."""
    return f"strand{index + 1}"

def find_input_problems(strands: list[StrandInput], material: str) -> list[str]:
    """Check user input before it is sent to NUPACK.
    Returns a list of problems; an empty list means the input is valid.
    Empty sequences must be caught here, as they crash NUPACK (segmentation fault).
    """
    if not strands:
        return ["Add at least one strand."]
    problems = []
    seen_names = set()
    for strand in strands:
        if not strand.sequence:
            problems.append(f"'{strand.name}' has no sequence.")
        elif not has_only_valid_bases(strand.sequence, material):
            problems.append(f"'{strand.name}' contains bases other than {', '.join(valid_bases(material))}.")
        if any(c in FORBIDDEN_NAME_CHARS for c in strand.name):
            problems.append(f"'{strand.name}': names cannot contain any of {' '.join(FORBIDDEN_NAME_CHARS)}")
        if strand.name in seen_names:
            problems.append(f"The name '{strand.name}' is used more than once.")
        seen_names.add(strand.name)
    return problems

def complex_label(complex_obj: object) -> str:
    """Turn a NUPACK complex (printed as '<Complex (a+b)>') into a clean label ('a+b')."""
    return str(complex_obj).replace("<Complex ", "").replace(">", "").replace("(", "").replace(")", "")

def complex_sort_key(label: str) -> tuple[int, str]:
    """Sort complexes by number of strands, then alphabetically."""
    return label.count("+"), label

def dot_bracket_pairs(structure: str) -> dict[int, int]:
    """Map each paired base index to its partner, from a dot-bracket structure.
    Strand breaks ('+') are removed first, so indices run continuously across strands,
    matching the NUPACK pair-probability matrix and the NucDraw coordinates.
    """
    stack = []
    pairs = {}
    for i, c in enumerate(structure.replace("+", "")):
        if c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            pairs[i] = j
            pairs[j] = i
    return pairs
