"""A–F grading for the nine scFM problem verdicts.

A is "no evidence of this problem", F is "the embedding *demonstrates*
the problem in the literature".
"""

from __future__ import annotations



def letter(score: float) -> str:
    """Map a 0..1 problem score to a letter grade.

    Parameters
    ----------
    score : float
        The verdict score, 0 (no problem) .. 1 (severe problem).
        ``None``, ``NaN`` and missing-data markers are mapped to "n/a".

    Returns
    -------
    str
        One of ``{"A", "B", "C", "D", "F", "n/a"}``.
    """
    try:
        s = float(score)
    except (TypeError, ValueError):
        return "n/a"
    if s != s:  # NaN
        return "n/a"
    if s < 0.15:
        return "A"
    if s < 0.35:
        return "B"
    if s < 0.60:
        return "C"
    if s < 0.80:
        return "D"
    return "F"


def grade_legend() -> str:
    """Return a markdown legend for the grade bands."""
    return (
        "# scFM problem grade legend\n\n"
        "| Grade | Score band | Plain language |\n"
        "|-------|------------|----------------|\n"
        "| A | 0.00 - 0.15 | No evidence of this problem |\n"
        "| B | 0.15 - 0.35 | Mild, monitor |\n"
        "| C | 0.35 - 0.60 | Moderate, address before publication |\n"
        "| D | 0.60 - 0.80 | Severe, paper claim likely overstated |\n"
        "| F | 0.80 - 1.00 | Embedding *demonstrates* the problem |\n"
        "| n/a | n/a | Problem not evaluable on this atlas |\n"
    )
