"""Log-formatting helpers for the atlas-level pipeline.

This module is the *only* place that knows how to render the
multi-line banners and the per-atlas timing summary table.  Every
other module in the pipeline calls :func:`log_box` and
:func:`print_atlas_timing_summary` so the visual style stays
consistent.

Design notes
------------

* The atlas-level pipeline emits a small set of *recurring* events
  (per-atlas preprocess header, per-task banner, per-task footer,
  per-atlas summary table).  Each event is rendered with a box so
  it is visually separated from the other INFO lines the user
  already gets from the existing ``|--- %(levelname)-8s %(message)s``
  format.  The single line that carries the ``|--- INFO`` prefix is
  the first line of the box; the rest is emitted as a continuation
  of the same log record so the whole box is a single block in a
  file / Nextflow log.

* The box uses Unicode box-drawing characters when
  :func:`sys.stdout.isatty` is True and falls back to pure ASCII
  (``+ - |``) otherwise.  This keeps file logs and TTY logs both
  readable and avoids `cat`-ing the log into a terminal that
  cannot render the Unicode glyphs.

* Time deltas are rendered with :func:`format_duration`, which
  always returns ``"Ns"``, ``"Nm Ms"`` or ``"Nh Mm Ss"`` — never a
  floating-point number.  This is friendlier for the per-atlas
  timing table the user wants.

Timings convention
------------------

The atlas-level timing summary is built from a flat dict the four
``create_metric_*`` functions and the four ``_precompute_*`` helpers
populate under a single private convention:

    args._checkatlas_timings = {
        "preprocess:cache_hit": float,        # seconds
        "preprocess:dimred":     float,
        "preprocess:annot":      float,
        "preprocess:batch_correction": float,
        "preprocess:cluster":    float,
        "metric:cluster":        (float, int),     # (seconds, n_rows)
        "metric:annot":          (float, int),
        "metric:batch_correction": (float, int),
        "metric:dimred":         (float, int),
    }

The ``__main__`` driver prints :func:`print_atlas_timing_summary`
once after the four ``create_metric_*`` calls.  Empty / missing
keys are skipped in the table so the same summary can be reused
when the user disables one or more tasks via ``--metric_* none``.
"""

from __future__ import annotations

import logging
import sys
from typing import Iterable, Mapping


# ──────────────────────────────────────────────────────────────────
# Terminal capability detection
# ──────────────────────────────────────────────────────────────────


def _box_chars() -> dict[str, str]:
    """Return the box-drawing characters for the current stdout.

    Falls back to pure ASCII when stdout is not a TTY (e.g. when
    checkatlas is invoked under Nextflow or piped to a file).
    """
    if sys.stdout.isatty():
        return {
            "tl": "╭", "tr": "╮", "bl": "╰", "br": "╯",
            "h": "─", "v": "│",
            "tee_l": "├", "tee_r": "┤", "tee_t": "┬", "tee_b": "┴",
            "cross": "┼",
        }
    return {
        "tl": "+", "tr": "+", "bl": "+", "br": "+",
        "h": "-", "v": "|",
        "tee_l": "+", "tee_r": "+", "tee_t": "+", "tee_b": "+",
        "cross": "+",
    }


# ──────────────────────────────────────────────────────────────────
# Duration formatting
# ──────────────────────────────────────────────────────────────────


def format_duration(seconds: float) -> str:
    """Render a wall-clock duration as a compact human string.

    Examples
    --------
    >>> format_duration(5.4)
    '5s'
    >>> format_duration(83.0)
    '1m 23s'
    >>> format_duration(3725.0)
    '1h 02m 05s'
    """
    seconds = max(0, int(round(float(seconds))))
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours}h {minutes:02d}m {secs:02d}s"
    if minutes:
        return f"{minutes}m {secs:02d}s"
    return f"{secs}s"


# ──────────────────────────────────────────────────────────────────
# Box rendering
# ──────────────────────────────────────────────────────────────────


def _wrap_box(title: str, body_lines: Iterable[str]) -> list[str]:
    """Build the lines of a single box, with no log prefix."""
    bc = _box_chars()
    title_str = f" {title} " if title else ""
    body = list(body_lines)
    width = max(
        [len(title_str) + 4]
        + [len(line) for line in body]
    )
    # Border: "╭─ TITLE ─────...─╮"
    border_chars = max(1, width - len(title_str) - 2)
    top = bc["tl"] + bc["h"] + title_str + bc["h"] * border_chars + bc["tr"]
    bottom = bc["bl"] + bc["h"] * (width - 2) + bc["br"]
    lines = [top]
    for line in body:
        # Pad to width-2 (between the two verticals)
        lines.append(bc["v"] + " " + line.ljust(width - 4) + " " + bc["v"])
    lines.append(bottom)
    return lines


def log_box(
    logger: logging.Logger,
    level: int,
    title: str,
    body_lines: Iterable[str],
) -> None:
    """Emit a box-style log record at the given level.

    The first line carries the standard ``|--- LEVEL ...`` prefix
    so the box is grep-able like any other log line.  Subsequent
    lines are emitted without the prefix (continuation of the same
    log record) so the box renders as one logical block.
    """
    lines = _wrap_box(title, body_lines)
    # Join into a single multi-line message; Python's logging module
    # preserves newlines when emitting the formatted record, so the
    # log handler renders the box as one block.
    logger.log(level, "\n".join(lines))


# ──────────────────────────────────────────────────────────────────
# Preprocess / per-task headers
# ──────────────────────────────────────────────────────────────────


# Stage labels used in the per-atlas summary table.  Kept short so
# the bordered table stays readable in a 100-column terminal.
STAGE_LABELS = {
    "preprocess:cache_hit":       "preprocess:cache_hit",
    "preprocess:dimred":          "preprocess:dimred",
    "preprocess:annot":           "preprocess:annot",
    "preprocess:batch_correction": "preprocess:batch_correction",
    "preprocess:cluster":          "preprocess:cluster",
    "metric:cluster":              "metric:cluster",
    "metric:annot":                "metric:annot",
    "metric:batch_correction":     "metric:batch_correction",
    "metric:dimred":               "metric:dimred",
}


def preprocess_header(
    logger: logging.Logger,
    atlas_name: str,
    n_obs: int,
    n_vars: int,
    *,
    run_dimred: bool,
    run_annot: bool,
    run_batch: bool,
    run_cluster: bool,
    cache_hit: bool,
) -> None:
    """Emit the per-atlas preprocess banner.

    Lists the four precompute tasks with a checkbox-style ``[x]`` /
    ``[ ]`` marker so the user can immediately see which
    precomputations are about to run (or whether the cache was hit
    and the per-task precompute block is being skipped).
    """
    def mark(b: bool) -> str:
        return "[x]" if b else "[ ]"

    body = [
        f"atlas: {atlas_name}   ({n_obs:,} cells × {n_vars:,} genes)",
        "precompute plan:",
        (
            f"  {mark(run_dimred)} dimred            "
            f"{mark(run_annot)} annot            "
            f"{mark(run_batch)} batch_correction "
            f"{mark(run_cluster)} cluster"
        ),
        (
            f"cache: {'HIT (skipping per-task precompute)'
                     if cache_hit else 'MISS (building per-task artefacts now)'}"
        ),
    ]
    log_box(logger, logging.INFO, f"checkatlas [preprocess] ─ {atlas_name}", body)


def task_header(
    logger: logging.Logger,
    *,
    atlas_name: str,
    n_obs: int,
    n_vars: int,
    task: str,         # "cluster" | "annot" | "batch_correction" | "dimred"
    keys: Mapping[str, list[str]],
) -> None:
    """Emit the per-task banner before the metric engine runs.

    ``keys`` maps a short label (e.g. ``"cluster labels"``,
    ``"ref keys"``, ``"batch keys"``, ``"embeddings"``) to the list
    of column names the column detector picked for this task.
    """
    body = [
        f"atlas: {atlas_name}   ({n_obs:,} cells × {n_vars:,} genes)",
    ]
    for label, names in keys.items():
        if not names:
            body.append(f"  {label}: (none detected)")
        else:
            body.append(f"  {label}: {', '.join(names)}")
    log_box(logger, logging.INFO, f"checkatlas [metric:{task}] ─ {atlas_name}", body)


def task_footer(
    logger: logging.Logger,
    *,
    task: str,
    atlas_name: str,
    elapsed: float,
    n_rows: int | None,
    output_path: str | None,
) -> None:
    """Emit the per-task footer after the metric engine returns."""
    bits = [f"elapsed: {format_duration(elapsed)}"]
    if n_rows is not None:
        bits.append(f"rows: {n_rows:,}")
    if output_path is not None:
        bits.append(f"saved: {output_path}")
    body = [f"{atlas_name}   |   " + "   |   ".join(bits)]
    log_box(logger, logging.INFO, f"checkatlas [metric:{task}] done", body)


def preprocess_footer(
    logger: logging.Logger,
    *,
    atlas_name: str,
    timings: Mapping[str, float],
    cache_hit: bool,
) -> None:
    """Emit the per-atlas preprocess footer after the precompute block."""
    lines = [f"{atlas_name}"]
    if cache_hit:
        lines.append(
            f"  cache hit: precompute skipped "
            f"({format_duration(timings.get('preprocess:cache_hit', 0.0))})"
        )
    else:
        for key in (
            "preprocess:dimred",
            "preprocess:annot",
            "preprocess:batch_correction",
            "preprocess:cluster",
        ):
            if key in timings:
                lines.append(
                    f"  {key.replace('preprocess:', ''):<19} "
                    f"{format_duration(timings[key])}"
                )
    log_box(logger, logging.INFO, "checkatlas [preprocess] done", lines)


# ──────────────────────────────────────────────────────────────────
# Per-atlas timing summary table
# ──────────────────────────────────────────────────────────────────


def print_atlas_timing_summary(
    logger: logging.Logger,
    atlas_name: str,
    timings: Mapping[str, object],
) -> None:
    """Print the bordered per-atlas summary table.

    ``timings`` is the per-atlas dict described in the module
    docstring.  The summary lists each stage (precompute +
    metric) in execution order with its elapsed time and (for the
    metric stages) the row count of the long-format result.
    """
    # Collect rows in execution order, skipping missing stages.
    rows: list[tuple[str, str, str]] = []  # (stage, time, extra)
    total = 0.0
    for stage in STAGE_LABELS:
        if stage not in timings:
            continue
        value = timings[stage]
        if isinstance(value, tuple):
            elapsed, n_rows = value
        else:
            elapsed, n_rows = float(value), None
        total += float(elapsed)
        extra = f"{n_rows:,} rows" if n_rows is not None else "—"
        rows.append((STAGE_LABELS[stage], format_duration(elapsed), extra))
    if not rows:
        return

    bc = _box_chars()
    stage_w = max(len("stage"), max(len(r[0]) for r in rows))
    time_w = max(len("time"), max(len(r[1]) for r in rows))
    extra_w = max(len("rows"), max(len(r[2]) for r in rows))
    width = stage_w + time_w + extra_w + 11  # 3 verticals + 2*3 spaces + 2 borders

    title = f" checkatlas [summary] ─ {atlas_name} "
    border = max(1, width - len(title) - 2)
    lines = [bc["tl"] + bc["h"] + title + bc["h"] * border + bc["tr"]]
    header = (
        f"{bc['v']}  "
        f"{'stage'.ljust(stage_w)}"
        f"   {bc['v']}   "
        f"{'time'.ljust(time_w)}"
        f"   {bc['v']}   "
        f"{'rows'.ljust(extra_w)}"
        f"   {bc['v']}"
    )
    lines.append(header)
    sep = (
        f"{bc['v']}  {bc['h'] * stage_w}   {bc['tee_l']}"
        f"   {bc['h'] * time_w}   {bc['tee_r']}"
        f"   {bc['h'] * extra_w}   {bc['v']}"
    )
    lines.append(sep)
    for stage, t, extra in rows:
        lines.append(
            f"{bc['v']}  {stage.ljust(stage_w)}"
            f"   {bc['v']}   {t.ljust(time_w)}"
            f"   {bc['v']}   {extra.ljust(extra_w)}"
            f"   {bc['v']}"
        )
    total_label = "TOTAL"
    total_time = format_duration(total)
    lines.append(sep)
    lines.append(
        f"{bc['v']}  {total_label.ljust(stage_w)}"
        f"   {bc['v']}   {total_time.ljust(time_w)}"
        f"   {bc['v']}   {''.ljust(extra_w)}"
        f"   {bc['v']}"
    )
    lines.append(bc["bl"] + bc["h"] * (width - 2) + bc["br"])

    logger.info("\n".join(lines))
