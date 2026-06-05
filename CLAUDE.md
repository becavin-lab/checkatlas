# CLAUDE.md

Guidance for Claude Code (claude.ai/code) when working in this repository.

## IMPORTANT

Never push to main. When asked to add a change or open a PR, check the
current branch first; if it's `main`, create a new branch and push to that.

Do NOT create a `Pipfile` — this project uses `pyproject.toml`.

## Style and conventions

- Python 3.9+. Use f-strings and other modern Python 3 syntax. Do not use
  `__future__` imports or `OrderedDict`.
- Double quotes for strings.
- No shebang lines on Python files unless under `scripts/`.
- Type hints are used extensively; mypy enforces.
- Code formatting is enforced with ruff.
- Tests live in `tests/` and use pytest.
- Helpers for distinct parts of the report — one method per section, one
  for parsing, one for general stats — are encouraged. They make
  `__init__` readable as a high-level outline. What to avoid is trivial
  wrappers: a helper that wraps one or two lines, just renames a
  one-liner, or is called once with no logical separation. Ask whether
  the function name is more meaningful than the code it hides; if not,
  inline it.
- Crash loudly on unexpected data. When parsing output from known
  bioinformatics tools, don't silently default known fields to empty
  dicts or zero values — that hides real format breakage behind a
  fake-looking report. Access documented fields directly
  (`parsed["total_reads"]`), not via `.get(key, 0)`. Catching a parse
  error to raise a friendlier message with the file path is fine;
  silently producing fake data is not. Reserve `.get(default)` for
  genuinely optional fields.
- Never use em-dashes (—) in any user-facing text: module descriptions,
  section titles, plot help text, docstrings, PR descriptions, commit
  messages. Use commas, semicolons, parentheses, or split into two
  sentences.