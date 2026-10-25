#!/usr/bin/env python3
# =============================================================================
# Module 12 Phase 3 — TSV -> formatted .xlsx packager
#
# Purpose : Assemble the supplementary workbook from sheets that R has already
#           built and asserted. This script does NO analysis and makes NO
#           decisions about content: it packages, types and formats.
#
# Usage   : tsv_to_xlsx.py <out.xlsx> <sheet name>=<file.tsv> [ ... ]
#           Sheet order on screen follows argument order.
#
# WHY PYTHON: the project conda R library has no Excel writer (openxlsx,
# openxlsx2, writexl and xlsx are all absent) and cannot gain one - that env's
# lib/R/etc/Makeconf has no CC, no CXX and no SHLIB_LIBADD, so every source
# package with compiled code fails to build. The system Python already carries
# openpyxl, which needs no compilation, so this route adds no dependency at all.
# =============================================================================
import csv
import sys

try:
    import openpyxl
    from openpyxl.styles import Alignment, Font
    from openpyxl.utils import get_column_letter
except ImportError:  # pragma: no cover - environment guard
    sys.exit("[FATAL] openpyxl not importable by %s" % sys.executable)

# Excel refuses to open a file with a cell longer than this, with an error that
# does not name the offending cell. Truncate defensively and say so in the cell.
CELL_MAX = 32767
SCI = "0.00E+00"          # p-values, FDR, q-values
DEC4 = "0.0000"           # effect sizes on the CLR scale
WIDTH_MIN, WIDTH_MAX = 8, 52


def is_sci_col(name):
    n = name.lower()
    return any(k in n for k in ("p-value", " p", "fdr", "q-value", "_p", "p (")) and "participants" not in n


def is_beta_col(name):
    n = name.lower()
    return any(k in n for k in ("beta", "se", "rho", "pearson", "prevalence", "odds"))


def typed(value):
    """Return a number where the source text is unambiguously numeric.

    Excel treats a numeric-looking string as text and quietly refuses to sort
    or filter it correctly, so this conversion is what makes the sheet usable
    rather than merely present. Empty and NA become a true blank cell.
    """
    if value is None:
        return None
    v = value.strip()
    if v in ("", "NA", "NaN", "nan", "-", "—"):
        return None
    if v in ("TRUE", "FALSE"):
        return v == "TRUE"
    try:
        return int(v)
    except ValueError:
        pass
    try:
        return float(v)
    except ValueError:
        return v[:CELL_MAX - 3] + "..." if len(v) > CELL_MAX else v


def add_sheet(wb, title, path, first):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))
    if not rows:
        sys.exit("[FATAL] empty sheet file: %s" % path)
    header, body = rows[0], rows[1:]

    ws = wb.active if first else wb.create_sheet()
    ws.title = title
    ws.append(header)
    for r in body:
        # A short row would silently shift values left; pad instead.
        if len(r) < len(header):
            r = r + [""] * (len(header) - len(r))
        elif len(r) > len(header):
            sys.exit("[FATAL] %s: row has %d fields, header has %d"
                     % (path, len(r), len(header)))
        ws.append([typed(c) for c in r])

    bold = Font(bold=True)
    for c in ws[1]:
        c.font = bold
        c.alignment = Alignment(vertical="top", wrap_text=True)

    for j, name in enumerate(header, start=1):
        letter = get_column_letter(j)
        if is_sci_col(name):
            fmt = SCI
        elif is_beta_col(name):
            fmt = DEC4
        else:
            fmt = None
        if fmt:
            for i in range(2, ws.max_row + 1):
                cell = ws.cell(row=i, column=j)
                if isinstance(cell.value, float):
                    cell.number_format = fmt
        longest = max([len(name)] + [len(str(ws.cell(row=i, column=j).value or ""))
                                     for i in range(2, ws.max_row + 1)])
        ws.column_dimensions[letter].width = max(WIDTH_MIN, min(WIDTH_MAX, longest + 2))

    # Frozen header + filter: on a 362-row sheet the header scrolls away
    # immediately without this, and the column names are the only thing that
    # makes the numbers interpretable.
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = "A1:%s%d" % (get_column_letter(len(header)), ws.max_row)
    return len(body), len(header)


def main(argv):
    if len(argv) < 3:
        sys.exit("usage: tsv_to_xlsx.py <out.xlsx> <sheet>=<tsv> [...]")
    out, specs = argv[1], argv[2:]
    wb = openpyxl.Workbook()
    for k, spec in enumerate(specs):
        if "=" not in spec:
            sys.exit("[FATAL] expected <sheet>=<tsv>, got: %s" % spec)
        title, path = spec.split("=", 1)
        if len(title) > 31:
            sys.exit("[FATAL] sheet name over Excel's 31-char limit: %s" % title)
        nrow, ncol = add_sheet(wb, title, path, first=(k == 0))
        print("  sheet '%s': %d rows x %d columns" % (title, nrow, ncol))
    wb.save(out)
    print("  wrote %s" % out)


if __name__ == "__main__":
    main(sys.argv)
