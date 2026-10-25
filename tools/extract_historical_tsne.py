#!/usr/bin/env python3
"""
Extract tSNE point coordinates from historical 5.4_combined_pca_tsne.pdf.

The PDF content stream is operator-based, but ggplot/cairographics packs
many operators per line (e.g. "Q q 25.62 210.54 61.03 56.33 re W n").
We tokenize the stream as a flat list of PDF operators + numeric/string
operands, then walk it.
"""
import re
import sys
import zlib
import csv
from pathlib import Path

PDF = Path(sys.argv[1] if len(sys.argv) > 1 else
           "/n/groups/patel/terry/nhanes_oral_mirco_cho/results/analyses_results/5_age_analyses_out/5.4_combined_pca_tsne.pdf")
OUT_CSV = Path(sys.argv[2] if len(sys.argv) > 2 else
               "/n/groups/patel/terry/nhanes_oral_mirco_cho/results/intermediate/historical_tsne_coords_k9.csv")

raw = PDF.read_bytes()

# ---- decompress single FlateDecode stream ----
m = re.search(rb"/Filter /FlateDecode[^>]*>>\s*stream\n", raw)
end = raw.find(b"\nendstream", m.end())
stream = zlib.decompress(raw[m.end():end]).decode("latin-1", errors="replace")
print(f"[*] Decompressed stream: {len(stream)} chars")

# ---- tokenize content stream ----
#
# Token types:
#   number   : matches re_num
#   name     : starts with "/" -- we just skip these as operands
#   string   : "(...)" -- preserve content for Tj
#   operator : alphabetic token
#
# Operators we care about:
#   q, Q                      : graphics state push/pop
#   re W n                    : "rectangle, clip, no-op end-path"
#   m, c, l, h, B, S, f, b, n : path construction / paint
#   rg, RG                    : non-stroke / stroke RGB color
#   Tm, Tj, TJ                : text matrix / show text
#   BT, ET                    : text begin / end

tok_re = re.compile(r"""
    \s+                            |
    \(([^)]*)\)                    |   # parenthesized string
    (/[A-Za-z][A-Za-z0-9_\-]*)     |   # name token
    (-?\d+\.\d+|-?\d+)             |   # number
    \[([^\]]*)\]                   |   # array (e.g. TJ argument)
    ([A-Za-z][A-Za-z*0-9]*)            # operator
""", re.VERBOSE)

operators = set("""
q Q re W n m c l h B S f b rg RG cm w J j d M Tm Tj TJ BT ET Tf gs Td TD T* Do
""".split())

tokens = []  # list of (kind, value)
pos = 0
while pos < len(stream):
    m = tok_re.match(stream, pos)
    if not m:
        pos += 1
        continue
    pos = m.end()
    if m.group(1) is not None:        # parenthesized string
        tokens.append(("STR", m.group(1)))
    elif m.group(2) is not None:      # name
        tokens.append(("NAME", m.group(2)))
    elif m.group(3) is not None:      # number
        tokens.append(("NUM", float(m.group(3))))
    elif m.group(4) is not None:      # array (TJ)
        tokens.append(("ARR", m.group(4)))
    elif m.group(5) is not None:      # operator
        tokens.append(("OP", m.group(5)))
    # else whitespace, skip
print(f"[*] Tokenized into {len(tokens)} tokens")

# ---- walk tokens ----
operand_stack = []
clip_stack = []          # list of dicts: {"kind": "q"} or {"kind": "clip", "rect": (x,y,w,h)}
points = []
texts  = []
current_rg = None
current_RG = None
path = []                # current path point list
current_Tm = None        # (tx, ty)
draw_idx = 0
n_circles_inside_clip = 0

for k, v in tokens:
    if k in ("NUM", "STR", "NAME", "ARR"):
        operand_stack.append((k, v))
        continue
    if k != "OP":
        continue
    op = v

    if op == "q":
        clip_stack.append({"kind": "q"})
        operand_stack.clear()
        continue
    if op == "Q":
        while clip_stack and clip_stack[-1]["kind"] != "q":
            clip_stack.pop()
        if clip_stack:
            clip_stack.pop()
        operand_stack.clear()
        continue
    if op == "re":
        # 4 numbers on the stack: x y w h
        if len(operand_stack) >= 4:
            nums = [operand_stack[-4][1], operand_stack[-3][1],
                    operand_stack[-2][1], operand_stack[-1][1]]
            operand_stack[:] = []
            # Stash the rect; next "W n" turns it into a clip.
            pending_rect = tuple(nums)
        else:
            pending_rect = None
        # Store pending_rect alongside the loop variable so W can find it
        re_pending = pending_rect
        # Don't push yet; wait for W
        # We use closures via list to keep state
        # Simplest: just attach to the loop via a sentinel
        _pending_re = pending_rect
        continue
    if op == "W":
        # marks the prior path as the clipping path; together with "n" closes
        # We assume the most recent "re" defined the clip rect, captured in
        # _pending_re. We push onto clip_stack now.
        if "_pending_re" in dir() and _pending_re is not None:
            clip_stack.append({"kind": "clip", "rect": _pending_re})
            _pending_re = None
        operand_stack.clear()
        continue
    if op == "n":
        # end path with no fill or stroke
        path = []
        operand_stack.clear()
        continue
    if op == "rg":
        if len(operand_stack) >= 3:
            current_rg = tuple(round(operand_stack[i][1], 6) for i in (-3, -2, -1))
        operand_stack.clear()
        continue
    if op == "RG":
        if len(operand_stack) >= 3:
            current_RG = tuple(round(operand_stack[i][1], 6) for i in (-3, -2, -1))
        operand_stack.clear()
        continue
    if op == "m":
        if len(operand_stack) >= 2:
            path = [(operand_stack[-2][1], operand_stack[-1][1])]
        operand_stack.clear()
        continue
    if op == "c":
        # 6 operands: x1 y1 x2 y2 x3 y3
        if len(operand_stack) >= 6:
            ops = [operand_stack[i][1] for i in range(-6, 0)]
            path += [(ops[0], ops[1]), (ops[2], ops[3]), (ops[4], ops[5])]
        elif len(operand_stack) >= 4:
            ops = [operand_stack[i][1] for i in range(-4, 0)]
            path += [(ops[0], ops[1]), (ops[2], ops[3])]
        operand_stack.clear()
        continue
    if op == "l":
        if len(operand_stack) >= 2:
            path.append((operand_stack[-2][1], operand_stack[-1][1]))
        operand_stack.clear()
        continue
    if op in ("B", "b", "S", "f", "F", "B*", "b*", "f*"):
        # close+fill+stroke (B) or stroke (S) -- treat as the end of a marker
        if path:
            xs = [p[0] for p in path]
            ys = [p[1] for p in path]
            cx = 0.5 * (min(xs) + max(xs))
            cy = 0.5 * (min(ys) + max(ys))
            radius = 0.5 * max(max(xs) - min(xs), max(ys) - min(ys))
            clip = next((e["rect"] for e in reversed(clip_stack) if e["kind"] == "clip"), None)
            points.append({
                "draw_idx": draw_idx,
                "clip": clip,
                "rgb_fill": current_rg,
                "rgb_stroke": current_RG,
                "cx": cx, "cy": cy,
                "radius": radius,
            })
            draw_idx += 1
            path = []
        operand_stack.clear()
        continue
    if op == "h":
        # close subpath
        operand_stack.clear()
        continue
    if op == "Tm":
        if len(operand_stack) >= 6:
            tx = operand_stack[-2][1]; ty = operand_stack[-1][1]
            current_Tm = (tx, ty)
        operand_stack.clear()
        continue
    if op == "Tj":
        # one string operand
        if operand_stack and operand_stack[-1][0] == "STR" and current_Tm is not None:
            clip = next((e["rect"] for e in reversed(clip_stack) if e["kind"] == "clip"), None)
            texts.append({
                "clip": clip,
                "x": current_Tm[0],
                "y": current_Tm[1],
                "text": operand_stack[-1][1],
            })
        operand_stack.clear()
        continue
    if op == "TJ":
        # array of strings; collect the strings out of the array
        if operand_stack and operand_stack[-1][0] == "ARR" and current_Tm is not None:
            arr = operand_stack[-1][1]
            # extract substrings inside () and any ASCII chars they contain
            txt = "".join(re.findall(r"\(([^)]*)\)", arr))
            clip = next((e["rect"] for e in reversed(clip_stack) if e["kind"] == "clip"), None)
            texts.append({
                "clip": clip,
                "x": current_Tm[0],
                "y": current_Tm[1],
                "text": txt,
            })
        operand_stack.clear()
        continue
    # any other operator just clears the stack
    operand_stack.clear()

print(f"[*] Parsed {len(points)} circle-like paths, {len(texts)} text strings")

# ---- identify panels ----
clip_freq = {}
for p in points:
    if p["clip"] is not None:
        clip_freq[p["clip"]] = clip_freq.get(p["clip"], 0) + 1
panel_clips = sorted(
    [c for c, n in clip_freq.items() if n >= 20 and c[2] > 30 and c[3] > 30],
    key=lambda r: (-r[1], r[0]),
)
print(f"[*] Found {len(panel_clips)} candidate panels:")
for c in panel_clips:
    print(f"    rect={c}  n_points={clip_freq[c]}")

# 2x2 grid: top-left, top-right, bottom-left, bottom-right by (y desc, x asc)
panels_sorted = sorted(panel_clips, key=lambda r: (-r[1], r[0]))
panel_A = panels_sorted[0] if len(panels_sorted) > 0 else None
panel_B = panels_sorted[1] if len(panels_sorted) > 1 else None
panel_C = panels_sorted[2] if len(panels_sorted) > 2 else None
panel_D = panels_sorted[3] if len(panels_sorted) > 3 else None
print(f"[*] panel A (PCA cluster) = {panel_A}")
print(f"[*] panel B (PCA phylum)  = {panel_B}")
print(f"[*] panel C (tSNE cluster)= {panel_C}")
print(f"[*] panel D (tSNE phylum) = {panel_D}")

# ---- find tick label texts adjacent to panel C ----
def parse_num(s):
    s2 = s.replace(" ", "").replace("−", "-").replace("–", "-").strip()
    # cairo sometimes uses an explicit minus-sign glyph rendered as space-prefixed
    if s2.startswith(" "): s2 = s2.lstrip()
    try:
        return float(s2)
    except ValueError:
        return None

def near_x_axis(t, panel):
    if panel is None: return False
    x0, y0, w, h = panel
    # x-tick labels sit just below the panel: y between y0-15 and y0+1
    return (y0 - 15 <= t["y"] <= y0 + 1) and (x0 - 5 <= t["x"] <= x0 + w + 5)

def near_y_axis(t, panel):
    if panel is None: return False
    x0, y0, w, h = panel
    # y-tick labels sit just left of the panel: x between x0-25 and x0
    return (y0 - 5 <= t["y"] <= y0 + h + 5) and (x0 - 25 <= t["x"] < x0 + 1)

c_x_ticks = sorted([t for t in texts if near_x_axis(t, panel_C)], key=lambda t: t["x"])
c_y_ticks = sorted([t for t in texts if near_y_axis(t, panel_C)], key=lambda t: t["y"])
print(f"[*] panel C x-tick candidates: {[(round(t['x'],2), t['text']) for t in c_x_ticks]}")
print(f"[*] panel C y-tick candidates: {[(round(t['y'],2), t['text']) for t in c_y_ticks]}")

x_pairs = [(t["x"], parse_num(t["text"])) for t in c_x_ticks if parse_num(t["text"]) is not None]
y_pairs = [(t["y"], parse_num(t["text"])) for t in c_y_ticks if parse_num(t["text"]) is not None]
print(f"[*] x_pairs: {x_pairs}")
print(f"[*] y_pairs: {y_pairs}")

def linfit(pairs):
    xs = [p[0] for p in pairs]; ys = [p[1] for p in pairs]
    n = len(xs)
    sx = sum(xs); sy = sum(ys)
    sxx = sum(x*x for x in xs); sxy = sum(xs[i]*ys[i] for i in range(n))
    if n*sxx - sx*sx == 0:
        return 1.0, 0.0
    slope = (n*sxy - sx*sy) / (n*sxx - sx*sx)
    intercept = (sy - slope*sx) / n
    return slope, intercept

if not x_pairs or not y_pairs:
    print("[!] Could not auto-calibrate axes; will fall back to panel-bbox heuristic.")
    # ggplot2 default expand on continuous scale: c(0.05, 0.05)
    # so the panel bbox spans (data_min - 0.05*range, data_max + 0.05*range).
    # The historical visible x ticks are -5 to 5 and y ticks are -2 to 4 (rough).
    # We will anchor on this assumption.
    x0, y0, w, h = panel_C
    # Assume the points fill the panel from edge to edge, expand=0.05.
    # Without true tick positions, set x: -5..5, y: -2..4. We'll revisit.
    mx = 10.0 / w; bx = -5.0 - mx * x0
    my = 6.0 / h;  by = -2.0 - my * y0
else:
    mx, bx = linfit(x_pairs)
    my, by = linfit(y_pairs)
print(f"[*] x transform: data = {mx:.6f} * pdf_x + {bx:.6f}")
print(f"[*] y transform: data = {my:.6f} * pdf_y + {by:.6f}")

# ---- extract panel C points in draw order ----
c_pts = sorted([p for p in points if p["clip"] == panel_C], key=lambda p: p["draw_idx"])
print(f"[*] panel C has {len(c_pts)} points")

rows = []
for i, p in enumerate(c_pts):
    dx = mx * p["cx"] + bx
    dy = my * p["cy"] + by
    rows.append({
        "row_index": i + 1,
        "pdf_x":     round(p["cx"], 4),
        "pdf_y":     round(p["cy"], 4),
        "radius_pt": round(p["radius"], 4),
        "rgb_r":     p["rgb_fill"][0] if p["rgb_fill"] else None,
        "rgb_g":     p["rgb_fill"][1] if p["rgb_fill"] else None,
        "rgb_b":     p["rgb_fill"][2] if p["rgb_fill"] else None,
        "tsne_dim1": round(dx, 6),
        "tsne_dim2": round(dy, 6),
    })

OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
with OUT_CSV.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["row_index","pdf_x","pdf_y","radius_pt","rgb_r","rgb_g","rgb_b","tsne_dim1","tsne_dim2"])
    w.writeheader()
    for r in rows:
        w.writerow(r)
print(f"[*] Wrote {OUT_CSV} ({len(rows)} rows)")
if rows:
    print(f"[*] tsne_dim1 range: {min(r['tsne_dim1'] for r in rows):.3f} .. {max(r['tsne_dim1'] for r in rows):.3f}")
    print(f"[*] tsne_dim2 range: {min(r['tsne_dim2'] for r in rows):.3f} .. {max(r['tsne_dim2'] for r in rows):.3f}")
