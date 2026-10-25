#!/usr/bin/env python3
"""Regenerate the author cards shown in README.md (assets/authors/*.svg).

Each card is one self-contained SVG (photo embedded) so the README grid can
never misalign or fail to load. Run this whenever a profile photo changes:

    python3 tools/make_author_cards.py
"""
import base64, html, io, os, sys, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "assets", "authors")
ORCID_ICON = "https://orcid.org/sites/default/files/images/orcid_16x16.png"
X_PHOTO = "https://unavatar.io/x/chiragjp"          # current X profile photo
X_PHOTO_FALLBACK = "https://pbs.twimg.com/profile_images/706108150315929600/cfsnuY7O_400x400.jpg"

PEOPLE = [  # github handle, display name, role, ORCID iD, photo source
    ("terry-b-cho", "Byeongyeon (Terry) Cho", "Lead developer",         "0000-0003-1370-2629", None),
    ("adkostic",    "Aleksandar D. Kostic",   "Principal investigator", "0000-0002-0837-4360", None),
    ("b-tierney",   "Braden Tierney",         "Supervisor",             "0000-0002-7533-8802", None),
    ("chiragjp",    "Chirag J. Patel",        "Principal investigator", "0000-0002-8756-8525", X_PHOTO),
]
W, H, R = 180, 140, 32               # card size, photo radius
NAME_PX, SMALL_PX = 12, 9.5
FONT = "-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif"
THEMES = {"light": {"name": "#1f2328", "muted": "#656d76"},
          "dark":  {"name": "#e6edf3", "muted": "#8d96a0"}}

def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "OralMicroNHANES-author-cards"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read(), r.headers.get_content_type()

def data_uri(url, fallback=None):
    try:
        b, t = fetch(url)
    except Exception:
        if not fallback: raise
        b, t = fetch(fallback)
    return f"data:{t};base64,{base64.b64encode(b).decode()}"

def text_width(t, px):   # average advances of a system UI semibold; textLength pins the real width to this
    w = 0.0
    for ch in t:
        if ch in "iljt.() ":              w += px * 0.45
        elif ch.isupper() or ch in "mw":  w += px * 0.75
        else:                             w += px * 0.55
    return round(w)

def main():
    os.makedirs(OUT, exist_ok=True)
    orcid = data_uri(ORCID_ICON)
    for gh, name, role, oid, photo in PEOPLE:
        src = photo or f"https://github.com/{gh}.png?size=128"
        img = data_uri(src, X_PHOTO_FALLBACK if photo else None)
        nw = text_width(name, NAME_PX); block = nw + 4 + 11; x0 = (W - block) / 2
        for theme, c in THEMES.items():
            svg = (
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="{FONT}">'
                f'<defs><clipPath id="c"><circle cx="{W/2}" cy="{14+6+R}" r="{R}"/></clipPath></defs>'
                f'<text x="{W/2}" y="14" text-anchor="middle" font-size="{SMALL_PX}" fill="{c["muted"]}">@{gh}</text>'
                f'<image href="{img}" x="{W/2-R}" y="{14+6}" width="{2*R}" height="{2*R}" clip-path="url(#c)" preserveAspectRatio="xMidYMid slice"/>'
                f'<text x="{x0:.1f}" y="{14+6+2*R+22}" font-size="{NAME_PX}" font-weight="600" fill="{c["name"]}" textLength="{nw}" lengthAdjust="spacing">{html.escape(name)}</text>'
                f'<image href="{orcid}" x="{x0+nw+4:.1f}" y="{14+6+2*R+12}" width="11" height="11"/>'
                f'<text x="{W/2}" y="{14+6+2*R+40}" text-anchor="middle" font-size="{SMALL_PX}" fill="{c["muted"]}">{role}</text>'
                f'</svg>')
            path = os.path.join(OUT, f"{gh}-{theme}.svg")
            io.open(path, "w", encoding="utf-8").write(svg)
        print(f"  {gh:12s} photo from {src.split('?')[0]}")
    print(f"wrote {2*len(PEOPLE)} cards to {OUT}")

if __name__ == "__main__":
    main()
