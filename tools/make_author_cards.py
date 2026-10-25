#!/usr/bin/env python3
"""Regenerate the author cards shown in README.md (assets/authors/*.svg).

Each author is drawn as two self-contained SVG images so each half can carry its
own link: the top half (handle + photo) links to GitHub, the bottom half (name,
ORCID mark, role) links to ORCID. Each half comes in a desktop and a phone size
and a light and a dark theme; README.md picks the right one with <picture>.

Photos are fetched and embedded as circular PNGs, so nothing in the README
depends on an external image service. Re-run whenever a photo changes:

    python3 tools/make_author_cards.py
"""
import base64, html, io, os, subprocess, tempfile, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "assets", "authors")
ORCID_ICON = "https://orcid.org/sites/default/files/images/orcid_16x16.png"
X_PHOTO = "https://unavatar.io/x/chiragjp"          # current X profile photo
X_PHOTO_FALLBACK = "https://pbs.twimg.com/profile_images/706108150315929600/cfsnuY7O_400x400.jpg"

PEOPLE = [  # github handle, display name, phone-size name lines, role, ORCID iD, photo source
    ("terry-b-cho", "Byeongyeon (Terry) Cho", ["Byeongyeon", "(Terry) Cho"], "Lead developer",         "0000-0003-1370-2629", None),
    ("adkostic",    "Aleksandar D. Kostic",   ["Aleksandar D.", "Kostic"],   "Principal investigator", "0000-0002-0837-4360", None),
    ("b-tierney",   "Braden Tierney",         ["Braden", "Tierney"],         "Supervisor",             "0000-0002-7533-8802", None),
    ("chiragjp",    "Chirag J. Patel",        ["Chirag J.", "Patel"],        "Principal investigator", "0000-0002-8756-8525", X_PHOTO),
]
FONT = "-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif"
THEMES = {"light": {"name": "#1f2328", "muted": "#656d76"},
          "dark":  {"name": "#e6edf3", "muted": "#8d96a0"}}
# size name -> card width, photo diameter, handle px, name px, role px
SIZES = {"desktop": dict(w=180, photo=64, handle=9.5, name=12,  role=9.5, orcid=11),
         "phone":   dict(w=82,  photo=44, handle=7.5, name=9,   role=7.5, orcid=9)}

def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "OralMicroNHANES-author-cards"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read(), r.headers.get_content_type()

def fetch_either(url, fallback):
    try:
        return fetch(url)[0]
    except Exception:
        if not fallback: raise
        return fetch(fallback)[0]

def circle_png(raw, size):
    """Crop a photo to a circle with transparent corners (Pillow if present, else ImageMagick)."""
    try:
        from PIL import Image, ImageDraw
        im = Image.open(io.BytesIO(raw)).convert("RGBA"); w, h = im.size; m = min(w, h)
        im = im.crop(((w - m) // 2, (h - m) // 2, (w - m) // 2 + m, (h - m) // 2 + m)).resize((size, size), Image.LANCZOS)
        mask = Image.new("L", (size * 4, size * 4), 0)
        ImageDraw.Draw(mask).ellipse((0, 0, size * 4 - 1, size * 4 - 1), fill=255)
        im.putalpha(mask.resize((size, size), Image.LANCZOS))
        out = io.BytesIO(); im.save(out, "PNG", optimize=True); return out.getvalue()
    except ImportError:
        with tempfile.TemporaryDirectory() as d:
            src, dst = os.path.join(d, "in"), os.path.join(d, "out.png"); open(src, "wb").write(raw)
            subprocess.run(["magick", src, "-resize", f"{size}x{size}^", "-gravity", "center", "-extent", f"{size}x{size}",
                            "(", "+clone", "-threshold", "101%", "-fill", "white", "-draw", f"circle {size//2},{size//2} {size//2},0", ")",
                            "-alpha", "off", "-compose", "copy_opacity", "-composite", dst], check=True)
            return open(dst, "rb").read()

def text_width(t, px):   # average advances of a system UI semibold; textLength pins the drawn width to this
    w = 0.0
    for ch in t:
        if ch in "iljt.() ":              w += px * 0.45
        elif ch.isupper() or ch in "mw":  w += px * 0.75
        else:                             w += px * 0.55
    return round(w)

def svg(w, h, body):
    return f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" font-family="{FONT}">{body}</svg>'

def main():
    os.makedirs(OUT, exist_ok=True)
    orcid_uri = "data:image/png;base64," + base64.b64encode(fetch(ORCID_ICON)[0]).decode()
    for gh, name, name_lines, role, oid, photo in PEOPLE:
        raw = fetch_either(photo or f"https://github.com/{gh}.png?size=256", X_PHOTO_FALLBACK if photo else None)
        for sz, S in SIZES.items():
            W, P = S["w"], S["photo"]
            img = "data:image/png;base64," + base64.b64encode(circle_png(raw, P * 2)).decode()   # 2x for retina
            for theme, c in THEMES.items():
                # top half: handle above the photo
                top_h = 14 + 4 + P + 2
                top = svg(W, top_h,
                    f'<text x="{W/2}" y="12" text-anchor="middle" font-size="{S["handle"]}" fill="{c["muted"]}">@{gh}</text>'
                    f'<image href="{img}" x="{(W-P)/2}" y="18" width="{P}" height="{P}"/>')
                # bottom half: name with the ORCID mark right after it, then the role
                if sz == "desktop":
                    nw = text_width(name, S["name"]); x0 = (W - (nw + 4 + S["orcid"])) / 2
                    body = (f'<text x="{x0:.1f}" y="14" font-size="{S["name"]}" font-weight="600" fill="{c["name"]}" textLength="{nw}" lengthAdjust="spacing">{html.escape(name)}</text>'
                            f'<image href="{orcid_uri}" x="{x0+nw+4:.1f}" y="{14-S["orcid"]+1}" width="{S["orcid"]}" height="{S["orcid"]}"/>'
                            f'<text x="{W/2}" y="32" text-anchor="middle" font-size="{S["role"]}" fill="{c["muted"]}">{role}</text>')
                    bot_h = 38
                else:
                    l1, l2 = name_lines; nw = text_width(l2, S["name"]); x0 = (W - (nw + 3 + S["orcid"])) / 2
                    body = (f'<text x="{W/2}" y="11" text-anchor="middle" font-size="{S["name"]}" font-weight="600" fill="{c["name"]}">{html.escape(l1)}</text>'
                            f'<text x="{x0:.1f}" y="23" font-size="{S["name"]}" font-weight="600" fill="{c["name"]}" textLength="{nw}" lengthAdjust="spacing">{html.escape(l2)}</text>'
                            f'<image href="{orcid_uri}" x="{x0+nw+3:.1f}" y="{23-S["orcid"]+1}" width="{S["orcid"]}" height="{S["orcid"]}"/>')
                    words = role.split(); lines = [role] if text_width(role, S["role"]) <= W - 4 else [words[0], " ".join(words[1:])]
                    for i, ln in enumerate(lines):
                        body += f'<text x="{W/2}" y="{36+i*10}" text-anchor="middle" font-size="{S["role"]}" fill="{c["muted"]}">{ln}</text>'
                    bot_h = 36 + (len(lines) - 1) * 10 + 4
                bot = svg(W, bot_h, body)
                io.open(os.path.join(OUT, f"{gh}-{sz}-{theme}-top.svg"), "w", encoding="utf-8").write(top)
                io.open(os.path.join(OUT, f"{gh}-{sz}-{theme}-bottom.svg"), "w", encoding="utf-8").write(bot)
        print(f"  {gh}")
    print(f"wrote {4*len(PEOPLE)*2} card halves to {OUT}")

if __name__ == "__main__":
    main()
