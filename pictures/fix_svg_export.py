#!/usr/bin/env python3
"""
Repair an Inkscape SVG so that PDF/EPS export matches the PNG export.

Figures built by importing a PDF into Inkscape pick up two artefacts that
break PDF export while leaving PNG export untouched:

  1. Every imported bitmap gets a redundant fully-white luminance <mask>.
     Cairo turns these into PDF soft masks; Apple's PDF engine (Preview,
     Acrobat, and anything using CoreGraphics) drops the masked image
     entirely, so the photographs vanish.  Ghostscript renders them, which
     is why the bug looks intermittent.
  2. The same photograph is embedded once per placement (27x here), and each
     copy is clipped down to ~14% of its area, so the export carries ~7x more
     image data than it draws.

This script strips (1) and, optionally, crops each placement to its clip
rectangle for (2).  Both are lossless: output renders pixel-identical to the
original (max channel difference of 1 over the whole page).

Usage:
    python3 fix_svg_export.py in.svg out.svg [--crop] [--margin=MM] [--transparent-margin]

Then export with:
    inkscape --export-area-drawing --export-type=pdf --export-filename=out.pdf out.svg
    inkscape --export-area-drawing --export-type=eps --export-filename=out.eps out.svg

Do NOT post-process the PDF with Ghostscript: re-distilling rebuilds the
mix-blend-mode groups and produces hard-edged banding.
"""
import re, sys, io, math, base64
from PIL import Image


def strip_opaque_masks(s):
    """Delete <mask> blocks whose image is entirely white, and their references."""
    opaque, out, pos = set(), [], 0
    for m in re.finditer(r'<mask\b[^<]*?>', s):
        st = m.start()
        en = s.index('</mask>', st) + len('</mask>')
        block = s[st:en]
        mid = re.search(r'id="([^"]+)"', block)
        img = re.search(r'<image\b[^<]*?base64,([^"]+)"', block)
        if mid and img:
            g = Image.open(io.BytesIO(base64.b64decode(img.group(1)))).convert('L')
            if g.getextrema()[0] == 255:
                opaque.add(mid.group(1))
                out.append(s[pos:st])
                pos = en
    out.append(s[pos:])
    s = ''.join(out)
    s = re.sub(r'\s*mask="url\(#([^)]+)\)"',
               lambda mo: '' if mo.group(1) in opaque else mo.group(0), s)
    return s, len(opaque)


def dedupe_images(s):
    """Collapse identical embedded images into one <image> in <defs> + <use> refs."""
    href, out, pos, n = None, [], 0, 0
    for m in re.finditer(r'<image\b[^<>]*?/>', s):
        t = m.group(0)
        h = re.search(r'xlink:href="(data:[^"]+)"', t)
        if not h:
            continue
        if href is None:
            href = h.group(1)
        elif h.group(1) != href:
            return s, 0          # more than one distinct image: leave alone
        keep = [f'{a}="{v.group(1)}"' for a in ('transform', 'clip-path', 'id')
                for v in [re.search(a + r'="([^"]*)"', t)] if v]
        out.append(s[pos:m.start()])
        out.append('<use xlink:href="#sharedPhoto" ' + ' '.join(keep) + ' />')
        pos, n = m.end(), n + 1
    if not n:
        return s, 0
    out.append(s[pos:])
    shared = ('<image id="sharedPhoto" width="1" height="1" preserveAspectRatio="none" '
              f'style="image-rendering:optimizeSpeed" xlink:href="{href}" />')
    return re.sub(r'(<defs\b[^<>]*id="defs1">)', r'\1' + shared.replace('\\', '\\\\'),
                  ''.join(out), count=1), n


def crop_to_clips(s):
    """Crop each <use> of the shared photo down to its own clip rectangle."""
    h = re.search(r'<image id="sharedPhoto"[^>]*xlink:href="data:image/[a-z]+;base64,([^"]+)"', s)
    if not h:
        return s, 0
    photo = Image.open(io.BytesIO(base64.b64decode(h.group(1))))
    if photo.mode == 'RGBA' and photo.getchannel('A').getextrema()[0] == 255:
        photo = photo.convert('RGB')
    W, H = photo.size

    clips = {}
    for m in re.finditer(r'<clipPath\b[^<>]*id="([^"]+)"[^<>]*>(.*?)</clipPath>', s, re.S):
        body = m.group(2)
        if body.count('<rect') == 1 and '<path' not in body:
            try:
                g = lambda k: float(re.search(k + r'="([-0-9.eE]+)"', body).group(1))
                clips[m.group(1)] = (g('x'), g('y'), g('width'), g('height'))
            except AttributeError:
                pass

    cache, newclips, out, pos, done = {}, [], [], 0, 0
    for m in re.finditer(r'<use\b[^<>]*?/>', s):
        t = m.group(0)
        cid = re.search(r'clip-path="url\(#([^)]+)\)"', t)
        tr = re.search(r'transform="matrix\(([^)]*)\)"', t)
        eid = re.search(r'id="([^"]+)"', t)
        if not cid or cid.group(1) not in clips or not tr:
            continue
        a, b, c, d, e, f = [float(v) for v in re.split(r'[,\s]+', tr.group(1).strip())]
        x, y, w, hh = clips[cid.group(1)]
        x0, y0 = max(0.0, x), max(0.0, y)
        x1, y1 = min(1.0, x + w), min(1.0, y + hh)
        if x1 <= x0 or y1 <= y0:
            continue
        box = (math.floor(x0 * W), math.floor(y0 * H),
               min(W, math.ceil(x1 * W)), min(H, math.ceil(y1 * H)))
        u0, v0, u1, v1 = box[0] / W, box[1] / H, box[2] / W, box[3] / H
        su, sv = u1 - u0, v1 - v0
        if box not in cache:
            buf = io.BytesIO()
            photo.crop(box).save(buf, 'PNG', optimize=True)
            cache[box] = 'data:image/png;base64,' + base64.b64encode(buf.getvalue()).decode()
        ncid = f'cropclip{done}'
        newclips.append(
            f'<clipPath clipPathUnits="userSpaceOnUse" id="{ncid}"><rect '
            f'x="{(x-u0)/su:.9f}" y="{(y-v0)/sv:.9f}" '
            f'width="{w/su:.9f}" height="{hh/sv:.9f}" /></clipPath>')
        out.append(s[pos:m.start()])
        out.append(
            f'<image width="1" height="1" preserveAspectRatio="none" '
            f'style="image-rendering:optimizeSpeed" transform="matrix('
            f'{a*su:.9f},{b*su:.9f},{c*sv:.9f},{d*sv:.9f},'
            f'{a*u0+c*v0+e:.9f},{b*u0+d*v0+f:.9f})" '
            f'clip-path="url(#{ncid})" id="{eid.group(1) if eid else "cropimg%d" % done}" '
            f'xlink:href="{cache[box]}" />')
        pos, done = m.end(), done + 1
    if not done:
        return s, 0
    out.append(s[pos:])
    s = ''.join(out)
    s = re.sub(r'<image id="sharedPhoto".*?/>', '', s, count=1, flags=re.S)
    return re.sub(r'(<defs\b[^<>]*id="defs1">)',
                  lambda mo: mo.group(1) + ''.join(newclips), s, count=1), done


def add_backdrop(path, margin_mm, white=True):
    """Insert a rect margin_mm larger than the drawing, behind everything.

    Gives the exported PDF/EPS a margin, since Inkscape's --export-margin is
    ignored when combined with --export-area-drawing.  Coordinates come from
    Inkscape's own bbox query, converted from CSS px to the document's mm.
    """
    import subprocess
    q = subprocess.run(['inkscape', '--query-x', '--query-y',
                        '--query-width', '--query-height', path],
                       capture_output=True, text=True)
    vals = [float(v) for v in q.stdout.split()]
    if len(vals) != 4:
        print('  bbox query failed; skipping backdrop'); return
    x, y, w, h = [v * 25.4 / 96 for v in vals]          # px -> mm
    m = margin_mm
    fill = '#ffffff' if white else 'none'
    op = '1' if white else '0'
    rect = (f'<rect id="exportBackdrop" x="{x-m:.6f}" y="{y-m:.6f}" '
            f'width="{w+2*m:.6f}" height="{h+2*m:.6f}" '
            f'style="fill:{fill};fill-opacity:{op};stroke:none" />')
    s = open(path).read()
    s = re.sub(r'<rect id="exportBackdrop".*?/>', '', s, flags=re.S)   # idempotent
    m2 = re.search(r'<g[^<>]*inkscape:groupmode="layer"[^<>]*>', s)
    if not m2:
        print('  no layer group found; skipping backdrop'); return
    s = s[:m2.end()] + rect + s[m2.end():]
    open(path, 'w').write(s)
    print(f'  added {margin_mm} mm {"white" if white else "transparent"} margin '
          f'({w:.1f} x {h:.1f} mm -> {w+2*m:.1f} x {h+2*m:.1f} mm)')


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if len(args) != 2:
        sys.exit(__doc__)
    src, dst = args
    s = open(src).read()
    s, n = strip_opaque_masks(s);  print(f'removed {n} redundant opaque masks')
    s, n = dedupe_images(s);       print(f'deduplicated {n} copies of the embedded photo')
    if '--crop' in sys.argv:
        s, n = crop_to_clips(s);   print(f'cropped {n} placements to their clip rects')
    open(dst, 'w').write(s)
    mg = [a for a in sys.argv[1:] if a.startswith('--margin')]
    if mg:
        val = float(mg[0].split('=')[1]) if '=' in mg[0] else 3.0
        add_backdrop(dst, val, white='--transparent-margin' not in sys.argv)
    print(f'{src}: {len(open(src).read())/1e6:.1f} MB  ->  {dst}: {len(s)/1e6:.1f} MB')


if __name__ == '__main__':
    main()
