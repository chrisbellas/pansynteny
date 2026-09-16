"""Render docs/pansynteny-banner.png -- the README banner.

The banner is a screenshot of the viewer with the wordmark drawn into the
empty region on the right, as a miniature synteny plot of the word itself:

    Pansynteny      s, the middle y and t are tracked as genes; the a holds
    Panysnteny      an anchor column down every row, white like the anchor
    Pannsteny       gene in the viewer. Row three has no y at all, so it
    Panyntseny      closes up, ends a letter short, and a blue skip-line
                    carries across to the copy two rows below -- the same
                    device the viewer uses for a gene absent in one genome.

Dropping the y also shifts t one slot left, which is the point: a deletion
moves everything downstream of it.

Usage:  python3 docs/make_banner.py [source.png]

Needs Pillow, and the source screenshot (default: ./pansynteny.png, a
1470x904 capture of the viewer in dark mode). Re-run after replacing the
screenshot; the crop and the blank region below may need adjusting to suit.
"""
from PIL import Image, ImageDraw, ImageFont
import sys

SRC  = sys.argv[1] if len(sys.argv) > 1 else "pansynteny.png"
CROP = (0, 0, 1470, 600)
RX0, RX1, RY0, RY1 = 850, 1285, 88, 458   # the empty region of the screenshot
CX = (RX0 + RX1) // 2

RED, TEAL, YEL = (255,103,112), (107,219,200), (210,204,75)
BLUE = (114,175,211)
WHITE, GREY = (255,255,255), (122,122,122)

WORD = list("Pansynteny")          # P0 a1 n2 s3 y4 n5 t6 e7 n8 y9
HL  = {3: RED, 4: TEAL, 6: YEL}    # s, the middle y, t
ANCHOR = 1                         # the 'a' -- same slot in every row, like the anchor gene

# Row 3 has no middle y at all: the letters after it close up, so the two n's
# become neighbours and t is carried one slot left.
ROWS = [
    [0,1,2,3,4,5,6,7,8,9],   # Pansynteny
    [0,1,2,4,3,5,6,7,8,9],   # Panysnteny
    [0,1,2,5,3,6,7,8,9],     # Pannsteny    <- y absent, nothing inserted
    [0,1,2,4,5,6,3,7,8,9],   # Panyntseny
]
MONO = "/usr/share/fonts/truetype/dejavu/DejaVuSansMono-Bold.ttf"
SANS = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"


def build(dst, alpha=95, size=44, lead=1.55, plain_first_row=True):
    im = Image.open(SRC).convert("RGB").crop(CROP)
    font, sub = ImageFont.truetype(MONO, size), ImageFont.truetype(SANS, 16)
    pitch, lh = font.getlength("M"), int(size * lead)
    ox = CX - pitch * len(WORD) / 2          # every row starts at the same left edge
    block_h = lh * len(ROWS)
    oy = RY0 + (RY1 - RY0 - block_h - 34) / 2
    xy = lambda r, s: (ox + s * pitch, oy + r * lh)

    ov = Image.new("RGBA", im.size, (0,0,0,0)); od = ImageDraw.Draw(ov)
    for gid, col in HL.items():
        rows_with = [r for r, row in enumerate(ROWS) if gid in row]
        for a, b in zip(rows_with, rows_with[1:]):
            x0, y0 = xy(a, ROWS[a].index(gid))
            x1, y1 = xy(b, ROWS[b].index(gid))
            if b - a == 1:
                od.polygon([(x0+2, y0+size*0.95), (x0+pitch-2, y0+size*0.95),
                            (x1+pitch-2, y1+size*0.10), (x1+2, y1+size*0.10)],
                           fill=col + (alpha,))
            else:                                    # absent in between: skip-line
                od.line([(x0+pitch/2, y0+size*0.95), (x1+pitch/2, y1+size*0.10)],
                        fill=BLUE + (150,), width=3)
    im = Image.alpha_composite(im.convert("RGBA"), ov).convert("RGB")

    d = ImageDraw.Draw(im)
    for r, row in enumerate(ROWS):
        for slot, gid in enumerate(row):
            if gid in HL:            fill = HL[gid]
            elif gid == ANCHOR:      fill = WHITE          # the anchor, every row
            elif r == 0 and plain_first_row: fill = WHITE
            else:                    fill = GREY
            d.text(xy(r, slot), WORD[gid], font=font, fill=fill)

    tag = "Gene synteny across a Panaroo pangenome"
    x0, y0, x1, y1 = d.textbbox((0,0), tag, font=sub)
    d.text((CX - (x1-x0)/2 - x0, oy + block_h + 14), tag, font=sub, fill=(138,138,138))
    im.save(dst, optimize=True); print(dst)


build("docs/pansynteny-banner.png")
