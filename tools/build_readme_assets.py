#!/usr/bin/env python3
from pathlib import Path
from math import sin, cos, pi
import argparse
import html
import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageFilter

WIDTH, HEIGHT, SCALE = 1280, 390, 2
INK = (244, 245, 255)
MUTED = (166, 179, 206)
VIOLET = (173, 147, 255)
PINK = (240, 158, 230)
ORANGE = (255, 182, 91)


def background():
    y, x = np.mgrid[0:HEIGHT*SCALE, 0:WIDTH*SCALE]
    xx, yy = x / SCALE, y / SCALE
    blend = np.clip(xx / WIDTH, 0, 1)[..., None]
    left = np.array([8., 18., 35.])
    right = np.array([24., 20., 48.])
    rgb = left * (1-blend) + right * blend
    glow = np.exp(-(((xx-1050)/250)**2 + ((yy-195)/200)**2))
    rgb += glow[..., None] * np.array([14., 6., 27.])
    return Image.fromarray(np.clip(rgb, 0, 255).astype(np.uint8), "RGB").convert("RGBA")


def font(directory, size, bold=False):
    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    return ImageFont.truetype(str(directory / name), round(size*SCALE))


def label(draw, position, text, face, colour, tracking=0):
    x, y = position
    if not tracking:
        draw.text((round(x*SCALE), round(y*SCALE)), text, font=face, fill=colour)
        return
    for letter in text:
        draw.text((round(x*SCALE), round(y*SCALE)), letter, font=face, fill=colour)
        x += draw.textlength(letter, font=face)/SCALE + tracking


def line(draw, points, colour, width):
    draw.line([(round(x*SCALE), round(y*SCALE)) for x,y in points],
              fill=colour, width=round(width*SCALE), joint="curve")


def dot(draw, x, y, radius, colour, outline=None):
    draw.ellipse(((x-radius)*SCALE, (y-radius)*SCALE,
                  (x+radius)*SCALE, (y+radius)*SCALE),
                 fill=colour, outline=outline, width=SCALE)


def frame(base, fonts, t):
    image = base.copy()
    d = ImageDraw.Draw(image)
    # Permanent content: every frame, including a paused GIF, is readable.
    label(d,(58,44),"NNQQ  /  MOLECULAR SIMULATION",fonts["small"],MUTED,1.4)
    label(d,(54,94),"From fibril",fonts["title"],INK)
    label(d,(54,166),"to solution.",fonts["title"],INK)
    label(d,(58,263),"Resolving the forces behind amyloid stability.",fonts["body"],MUTED)
    line(d,[(59,315),(87,315)],VIOLET,3)
    label(d,(100,302),"SMD",fonts["tag"],INK)
    label(d,(165,302),"/",fonts["tag"],(91,107,139))
    label(d,(189,302),"FREE ENERGY",fonts["tag"],INK)
    label(d,(345,302),"/",fonts["tag"],(91,107,139))
    label(d,(369,302),"SAPT",fonts["tag"],INK)
    # Faint background curves create depth without encoding any measured data.
    layer = Image.new("RGBA", image.size)
    ld = ImageDraw.Draw(layer)
    for i in range(5):
        points=[(760+j*5, 75+i*57+15*sin(j*.07+i*.5)) for j in range(95)]
        line(ld,points,(154,131,221,15),.7)
    image = Image.alpha_composite(image,layer)
    d=ImageDraw.Draw(image)
    # Abstract, stacked ribbons: intentionally not atom-resolved geometry.
    for sheet in range(2):
        for row in range(6):
            points=[]
            for j in range(64):
                s=j/63
                x=790+sheet*56+s*202
                y=104+row*31+sheet*14+7*sin(2*pi*s+row*.3)+1.2*sin(t+row*.35)
                points.append((x,y))
            colour=VIOLET if sheet else (112,126,176)
            line(d,points,tuple(int(c*.38) for c in colour),9)
            line(d,points,colour,3.2)
            for j in (4,23,42,60):
                x,y=points[j]
                dot(d,x,y,3.7,colour)
                dot(d,x-.7,y-.8,.9,(232,231,255))
    # A separate floating motif suggests solution without depicting an MD path.
    move=6*sin(t)
    points=[(1095+j*2+move,205+8*sin(j*.2)+2*cos(t)) for j in range(43)]
    line(d,points,(88,62,112),9)
    line(d,points,PINK,3.2)
    for j in (2,15,28,40):
        x,y=points[j]
        dot(d,x,y,4,PINK)
    for i,(x,y,r,c) in enumerate([
        (1083,79,6,VIOLET),(1188,112,8,ORANGE),(1124,145,5,PINK),
        (1200,240,6,VIOLET),(1100,304,7,ORANGE),(1173,323,5,PINK),
        (760,275,4,VIOLET),(970,65,4,PINK),
    ]):
        x+=5*cos(t+i*.9);y+=4*sin(t+i*.9)
        dot(d,x,y,r+5,tuple(int(v*.19) for v in c))
        dot(d,x,y,r,c)
        dot(d,x-r*.23,y-r*.25,r*.23,(248,244,255))
    label(d,(818,354),"ILLUSTRATION  /  NOT SIMULATION DATA",fonts["note"],(132,139,168),.35)
    return image.convert("RGB").resize((WIDTH,HEIGHT),Image.Resampling.LANCZOS)


def badge(path, left, right, accent, width):
    split=round(width*.47)
    # Plain SVG; no scripts, foreignObject, CSS, animation, or external resources.
    svg=f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="28" role="img" aria-label="{html.escape(left+' '+right)}">
<title>{html.escape(left+' '+right)}</title>
<rect x=".5" y=".5" width="{width-1}" height="27" rx="6" fill="#131d30" stroke="#37435b"/>
<path d="M {split} 6 V 22" stroke="#37435b"/>
<g font-family="DejaVu Sans,Arial,sans-serif" font-size="11" text-anchor="middle">
<text x="{split/2}" y="18" fill="#d6deee">{html.escape(left)}</text>
<text x="{(split+width)/2}" y="18" fill="{accent}" font-weight="bold">{html.escape(right)}</text>
</g></svg>'''
    path.write_text(svg,encoding="utf-8")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--font-dir",type=Path,default=Path("/usr/share/fonts/truetype/dejavu"))
    args=parser.parse_args()
    out=Path(__file__).resolve().parents[1]/"docs"/"assets"
    out.mkdir(parents=True,exist_ok=True)
    fonts={"small":font(args.font_dir,12,True),"title":font(args.font_dir,64,True),
           "body":font(args.font_dir,20),"tag":font(args.font_dir,15,True),
           "note":font(args.font_dir,10)}
    base=background()
    frames=[frame(base,fonts,2*pi*i/60) for i in range(60)]
    frames[0].save(out/"amyloid-static.png",optimize=True)
    # One palette for every frame prevents colour flicker during looping.
    contact=Image.new("RGB",(WIDTH,HEIGHT*4))
    for i,j in enumerate((0,15,30,45)):
        contact.paste(frames[j],(0,HEIGHT*i))
    palette=contact.quantize(colors=256,method=Image.Quantize.MEDIANCUT)
    indexed=[f.quantize(palette=palette,dither=Image.Dither.NONE) for f in frames]
    indexed[0].save(out/"amyloid-motion.gif",save_all=True,append_images=indexed[1:],
                    duration=100,loop=0,optimize=True,disposal=1)
    badge(out/"badge-chemrxiv.svg","ChemRxiv","PREPRINT","#c9b8ff",166)
    badge(out/"badge-gromacs.svg","GROMACS","2022.2","#b9d4ff",154)
    badge(out/"badge-python.svg","Python","3.11.8","#ffd49f",138)
    badge(out/"badge-nnqq.svg","Peptide","NNQQ","#f4b9ec",130)
    print(f"Generated {len(frames)} frames; GIF: {(out/'amyloid-motion.gif').stat().st_size/1024:.0f} KiB")


if __name__=="__main__":
    main()
