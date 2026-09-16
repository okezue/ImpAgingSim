"""Rebuild the geometry and sequence artwork; deterministic, explicitly schematic.

SVGs retain editable text. Figma imports use the same geometry plus native text
layers via temporary JSON specifications (pass --figma-json-dir).
"""
from pathlib import Path
import argparse, json, html, math, subprocess
import numpy as np
from scipy.interpolate import make_interp_spline

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'docs/figures'
A='#007F86'; B='#D97943'; INK='#243345'; MUTED='#657385'; PALE='#D9E1E7'; PURPLE='#7556A5'

class SVG:
 def __init__(self,w,h,title):
  self.w=w;self.h=h;self.title=title;self.shapes=[];self.labels=[]
  self.rect(0,0,w,h,'#FFFFFF')
 def rect(self,x,y,w,h,fill,stroke='none',sw=1,rx=0):
  self.shapes.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')
 def path(self,d,color=INK,sw=2,fill='none',dash=None,opacity=1):
  extra=f' stroke-dasharray="{dash}"' if dash else ''
  self.shapes.append(f'<path d="{d}" fill="{fill}" stroke="{color}" stroke-width="{sw}" stroke-linecap="round" stroke-linejoin="round" opacity="{opacity}"{extra}/>')
 def line(self,x1,y1,x2,y2,color=INK,sw=2,dash=None):
  self.path(f'M{x1:.3f},{y1:.3f} L{x2:.3f},{y2:.3f}',color,sw,dash=dash)
 def circle(self,x,y,r,fill,stroke='none',sw=1,opacity=1):
  self.shapes.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{r:.3f}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" opacity="{opacity}"/>')
 def bead(self,x,y,r,v,opacity=1):
  self.circle(x,y,r,A if v else B,'#FFFFFF',1.2,opacity)
 def text(self,x,y,s,size=20,color=INK,weight='Regular',anchor='start'):
  self.labels.append(dict(x=x,y=y,text=s,size=size,color=color,weight=weight,anchor=anchor))
 def arrow(self,x1,y1,x2,y2,color=MUTED,sw=2,head=7):
  self.line(x1,y1,x2,y2,color,sw)
  a=math.atan2(y2-y1,x2-x1)
  p=[(x2-head*math.cos(a-v),y2-head*math.sin(a-v)) for v in [-.45,.45]]
  self.path(f'M{p[0][0]},{p[0][1]} L{x2},{y2} L{p[1][0]},{p[1][1]}',color,sw)
 def panel(self,x,letter):
  self.text(x,39,letter,24,weight='Bold')
 def svg(self,text=True):
  labels=''
  if text:
   for t in self.labels:
    labels+=f'<text x="{t["x"]}" y="{t["y"]}" fill="{t["color"]}" font-family="DejaVu Sans, Arial, sans-serif" font-size="{t["size"]}" font-weight="{700 if t["weight"]=="Bold" else 600 if t["weight"]=="Semi Bold" else 400}" text-anchor="{t["anchor"]}">{html.escape(t["text"])}</text>'
  return f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.w}" height="{self.h}" viewBox="0 0 {self.w} {self.h}" role="img"><title>{html.escape(self.title)}</title>{"".join(self.shapes)}{labels}</svg>'
 def save(self,name,tmp=None):
  OUT.mkdir(parents=True,exist_ok=True);(OUT/f'{name}.svg').write_text(self.svg())
  if tmp:
   tmp.mkdir(parents=True,exist_ok=True)
   (tmp/f'{name}.json').write_text(json.dumps(dict(name=name,width=self.w,height=self.h,svg=self.svg(False),labels=self.labels),default=lambda v:v.item()))

def model(tmp):
 s=SVG(1200,470,'Bead-spring A/B chain and a periodic multichain melt; schematic geometry')
 s.panel(34,'A')
 s.panel(604,'B')
 # Forty beads, equally spaced along a deliberate coiled path (schematic).
 ctrl=np.array([[75,184],[103,109],[180,109],[201,174],[143,222],[157,302],[250,336],[327,300],[289,234],[255,158],[327,117],[423,141],[449,223],[420,298]])
 raw=make_interp_spline(np.linspace(0,1,len(ctrl)),ctrl,k=3)(np.linspace(0,1,2000))
 arclen=np.r_[0,np.cumsum(np.linalg.norm(np.diff(raw,axis=0),axis=1))]
 p=np.column_stack([np.interp(np.linspace(0,arclen[-1],40),arclen,raw[:,d]) for d in (0,1)])
 seq=np.array([1]*8+[0]*7+[1]*4+[0]*8+[1]*6+[0]*7)
 s.path('M'+' L'.join(f'{x:.2f},{y:.2f}' for x,y in p),MUTED,2.6)
 for (x,y),v in zip(p,seq):s.bead(x,y,8.4,v)
 s.text(51,82,'A',18,A,weight='Semi Bold');s.bead(85,76,7,1)
 s.text(123,82,'B',18,B,weight='Semi Bold');s.bead(157,76,7,0)
 # Orthographic box and many short chain excerpts. It is not a trajectory.
 ox,oy=686,159;W=345;H=237;dx=99;dy=-74
 back=[(ox+dx,oy+dy),(ox+dx+W,oy+dy),(ox+dx+W,oy+dy+H),(ox+dx,oy+dy+H)]
 front=[(ox,oy),(ox+W,oy),(ox+W,oy+H),(ox,oy+H)]
 for i in range(4):
  s.line(*back[i],*back[(i+1)%4],PALE,1.5,'5 6')
  s.line(*front[i],*back[i],PALE,1.5,'5 6' if i==3 else None)
 rng=np.random.default_rng(2718)
 # Project mini chains at a range of depths, ensuring all beads are inside box.
 for j in range(17):
  u=rng.uniform(27,W-82);v=rng.uniform(32,H-35);z=rng.uniform(0,1)
  phase=rng.uniform(0,2*np.pi);step=rng.uniform(4,6);t=np.arange(14)
  xs=ox+u+z*dx+step*t;ys=oy+v+z*dy+12*np.sin(t*.43+phase)
  vals=((t+j)%14<7).astype(int) if j%3 else rng.integers(0,2,14)
  s.path('M'+' L'.join(f'{x:.2f},{y:.2f}' for x,y in zip(xs,ys)),MUTED,1.0,opacity=.54)
  for x,y,b in zip(xs,ys,vals):s.bead(x,y,3.8,b,.74)
 for i in range(4):s.line(*front[i],*front[(i+1)%4],INK,1.6)
 # A chain crossing one face, reappearing at the translated opposite face.
 yy=298
 left=[(ox-19,yy+5),(ox+3,yy),(ox+25,yy-7),(ox+45,yy-4)]
 right=[(x+W,y) for x,y in left[:2]]
 for coords in (left,right):
  s.path('M'+' L'.join(f'{x},{y}' for x,y in coords),PURPLE,2.8)
  for x,y in coords:s.circle(x,y,5.8,A,PURPLE,1.5)
 s.line(ox+3,yy,ox+3,423,PURPLE,1,'4 5')
 s.line(ox+W+3,yy,ox+W+3,423,PURPLE,1,'4 5')
 s.arrow(ox+8,423,ox+W-2,423,PURPLE,1.6)
 s.text(858,450,'L',18,PURPLE,anchor='middle')
 s.save('polymer-model',tmp)

def sequence(tmp):
 s=SVG(1200,555,'Bernoulli inheritance gates control correlation amplitude independently of Markov range')
 s.panel(34,'A')
 s.panel(755,'B')
 z=np.array([1,1,1,1,1,0,0,0,0,1,1,1]); mask=np.array([1,1,0,1,1,0,1,0,1,1,0,1]);u=np.array([0,0,0,0,0,1,0,0,0,0,0,0]);out=np.where(mask,z,u)
 xs=224+np.arange(len(z))*40
 s.text(40,133,'Backbone z',19,weight='Semi Bold')
 s.text(40,227,'Mask b',19,weight='Semi Bold')
 s.text(40,331,'Fresh draw u',19,weight='Semi Bold')
 s.text(40,441,'Expressed ψ',19,weight='Semi Bold')
 s.line(xs[0],135,xs[-1],135,PALE,2)
 for x,zz,m,uu,oo in zip(xs,z,mask,u,out):
  s.bead(x,135,11.5,zz)
  s.rect(x-12,218,24,24,'#F0F3F6' if m else '#F0EAF7','none',rx=4)
  s.text(x,237,str(m),17,INK if m else PURPLE,anchor='middle')
  if m:
   s.line(x,149,x,210,PALE,1.4)
   s.arrow(x,252,x,422,PALE,1.4,6)
  else:
   s.line(x,150,x,203,PALE,1.4,'3 5')
   s.line(x-5,192,x+5,202,PURPLE,1.7);s.line(x-5,202,x+5,192,PURPLE,1.7)
   s.bead(x,334,10.5,uu)
   s.circle(x,334,16,'none',PURPLE,1.3)
   s.arrow(x,356,x,422,PURPLE,1.7,6)
 # output bead row has no backbone connects until the final stage
 s.line(xs[0],438,xs[-1],438,MUTED,2)
 for x,oo in zip(xs,out):s.bead(x,438,11.5,oo)
 # Discrete correlation panel. No misleading line connects ell=0 to ell=1.
 x0,y0,w,h=817,400,321,266
 def xy(l,g):return x0+w*l/20,y0-h*g
 s.line(x0,y0,x0+w+10,y0,INK,1.4);s.line(x0,y0,x0,116,INK,1.4)
 for gg in [0,.5,1]:
  x,y=xy(0,gg)
  if gg>0:s.line(x0,y,x0+w,y,'#ECF0F3',1)
  s.text(x0-16,y+6,f'{gg:g}',17,MUTED,anchor='end')
 for l in [0,10,20]:
  x,y=xy(l,0);s.line(x,y,x,y+5,INK,1.2);s.text(x,y+27,str(l),17,MUTED,anchor='middle')
 s.text(x0,87,'Γ(ℓ)',18)
 s.text(x0+w/2,457,'Lag ℓ',18,anchor='middle')
 curves=[(1,.97,INK,None),(.5,.97,PURPLE,None),(1,.85,MUTED,'6 6')]
 for k,pi,col,dash in curves:
  xx=np.arange(1,21);yy=k*k*(2*pi-1)**xx
  s.path('M'+' L'.join(f'{xy(l,g)[0]},{xy(l,g)[1]}' for l,g in zip(xx,yy)),col,2.4,dash=dash)
  for l,g in zip(xx,yy):s.circle(*xy(l,g),2.7,col)
 s.circle(*xy(0,1),5.8,'white',INK,2)
 # concise line legends, no information boxes
 for j,(label,col,dash) in enumerate([('κ = 1,  π = 0.97',INK,None),('κ = ½, π = 0.97',PURPLE,None),('κ = 1,  π = 0.85',MUTED,'6 6')]):
  y=493+j*31;s.line(797,y-5,829,y-5,col,2.4,dash);s.text(844,y,label,18)
 # Reclaim the space formerly occupied by panel headings.
 s.shapes=[s.shapes[0],'<g transform="translate(0,-35)">',*s.shapes[1:],'</g>']
 for label in s.labels:
  if label['size']!=24: label['y']-=35
 s.save('sequence-construction',tmp)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--figma-json-dir',type=Path)
 ap.add_argument('--export',action='store_true',help='Use Inkscape to outline text and export PNGs')
 args=ap.parse_args()
 model(args.figma_json_dir);sequence(args.figma_json_dir)

 if args.export:
  for name in ['polymer-model','sequence-construction']:
   src=OUT/f'{name}.svg'; outlined=OUT/f'{name}.outlined.svg'
   subprocess.run(['inkscape',str(src),'--export-text-to-path','--export-plain-svg',f'--export-filename={outlined}'],check=True,capture_output=True)
   outlined.replace(src)
   subprocess.run(['inkscape',str(src),f'--export-filename={src.with_suffix(".png")}','--export-width=2160'],check=True,capture_output=True)
