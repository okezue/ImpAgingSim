"""Exact editable schematic for effective chromatin marks; not simulation snapshots.

BioRender drafts informed the visual organization. All final geometry is explicit
here so state identity, chain connectivity and the local neighborhood can be audited.
"""
from pathlib import Path
import subprocess
from render_schematic_figures import SVG, A, B, INK, MUTED, PURPLE
import numpy as np

OUT=Path(__file__).resolve().parents[1]/'docs/figures'
s=SVG(1200,500,'Chromatin-inspired effective mark dynamics, spatial feedback and density versus site identity; schematic')
for x,p in [(28,'A'),(364,'B'),(906,'C')]:s.panel(x,p)

def nuc(x,y,r,marked):
    color=B if marked else A
    s.circle(x,y,r,color,'#FFFFFF',1.5)
    # The same DNA-wrap motif on both states: changing the mark does not create a nucleosome.
    for k in [-.3,.3]:
        s.path(f'M{x-r*.78},{y+r*(.55+k)} Q{x+r*.05},{y+r*(.8+k)} {x+r*.74},{y+r*(-.5+k)}', '#FFFFFF',max(.8,r*.10),opacity=.68)
    if marked:
        s.path(f'M{x+r*.70},{y-r*.7} Q{x+r*1.15},{y-r*1.35} {x+r*1.43},{y-r*1.14}',INK,max(1,r*.065))
        s.circle(x+r*1.45,y-r*1.13,max(1.5,r*.115),B)

def chain(points,marks,r=10.5,motif=True):
    s.path('M'+' L'.join(f'{x:.2f},{y:.2f}' for x,y in points),MUTED,1.9)
    for (x,y),m in zip(points,marks):
        if motif:nuc(x,y,r,m)
        else:s.bead(x,y,r,not m)

def curved_arrow(d,tip,prev,color=MUTED):
    s.path(d,color,2)
    import math
    angle=math.atan2(tip[1]-prev[1],tip[0]-prev[0]);head=8
    for a in [-.45,.45]:
        s.line(tip[0],tip[1],tip[0]-head*math.cos(angle+a),tip[1]-head*math.sin(angle+a),color,2)

# A: reversible conversion of the same coarse-grained chromatin site.
nuc(96,241,29,False);nuc(266,241,29,True)
s.text(96,314,'Unmarked',18,anchor='middle')
s.text(266,314,'Marked',18,anchor='middle')
curved_arrow('M108,192 C143,148 220,148 252,192',(252,192),(235,167),PURPLE)
curved_arrow('M252,351 C220,390 143,390 108,351',(108,351),(124,371))
s.text(179,147,'Writing',18,PURPLE,anchor='middle')
s.text(179,411,'Turnover',18,MUTED,anchor='middle')

# B: one continuous folded strand. Six marked centers lie within r_c of
# the central unmarked target, including contour-distant segments.
base=np.array([[-150,60],[-130,33],[-100,27],[-65,17],[-28,0],[-65,-35],[-90,-66],[-90,-101],[-60,-128],[-24,-118],[-20,-80],[-14,-24],[0,0],[14,-24],[31,-67],[52,-105],[90,-121],[123,-102],[136,-64],[118,-27],[70,-11],[28,0],[68,30],[99,62],[104,99],[75,130],[43,119],[28,82],[14,24],[-14,24],[-38,61],[-65,98],[-100,126],[-134,114]],float)
center=np.array([617,243]);points=base+center
neighbor=[4,11,13,21,28,29];marks=[i in neighbor for i in range(len(points))]
# Dashed circle is the local interaction neighborhood, not a condensate surface.
s.circle(*center,36,'#F4EFF9',PURPLE,1.2)
# Replace solid circular stroke with an explicit dashed circle.
s.shapes[-1]=s.shapes[-1].replace('stroke-width="1.2"','stroke-width="1.2" stroke-dasharray="4 5"')
for i in neighbor:s.line(*center,*points[i],PURPLE,1.1,'3 4')
chain(points,marks)
s.circle(*center,13.2,'none',PURPLE,2)

# Local writing: isolate the target conversion below the neighborhood.
nuc(574,418,15,False);s.arrow(603,418,640,418,PURPLE,2);nuc(669,418,15,True)

# C: same coarse spatial organization, different site labels. These are
# illustrative matched states, not time samples or a measured trajectory.
mini=np.array([[-61,-6],[-49,-38],[-19,-53],[10,-43],[35,-22],[17,4],[-11,0],[-38,15],[-27,42],[3,57],[34,44],[55,18],[75,-8],[91,13]],float)
first={1,2,3,4,5,6,7,8,9,10};second={0,2,3,4,5,6,7,8,10,11}
for y,state in [(174,first),(362,second)]:
    p=mini+np.array([1028,y]);chain(p,[i in state for i in range(len(p))],8.7,False)
# Sparse rings explicitly identify examples of changed site identity.
for i in sorted(first.symmetric_difference(second)):
    for y in [174,362]:
        x0,y0=mini[i]+[1028,y];s.circle(x0,y0,11.5,'none',PURPLE,1.3)
s.arrow(1040,251,1040,277,MUTED,1.8)

s.text(959,95,'t',19,MUTED)
s.text(959,283,'t + Δt',19,MUTED)

OUT.mkdir(parents=True,exist_ok=True)
svg=OUT/'chromatin-feedback-concept.svg';svg.write_text(s.svg())
subprocess.run(['inkscape',str(svg),'--export-text-to-path','--export-plain-svg','--export-type=svg',f'--export-filename={svg}'],check=True,capture_output=True)
subprocess.run(['inkscape',str(svg),'--export-width=2400',f'--export-filename={OUT/"chromatin-feedback-concept.png"}'],check=True,capture_output=True)
print(svg)
