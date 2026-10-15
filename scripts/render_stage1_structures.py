#!/usr/bin/env python3
"""Geometry-based structure comparison, not a schematic or generated molecule."""
import json
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from gaussian_io import ROOT, BONDS
from figure_fonts import font

data = json.loads((ROOT/"analysis/stage1/audit.json").read_text(encoding="utf-8"))
systems = data["systems"]
reference = np.array(systems[1]["coordinates_angstrom"])
base = reference[:6] - reference[:6].mean(axis=0)
_, _, vt = np.linalg.svd(base)
normal = vt[-1]
along = reference[5] - reference[0]
along /= np.linalg.norm(along)
across = np.cross(normal, along)
across /= np.linalg.norm(across)
if np.dot(reference[6]-reference[5], across) < 0:
    across *= -1
# Oblique view separates the two out-of-plane methyl hydrogens.
across = np.cos(np.pi/6)*across + np.sin(np.pi/6)*normal
projections=[]
for system in systems:
    xyz = np.array(system["coordinates_angstrom"])
    first = xyz[:6]-xyz[:6].mean(axis=0)
    u,_,vt = np.linalg.svd(first.T@base)
    rotation = u@np.diag([1,1,np.linalg.det(u@vt)])@vt
    aligned = (xyz-xyz[:6].mean(axis=0))@rotation
    projections.append(np.column_stack([aligned@along,aligned@across]))

image=Image.new("RGB",(2400,1030),"#f5f7fa")
draw=ImageDraw.Draw(image)
draw.text((65,45),"Optimized sorbic acid: three distinct local minima",font=font(48,True),fill="#15283e")
draw.text((65,114),"Gas-phase PBE0/6-31G(d)  |  Optimized (2E,4E)-sorbic acid conformers",font=font(29),fill="#506278")
all_points=np.concatenate(projections)
span=np.ptp(all_points,axis=0)
scale=min(630/span[0],330/span[1])
palette=["#3974a2","#ce8133","#788744"]
for i,(system,coords) in enumerate(zip(systems,projections)):
    left=50+i*780
    draw.rounded_rectangle((left,190,left+750,910),radius=22,fill="white",outline="#dce2e9",width=2)
    draw.rounded_rectangle((left+24,217,left+37,278),radius=6,fill=palette[i])
    draw.text((left+57,220),system["id"],font=font(37,True),fill="#15283e")
    if system["id"]=="conf-02":
        draw.text((left+310,233),"LOWEST OF THREE",font=font(21,True),fill=palette[i])
    center=(coords.max(axis=0)+coords.min(axis=0))/2
    points=(coords-center)*[scale,-scale]+[left+375,490]
    for a,b,order in BONDS:
        first,second=points[a-1],points[b-1]
        direction=second-first
        sideways=np.array([-direction[1],direction[0]])/np.linalg.norm(direction)*4
        for shift in ([sideways,-sideways] if order==2 else [np.zeros(2)]):
            draw.line([tuple(first+shift),tuple(second+shift)],fill="#697581",width=5)
    for j,(point,number) in enumerate(zip(points,system["atom_numbers"])):
        radius=18 if number==1 else 25
        color={1:"#e9edf1",6:"#354658",8:"#c9484a"}[number]
        draw.ellipse((point[0]-radius,point[1]-radius,point[0]+radius,point[1]+radius),fill=color,outline="#a7afb6",width=1)
        if number!=1 or j==15:
            label=("H" if number==1 else "O" if number==8 else "C")+str(j+1)
            draw.text((point[0],point[1]),label,font=font(16,True),anchor="mm",fill="#142537" if number==1 else "white")
    draw.text((left+37,710),f"Relative E:  {system['relative_E_kJmol']:.2f} kJ/mol",font=font(29,True),fill="#15283e")
    draw.text((left+37,762),f"Lowest mode:  {system['min_frequency_cm1']:.2f} cm-1",font=font(26),fill="#506278")
    draw.text((left+37,811),"E/E retained  |  42 real modes",font=font(25),fill="#506278")
    draw.text((left+37,855),"O-H orientation: "+("anti" if i==2 else "syn"),font=font(25),fill="#506278")
draw.text((65,965),"C1: methyl carbon; C6: carboxyl carbon; H16: hydroxyl proton. Same oblique view and scale for all structures.",font=font(27),fill="#506278")
path=ROOT/"analysis/stage1/optimized-structures.png"
image.save(path)
print(path)
