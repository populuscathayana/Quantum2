#!/usr/bin/env python3
"""Reconstruct s/sp/p/6d Gaussian MOs from fchk; verify analytic overlap.

No SCF or quantum-chemistry engine is run. Cartesian primitive normalization
and Gaussian shell ordering follow the formatted-checkpoint specification.
Isosurfaces use linear marching tetrahedra and an orthographic projection.
"""
import itertools
import json
import math

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from gaussian_io import ROOT, read_fchk, BONDS
from figure_fonts import font

SOURCE = ROOT / "data/calculations/conf-02/uv.fchk"
DEST = ROOT / "analysis/stage2"


def dfact(n):
    return math.prod(range(n, 0, -2))


def basis_functions(f):
    xyz = np.array(f["Current cartesian coordinates"]).reshape(-1,3)
    exp = np.array(f["Primitive exponents"])
    cs = np.array(f["Contraction coefficients"])
    cp = np.array(f["P(S=P) Contraction coefficients"])
    angular = {0:[(0,0,0)], 1:[(1,0,0),(0,1,0),(0,0,1)],
               2:[(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)]}
    basis, offset = [], 0
    for kind, count, atom in zip(f["Shell types"], f["Number of primitives per shell"], f["Shell to atom map"]):
        assert kind in (-1,0,1,2)
        sl = slice(offset, offset+count)
        powers = angular[0]+angular[1] if kind == -1 else angular[kind]
        for power in powers:
            a = exp[sl]
            coeff = cp[sl] if kind == -1 and sum(power)==1 else cs[sl]
            norm = (2*a/np.pi)**0.75 * np.sqrt((4*a)**sum(power)/math.prod(dfact(2*l-1) for l in power))
            basis.append((xyz[atom-1], power, a, coeff*norm))
        offset += count
    assert len(basis) == f["Number of basis functions"]
    return basis


def overlap_one(b1,b2):
    A, la, a, ca = b1
    B, lb, b, cb = b2
    a, b = a[:,None], b[None,:]
    p = a+b
    P = (a[...,None]*A + b[...,None]*B)/p[...,None]
    val = ca[:,None]*cb[None,:]*np.exp(-a*b/p*np.sum((A-B)**2))
    for axis in range(3):
        integral = np.zeros_like(p)
        for i in range(la[axis]+1):
            for j in range(lb[axis]+1):
                if (i+j)%2:
                    continue
                k = (i+j)//2
                integral += (math.comb(la[axis],i)*math.comb(lb[axis],j)*
                    (P[...,axis]-A[axis])**(la[axis]-i)*(P[...,axis]-B[axis])**(lb[axis]-j)*
                    dfact(2*k-1)/(2*p)**k*np.sqrt(np.pi/p))
        val *= integral
    return val.sum()


def evaluate(points, basis, coeff):
    values = np.zeros((len(points),len(coeff)))
    for k,(center,power,exps,weights) in enumerate(basis):
        r = points-center
        r2 = np.sum(r*r,axis=1)
        radial = np.exp(-r2[:,None]*exps) @ weights
        ao = radial*np.prod(r**np.array(power),axis=1)
        values += ao[:,None]*coeff[:,k]
    return values


def surface(grid, axes, iso):
    # The six tetrahedra consistently share the long 000-111 diagonal.
    corners = np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]])
    tets = [[0,1,3,7],[0,3,2,7],[0,2,6,7],[0,6,4,7],[0,4,5,7],[0,5,1,7]]
    shape = np.array(grid.shape)-1
    indices = np.array(np.meshgrid(*(np.arange(s) for s in shape),indexing="ij")).reshape(3,-1).T
    cubev = np.array([grid[c[0]:c[0]+shape[0],c[1]:c[1]+shape[1],c[2]:c[2]+shape[2]].ravel() for c in corners]).T
    keep = (cubev.min(axis=1)<iso)&(cubev.max(axis=1)>iso)
    indices, cubev = indices[keep],cubev[keep]
    step = np.array([x[1]-x[0] for x in axes])
    origin = np.array([x[0] for x in axes])
    triangles = []
    for tet in tets:
        vals = cubev[:,tet]
        loc = (indices[:,None,:]+corners[tet][None,:,:])*step+origin
        masks = np.sum((vals>iso)*np.array([1,2,4,8]),axis=1)
        for mask in range(1,15):
            rows = masks==mask
            if not rows.any():
                continue
            inside = [k for k in range(4) if mask&(1<<k)]
            outside = [k for k in range(4) if not mask&(1<<k)]
            def edge(i,j):
                frac = (iso-vals[rows,i])/(vals[rows,j]-vals[rows,i])
                return loc[rows,i]+frac[:,None]*(loc[rows,j]-loc[rows,i])
            if len(inside)==1:
                triangles.append(np.stack([edge(inside[0],j) for j in outside],axis=1))
            elif len(outside)==1:
                triangles.append(np.stack([edge(outside[0],j) for j in inside],axis=1))
            else:
                a,b = inside; c,d = outside
                ac,ad,bc,bd = edge(a,c),edge(a,d),edge(b,c),edge(b,d)
                triangles.extend([np.stack([ac,ad,bc],axis=1),np.stack([ad,bd,bc],axis=1)])
    return np.concatenate(triangles)


def main():
    f = read_fchk(SOURCE)
    basis = basis_functions(f)
    coeff = np.array(f["Alpha MO coefficients"]).reshape(136,136)
    overlap = np.empty((136,136))
    for i in range(136):
        for j in range(i+1):
            overlap[i,j] = overlap[j,i] = overlap_one(basis[i],basis[j])
    norm_error = float(np.max(np.abs(coeff @ overlap @ coeff.T - np.eye(136))))
    assert norm_error < 1e-6, norm_error
    xyz = np.array(f["Current cartesian coordinates"]).reshape(-1,3)
    center = xyz[:8].mean(axis=0)
    _,_,vt = np.linalg.svd(xyz[:8]-center)
    frame = vt.T
    if ((xyz[5]-xyz[0]) @ frame)[0] < 0:
        frame[:,0] *= -1
    frame[:,2] = np.cross(frame[:,0],frame[:,1])
    xyz_local = (xyz-center) @ frame
    bounds = np.array([xyz_local.min(axis=0)-3.5,xyz_local.max(axis=0)+3.5])
    axes = [np.arange(lo,hi+0.18,0.18) for lo,hi in bounds.T]
    gridpts = np.stack(np.meshgrid(*axes,indexing="ij"),axis=-1)
    flat = gridpts.reshape(-1,3) @ frame.T + center
    selected = coeff[[28,29,30]]
    values = np.concatenate([evaluate(flat[i:i+15000],basis,selected) for i in range(0,len(flat),15000)])
    volumes = values.reshape(*gridpts.shape[:-1],3)
    # A pi-type MO is approximately odd under reflection in the heavy-atom plane.
    probe = gridpts.reshape(-1,3)[::200]
    reflected = probe.copy(); reflected[:,2] *= -1
    a = evaluate(probe @ frame.T+center,basis,selected)
    b = evaluate(reflected @ frame.T+center,basis,selected)
    parity = np.sqrt(np.sum((a+b)**2,axis=0)/np.sum(a*a,axis=0))
    angle = np.deg2rad(32)
    view = np.array([[1,0,0],[0,np.cos(angle),-np.sin(angle)],[0,np.sin(angle),np.cos(angle)]])
    atomview = xyz_local @ view.T
    image = Image.new("RGB", (2400,1830),"white")
    draw = ImageDraw.Draw(image)
    draw.text((75,45),"Sorbic acid | orbitals involved in the two lowest singlet excitations",font=font(46),fill="#26333e")
    draw.text((75,113),"conf-02   /   PBE0/6-31G(d), SMD(diethyl ether)   /   isovalue = +/-0.035 bohr^(-3/2)",font=font(29),fill="#54616a")
    colors = [np.array([57,116,162]),np.array([206,129,51])]
    labels = ["MO 29  /  HOMO-1", "MO 30  /  HOMO", "MO 31  /  LUMO"]
    triangles_count = []
    for index in range(3):
        all_tri, phases = [], []
        for phase, iso in enumerate([0.035,-0.035]):
            tri = surface(volumes[...,index],axes,iso) @ view.T
            all_tri.append(tri); phases.extend([phase]*len(tri))
        triangles = np.concatenate(all_tri); phases=np.array(phases)
        triangles_count.append(len(triangles))
        scale = 62
        ycenter = 440+index*510
        def project(p): return np.stack([p[...,0]*scale+1480, -p[...,1]*scale+ycenter],axis=-1)
        draw.text((75,ycenter-190),labels[index],font=font(48),fill="#26333e")
        ev = f["Alpha Orbital Energies"][28+index]*27.211386245988
        draw.text((75,ycenter-129),f"{ev:.3f} eV",font=font(42),fill="#54616a")
        primitives=[]
        for k in range(len(triangles)):
            t = triangles[k]
            normal = np.cross(t[1]-t[0],t[2]-t[0]); normal /= max(np.linalg.norm(normal),1e-20)
            light = 0.7+0.3*abs(normal @ np.array([0.2,-0.3,0.9327379]))
            color = tuple(np.clip(colors[phases[k]]*light+255*(1-light)*0.12,0,255).astype(int))
            primitives.append((t[:,2].mean(),"tri",project(t),color))
        for i,j,order in BONDS:
            start,end = atomview[i-1],atomview[j-1]
            for frac in np.linspace(0,1,25,endpoint=False):
                p,q = start+(end-start)*frac,start+(end-start)*(frac+1/25)
                primitives.append(((p[2]+q[2])/2,"bond",project(np.array([p,q])),None))
        for k, p in enumerate(atomview):
            primitives.append((p[2],"atom",project(p),(k,f["Atomic numbers"][k])))
        for _, kind, p, color in sorted(primitives,key=lambda x:x[0]):
            if kind=="tri":
                draw.polygon([tuple(pt) for pt in p],fill=color)
            elif kind=="bond":
                draw.line([tuple(pt) for pt in p],fill="#737b81",width=10)
            else:
                k,z=color; radius=12 if z==1 else 20
                fill={1:"#edf0f2",6:"#8a9298",8:"#38424b"}[z]
                draw.ellipse((p[0]-radius,p[1]-radius,p[0]+radius,p[1]+radius),fill=fill,outline="#505a63",width=2)
        if index<2:
            draw.line((75,ycenter+278,2325,ycenter+278),fill="#e5e9ed",width=2)
    draw.text((75,1770),"Blue / orange: opposite orbital phases (not charge). All panels use the same orientation, scale and isovalue.",font=font(28),fill="#54616a")
    image.save(DEST / "frontier-orbitals.png")
    metadata={"source":str(SOURCE.relative_to(ROOT)), "mo_indices":[29,30,31],
        "max_all_MO_orthonormality_error":norm_error,"pi_odd_reflection_relative_residual":parity.tolist(),
        "grid_spacing_bohr":0.18,"isovalue_au":0.035,"grid_shape":list(volumes.shape[:3]),
        "triangles_per_orbital":triangles_count,"renderer":"normalized Cartesian GTO; linear marching tetrahedra; depth-sorted orthographic projection"}
    (DEST/"orbital-validation.json").write_text(json.dumps(metadata,indent=2)+"\n", encoding="utf-8")
    print(json.dumps(metadata,indent=2))


if __name__ == "__main__":
    main()
