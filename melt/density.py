from __future__ import annotations
import numpy as np

def grid_density_cic(pos,weights,box_size,grid_size):
    G=int(grid_size);L=float(box_size)
    h=L/G
    grid=np.zeros((G,G,G),dtype=np.float64)
    p=(pos%L)/h
    i0=np.floor(p).astype(np.int64)
    f=p-i0
    i0=i0%G
    i1=(i0+1)%G
    wx0,wy0,wz0=1.0-f[:,0],1.0-f[:,1],1.0-f[:,2]
    wx1,wy1,wz1=f[:,0],f[:,1],f[:,2]
    w=weights.astype(np.float64)
    for ix,iy,iz,wxa,wya,wza in [
        (i0[:,0],i0[:,1],i0[:,2],wx0,wy0,wz0),
        (i1[:,0],i0[:,1],i0[:,2],wx1,wy0,wz0),
        (i0[:,0],i1[:,1],i0[:,2],wx0,wy1,wz0),
        (i0[:,0],i0[:,1],i1[:,2],wx0,wy0,wz1),
        (i1[:,0],i1[:,1],i0[:,2],wx1,wy1,wz0),
        (i1[:,0],i0[:,1],i1[:,2],wx1,wy0,wz1),
        (i0[:,0],i1[:,1],i1[:,2],wx0,wy1,wz1),
        (i1[:,0],i1[:,1],i1[:,2],wx1,wy1,wz1),
    ]:
        np.add.at(grid,(ix,iy,iz),w*wxa*wya*wza)
    return grid

def density_field_A(pos,types,box_size,grid_size):
    return grid_density_cic(pos,types.astype(np.float64),box_size,grid_size)

def density_field_B(pos,types,box_size,grid_size):
    return grid_density_cic(pos,1.0-types.astype(np.float64),box_size,grid_size)

def density_contrast(pos,types,box_size,grid_size):
    w=2.0*types.astype(np.float64)-1.0
    return grid_density_cic(pos,w,box_size,grid_size)
