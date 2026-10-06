#!/usr/bin/env python3
"""Finite interval triangle witness for the nonconvex exact CSS domain.

The coefficients are those in the appendix of Colbrook--Sadeghi--Stepaniants.
The pointwise map error is their displayed bound (map-centre-error).
Floating point selects a triangle; Arb certifies the three determinant signs.
"""
import json
import numpy as np
from scipy.spatial import ConvexHull
from flint import arb, acb, ctx

COEFFICIENTS = [
 '1.50795385361929982e-2', '5.76689422728655252e-3',
 '1.01557038850826134e-3', '2.71368762415077922e-4',
 '5.94916615455076681e-5', '1.40376217202844726e-5',
 '3.19906904482506334e-6', '7.39433442683568780e-7',
 '1.69810772559526110e-7', '3.90984465868492322e-8',
 '8.99107926650920794e-9', '2.06853696471111575e-9',
 '4.75778692314195368e-10', '1.09440420794361152e-10',
 '2.51724454756957816e-11', '5.78995922871089628e-12',
 '1.33169815228079109e-12', '3.06233427501141624e-13',
 '7.09000044668921563e-14', '1.71384604459278807e-14']

def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def interval_point(k,N):
    w=acb(0,2*arb.pi()*k/N).exp()
    z=w
    for j,c in enumerate(COEFFICIENTS,1):z+=arb(c)*w**(13*j+1)
    return z.real,z.imag
def sub(a,b):return a[0]-b[0],a[1]-b[1]
def length(a):return (a[0]**2+a[1]**2).sqrt()

N=4096
w=np.exp(2j*np.pi*np.arange(N)/N)
z=w.copy()
for j,c in enumerate(COEFFICIENTS,1):z+=float(c)*w**(13*j+1)
points=np.column_stack([z.real,z.imag])
hull=ConvexHull(points)
facets=hull.equations
distances=-(points @ facets[:,:2].T+facets[:,2]).max(axis=1)
point_index=int(distances.argmax());point=points[point_index]
nearest=int((point @ facets[:,:2].T+facets[:,2]).argmax())
ends=hull.simplices[nearest]
best=None
for third in hull.vertices:
    indices=[int(ends[0]),int(ends[1]),int(third)]
    if len(set(indices))!=3:continue
    triangle=points[indices]
    if cross(triangle[1]-triangle[0],triangle[2]-triangle[0])<0:
        indices[0],indices[1]=indices[1],indices[0];triangle=points[indices]
    margins=[cross(triangle[(i+1)%3]-triangle[i],point-triangle[i])/
             np.linalg.norm(triangle[(i+1)%3]-triangle[i]) for i in range(3)]
    if min(margins)>0 and (best is None or min(margins)>best[0]):
        best=(min(margins),indices)
if best is None:raise ArithmeticError('no triangle witness')

ctx.prec=512
p=interval_point(point_index,N)
vertices=[interval_point(k,N) for k in best[1]]
# The additional allowance covers rounding of all twenty displayed decimal
# coefficients relative to the source's exact hexadecimal centre.
epsilon=arb('4.6685340803e-10')+arb('1e-18')
certified=[]
for i in range(3):
    edge=sub(vertices[(i+1)%3],vertices[i]);rel=sub(p,vertices[i])
    det=cross(edge,rel)
    # Each difference vector changes by at most 2 epsilon in Euclidean norm.
    perturbation=2*epsilon*(length(edge)+length(rel))+4*epsilon**2
    if not det>perturbation:raise ArithmeticError('uncertain perturbed triangle sign')
    margin=det/length(edge)
    if not margin>arb('.01'):raise ArithmeticError('insufficient centre margin')
    certified.append(dict(centre_inset_lower=margin.lower().str(30),
                          perturbed_determinant_lower=(det-perturbation).lower().str(30)))
# Under convexity of the exact domain, the triangle of its three boundary
# points is in its closure. All three strict signs put its fourth boundary
# point in the interior of that triangle, a contradiction.
print(json.dumps(dict(status='EXACT_DOMAIN_NONCONVEX',nodes=N,point=point_index,
                     triangle=best[1],map_error=str(epsilon),edges=certified),indent=2))
