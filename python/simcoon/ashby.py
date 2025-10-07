# -*- coding: utf-8 -*-

import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib import rcParams
import math as m

def unit_vector(pt_a, pt_b):
    b_a = [pt_b[0]-pt_a[0],pt_b[1]-pt_a[1]]
    distance = m.sqrt(b_a[0]**2+b_a[1]**2)
    return np.array([b_a[0]*(1.0/distance),b_a[1]*(1.0/distance)])

def length(pt_a, pt_b):
    b_a = [pt_b[0]-pt_a[0],pt_b[1]-pt_a[1]]
    distance = m.sqrt(b_a[0]**2+b_a[1]**2)
    return distance      
           
def uv_2(uv1, uv2):
    b_a = [uv1[0]-uv2[0],uv1[1]-uv2[1]]
    distance = m.sqrt(b_a[0]**2+b_a[1]**2)
    return np.array([b_a[0]*(1.0/distance),b_a[1]*(1.0/distance)])

def poly_convexHull(points, color, coef_multi=0.1, rad=0.3, lw=2):
    """
    Plot the convex hull around a set of points as a 
    shaded polygon.
    """
    pt_env = points
    for i in range(0,len(points)):
        pt_env = np.append(pt_env, points[i] + (2.0*coef_multi*np.random.rand(10,2)-coef_multi)*points[i], axis=0)

    hull_env = ConvexHull(pt_env)
    hull_indices_env = hull_env.vertices

    u_v = np.zeros((len(hull_indices_env),2))
    
    verts = np.zeros((len(hull_indices_env),2))
    dist = np.zeros((len(hull_indices_env),2))
    
    for i in range(0,len(hull_indices_env)):
        verts[i] = pt_env[hull_indices_env[i]]
        u_v[i] = unit_vector(pt_env[hull_indices_env[i-1]], pt_env[hull_indices_env[i]]) 
        dist[i] = length(pt_env[hull_indices_env[i-1]], pt_env[hull_indices_env[i]])         

    verts2 = np.zeros((3*len(verts)+1,2))
    for i in range(0,len(hull_indices_env)-1):
        verts2[i*3] = verts[i] - u_v[i]*dist[i]*rad
        verts2[i*3+1] = verts[i]
        verts2[i*3+2] = verts[i] + u_v[i+1]*dist[i+1]*rad
    verts2[-4] = verts[-1] - u_v[-1]*dist[-1]*rad
    verts2[-3] = verts[-1]
    verts2[-2] = verts[-1] + u_v[0]*dist[0]*rad
    verts2[-1] = verts2[0]

#    for pt in pt_env:
#        plt.plot(pt[0], pt[1], 'ko')

#    for pt in points:
#        plt.plot(pt[0], pt[1], 'ro')

    codes = [Path.MOVETO,]
    for j in range(len(verts)):
        codes.extend([Path.CURVE3, Path.CURVE3, Path.LINETO,])
    #codes.append(Path.CURVE3)

    path = Path(verts2, codes)
    patch = patches.PathPatch(path, facecolor=color, lw=0, alpha=0.2)
    edge = patches.PathPatch(path, edgecolor=color, facecolor='none', lw=lw)
    plt.gca().add_patch(patch)
    plt.gca().add_patch(edge)


def poly_enclose(points, color, inc=1.2, rad=0.3, lw=2):
    """
        Plot the convex hull around a set of points as a
        shaded polygon.
        """
    hull = ConvexHull(points)
    
    
    cent = np.mean(points, 0)
    pts = []
    for pt in points[hull.vertices]:
        pts.append(pt.tolist())
    #        pts.append(pt.tolist())
    
    #    pts.sort(key=lambda p: np.arctan2(p[1] - cent[1],
    #                                    p[0] - cent[0]))
    #    pts = pts[0::2]  # Deleting duplicates
    pts.insert(len(pts), pts[0])


    verts = inc*(np.array(pts)- cent) + cent
    verts2 = np.zeros((3*verts.shape[0]-2,2))
    verts2[0::3] = verts
    verts2[1::3,:] = (1-rad)*verts[0:-1,:] + rad*verts[1:,:]
    verts2[2::3,:] = rad*verts[0:-1,:] + (1-rad)*verts[1:,:]
    verts2[0:-1] = verts2[1:]
    verts2[-1] = verts2[0]
    
    codes = [Path.MOVETO, Path.LINETO, Path.CURVE3,]
    for j in range(len(pts)-2):
        codes.extend([Path.CURVE3, Path.LINETO, Path.CURVE3,])
    codes.append(Path.CURVE3)


    path = Path(verts2, codes)
    patch = patches.PathPatch(path, facecolor=color, lw=0, alpha=0.2)
    edge = patches.PathPatch(path, edgecolor=color, facecolor='none', lw=lw)
    plt.gca().add_patch(patch)
    plt.gca().add_patch(edge)

def ellip_enclose(points, color, inc=1.2, lw=2, nst=2):
    """
        Plot the minimum ellipse around a set of points.
        
        Based on:
        https://github.com/joferkington/oost_paper_code/blob/master/error_ellipse.py
        """
    
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]
    
    x = points[:,0]
    y = points[:,1]
    cov = np.cov(x, y)
    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    w, h = 2 * nst * np.sqrt(vals)
    center = np.mean(points, 0)
    ell = patches.Ellipse(center, width=inc*w, height=inc*h, angle=theta,
                          facecolor=color, alpha=0.2, lw=0)
    edge = patches.Ellipse(center, width=inc*w, height=inc*h, angle=theta,
                                                 facecolor='none', edgecolor=color, lw=lw)
    plt.gca().add_artist(ell)
    plt.gca().add_artist(edge)
