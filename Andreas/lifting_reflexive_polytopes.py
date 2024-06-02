from cytools import Polytope
import numpy as np


def trivial_lift(poly):
    
    starting_dim = poly.dim()
    
    pts = poly.points()
    
    origin = np.zeros(starting_dim)
    
    new_pts = np.array([np.append(origin,1),np.append(origin,-1)])
    
    pts_lifted = np.append(pts,np.zeros((len(pts),1)),axis=1)
    
    pts_lifted = np.append(pts_lifted,new_pts,axis=0).astype(int)
    
    p_lifted = Polytope(pts_lifted)
    
    return p_lifted
