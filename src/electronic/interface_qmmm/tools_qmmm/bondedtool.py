#!/usr/bin/env python3 
import numpy as np 

def cal_distance(coor_list):
    assert len(coor_list) == 2, coor_list
    a = np.array(coor_list[0])
    b = np.array(coor_list[1])
    dis = ((a - b) ** 2).sum() ** 0.5
    return dis 


def cal_angle(coor_list):
    assert len(coor_list) == 3, coor_list
    a = cal_distance(coor_list=coor_list[:2])
    b = cal_distance(coor_list=coor_list[1:])
    c = cal_distance(coor_list=[coor_list[0],coor_list[-1]])
    ang = np.degrees(np.arccos((a ** 2 + b ** 2 - c ** 2) / (2 * a * b )))

    return ang


def cal_dihedral(coor_list):
    assert len(coor_list) == 4, coor_list
    a = np.array(coor_list[0])
    b = np.array(coor_list[1])
    c = np.array(coor_list[2])
    d = np.array(coor_list[3])

    # calculate bond vectors
    v_ab = b - a 
    v_bc = c - b 
    v_cd = d - c

    # calculate vectors perpendicular to the bonds
    v1 = np.cross(v_ab, v_bc) 
    # normalize to length 1.0
    v1 /= np.linalg.norm(v1)

    v2 = np.cross(v_bc, v_cd)
    v2 /= np.linalg.norm(v2)

    # calculate the torsion angle between two vectors, and convert to degrees
    dih = np.degrees(np.arccos(np.dot(v1,v2)))

    if np.dot(v1, v_cd) < 0:
        dih *= -1

    return dih

