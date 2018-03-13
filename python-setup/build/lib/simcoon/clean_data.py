# -*- coding: utf-8 -*-
import numpy as np
from scipy import interpolate
from math import factorial
import os

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def split_monotonic(data):

    np_split_indexes= np.where(np.diff(np.sign(np.diff(data))) != 0)[0]+2
    split_indexes= [1] + list(np_split_indexes) + [len(data)]
    return [data[a-1:b] for a,b in zip(split_indexes[:-1],split_indexes[1:])]

def split_monotonic_func(x,y):

    d = []
    np_split_indexes= np.where(np.diff(np.sign(np.diff(y))) != 0)[0]+2
    split_indexes= [1] + list(np_split_indexes) + [len(y)]
    d.append([x[a-1:b] for a,b in zip(split_indexes[:-1],split_indexes[1:])])
    d.append([y[a-1:b] for a,b in zip(split_indexes[:-1],split_indexes[1:])])
    return d;

def interpolate_multi_monotonic(x,y, num):
    
    d = split_monotonic_func(x,y)
    d2 = []
    x2_list = []
    y2_list = []
    for xi, yi in zip(d[0], d[1]):
        f_x2y2 = interpolate.interp1d(xi, yi)
        x2i = np.linspace(xi[0],xi[-1],num)
        y2i = f_x2y2(x2i)
        x2_list.append(x2i)
        y2_list.append(y2i)

    x_tot = np.array(x2_list[0])
    for a in x2_list[1:]:
        x_tot = np.append(x_tot,a[1:])

    y_tot = np.array(y2_list[0])
    for a in y2_list[1:]:
        y_tot = np.append(y_tot,a[1:])

    d2.append(x_tot)
    d2.append(y_tot)
    return d2;

def interpolatespline_monotonic(xi,yi,x2i,num):

    f_x2y2_in = interpolate.interp1d(xi,yi)
    xi_inter = np.linspace(np.amin(xi),np.amax(xi),num)
    yi_inter = f_x2y2_in(xi_inter)
    f_x2y2_out = interpolate.InterpolatedUnivariateSpline(xi_inter,yi_inter)

    xi_min = np.amin(xi)
    xi_max = np.amax(xi)

    y2i = np.zeros(len(x2i))
    
    for i in range(0,len(x2i)):
    
        if(x2i[i] < xi_min):
            y2i[i] = f_x2y2_out(x2i[i])
        elif(x2i[i] > xi_max):
            y2i[i] = f_x2y2_out(x2i[i])
        else:
            y2i[i] = f_x2y2_in(x2i[i])

    return y2i

def interpolatespline_multi_monotonic(x,y,x_res,num):
    
    d = split_monotonic_func(x,y)
    d_res = split_monotonic(x_res)
    
    assert(len(d[0])<=len(d_res))
    
    d2 = []
    x2_list = []
    y2_list = []
    for xi, yi, x2i in zip(d[0], d[1], d_res):
        
        y2i = interpolatespline_monotonic(xi,yi,x2i,num)
        x2_list.append(x2i)
        y2_list.append(y2i)
    
    x_tot = np.array(x2_list[0])
    for a in x2_list[1:]:
        x_tot = np.append(x_tot,a[1:])
    
    y_tot = np.array(y2_list[0])
    for a in y2_list[1:]:
        y_tot = np.append(y_tot,a[1:])
    
    d2.append(x_tot)
    d2.append(y_tot)
    return d2;



