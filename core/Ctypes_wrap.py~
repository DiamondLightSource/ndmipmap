#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
using Ctypes to access C mipmap function
"""
import ctypes
from numpy.ctypeslib import ndpointer

def mipmapC():
    dll = ctypes.CDLL('/c_function/mipmap_core.so', mode=ctypes.RTLD_GLOBAL)
    func = dll.mipmap_core
    func.restype = None
    
    func.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),#input
                ctypes.c_float, # lambda
                ctypes.c_int, # iterations 
                ctypes.c_float, # epsil - tolerance
                ctypes.c_int, # methTV = 0 - isotropic, 1- anisotropic
                ctypes.c_int, # nonegativity = 0 (off)
                ctypes.c_int, # printing =0 (off)
                ctypes.c_int, # dimX
                ctypes.c_int, # dimY
                ctypes.c_int, # dimZ
                ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")] # output
    return func

