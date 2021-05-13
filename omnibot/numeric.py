# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 18:21:32 2020

@author: Siro Moreno
"""

import numpy as np

l = 0.2096
L = 0.2096
L_2 = (L + l) / (2 ** 0.5)
l_2 = 0.0
r = 0.0667
m = 15.75
I_w = 0.00266
I_z = 0.461

phi_dot_max = 2 * np.pi * 7000 / (49 * 60)
psi_dot_max = phi_dot_max * 0.0667 / (2 ** 0.5 * 0.29642)


def max_speed_axes_2(psi_dot):
    L_2 = 0.29642
    r = 0.0667
    phi_dot_max = 2 * np.pi * 7000 / (49 * 60)
    return phi_dot_max * r / (2 ** 0.5) - np.abs(L_2 * psi_dot)
