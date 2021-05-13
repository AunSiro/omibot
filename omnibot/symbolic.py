# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:05:38 2020

@author: Siro Moreno
"""
from sympy import (
    symbols,
    pi,
    cos,
    sin,
    simplify,
    integrate,
    Eq,
    solve,
    dsolve,
    Matrix,
    preorder_traversal,
    Float,
    solve_linear_system,
    eye,
    zeros,
    lambdify,
    sqrt,
)
from sympy.physics.mechanics import dynamicsymbols
from sympy.functions import sign


def integerize(expr):
    expr2 = expr
    for a in preorder_traversal(expr):
        if isinstance(a, Float):
            expr2 = expr2.subs(a, round(a))
    return expr2


def roundize(expr, n=4):
    expr2 = expr
    for a in preorder_traversal(expr):
        if isinstance(a, Float):
            expr2 = expr2.subs(a, round(a, n))
    return expr2


def generic_omnibot_mats(n=4, null_beta=True, equal_r=True):
    t, r, d, s = symbols("t r d s")
    alpha, beta = dynamicsymbols("alpha beta")

    W = Matrix([[r, cos(alpha)], [0, sin(alpha)]])
    T = Matrix([[cos(beta), -sin(beta)], [sin(beta), cos(beta)]])
    A = Matrix([[1, 0, -d], [0, 1, s]])
    WTA = W.inv() @ T.inv() @ A
    WTA.simplify()
    r_n = WTA[0, :]
    s_n = WTA[1, :]

    R_list = []
    S_list = []

    for ii in range(n):
        r_ii = r_n.subs(alpha, symbols("alpha_" + str(ii + 1)))
        r_ii = r_ii.subs(d, symbols("d_" + str(ii + 1)))
        r_ii = r_ii.subs(s, symbols("s_" + str(ii + 1)))
        s_ii = s_n.subs(alpha, symbols("alpha_" + str(ii + 1)))
        s_ii = s_ii.subs(d, symbols("d_" + str(ii + 1)))
        s_ii = s_ii.subs(s, symbols("s_" + str(ii + 1)))
        if not equal_r:
            r_ii = r_ii.subs(r, symbols("r_" + str(ii + 1)))
            s_ii = s_ii.subs(r, symbols("r_" + str(ii + 1)))
        if null_beta:
            r_ii = r_ii.subs(beta, 0)
            s_ii = s_ii.subs(beta, 0)
        else:
            r_ii = r_ii.subs(beta, symbols("beta_" + str(ii + 1)))
            s_ii = s_ii.subs(beta, symbols("beta_" + str(ii + 1)))

        S_list.append(s_ii)
        R_list.append(r_ii)

    R = Matrix(R_list)
    S = Matrix(S_list)
    return R, S


t, r, d, s, m, I_z, I_w, l, L, l_2, L_2 = symbols("t r d s m I_z I_w l L l_2 L_2")
x, y, alpha, beta, sigma, psi, theta = dynamicsymbols("x y alpha beta sigma psi theta")
psi_dot = psi.diff()
psi_dot_dot = psi.diff(t, 2)

q = Matrix([x, y, psi] + [dynamicsymbols(f"phi_{i+1}") for i in range(4)])
q_dot = q.diff(t)
q_dot_dot = q.diff(t, 2)
q_r = Matrix([x, y, psi])
q_r_dot = q_r.diff(t)
q_r_dot_dot = q_r.diff(t, 2)
Gamma = Matrix([dynamicsymbols(f"tau_{i+1}") for i in range(4)])
M_w = eye(4) * I_w
M_r = Matrix([[m, 0, 0], [0, m, 0], [0, 0, I_z]])
R_psi = Matrix([[cos(psi), -sin(psi), 0], [sin(psi), cos(psi), 0], [0, 0, 1]])


def dejabot_mats_body():
    R, S = generic_omnibot_mats()
    for ii in range(4):
        alpha = pi / 4 * (1 - 2 * ((int((ii + 1) / 2)) % 2))
        s = L * (1 - 2 * ((int((ii) / 2)) % 2))
        d = l * (1 - 2 * (ii % 2))
        R[ii, :] = R[ii, :].subs(symbols("alpha_" + str(ii + 1)), alpha)
        R[ii, :] = R[ii, :].subs(symbols("s_" + str(ii + 1)), s)
        R[ii, :] = R[ii, :].subs(symbols("d_" + str(ii + 1)), d)
        S[ii, :] = S[ii, :].subs(symbols("alpha_" + str(ii + 1)), alpha)
        S[ii, :] = S[ii, :].subs(symbols("s_" + str(ii + 1)), s)
        S[ii, :] = S[ii, :].subs(symbols("d_" + str(ii + 1)), d)
    return integerize(R), integerize(S)


def dejabot_mats_2():
    R, S = generic_omnibot_mats(null_beta=False)
    s_list = [l_2, L_2, -L_2, -l_2]
    d_list = [L_2, l_2, -l_2, -L_2]
    for ii in range(4):
        alpha = pi / 4 * (1 - 2 * ((int((ii + 1) / 2)) % 2))
        s = s_list[ii]
        d = d_list[ii]
        R[ii, :] = R[ii, :].subs(symbols("alpha_" + str(ii + 1)), alpha)
        R[ii, :] = R[ii, :].subs(symbols("s_" + str(ii + 1)), s)
        R[ii, :] = R[ii, :].subs(symbols("d_" + str(ii + 1)), d)
        R[ii, :] = R[ii, :].subs(symbols("beta_" + str(ii + 1)), pi / 4)
        S[ii, :] = S[ii, :].subs(symbols("alpha_" + str(ii + 1)), alpha)
        S[ii, :] = S[ii, :].subs(symbols("s_" + str(ii + 1)), s)
        S[ii, :] = S[ii, :].subs(symbols("d_" + str(ii + 1)), d)
        S[ii, :] = S[ii, :].subs(symbols("beta_" + str(ii + 1)), pi / 4)
    return integerize(R), integerize(S)


def calc_matrices(axes="body"):
    if axes == "body":
        R, S = dejabot_mats_body()
    elif str(axes) == "2":
        R, S = dejabot_mats_2()
    else:
        raise NameError("Unrecognized axes name:" + str(axes))

    H = simplify(M_r + R_psi @ R.T @ M_w @ R @ R_psi.T)
    K = simplify(R_psi @ R.T @ M_w @ R @ R_psi.diff().T)
    F_a = R_psi @ R.T @ Gamma
    F_a = Matrix(
        [
            F_a[0].factor(sin(psi), cos(psi)),
            F_a[1].factor(sin(psi), cos(psi)),
            F_a[2].factor(),
        ]
    )
    # A = R_psi@R.T COMPROBAR
    H_inv = simplify(H.inv())
    R_inv = simplify(R.pinv())
    return R, S, H, K, F_a, H_inv, R_inv


def calc_vectors(axes="body"):
    if axes == "body":
        R, S, H, K, F_a, H_inv, R_inv = calc_matrices(axes="body")
        q_w_dot = simplify(R @ R_psi.T @ q_r_dot)
    elif str(axes) == "2":
        R, S, H, K, F_a, H_inv, R_inv = calc_matrices(axes="2")
        w, w_1, w_2, a, a_1, a_2 = symbols("w w_1 w_2 a a_1 a_2")
        q_2_dot = Matrix([w_1, w_2, psi_dot])
        q_2_dot_dot = Matrix([a_1, a_2, psi_dot_dot])
        q_w_dot = simplify(R @ q_2_dot)
    else:
        raise NameError("Unrecognized axes name:" + str(axes))


# -----------  Electric variables and motor models ------------------------


a, b, k_e, k_m, n, r_e, tau_m, mu, v_max = symbols("a b K_e K_m N R_e tau_m mu V_{max}")
phi = dynamicsymbols("phi")
phi_dot = phi.diff(t)
