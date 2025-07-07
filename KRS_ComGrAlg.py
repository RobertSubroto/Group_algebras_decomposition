# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 14:59:36 2024

@author: Bobby Subroto
"""

import math
import sympy 
from sympy.ntheory.factor_ import multiplicity
from sympy.ntheory.factor_ import divisors
from sympy.ntheory.residue_ntheory import n_order
import itertools
import pandas as pd
from functools import reduce 



def Div_list(p, m_vec):
    return [list(item) for item in itertools.product(*[divisors(n) for n in m_vec])]


def var_nu(d_vec, q):
    return reduce(math.lcm, [n_order(q, n) for n in d_vec])


def var_eta(d_vec, q):
    nominator = math.prod([sympy.totient(n) for n in d_vec])
    denominator = var_nu(d_vec, q)
    return int(nominator / denominator)


def KRS_decomposition(p, t, m_vec):
    q = p ** t
    p_adic_list = [multiplicity(p, n) for n in m_vec]
    p_rad_list = [int(n / (p ** multiplicity(p, n))) for n in m_vec]
    
    div_list = Div_list(p, p_rad_list)
    temp_list = []
    for n in div_list:
        temp_list.extend([[p, var_nu(n, q), p_adic_list]] * var_eta(n, q))
    index = pd.Index(temp_list)
    
    df = index.value_counts().reset_index()
    df.columns = ['Irreducible components', 'Count']
    return df


z = KRS_decomposition(2, 1, [5,5,32])








