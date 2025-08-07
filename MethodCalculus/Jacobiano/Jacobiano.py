
def jacobiano(bus_data, p_spec, q_spec, v_st, p_st, abs_ybus, angle_ybus, pv_buses, pq_buses, n):

    import sympy as sp
    import numpy as np

    var_abs_name, var_angle_name = [], []
    for r in range(n):
        var_abs_name.append('p_'+str(r))
        var_angle_name.append('v_'+str(r))
    abs_vars, angle_vars = sp.symbols(var_abs_name), sp.symbols(var_angle_name)

    bus_data_slack = bus_data.loc[(bus_data[2] == 'SLACK')].index
    bus_data_pv = bus_data.loc[(bus_data[2] == 'PV')].index
    bus_data_pq = bus_data.loc[(bus_data[2] == 'PQ')].index
    for s in bus_data_slack:
        abs_vars[s], angle_vars[s] = bus_data.iloc[s, 3], bus_data.iloc[s, 4]
    for t in bus_data_pv:
        abs_vars[t] = bus_data.iloc[t, 3]

    delta_p, delta_q = [], []
    for u in (list(bus_data_pv) + list(bus_data_pq)):
        aux = 0
        for v in range(n):
            if u != v:
                aux += abs_vars[v] * abs_ybus[u][v] * sp.cos(angle_vars[u] - angle_vars[v] - angle_ybus[u][v])
        aux1 = p_spec[u] - ((abs_vars[u]**2) * abs_ybus[u][u] * sp.cos(angle_ybus[u][u]) + abs_vars[u] * aux)
        delta_p.append(aux1)

    for w in bus_data_pq:
        aux = 0
        for x in range(n):
            if w != x:
                aux += abs_vars[x] * abs_ybus[w][x] * sp.sin(angle_vars[w] - angle_vars[x] - angle_ybus[w][x])
        aux1 = q_spec[w] - (-(abs_vars[w] ** 2) * abs_ybus[w][w] * sp.sin(angle_ybus[w][w]) + abs_vars[w] * aux)
        delta_q.append(aux1)

    funct_delta = delta_p + delta_q

    # Now discard the values that we already know
    for y in bus_data_slack:
        abs_vars[y], angle_vars[y] = 0, 0
    for z in bus_data_pv:
        abs_vars[z] = 0

    vars_of_diff = angle_vars + abs_vars
    n_j = len(pv_buses) + 2*len(pq_buses)
    j = np.zeros((n_j, n_j)) * sp.Symbol('x')

    k_col = 0
    for a in range(len(vars_of_diff)):

        k_row = 0
        if vars_of_diff[a] != 0:

            for b in range(n_j):
                j[k_row, k_col] = sp.Derivative(funct_delta[b], vars_of_diff[a]).doit()
                k_row += 1

            k_col += 1

    # Now get the values to eval the jacobian
    sym_values = []
    num_values = []

    for c in range(len(vars_of_diff)):
        if vars_of_diff[c] != 0:
            sym_values.append(vars_of_diff[c])

    for d in range(len(angle_vars)):
        if angle_vars[d] != 0:
            angle_vars[d] = p_st[d]
            num_values.append(angle_vars[d])

    for e in range(len(abs_vars)):
        if abs_vars[e] != 0:
            abs_vars[e] = v_st[e]
            num_values.append(abs_vars[e])

    values_to_sub = dict(zip(sym_values, num_values))
    j = np.array(sp.Array(j).subs(values_to_sub).simplify(), dtype=float)

    return j
