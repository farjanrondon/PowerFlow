import numpy as np


def gauss_seidel(ybus, iter_max, err, bus_data):

    import numpy as np
    import time
    import copy

    # For the numeric method, and for this specific use, is important
    # the model of the loads, we will use the ZIP model (impedance, current
    # power)
    p_load = bus_data[7] * (bus_data[9]*(bus_data[3])**2 + bus_data[10]*bus_data[3] + bus_data[11])
    p_spec = bus_data[5] - p_load
    q_load = bus_data[8] * (bus_data[9]*(bus_data[3])**2 + bus_data[10]*bus_data[3] + bus_data[11])
    q_spec = bus_data[6] - q_load

    s_load = p_load + q_load*1j

    # The Gauss-Seidel method can work with complex numbers
    start_values = bus_data[3] * np.exp(bus_data[4]*1j)
    start_values_aux = copy.deepcopy(start_values)

    # To count the iterations
    k = 1

    time_o = time.time()
    n = len(bus_data[1])

    while True:

        for d in range(n):

            if bus_data.iloc[d, 2] == 'PV':

                # Estimate of the reactive power
                sum_q = 0
                for e in range(n):

                    if (d != e) and (e < d):
                        sum_q += start_values_aux[e] * ybus[d][e]
                    else:
                        sum_q += start_values_aux[e] * ybus[d][e]
                q_estimate = - (start_values_aux[d].conjugate() * sum_q).imag

                # Estimate the bus voltage
                sum_v = 0
                for f in range(n):
                    if d != f:
                        sum_v += start_values_aux[f] * ybus[d][f]
                v_estimate = (1/ybus[d][d]) * (((p_spec[d] + q_estimate*1j)/(start_values_aux[d])).conjugate() - sum_v)

                # Bus voltage correction
                v_bus = (v_estimate / abs(v_estimate)) * abs(start_values_aux[d])

                # Finally for this kind of bus, save the result
                start_values_aux[d] = v_bus

            elif bus_data.iloc[d, 2] == 'PQ':

                # Only calculate the estimate bus voltage
                sum_v = 0
                for g in range(n):
                    if d != g:
                        sum_v += start_values_aux[g] * ybus[d][g]
                s_spec = p_spec[d] + q_spec[d] * 1j
                v_bus = (1 / ybus[d][d]) * ((s_spec / (start_values_aux[d])).conjugate() - sum_v)

                # Finally for this kind of bus, save the result
                start_values_aux[d] = v_bus

            else:
                continue

        abs_voltage_err = abs(start_values - start_values_aux)
        err_gs = max(abs_voltage_err)
        if (err_gs <= err) or (k == iter_max):
            break

        start_values = copy.deepcopy(start_values_aux)
        k += 1

    time_f = time.time()
    time_elapsed = time_f - time_o
    buses_voltages = start_values_aux

    return buses_voltages, time_elapsed, k, s_load, err_gs


def newton_raphson(ybus, iter_max, err, bus_data):

    import numpy as np
    import time
    import copy
    from pandas import DataFrame as Df
    from MethodCalculus.Jacobiano.Jacobiano import jacobiano
    from pandas import Series as Sr

    p_load = bus_data[7] * (bus_data[9] * (bus_data[3]) ** 2 + bus_data[10] * bus_data[3] + bus_data[11])
    q_load = bus_data[8] * (bus_data[9] * (bus_data[3]) ** 2 + bus_data[10] * bus_data[3] + bus_data[11])

    p_spec = bus_data[5] - p_load
    q_spec = bus_data[6] - q_load
    s_load = p_load + q_load * 1j

    abs_ybus, angle_ybus = abs(ybus), np.angle(ybus)

    # Start values, is important separate the modules from the angles; the NR method isn't work
    # with complex numbers
    bus_data_aux = copy.deepcopy(bus_data)
    v_st, p_st = bus_data_aux[3], bus_data_aux[4]
    # v_st: volts_start_values
    # p_st: phase_start_values

    # Iterations counter
    k = 1
    n = len(bus_data_aux[1])

    # Run time count
    time_o = time.time()
    while True:

        pv_buses = bus_data_aux.loc[(bus_data_aux[2] == 'PV')]
        pq_buses = bus_data_aux.loc[(bus_data_aux[2] == 'PQ')]
        pv_buses_index = pv_buses.index
        pq_buses_index = pq_buses.index

        delta_p, delta_q = [], []

        x1_old = pv_buses[4]
        x2_old, x3_old = pq_buses[4], pq_buses[3]
        x_old = Df(np.r_[x1_old, x2_old, x3_old])

        for o in range(n):

            p_aux = 0
            for p in range(n):
                if o != p:
                    p_aux += v_st[p] * abs_ybus[o][p] * np.cos(p_st[o] - p_st[p] - angle_ybus[o][p])

            if bus_data.iloc[o, 2] == 'PV':
                aux = p_spec[o] - (v_st[o] ** 2) * abs_ybus[o][o] * np.cos(angle_ybus[o][o]) - v_st[o] * p_aux
                delta_p.append(aux)
            elif bus_data.iloc[o, 2] == 'PQ':
                aux = p_spec[o] - (v_st[o] ** 2) * abs_ybus[o][o] * np.cos(angle_ybus[o][o]) - v_st[o] * p_aux
                delta_p.append(aux)

                q_aux = 0
                for q in range(n):
                    if o != q:
                        q_aux += v_st[q] * abs_ybus[o][q] * np.sin(p_st[o] - p_st[q] - angle_ybus[o][q])
                aux = q_spec[o] + (v_st[o] ** 2) * abs_ybus[o][o] * np.sin(angle_ybus[o][o]) - v_st[o] * q_aux
                delta_q.append(aux)

        funct_old = Df(np.r_[delta_p, delta_q])
        j = jacobiano(bus_data, p_spec, q_spec, v_st, p_st, abs_ybus, angle_ybus, pv_buses, pq_buses, n)

        x_new = x_old - (np.linalg.inv(j) @ funct_old)

        abs_voltage_err = abs(x_old - x_new)
        err_nr = max(abs_voltage_err[0])
        if err_nr <= err or k == iter_max:
            break

        p_to_new_iter = list(p_st)
        v_to_new_iter = list(v_st)
        count_aux = 0
        for f in list(pv_buses_index)+list(pq_buses_index):
            p_to_new_iter[f] = x_new.iloc[count_aux, 0]
            count_aux += 1
        for g in list(pq_buses_index):
            v_to_new_iter[g] = x_new.iloc[count_aux, 0]
            count_aux += 1

        p_st = copy.deepcopy(Sr(p_to_new_iter))
        v_st = copy.deepcopy(Sr(v_to_new_iter))

        bus_data_aux[3], bus_data_aux[4] = v_st, p_st

        k += 1

    time_f = time.time()
    time_elapsed = time_f - time_o
    buses_voltages = v_st * np.exp(p_st * 1j)

    return buses_voltages, time_elapsed, k, s_load, err_nr


def fast_descoupled(ybus, iter_max, err, bus_data):

    import numpy as np
    import time
    import copy
    from pandas import DataFrame as Df
    from pandas import Series as Sr

    p_load = bus_data[7] * (bus_data[9] * (bus_data[3]) ** 2 + bus_data[10] * bus_data[3] + bus_data[11])
    q_load = bus_data[8] * (bus_data[9] * (bus_data[3]) ** 2 + bus_data[10] * bus_data[3] + bus_data[11])

    p_spec = bus_data[5] - p_load
    q_spec = bus_data[6] - q_load
    s_load = p_load + q_load * 1j

    abs_ybus, angle_ybus, B = abs(ybus), np.angle(ybus), ybus.imag

    # Start values, is important separate the modules from the angles; the NR method isn't work
    # with complex numbers
    bus_data_aux = copy.deepcopy(bus_data)
    v_st, p_st = bus_data_aux[3], bus_data_aux[4]
    # v_st: volts_start_values
    # p_st: phase_start_values

    # Iterations counter
    k = 1
    n = len(bus_data_aux[1])

    # Imag matrix
    B1, B2 = copy.deepcopy(B), copy.deepcopy(B)
    slack_buses_index = sorted(list(bus_data_aux.loc[(bus_data_aux[2] == 'SLACK')].index), reverse=True)
    pv_buses_index = sorted(list(bus_data_aux.loc[(bus_data_aux[2] == 'PV')].index), reverse=True)

    to_delete_in_b1 = slack_buses_index
    B1 = np.delete(B1, to_delete_in_b1, axis=0)
    B1 = np.delete(B1, to_delete_in_b1, axis=1)

    to_delete_in_b2 = slack_buses_index + pv_buses_index
    B2 = np.delete(B2, to_delete_in_b2, axis=0)
    B2 = np.delete(B2, to_delete_in_b2, axis=1)

    # Run time count
    time_o = time.time()
    while True:

        pv_buses = bus_data_aux.loc[(bus_data_aux[2] == 'PV')]
        pq_buses = bus_data_aux.loc[(bus_data_aux[2] == 'PQ')]
        pv_buses_index = pv_buses.index
        pq_buses_index = pq_buses.index

        delta_p, delta_q = [], []

        x1_old = pv_buses[4]
        x2_old, x3_old = pq_buses[4], pq_buses[3]
        t_old = Sr(np.r_[x1_old, x2_old])
        v_old = Sr(np.array(x3_old))

        v_for_delta_p = []
        v_for_delta_q = []
        for o in range(n):

            p_aux = 0
            for p in range(n):
                if o != p:
                    p_aux += v_st[p] * abs_ybus[o][p] * np.cos(p_st[o] - p_st[p] - angle_ybus[o][p])

            if bus_data.iloc[o, 2] == 'PV':
                aux = p_spec[o] - (v_st[o] ** 2) * abs_ybus[o][o] * np.cos(angle_ybus[o][o]) - v_st[o] * p_aux
                delta_p.append(aux)
                v_for_delta_p.append(bus_data_aux.iloc[o, 3])
            elif bus_data.iloc[o, 2] == 'PQ':
                aux = p_spec[o] - (v_st[o] ** 2) * abs_ybus[o][o] * np.cos(angle_ybus[o][o]) - v_st[o] * p_aux
                delta_p.append(aux)
                v_for_delta_p.append((bus_data_aux.loc[o, 3]))
                v_for_delta_q.append((bus_data_aux.loc[o, 3]))

                q_aux = 0
                for q in range(n):
                    if o != q:
                        q_aux += v_st[q] * abs_ybus[o][q] * np.sin(p_st[o] - p_st[q] - angle_ybus[o][q])
                aux = q_spec[o] + (v_st[o] ** 2) * abs_ybus[o][o] * np.sin(angle_ybus[o][o]) - v_st[o] * q_aux
                delta_q.append(aux)

        funct_p = np.array(delta_p) / np.array(v_for_delta_p)
        funct_q = np.array(delta_q) / np.array(v_for_delta_q)

        t_new = t_old - (np.linalg.inv(B1) @ funct_p)
        v_new = v_old - (np.linalg.inv(B2) @ funct_q)

        abs_voltage_err = abs(v_old - v_new)
        phase_voltage_err = abs(t_old - t_new)
        err_fd = max(max(abs_voltage_err), max(phase_voltage_err))
        if err_fd <= err or k == iter_max:
            break

        p_to_new_iter = list(p_st)
        v_to_new_iter = list(v_st)
        count_aux = 0
        for f in list(pv_buses_index) + list(pq_buses_index):
            p_to_new_iter[f] = t_new[count_aux]
            count_aux += 1

        count_aux = 0
        for g in list(pq_buses_index):
            v_to_new_iter[g] = v_new[count_aux]
            count_aux += 1

        p_st = copy.deepcopy(Sr(p_to_new_iter))
        v_st = copy.deepcopy(Sr(v_to_new_iter))

        bus_data_aux[3], bus_data_aux[4] = v_st, p_st

        k += 1

    time_f = time.time()
    time_elapsed = time_f - time_o
    buses_voltages = v_st * np.exp(p_st * 1j)

    return buses_voltages, time_elapsed, k, s_load, err_fd


def dc_power_flow(bus_data, line_data, trx_data):

    from pandas import DataFrame as Df
    from pandas import Series as Sr
    import numpy as np

    p_load = bus_data[7] * (bus_data[9] * (bus_data[3]) ** 2 + bus_data[10] * bus_data[3] + bus_data[11])
    p_spec = bus_data[5] - p_load

    bus_data_aux = Df(np.c_[bus_data[1], bus_data[4], p_spec])
    line_data_aux = Df(np.c_[np.r_[line_data[1], trx_data[1]], np.r_[line_data[2], trx_data[2]], np.r_[Sr(1/(1/np.array(line_data[3])).imag), (1/trx_data[6])]])

    n = len(bus_data[1])
    X = np.zeros((n, n))

    for h in range(n):
        for k in range(n):

            # First the elements of the diagonal
            if h == k:
                tmp_data = line_data_aux.loc[(line_data_aux[0] == h+1) | (line_data_aux[1] == k+1)]
                if not tmp_data.empty:
                    X[h][k] = sum(tmp_data[2])

            # For the elements out of the diagonal
            if h != k:
                tmp_data = line_data_aux.loc[(line_data_aux[0] == h+1) & (line_data_aux[1] == k+1)]
                if not tmp_data.empty:
                    X[h][k] = -sum(tmp_data[2])
                    X[k][h] = -sum(tmp_data[2])

    slack_buses_index = list(bus_data.loc[(bus_data[2] == 'SLACK')].index)
    X = np.delete(X, slack_buses_index, axis=0)
    X = np.delete(X, slack_buses_index, axis=1)

    p_spec_arr = np.array(p_spec)
    p_spec_arr = np.delete(p_spec_arr, slack_buses_index, axis=0)

    dc_result = np.linalg.inv(X) @ p_spec_arr
    dc_result = list(dc_result)

    for m in slack_buses_index:
        dc_result.insert(m, bus_data.iloc[m, 4])

    dc_result = Df(np.c_[bus_data[1], np.array(dc_result)*180/np.pi])

    return dc_result


def get_result_table(bus_data, vbus, p_calc, q_calc, p_gen, q_gen, s_load):

    import numpy as np
    from pandas import DataFrame as Df

    result_table = bus_data.drop(columns=[11, 10])
    result_table[2], result_table[3] = abs(vbus), Df(np.angle(vbus))*180/np.pi
    result_table[4], result_table[5] = Df(p_calc), Df(q_calc)
    result_table[6], result_table[7] = Df(p_gen), Df(q_gen)
    result_table[8], result_table[9] = Df(np.array(s_load).real), Df(np.array(s_load).imag)

    return result_table
