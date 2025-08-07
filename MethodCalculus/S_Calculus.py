
def sflow(bus_data, line_data, trx_data, vbus):

    import pandas as pd
    import numpy as np
    from pandas import DataFrame as DF

    line_data[4] = 0.5 * line_data[4]
    line_data[5] = line_data[4]
    if isinstance(trx_data, pd.DataFrame):
        trx_data = trx_data.drop(columns=[6])
        system_data = np.r_[line_data, trx_data]
        system_data = pd.DataFrame(system_data)
    else:
        system_data = line_data

    s_ij = []
    s_ji = []
    for h in range(len(system_data[1])):

        i = system_data.iloc[h, 1]
        j = system_data.iloc[h, 2]

        # Power flow in one way (from i to j)
        aux1 = abs(vbus[i-1])**2 * system_data.iloc[h, 4].conjugate() + vbus[i-1]*((vbus[i-1]-vbus[j-1])*system_data.iloc[h, 3]).conjugate()
        s_ij.append(aux1)

        # Power flow in the other way (from j to i)
        aux2 = abs(vbus[j-1])**2 * system_data.iloc[h, 5].conjugate() + vbus[j-1]*((vbus[j-1]-vbus[i-1])*system_data.iloc[h, 3]).conjugate()
        s_ji.append(aux2)

    s_ij, s_ji = np.array(s_ij), np.array(s_ji)
    s_ij_real, s_ij_imag, s_ji_real, s_ji_imag = s_ij.real, s_ij.imag, s_ji.real, s_ji.imag
    s_ij_real, s_ij_imag, s_ji_real, s_ji_imag = DF(s_ij_real), DF(s_ij_imag), DF(s_ji_real), DF(s_ji_imag)

    system_data = system_data.drop(columns=[5, 4, 3])
    system_data[3], system_data[4] = s_ij_real, s_ij_imag
    system_data[5], system_data[6] = s_ji_real, s_ji_imag
    system_data[7], system_data[8] = s_ij_real + s_ji_real, s_ij_imag + s_ji_imag

    return system_data


def scalc(bus_data, ybus, vbus):

    import numpy as np

    s_calc = []

    for k in range(len(bus_data[1])):
        aux3 = 0
        for m in range(len(bus_data[1])):
            aux3 += vbus[k] * (vbus[m] * ybus[k][m]).conjugate()

        s_calc.append(aux3)

    s_calc = np.array(s_calc)
    p_calc, q_calc = s_calc.real, s_calc.imag

    return s_calc, p_calc, q_calc


def sgen(bus_data, s_calc, s_load):

    import numpy as np

    s_gen = []

    for n in range(len(bus_data[1])):

        if bus_data.iloc[n, 2] == 'SLACK':
            s_gen.append(s_calc[n] + s_load[n])
        elif bus_data.iloc[n, 2] == 'PV':
            s_gen.append(bus_data.iloc[n, 5] + (s_calc[n].imag + s_load[n].imag)*1j)
        else:
            s_gen.append(0+0j)

    s_gen = np.array(s_gen)
    p_gen, q_gen = s_gen.real, s_gen.imag

    return s_gen, p_gen, q_gen


def sgbal(s_gen, s_load, p_loss, q_loss):

    import pandas as pd

    s_gbal = pd.DataFrame({0: [0], 1: [0], 2: [0], 3: [0], 4: [0], 5: [0], 6: [0], 7: [0]})

    sum1, sum2, sum3, sum4 = sum(s_gen), sum(s_load), sum(p_loss), sum(q_loss)

    delta_p = sum1.real - sum2.real - sum3
    delta_q = sum1.imag - sum2.imag - sum4

    s_gbal.iloc[0] = [sum1.real, sum1.imag, sum2.real, sum2.imag, sum3, sum4, delta_p, delta_q]

    return s_gbal
