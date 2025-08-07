import pandas as pd


def y_bus(bus_data, line_data, trx_data, shunt_elements_data):

    import numpy as np

    # To use the nodal incidence matrix method is important to first create the
    # incidence matrix; the best option consist in this easy algorithmic method:

    # It's known for a graph of 'n' nodes, that the incidence matrix has a
    # number of rows equal to the sum of factors of 'n', in other words:
    # number of rows = n+(n-1)+(n-2)+(n-3)+ ... + 1
    # On the other hand, this matrix has a number of columns equal to 'n'. For this
    # program, 'n' is just the number of buses of the electrical system.

    # Matrix of zeros with the dimension of the incidence matrix
    n = len(bus_data[1])
    n_rows = int(sum(np.linspace(1, n, n)))
    incidence_matrix = []

    for a in range(1, n+1):

        while n-a > 0:

            if a == 1:
                aux_a = np.c_[np.ones((n-a, 1)), (-1) * np.eye(n-a)]
                incidence_matrix = aux_a
            else:
                aux_b = np.c_[np.zeros((n-a, a-1)), np.ones((n-a, 1)), (-1) * np.eye(n-a)]
                incidence_matrix = np.r_[incidence_matrix, aux_b]
            break

        if n-a == 0:
            incidence_matrix = np.r_[incidence_matrix, np.eye(n)]

    # Lets to create the primitive matrix of admittances, this is a square matrix who
    # shape is the same number of rows of the incidence matrix; additionally this matrix
    # is diagonal
    y_primitive = (0+0j) * np.zeros((n_rows, n_rows))

    # Is important to create a two group of data, one of them need to be only with lines
    # admittances, this includes lines and trx, the other group is only for shunt elements
    # and capacitive effects from lines

    if isinstance(trx_data, pd.DataFrame):
        df_lines = pd.DataFrame(np.r_[line_data[0], trx_data[0]])
        lines_connection = (line_data[1].map(str) + line_data[2].map(str)).map(int)
        trx_connection = (trx_data[1].map(str) + trx_data[2].map(str)).map(int)
        df_lines[1] = np.r_[lines_connection, trx_connection]
        df_lines[2] = np.r_[line_data[3], trx_data[3]]
    else:
        df_lines = line_data
        lines_connection = (line_data[1].map(str) + line_data[2].map(str)).map(int)
        df_lines[1] = lines_connection
        df_lines = df_lines.drop(columns=[4, 2])

    if isinstance(shunt_elements_data, pd.DataFrame):
        df_shunt = pd.DataFrame(np.r_[line_data[1], trx_data[1], line_data[2], trx_data[2], shunt_elements_data[1]])
        df_shunt[1] = np.r_[0.5 * line_data[4], trx_data[4], 0.5 * line_data[4], trx_data[5], shunt_elements_data[2]]
    else:
        df_shunt = pd.DataFrame(np.r_[line_data[1], trx_data[1], line_data[2], trx_data[2]])
        df_shunt[1] = np.r_[0.5 * line_data[4], trx_data[4], 0.5 * line_data[4], trx_data[5]]

    # First the lines elements
    k = 1
    y_index = -1
    while True:

        for b in range(k+1, n+1):

            # Is important compare the index of every connection between buses
            ij = int(str(k) + str(b))
            ji = int(str(b) + str(k))

            tmp_data = df_lines.loc[(df_lines[1] == ij) | (df_lines[1] == ji)]

            y_index += 1
            if not tmp_data.empty:
                y_primitive[y_index][y_index] = sum(tmp_data[2])

        k += 1
        if k == n:
            break

    # Now the shunt elements
    for c in range(1, n+1):

        tmp_data = df_shunt.loc[df_shunt[0] == c]

        y_index += 1
        if not tmp_data.empty:
            y_primitive[y_index][y_index] = sum(tmp_data[1])

    # Finally the nodal admittance matriz is the product of 'y_primitive' with 'incidence_matrix'
    ybus_matrix = incidence_matrix.transpose() @ y_primitive @ incidence_matrix

    return ybus_matrix
