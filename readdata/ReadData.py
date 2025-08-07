
def read_data(archive_name):

    import pandas as pd
    import numpy as np

    # CONFIG sheet
    sh_config = pd.read_excel(archive_name, 'CONFIG')
    # methods to apply. 'Y' == apply. 'N' == no apply :
    apply = list(sh_config.iloc[0:4, 1])
    # max iterations :
    iter_max = sh_config.iloc[6, 1]
    # error :
    err = sh_config.iloc[5, 1]
    # name :
    name = sh_config.iloc[7, 1]

    # BUS sheet
    sh_bus = pd.read_excel(archive_name, 'BUS', header=None, skiprows=1)
    sh_bus[4] = sh_bus[4] * np.pi / 180

    # LINES sheet
    sh_lines = pd.read_excel(archive_name, 'LINES', header=None, skiprows=1)
    # For lines is important to have the expression of the impedance of evey line, for that reason
    # the next step is to build these impedance, as well as the derivation elements to ground and modify
    # the SH_LINES variable.

    # Impedance lines
    z = sh_lines[3] + sh_lines[4]*1j
    # Admittance lines
    y_lines = 1 / z

    # Admittances shunt
    y_shunt = sh_lines[5]*1j

    sh_lines = sh_lines.drop(columns=[3, 4, 5])
    sh_lines[3] = y_lines
    sh_lines[4] = y_shunt

    # TRX sheet
    sh_trx = pd.read_excel(archive_name, 'TRX', header=None, skiprows=1)
    # A possible user could test this program with a system without transformers, for that reason
    # this case would be considered.

    if sh_trx.empty:
        sh_trx = 'no trx'
    else:
        # TAP value
        tap = sh_trx[5]

        z_cc = sh_trx[3] + sh_trx[4] * 1j
        y_cc = 1 / z_cc

        # 'pi' model values
        y_ij = y_cc * tap
        y_io = (1 - tap) * y_cc
        y_jo = ((tap * tap) - tap) * y_cc

        # For the DC method is important the reactance of shortcircuit, the 'Xcc' row (index 4)
        # isn't dropped.
        sh_trx[6] = sh_trx[4]
        sh_trx[3] = y_ij
        sh_trx[4] = y_io
        sh_trx[5] = y_jo

    # 'SHUNT_ELEMENTS' sheet
    sh_shunt_elements = pd.read_excel(archive_name, 'SHUNT_ELEMENTS', header=None, skiprows=1)
    # Let's contemplate the possibility of a system without shunt elements.

    if sh_shunt_elements.empty:
        sh_shunt_elements = 'no shunt elements'
    else:
        z_shunt_element = sh_shunt_elements[2] + sh_shunt_elements[3] * 1j
        y_shunt_element = 1 / z_shunt_element

        sh_shunt_elements = sh_shunt_elements.drop(columns=3)
        sh_shunt_elements[2] = y_shunt_element

    return apply, iter_max, err, name, sh_bus, sh_lines, sh_trx, sh_shunt_elements
