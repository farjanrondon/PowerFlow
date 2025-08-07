
from readdata.ReadData import read_data as rd
from Ybus.Ybus import y_bus
from OutputCopy.OutputCopy import output_copy
from MethodCalculus.Method import *
from MethodCalculus.S_Calculus import *
from MethodCalculus.WriteData import *

# First lets read the archive with the input data
archive_input = 'data_io.xlsx'
data_io = rd(archive_input)

# Structure of 'data_io' (a tuple)
'''
index 0: numeric methods to apply
index 1: max number of iterations
index 2: error of every numeric method to apply
index 3: name of the archive with the results
index 4: BUS sheet, this is a DataFrame
index 5: LINES sheet, this is a DataFrame too
index 6: TRX sheet, this is a DataFrame and is important to check if it's empty or not
index 7: SHUNT ELEMENTS sheet, a DataFrame, needs to be checked too
'''

to_apply = data_io[0]
iter_max = data_io[1]
err = data_io[2]
name = data_io[3]
bus_data = data_io[4]
line_data = data_io[5]
trx_data = data_io[6]
shunt_elements_data = data_io[7]

# Now we need to create the nodal admittance matrix; between the inspection method
# and the nodal incidence matrix method is best to use the fastest method
YBUS_matrix = y_bus(bus_data, line_data, trx_data, shunt_elements_data)

# Is important create a copy of the input file with the results, this is so
# relevant because the idea of the program is run multiple cases for study
# the power flow methodology
output_copy(name)

# With the 'to_apply' info, we proceed to verify witch numeric method needs to
# apply this program

# Gauss-Seidel
if to_apply[0] == 'Y':
    result_gs = gauss_seidel(YBUS_matrix, iter_max, err, bus_data)
    s_flow = sflow(bus_data, line_data, trx_data, result_gs[0])
    s_calc = scalc(bus_data, YBUS_matrix, result_gs[0])
    s_gen = sgen(bus_data, s_calc[0], result_gs[3])
    s_gbal = sgbal(s_gen[0], result_gs[3], s_flow[7], s_flow[8])
    result_table = get_result_table(bus_data, result_gs[0], s_calc[1], s_calc[2], s_gen[1], s_gen[2], result_gs[3])
    write_data('GS', result_table, s_flow, s_gbal, result_gs[4], result_gs[2], result_gs[1], name)

if to_apply[1] == 'Y':
    result_nr = newton_raphson(YBUS_matrix, iter_max, err, bus_data)
    s_flow = sflow(bus_data, line_data, trx_data, result_nr[0])
    s_calc = scalc(bus_data, YBUS_matrix, result_nr[0])
    s_gen = sgen(bus_data, s_calc[0], result_nr[3])
    s_gbal = sgbal(s_gen[0], result_nr[3], s_flow[7], s_flow[8])
    result_table = get_result_table(bus_data, result_nr[0], s_calc[1], s_calc[2], s_gen[1], s_gen[2], result_nr[3])
    write_data('NR', result_table, s_flow, s_gbal, result_nr[4], result_nr[2], result_nr[1], name)

if to_apply[2] == 'Y':
    result_fd = fast_descoupled(YBUS_matrix, iter_max, err, bus_data)
    s_flow = sflow(bus_data, line_data, trx_data, result_fd[0])
    s_calc = scalc(bus_data, YBUS_matrix, result_fd[0])
    s_gen = sgen(bus_data, s_calc[0], result_fd[3])
    s_gbal = sgbal(s_gen[0], result_fd[3], s_flow[7], s_flow[8])
    result_table = get_result_table(bus_data, result_fd[0], s_calc[1], s_calc[2], s_gen[1], s_gen[2], result_fd[3])
    write_data('FD', result_table, s_flow, s_gbal, result_fd[4], result_fd[2], result_fd[1], name)

if to_apply[3] == 'Y':
    result_dc = dc_power_flow(bus_data, line_data, trx_data)
    write_data_dc(name, result_dc)
