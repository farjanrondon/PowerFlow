
def write_data(id_table, result_table, power_flow_table, power_balance_table, error, iterations, time, name):

    from pandas import ExcelWriter as EWriter
    from openpyxl import load_workbook

    with EWriter(path=name, mode='a', engine='openpyxl', if_sheet_exists='overlay') as writer:
        sh_name1 = 'RESULTS ' + id_table
        result_table.to_excel(writer, sheet_name=sh_name1, header=False, index=False, startrow=2)

        sh_name2 = 'POWER FLOW ' + id_table
        power_flow_table.to_excel(writer, sheet_name=sh_name2, header=False, index=False, startrow=1)

        sh_name3 = 'POWER BALANCE ' + id_table
        power_balance_table.to_excel(writer, sheet_name=sh_name3, header=False, index=False, startrow=2)

    workbook_case = load_workbook(name)
    sh_result = workbook_case[sh_name1]
    sh_result['C1'], sh_result['E1'], sh_result['G1'] = iterations, error, time
    workbook_case.save(name)

    return 0


def write_data_dc(name, bus_angle):

    from pandas import ExcelWriter as EWriter

    with EWriter(path=name, mode='a', engine='openpyxl', if_sheet_exists='overlay') as writer:
        bus_angle.to_excel(writer, sheet_name='RESULTS DC', header=False, index=False, startrow=1)

    return 0
