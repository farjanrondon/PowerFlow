
def output_copy(name):

    import shutil
    from openpyxl import load_workbook

    # First create the copy of the archive
    shutil.copyfile('data_io.xlsx', f'{name}')

    # We need to consider possibility that the user execute the program many
    # times, so lets to increase the number of the 'name' variable, and
    # this new data will be write in the original data input archive
    number_name = int(name[4])
    number_name += 1

    new_name = 'case' + str(number_name) + '.xlsx'

    workbook_dataio = load_workbook('data_io.xlsx')
    sh_config = workbook_dataio['CONFIG']
    sh_config['B9'] = new_name
    workbook_dataio.save('data_io.xlsx')

    return 0
