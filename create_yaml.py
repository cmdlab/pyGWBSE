import yaml

def get_true_false_input(input_msg):
    while True:
        user_input = input(input_msg)

        if user_input == '' or user_input is None:
            return True
        elif user_input == 'False' or user_input == 'false' or user_input == 'f' or user_input == 'F':
            return False
        else:
            print('Please enter either f or F or false or leave blank for the default value')


def get_integer_input(input_msg):
    while True:
        user_input = input(input_msg)

        try:
            if user_input == '' or user_input is None:
                break

            int(user_input)
            is_integer = True
        except ValueError:
            is_integer = False

        if is_integer:
            break
        else:
            print('Please Enter Integer Value or leave blank for default value')

    return user_input

def switch_operation(dict_file, num):
    if num == 1:
        dict_file["PARAMS"]["kapr"] = int(get_integer_input('KPOINT, parallelization, see VASP documentation, Default value is 4: ') or "4")
    elif num == 2:
        dict_file["PARAMS"]["ppn"] = int(get_integer_input('processors, available per kpoint (NCPU/KPAR), Default value is 13: ') or "13")
    elif num == 3:
        dict_file["PARAMS"]["reciprocal_density"] = int(get_integer_input('reciprocal density, Default value is 10: ') or "10")
    elif num == 4:
        dict_file["PARAMS"]["encutgw"] = int(get_integer_input('Input for encutgw, Default value is 100: ') or "100")
    elif num == 5:
        dict_file["PARAMS"]["nbgwfactor"] = int(get_integer_input('Input for nbgwfactor, Default value is 5: ') or "5")
    elif num == 6:
        dict_file["PARAMS"]["nomegagw"] = int(get_integer_input('Input for nomegagw, Default value is 50: ') or "50")
    elif num == 7:
        dict_file["PARAMS"]["convparam"] = input('Input for convparam, Default value is NOMEGA: ') or "NOMEGA"
    elif num == 8:
        dict_file["PARAMS"]["convsteps"] = int(get_integer_input('Input for convsteps, Default value is 10: ') or "10")
    elif num == 9:
        dict_file["PARAMS"]["conviter"] = int(get_integer_input('Input for conviter, Default value is 5: ') or "5")
    elif num == 10:
        dict_file["STRUCTURE"]["source"] = input(
            'MID(get structure from MP database with material id)/POSCAR (provide a POSCAR file), Default value is '
            'MID: ') or "MID"
    elif num == 11:
        dict_file["STRUCTURE"]["mat_name"] = input(
            'identifier for database entry (material id if source=MID), Default value is NEW_MAT: ') or "NEW_MAT"
    elif num == 12:
        dict_file["STRUCTURE"]["material_id"] = input(
            'material_id of the input structure in MP database, Default value is mp-149: ') or "mp-149"
    elif num == 13:
        dict_file["WFLOW_DESIGN"]["emc_fw"] = get_true_false_input(
            'set true for Effective mass calculation, Default value is true[leave blank], to input false -- type('
            'f/F/false/Flase): ')
    elif num == 14:
        dict_file["WFLOW_DESIGN"]["wannier_fw"] = get_true_false_input(
            'set true for wannier interpolation, Default value is true[leave blank], to input false -- type('
            'f/F/false/Flase): ')
    elif num == 15:
        dict_file["WFLOW_DESIGN"]["conv_fw"] = get_true_false_input(
            'true for performing NBAND convergence calculation. If true GW will be performed by using '
            'nbnd_conv_start*ppn, Default value is true[leave blank], to input false -- type(f/F/false/Flase): ')
    elif num == 16:
        dict_file["WFLOW_DESIGN"]["gw_fw"] = get_true_false_input(
            'set true for performing GW calculation, Default value is true[leave blank], to input false -- type('
            'f/F/false/Flase): ')
    elif num == 17:
        dict_file["WFLOW_DESIGN"]["bse_fw"] = get_true_false_input(
            'set true for performing BSE calculation, Default value is true[leave blank], to input false -- type('
            'f/F/false/Flase): ')

def create_yaml_file():
    dict_file = {
        "PARAMS": {
            "kpar": int(get_integer_input('KPOINT, parallelization, see VASP documentation, Default value is 4: ') or "4"),
            "ppn": int(get_integer_input('processors, available per kpoint (NCPU/KPAR), Default value is 13: ') or "13"),
            "reciprocal_density": int(get_integer_input('reciprocal density, Default value is 10: ') or "10"),
            "encutgw": int(get_integer_input('Input for encutgw, Default value is 100: ') or "100"),
            "nbgwfactor": int(get_integer_input('Input for nbgwfactor, Default value is 5: ') or "5"),
            "nomegagw": int(get_integer_input('Input for nomegagw, Default value is 50: ') or "50"),
            "convparam": input('Input for convparam, Default value is NOMEGA: ') or "NOMEGA",
            "convsteps": int(get_integer_input('Input for convsteps, Default value is 10: ') or "10"),
            "conviter": int(get_integer_input('Input for conviter, Default value is 5: ') or "5")
        },
        "STRUCTURE": {
            "source": input(
                'MID(get structure from MP database with material id)/POSCAR (provide a POSCAR file), Default value '
                'is MID: ') or "MID",
            "mat_name": input(
                'identifier for database entry (material id if source=MID), Default value is NEW_MAT: ') or "NEW_MAT",
            "material_id": input(
                'material_id of the input structure in MP database, Default value is mp-149: ') or "mp-149"
        },
        "WFLOW_DESIGN": {
            "emc_fw": get_true_false_input(
                'set true for Effective mass calculation, Default value is true[leave blank], to input false -- type('
                'f/F/false/Flase): '),
            "wannier_fw": get_true_false_input(
                'set true for wannier interpolation, Default value is true[leave blank], to input false -- type('
                'f/F/false/Flase): '),
            "conv_fw": get_true_false_input(
                'true for performing NBAND convergence calculation. If true GW will be performed by using '
                'nbnd_conv_start*ppn, Default value is true[leave blank], to input false -- type(f/F/false/Flase): '),
            "gw_fw": get_true_false_input(
                'set true for performing GW calculation, Default value is true[leave blank], to input false -- type('
                'f/F/false/Flase): '),
            "bse_fw": get_true_false_input(
                'set true for performing BSE calculation, Default value is true[leave blank], to input false -- type('
                'f/F/false/Flase): ')
        }
    }

    print('Your Inputs are: \n')
    for key, value in dict_file.items():
        print(key, '--')
        for k, v in value.items():
            print('\t' + k, ' : ', v)
        print('-----------------')

    while True:
        change_input = input('Do you want to change your input, type[y/Y/n/N]: ')

        if change_input == '' or change_input is None or change_input == 'N' or change_input == 'n':
            break

        print('Please select Number to change your input')
        print('1. kpar')
        print('2. ppn')
        print('3. reciprocal_density')
        print('4. encutgw')
        print('5. nbgwfactor')
        print('6. nomegagw')
        print('7. convparam')
        print('8. convsteps')
        print('9. conviter')
        print('10. source')
        print('11. mat_name')
        print('12. material_id')
        print('13. emc_fw')
        print('14. wannier_fw')
        print('15. conv_fw')
        print('16. gw_fw')
        print('17. bse_fw')

        get_input = int(get_integer_input('Please select Number to change your input, or leave blank for no change: ') or "0")

        if get_input == 0 or get_input is None or get_input == '':
            print('bye')
            break

        switch_operation(dict_file=dict_file, num=get_input)

    print('Your Inputs after updates are: \n')
    for key, value in dict_file.items():
        print(key, '--')
        for k, v in value.items():
            print('\t' + k, ' : ', v)
        print('-----------------')

    with open('input.yaml', 'w') as f:
        yaml.safe_dump(dict_file, f, sort_keys=False, default_flow_style=False)

    return
