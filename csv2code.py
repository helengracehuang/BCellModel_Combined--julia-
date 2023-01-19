import sys
#import deeptools
import numpy as np
from itertools import chain
import argparse
import os
import csv

IKB_NFKB_decay = 0.0144
IKB_NFKB_basaldecay = 0.0144

species_list = {}
process_list = {
    "initialConc": 0,
    "basalsynthesis": 1,
    "synthesis": 2,
    "translation": 3,
    "basaldecay": 4,
    "decay": 5,
    "association": 6,
    "dissociation": 7,
    "activation": 8,
    "deactivation": 9,
    "transport_in": 10,
    "transport_out": 11
}
species_params = {}

def csv2code_index(input, output):
    """
    reads in each line in the input files (.csv), convert to code
    :param name of files to be read and write
    """
    print("Converting to code (indexing constants)...")
    f = open(input,'r')
    out = open(output,'w')
    next(f) # skip first line (header)

    for line in f:
        line = line.strip()
        array = line.split(',')
        index = int(array[1]) + 2
        species = array[0]
        species_name = array[2]
        species_list[species.upper()] = species_name
        out.write('# ' + str(index) + ' : ' + species_name + '\n') # comment
        out.write('const ' + species.upper() + ' = ' + str(index) + ';' + '\n') # declare constant
    # print(species_list)
    out.write('\n')
    out.close()
    print("Fininshed converting csv to constant indices.")

def csv2code_name(input, output):
    """
    reads in each line in the input files (.csv), convert to code
    :param name of files to be read and write
    """
    print("Converting to code (species names)...")
    f = open(input,'r')
    out = open(output,'a')
    next(f) # skip first line (header)

    for line in f:
        line = line.strip()
        array = line.split(',')
        # index = int(array[1]) + 2
        index = int(array[1])
        species = array[0]
        species_name = array[2]
        # print(species_name)
        species_list[species.upper()] = species_name
        if species.startswith('t'):
            species_name = species_name[1:] + " transcript"
        elif species.startswith('n'):
            species_name =  "nuclear " + species_name[:-1]
        else:
            species_name = species_name + " protein"
        out.write('# ' + str(index) + ' : ' + species_name + '\n') # comment
        out.write('speciesNames[' + species.upper() + '] = \"' + species_name + '\";' + '\n') # declare constant
    # print(species_list)
    out.write('\n')
    out.close()
    print("Fininshed converting csv to species name list.")

def csv2code_param(input, output):
    """
    reads in each line in the input files (.csv), convert to code
    :param name of files to be read and write
    """
    print("Converting to code (parameters)...")
    f = open(input,'r')
    out = open(output,'a')
    next(f) # skip first line (header)

    for line in f:
        line = line.strip()
        array = line.split(',')
        species = array[6]
        process = array[7]
        param = float(array[3])

        IkB_bound_NFkB = False
        if ":" in species:
            IkB_bound_NFkB = True
            if species.startswith('n'):
                species_array = species[1:].split(':')
                species = 'n' + species_array[1] + species_array[0]
            else:
                species_array = species.split(':')
                species = species_array[1] + species_array[0]
        species = species.upper()
        if species in species_params:
            species_params[species][process_list[process]] = param
        else:
            species_params[species] = [0] * 12
            species_params[species][process_list[process]] = param
            if IkB_bound_NFkB:
                species_params[species][process_list["basaldecay"]] = IKB_NFKB_basaldecay
                species_params[species][process_list["decay"]] = IKB_NFKB_decay

    # idx = 3
    idx = 1
    for species in species_params:
        species_name = species_list[species]
        out.write('# ' + str(idx) + ' : ' + species_name + '\n') # comment
        out.write('setParams!((@view rates[' + species + ', 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=' + str(species_params[species][1]) +
        ', synthesis=' + str(species_params[species][2]) +
        ', translation=' + str(species_params[species][3]) +
        ', basaldecay=' + str(species_params[species][4]) +
        ', decay=' + str(species_params[species][5]) +
        ', association=' + str(species_params[species][6]) +
        ', dissociation=' + str(species_params[species][7]) +
        ', activation=' + str(species_params[species][8]) +
        ', deactivation=' + str(species_params[species][9]) +
        ', transport_in=' + str(species_params[species][10]) +
        ', transport_out=' + str(species_params[species][11]) + ');\n') # parameters
        idx += 1
    out.close()
    print("Fininshed converting csv to parameters.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Convert csv file from spreadsheet to code.')
    parser.add_argument('-i', '--input', required=True, dest='input',
                        help='*.csv file from species spreadsheet')
    parser.add_argument('-n', '--input2', required=True, dest='input2',
                        help='*.csv file from reaction spreadsheet')
    parser.add_argument('-o', '--output', required=True, dest='output',
                        help='output file name')


    args = parser.parse_args()
    input_species = args.input
    input_params = args.input2
    output = args.output

    csv2code_index(input_species, output)
    csv2code_name(input_species, output)
    csv2code_param(input_params, output)
