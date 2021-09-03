###############################
#OBJECTIVE
###############################

import pandas as pd
import subprocess
import time

data_file = 'data/c12_1797_err.dat'
POINT = pd.read_table(data_file, sep=' ', header=None)
POINT.columns = ['theta', 'sigma', 'err']
# POINT['sigma'] = POINT['sigma']

def runcommand(bashCommand):
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return [output, error]

def obj(parameters):
    score = 0
    for i, theta in enumerate(POINT['theta']):
        command = './../../bin/antip_scan 608.0 12.0 6.0 ' + \
                  str(parameters[0]) + ' ' +\
                  str(parameters[1]) + ' ' + \
                  str(parameters[2]) + ' ' + \
                  str(parameters[3]) + \
                  ' 1.25 ' + \
                  str(parameters[4]) + ' ' + \
                  str(parameters[5]) + ' ' + \
                  'ang ' + \
                  str(theta) + \
                  ' -1'
        output = runcommand(command)
        s = float(output[0].split()[1])
        # print(output,' theta ',theta,  ' and ', POINT['sigma'][i])
        score = score + ((s - POINT['sigma'][i])**2)/POINT['sigma'][i]

    return score


def obj_errweighted(parameters):
    score = 0
    for i, theta in enumerate(POINT['theta']):
        command = './../../../../mnt/d/progetti/elastic-scattering/bin/antip_scan 608.0 12.0 6.0 ' + \
                  str(parameters[0]) + ' ' + \
                  str(parameters[1]) + ' ' + \
                  str(parameters[2]) + ' ' + \
                  str(parameters[3]) + \
                  ' 1.25 ' + \
                  str(parameters[4]) + ' ' + \
                  str(parameters[5]) + ' ' + \
                  'ang ' + \
                  str(theta) + \
                  ' -1'
        output = runcommand(command)
        s = float(output[0].split()[1])
        # print(output,' theta ',theta,  ' and ', POINT['sigma'][i])
        partial_score = ((s - POINT['sigma'][i]) ** 2) / POINT['sigma'][i]
        score = score + partial_score*(1/POINT['err'][i])

    return score

# start = time.time()
# print(obj([79.2,50.8,1.066,1.260,0.361,0.465]))
# stop = time.time()
#
# print(stop-start)
