#!/usr/bin/env python3
import os
#from pathlib import path
import copy
import logging
import math
from collections import defaultdict
import argparse
import numpy as np
import subprocess
from datetime import datetime
import glob
import logging
log=logging.getLogger(__name__)
log.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
formatter = logging.Formatter('%(filename)s:%(lineno)04d - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
log.handlers=[console_handler]

def execute_mathematica(program,args,dir):
    log.info(f'executing mathematica program {program} with arguments {args} in {dir}')
    for i in range(10):#arbitrary max number of tries
        process=subprocess.Popen([f"wolframscript",program]+[a for a in args], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True,cwd=dir)
        outputs=[]
        for line in process.stdout:
            if (len(line.strip())>0):
                log.info(line.strip())
                if 'The product exited because of a license error' in line:
                    process.wait()
                    break
                outputs.append(line.strip())
        else:
            process.wait()
            break
    retcode=process.returncode
    if retcode !=0:
        raise Exception(f"error executing {program} with args {args} in {dir}")
    return outputs

def parse_cmd_line():
    parser = argparse.ArgumentParser(description='plot fourier transform of given pattern')
    parser.add_argument('--file', type=str,help='file containing studied grayscale image')
    parser.add_argument('--interval',help='Domain of analysed frequencies, square region [-intervall,intervall]^2' ,type=int)
    parser.add_argument('--threshold',help='threshfold of amplitude of a peak in order to be shown in a diffraction diagramm' ,type=float,default=0.01)
    parser.add_argument('--thresholdgauge',help='threshfold of amplitude of a peak in order to be considered for gauge calculations' ,type=float,default=0.01)
    parser.add_argument('--multiples',help='multiples of lattice basis vectors up to which we verify gauge linearity' ,type=int,default=4)
    parser.add_argument('--symmetries',help='Symmetries to verify under criterion of indistinguishability, indicated as element of GL(naturals,mu) with regards to the mu base vectors' ,type=str,nargs="*",default=[])
    args = parser.parse_args()
    return args
def main():
    args=parse_cmd_line()

    program=os.path.join(os.path.dirname(os.path.realpath(__file__)),'do_fourier_calculation.wls')
    execute_mathematica(program,[args.file,str(args.interval),str(args.threshold),str(args.multiples),str(args.thresholdgauge)]+args.symmetries,'.')
if __name__=="__main__":
    main()
