#!/usr/bin/env python3
import os
#from pathlib import path
import miscellaneous
import copy
import math
from collections import defaultdict
try:
    import userhomog_settings as homog_settings
except ImportError:
    import homog_settings
import argparse
import pandas as pd
from collections import defaultdict
from itertools import combinations
import numpy as np
import subprocess
from datetime import datetime
import homogenise_pattern
import plot_pattern
import glob

def parse_cmd_line():
    parser = argparse.ArgumentParser(description='plot fourier transform of given pattern')
    parser.add_argument('--file', type=str,help='file containing studied grayscale image')
    parser.add_argument('--intervall',help='how many pixels to include in diffraction diagram: square region from -intervall to +intervall' ,type=int)
    parser.add_argument('--threshold',help='threshfold of amplitude of a peak in order to be shown in a diffraction diagramm' ,type=float,default=0.01)
    parser.add_argument('--thresholdgauge',help='threshfold of amplitude of a peak in order to be considered for gauge calculations' ,type=float,default=0.01)
    parser.add_argument('--multiples',help='multiples of lattice basis vectors up to which we verify gauge linearity' ,type=int,default=4)
    parser.add_argument('--symmetries',help='Symmetries to verify under criterion of indistinguishability, indicated as element of GL(naturals,mu) with regards to the mu base vectors' ,type=str,nargs="*")
    args = parser.parse_args()
    return args
def main():
    args=parse_cmd_line()
    print(f'args={args}')

    program=os.path.join(os.path.dirname(os.path.realpath(__file__)),'do_fft_calculation.wls')
    miscellaneous.execute_mathematica(program,[args.file,str(args.intervall),str(args.threshold),str(args.multiples),str(args.thresholdgauge)]+args.symmetries,'.')
if __name__=="__main__":
    main()
