#!/usr/bin/python3
import argparse
import numpy as np
from PIL import Image
import logging_remote

def parse_cmd_line():
    parser = argparse.ArgumentParser(description='calculates numerically the fourier transform of a given array of wavectors for an image')
    parser.add_argument('image', type=str,default='image file to analyse')
    parser.add_argument('wavevectorfile',help=',-seperated csv file of 2 columns. This function adds two columns corresponding to absolute value and phase of the fourier transform', type=str)
    args = parser.parse_args()
    return args
def integrate_wavevectors(grayscale_values,wavevectors):
    imagesize=grayscale_values.shape#height then width
    data=np.hstack((wavevectors,np.zeros((wavevectors.shape[0],2))))
    i, j = np.meshgrid(np.arange(imagesize[0]), np.arange(imagesize[1]), indexing='ij')#j numbers columns and i numbers rows

    for ivec,wavevec in enumerate(wavevectors):
        exponent = -2j * np.pi * (((j-(imagesize[1]/2-0.5)) / np.sqrt(np.prod(imagesize))) * wavevec[0] - ((i-(imagesize[0]/2-0.5)) / np.sqrt(np.prod(imagesize))) * wavevec[1])
        #-(imagesize[1]/2-0.5) in order to put the origin in the center of the image
        #normalisation with np.sqrt(np.prod(imagesize)) allows to treat rectangular images while maintaining a coherence with the fft, i has negative sign because rows are numbered downwards, we subtract half of the imagesize so that (i,j)=(0,0) lies in the center of the image and not the corner
        result = np.sum(grayscale_values * np.exp(exponent))/np.prod(imagesize)
        data[ivec,2]=np.abs(result)
        data[ivec,3]=np.angle(result)
    return data
def main():
    args=parse_cmd_line()
    image = Image.open(args.image).convert('L')  # 'L' mode converts image to grayscale
    # Convert image to numpy array
    grayscale_values = np.array(image)/255
    #wavevectors
    wavevectors=np.genfromtxt(args.wavevectorfile,delimiter=',')
    data=integrate_wavevectors(grayscale_values,wavevectors)
    to_save=args.wavevectorfile
    np.savetxt(to_save,data,delimiter=',',header='k_x,k_y,Abs(fourier(rho)),Arg(fourier(rho))',fmt='%.10f') 

if __name__=="__main__":
    main()
