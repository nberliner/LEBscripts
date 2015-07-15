#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert super-resolution localisation data into movieInfo.mat Matlab files
to be used with utrack tracking software (http://lccb.hms.harvard.edu/software.html)

Author: Niklas Berliner (niklas.berliner@gmail.com)

(C) 2015 Niklas Berliner

"""

from __future__ import print_function

import os
import sys
import scipy.io
import numpy as np

import lib.readLocalisations  as rL


def createUtrackMatlabLocalisationFile(posList):
    """ Takes a list of localisations and the amplitude of the fit and
    generates a .mat file that can be read in with Matlab and used with
    u-track for tracking
    
    Note: Make sure the frames are counted from 0 and not 1 !!
    
    input: posList          type: dict()
    
    output outMat           type: dict() (see scipy.io.savemat)
    """
    # Create the empty array that holds the data.
    # See the u-track manual for a idea of how u-track stores the
    # localisation data
    matrixSize = np.int(np.max(posList.keys())) + 1 # since the frames are counted from 0,
                                                    # we need to add +1

    data    = np.array(np.zeros((matrixSize,1)) , np.dtype([('xCoord', 'O'), ('yCoord', 'O'), ('amp', 'O')]))

    # Set the values of the frames for which a localisation was found
    for frame in posList:
        try:
            xCoord = [np.array([ (item[0], item[1]) for item in posList[frame] ])]
            yCoord = [np.array([ (item[2], item[3]) for item in posList[frame] ])]
            amp    = [np.array([ (item[4], item[5]) for item in posList[frame] ])]
        except:
            print("Error occured trying to convert the localisations!\nNot saving to disk. Sorry!\n")
            sys.exit(1)

        data[frame]['xCoord'] = xCoord
        data[frame]['yCoord'] = yCoord
        data[frame]['amp']    = amp
    
    # The remaining values need to be set the an empty array, zero doesn't work
    # with utrack
    initialized = posList.keys()
    for index in range(matrixSize):
        if index in initialized:
            continue
        else:
            data[index]['xCoord'][0] = np.array([], dtype=np.float64)
            data[index]['yCoord'][0] = np.array([], dtype=np.float64)
            data[index]['amp'][0]    = np.array([], dtype=np.float64)
    
    outMat = dict()
    outMat['movieInfo'] = data
    return outMat


def convert(dataFrame):

    # Check the indices of the x,y,frame columns
    indexX = list(pos.columns).index('x')
    indexY = list(pos.columns).index('y')
    indexT = list(pos.columns).index('frame')
    
    # For thunderstorm input
    try:
        indexUncertainty = list(pos.columns).index('uncertainty')
        indexAmplitude   = list(pos.columns).index('intensity [photon]')
    except:
        indexUncertainty = None
        indexAmplitude   = None
    
    # Convert the DataFrame to an numpy array
    data = dataFrame.values
    
    # Put the localisations into a dict()
    newData = dict()
    for row in range(np.shape(pos)[0]):
        X = data[row,:][indexX]
        Y = data[row,:][indexY]
        T = data[row,:][indexT]
        
        if indexUncertainty is None:
            newData.setdefault(T, list()).append( (X, 0, Y, 0, 1, 0) )
        else:
            U   = data[row,:][indexUncertainty]
            Amp = data[row,:][indexAmplitude]
            newData.setdefault(T, list()).append( (X, U, Y, U, Amp, 0) )

    return newData



if __name__ == '__main__':

    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input file with localisations")
    parser.add_argument("-o", "--output", help="Where the data should be written to")
    parser.add_argument("-t", "--type", help="Type of the localisation file (xyt (0 default), rapidSTORM (1), thunderstorm (2)")
    parser.add_argument("-p", "--pixelSize", help="Pixel size of the data. Used as conversion factor to convert into pixels. Set to 1 if in doubt, then no change will happen.")
    args = parser.parse_args()
    
    PWD = os.getcwd()
    
    if args.file is None:
        print('Input file not specified.\nExiting.\n')
        sys.exit()
    else:
        fileNameIn = os.path.join(PWD, args.file)
        
    if args.output is None:
        output = '.'.join(fileNameIn.split('.')[:-1]) + '_movieInfo.mat'
        print("Output file set to %s" %output)
    else:
        output = os.path.join(PWD, args.output)
    
    if args.pixelSize == '1' or args.pixelSize is None:
        pixelSize = 1.0
    else:
        pixelSize = float(args.pixelSize)
        
    if args.type == '0' or args.type == 'xyt' or args.type is None:
        pos = rL.readXYTLocalisations(fileNameIn, pixelSize)
    elif args.type == '1' or args.type == 'rapidSTORM':
        pos = rL.readRapidStormLocalisations(fileNameIn, pixelSize)
    elif args.type == '2' or args.type == 'tunderstorm':
        pos = rL.readThunderstormLocalisations(fileNameIn, pixelSize)
    else:
        print('Input type not understood.\nExiting.\n')
        sys.exit()
        
    print('Converting..')
    # Convert the pandas DataFrame to dict (needed for compatibility reasons
    # to createUtrackMatlabLocalisationFile() )
    pos = convert(pos)
    
    # Convert the data
    outMat = createUtrackMatlabLocalisationFile(pos)
    
    # Save to disc
    scipy.io.savemat(output, outMat)

    print('Done.')





















