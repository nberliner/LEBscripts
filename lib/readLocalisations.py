# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 12:47:19 2015

@author: berliner
"""

import xml.etree.ElementTree as ET
import numpy as np
from pandas import DataFrame

class indexGenerator():
    """ Generates an increasing index for each localisation based on the frame
    number. Assumes the frames to be ordered. For use with pandas multiindex
    DataFrames """
    def __init__(self):
        self.currentFrame  = 0
        self.previousFrame = 0
        self.idx           = -1

    def __call__(self, frame):
        if frame == self.currentFrame:
            self.idx += 1
            return self.idx
        else:
            # update the frames
            self.previousFrame = self.currentFrame
            self.currentFrame  = frame
            self.idx = 0
            return self.idx

def generateIndex(allData, imageNumber):
    idx = indexGenerator()
    # assemble the index structure for the DataFrame
    level1 = [ 'frame_' + str(int(frame)) for frame in allData[:,imageNumber] ]
    level2 = [ idx(frame) for frame in allData[:,imageNumber] ]
    index  = [ level1, level2 ]
    return index

def readRapidStormLocalisations(fname, photonConversion=1.0, pixelSize=1.0):
    """ Read rapidStorm localisations from text file.
    
        photonConversion should be set to convert the photon counts correctly
       
        with,
            x    = x position
            sx   = error of x position
            y    = y position
            sy   = error of y position
            amp  = amplitude
    
    """
    assert( isinstance(pixelSize, float) or isinstance(pixelSize, int) ) # int for backwards compatibility
    
    photonConversion = float(photonConversion)
    
#    idx       = indexGenerator()
    pixelSize = float(pixelSize)
    
    xPosition            = None
    xPositionUncertainty = None
    yPosition            = None
    yPositionUncertainty = None
    amplitude            = None
    imageNumber          = None
    PSFpositionX         = None
    PSFpositionY         = None
    fitResidues          = None
    
    with open(fname, 'r') as f:
        # check the file structure
        header = f.readline()
        # 
        root = ET.fromstring(header[2:])
        for index, child in enumerate(root):
            if child.attrib['identifier'] == 'Position-0-0':
                xPosition = index
            elif child.attrib['identifier'] == 'Position-0-0-uncertainty':
                xPositionUncertainty = index
            elif child.attrib['identifier'] == 'Position-1-0':
                yPosition = index
            elif child.attrib['identifier'] == 'Position-1-0-uncertainty':
                yPositionUncertainty = index
            elif child.attrib['identifier'] == 'Amplitude-0-0':
                amplitude  = index
            elif child.attrib['identifier'] == 'ImageNumber-0-0':
                imageNumber = index
            elif child.attrib['identifier'] == 'PSFWidth-0-0':
                PSFpositionX = index
            elif child.attrib['identifier'] == 'PSFWidth-1-0':
                PSFpositionY = index
            elif child.attrib['identifier'] == 'FitResidues-0-0':
                fitResidues = index
            elif child.attrib['identifier'] == 'LocalBackground-0-0':
                localBackground = index
            
        # make sure all necessary fields have been identified
        assert( xPosition            != None )
        assert( yPosition            != None )
        assert( amplitude            != None )
        assert( imageNumber          != None )
        
    ## read the full file and add a NaN column
    allData  = np.loadtxt(fname, skiprows=1)
    rowCount = np.shape(allData)[0]
    
    # calculate the SNR
    SNR      = np.zeros((rowCount,1))
    SNR[:,0] = allData[:,amplitude] / allData[:,localBackground]
    allData  = np.concatenate((allData,SNR),axis=1)
    SNRindex = np.shape(allData)[1] - 1
    
    # add groupID
#    groupID = np.reshape( np.arange(rowCount), (-1,1) )
    zeros   = np.zeros((rowCount,1))
    zeros.fill(np.NaN)
    allData = np.concatenate((allData,zeros),axis=1)
    
#    groupIDindex = -2
    
    # convert the amplitude to photon count
    allData[:,amplitude] /= photonConversion
    
    # convert the pixelSize
    allData[:,[xPosition,yPosition]] /= pixelSize

    
    # check which data is there and what needs to be added
    if xPositionUncertainty is None:
        xPositionUncertainty = -1 # select the last column of all zeros
        yPositionUncertainty = -1
    if PSFpositionX is None:
            PSFpositionX = -1
            PSFpositionY = -1
    if fitResidues is None:
        fitResidues = -1
        
    dataIndexes = [xPosition, yPosition, xPositionUncertainty, yPositionUncertainty, \
                   amplitude, imageNumber, fitResidues, SNRindex]
    
    # assemble the index structure for the DataFrame
    index = generateIndex(allData, imageNumber)
#    level1 = [ 'frame_' + str(int(frame)) for frame in allData[:,imageNumber] ]
#    level2 = [ idx(frame) for frame in allData[:,imageNumber] ]
#    index  = [ level1, level2 ]
    
    columns = ['x','y','Uncertainty x','Uncertainty y','Photon Count', 'frame', \
               'FitResidue', 'SNR']
    
    data = DataFrame(allData[:,dataIndexes], columns=columns, index=index)

    return data




def readXYTLocalisations(fname, pixelSize=1.0):
    """
    Read a generic xyt file. The first line is used as header information.
    The following columns must be present (case sensitive!),
    
            x   y   frame
    
    the remaining columns are read and can be used for filtering.
    """
    # Read the header
    with open(fname, 'r') as f:
        header = f.readline().strip().split('\t')
        # Check the required columns
        assert( 'x' in header )
        assert( 'y' in header )
        assert( 'frame' in header )
    
    # Read the data
    allData  = np.loadtxt(fname, skiprows=1)
    
    # Sort ascending frames, thanks to: http://stackoverflow.com/a/2828121
    frameIndex = header.index('frame')
    allData = allData[allData[:,frameIndex].argsort()]
    
    # Assemble the index structure for the DataFrame
    index = generateIndex(allData, frameIndex)
#    idx         = indexGenerator()
#    level1 = [ 'frame_' + str(int(frame)) for frame in allData[:,frameIndex] ]
#    level2 = [ idx(frame) for frame in allData[:,frameIndex] ]
#    index  = [ level1, level2 ]
    
    # Put the data together
    data = DataFrame(allData, index=index, columns=header)
    
    # Convert from nm to px    
    data[['x','y']] = data[['x','y']] / float(pixelSize)
    
    return data



def readLEBLocalisations(fname, skiprows, columns):
    # Read the data
    allData  = np.loadtxt(fname, skiprows=skiprows)

    # Sort ascending frames, thanks to: http://stackoverflow.com/a/2828121
    frameIndex = columns.index('frame')
    allData = allData[allData[:,frameIndex].argsort()]
    
    # Assemble the index structure for the DataFrame
    index = generateIndex(allData, frameIndex)
    
    # Put the data together
    data = DataFrame(allData, index=index, columns=columns)
    
    return data






