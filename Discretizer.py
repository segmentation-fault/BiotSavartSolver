'''
Simple Biot Savart Solver for arbitrarily shaped wires
Here the dicretization functions
Copyright (C) 2012  Antonio Franco (antonio_franco@live.it)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from numpy import array , append, linalg
from math import floor

class Discretizer(object):
    '''
    Discretizes a wire
    '''


    def __init__(self):
        '''
        Nothing
        '''
        return
    
    def plotme(self,ax,dSegment):
        X = dSegment[0]
        Y = dSegment[1]
        Z = dSegment[2]
        
        ax.plot(X,Y,Z,'x')
        return
        
    def Discretize_Uniform(self,myWire,dL,n):
        '''
        Discretizes the n-th segment of myWire in chuncks of length
        dL returning them in an array. n is zero based and represent the segment
        composed by the n-th vertex and the n+1-th vertex
        '''
        if len(myWire.coordz[0]) < n+1:
            return []
        
        Segment = [array([myWire.coordz[0][n],myWire.coordz[0][n+1]]),array([myWire.coordz[1][n],myWire.coordz[1][n+1]]),array([myWire.coordz[2][n],myWire.coordz[2][n+1]])]
        
        #Segment vector calculation
        ux = Segment[0][1]-Segment[0][0] 
        uy = Segment[1][1]-Segment[1][0]        
        uz = Segment[2][1]-Segment[2][0]
        
        u = array([ux,uy,uz])
        
        SegmentLength = linalg.norm(u)
        
        #Segment shorter than dL
        if SegmentLength<dL:
            return Segment 
        
        #How many segments?
        lengthRatio = SegmentLength/dL
        nSegments = floor(lengthRatio)
        
        dSegmentx = array([Segment[0][0]])
        dSegmenty = array([Segment[1][0]])
        dSegmentz = array([Segment[2][0]])
                    
        #length vector calculation
        v=dL*u/linalg.norm(u)        
            
        for i in range (1, int(nSegments)+1):
            dSegmentx=append(dSegmentx,dSegmentx[i-1]+v[0])
            dSegmenty=append(dSegmenty,dSegmenty[i-1]+v[1])
            dSegmentz=append(dSegmentz,dSegmentz[i-1]+v[2])
                 
        if nSegments != lengthRatio:
            dSegmentx=append(dSegmentx,Segment[0][1])
            dSegmenty=append(dSegmenty,Segment[1][1])
            dSegmentz=append(dSegmentz,Segment[2][1])
            
        dSegment = [dSegmentx,dSegmenty,dSegmentz]
                    
        return dSegment