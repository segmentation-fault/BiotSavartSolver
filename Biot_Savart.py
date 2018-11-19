'''
Simple Biot Savart Solver for arbitrarily shaped wires
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

from numpy import cross,zeros,pi,array,linalg
import Discretizer

class BSSolver(object):
    '''
    Calculates the B field at given points for a Wire arbitrarily shaped
    '''


    def __init__(self):
        '''
        Nothing is a static class
        '''
        return
    
    def Solve(self,myWire,points,dL):
        '''
        Solves the B field at points points for a given wire myWire discretized uniformly
        at a length dL
        '''
        B = zeros((len(points[0]),3))
        Discr = Discretizer.Discretizer()
        
        for i in range (0, int(len(myWire.coordz[0]))-1):
            dSegment=Discr.Discretize_Uniform(myWire, dL, i)
            for p in range (0,int(len(points[0]))):  
                point = array([points[0][p],points[1][p],points[2][p]])
                B[p] = B[p] + self.__solveSegment(dSegment, point, myWire.I)
        
        return B
    
    def __solveSegment(self,dSegment,point,I_complex):
        '''
        Solves B in a given point for a given discretized Segment
        '''
        B = zeros((1,3))
        
        for i in range (0, int(len(dSegment[0]))-1):
            dL = [array([dSegment[0][i],dSegment[1][i],dSegment[2][i]]),array([dSegment[0][i+1],dSegment[1][i+1],dSegment[2][i+1]])]
            B = B + self.__dB(dL, point, I_complex)
        
        return B
    
    def __dB(self,dL,point,I_complex):      
        '''
        Solves dB according to Biot Savart
        '''       
        #Magnetic permittivity of the vacuum
        mu = 4*pi*1e-7;
        
        #Wire vector calculation
        ux = dL[1][0]-dL[0][0] 
        uy = dL[1][1]-dL[0][1]        
        uz = dL[1][2]-dL[0][2]
        
        vdL = array([ux,uy,uz])
        
        #Reference point in the middle of the segment
        ux = (dL[1][0]+dL[0][0])/2 
        uy = (dL[1][1]+dL[0][1])/2        
        uz = (dL[1][2]+dL[0][2])/2
        
        rp = array([ux,uy,uz])
        
        #Displacement vector
        ux = point[0]-rp[0]
        uy = point[1]-rp[1]       
        uz = point[2]-rp[2]
        
        r = array([ux,uy,uz])
        
        #Integral calculation according to Biot Savart
        r3 = linalg.norm(r)**3
        
        Integrand = cross(vdL,r)/r3
        
        dB = mu/(4*pi)*I_complex*Integrand
        
        return dB
    
    def absVector(self,B):
        '''
        Returns the norm of the vector field B
        '''
        aB = zeros(len(B))
        
        for i in range (0,len(aB)):
            aB[i]=abs(linalg.norm(B[i]))
        
        return aB
    
    