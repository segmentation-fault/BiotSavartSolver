'''
Simple Biot Savart Solver for arbitrarily shaped wires
Here a simple Wire builder
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

from numpy import array,pi,cos,sin,r_,linspace,zeros,concatenate,append
import cmath

class Wire:
    '''
    Implements an arbitrary shaped wire
    '''    
    coordz = []
    '''Coordinates of the vertex of the wire in the form [X,Y,Z]'''
    I = complex(1,0)
    '''Complex current carried by the wire'''
    def __init__(self):
        '''By default initited as a toroidal coil with
            R1 = 10
            R2 = 1
            N = 100
            step = 0.001
            and current 1A with 0 phase
        '''
        R1 = 10
        R2 = 1
        N = 100
        step = 0.001
        self.Create_Toroidal_Coil(R1, R2, N, step)
        self.Set_Current(1, 0)
        return
    
    def Set_Current(self,modulus,angle):
        '''Sets current with absolute value modulus and phase angle (in radians)'''
        self.I = cmath.rect(modulus,angle)
        return
    
    def Create_Toroidal_Coil(self, R1 , R2 , N , step ):
        '''
        Create_Toroidal_Coil( R1 , R2 , N , step )
        Creates a toroidal coil of major radius R1, minor radius R2 with N
         turns and a step step
         Initiates coordz
        '''
        a=R1
        b=R2
        c=N
        
        t=r_[0:2*pi:step]
        
        X=(a+b*sin(c*t))*cos(t);
        Y=(a+b*sin(c*t))*sin(t);
        Z=b*cos(c*t);
        
        self.coordz = [X,Y,Z]
        
        return
    
    def Create_Solenoid(self, R , N , l , step ):
        '''
        Create_Solenoid(self, R , N , l , step )
        Creates a solenoid whose length is l with radius R, N turns with step
        step along the z axis
        '''
        a = R;
        b = l/(2*pi*N);
        T = l/b;
        
        t=r_[0:T:step]
        
        X = a*cos(t);
        Y = a*sin(t);
        Z = b*t;
        
        self.coordz = [X,Y,Z]
        return
    
    def Create_Loop(self,center,radius,NOP,Orientation='xy'):
        '''
        Create_Loop(self,center,radius,NOP)
        a circle with center defined as
        a vector CENTER, radius as a scaler RADIS. NOP is 
        the number of points on the circle.
        '''
        t=linspace(0,2*pi,NOP)
        
        if Orientation == 'xy':
            X = center[0]+radius*sin(t)
            Y = center[1]+radius*cos(t)
            Z = zeros(NOP)
        elif Orientation == 'xz':
            X = center[0]+radius*sin(t)
            Z = center[1]+radius*cos(t)
            Y = zeros(NOP)
        elif Orientation == 'yz':
            Y = center[0]+radius*sin(t)
            Z = center[1]+radius*cos(t)
            X = zeros(NOP)
        
        self.coordz = [X,Y,Z]
        return
    
    def AugmentWire(self,Theta,Phi,length,Origin=None):
        '''
        AugmentWire(self,Theta,Phi,length,Origin=None)
        augments the existing wire by a segment lenght long, starting from point
        Origin, with inclination Theta and Azimuth Phi. If origin = None then the last
        calculated vertex is used
        '''
        
        #If an origin is not specified, the last vertex is assumed as origin
        if not Origin is None:
            newWire = self.__Create_Wire(Origin, Theta, Phi, length)
        else:
            temp = array(self.coordz)
            newOrigin = temp[:,1]
            newWire = self.__Create_Wire(newOrigin, Theta, Phi, length)
        
        #If no coordinates are present, then we simply put the new vertices in the list
        if len(self.coordz)==0:
            self.coordz = newWire
        elif Origin != None:
            X = concatenate((self.coordz[0],newWire[0]),axis = 1)
            Y = concatenate((self.coordz[1],newWire[1]),axis = 1)
            Z = concatenate((self.coordz[2],newWire[2]),axis = 1) 
            self.coordz = [X,Y,Z]
        else:
            X = append(self.coordz[0],newWire[0][1])
            Y = append(self.coordz[1],newWire[1][1])
            Z = append(self.coordz[2],newWire[2][1])
            self.coordz = [X,Y,Z]
            
        return
    
    def __Create_Wire(self,Origin,Theta,Phi,length):
        '''
        Create_Wire(self,Origin,Theta,Phi,length)
        creates a single wire lenght long, starting from point
        Origin, with inclination Theta and Azimuth Phi and
        returns its coordinates in the form [X,Y,Z]
        '''
        
        #Computes the unit vector
        ux = cos(Phi)*sin(Theta)
        uy = sin(Phi)*sin(Theta)
        uz = cos(Theta)
        
        u = array([ux,uy,uz])
        
        #Computes the second vertex
        P2 = Origin + length * u
        
        X = array([Origin[0],P2[0]])
        Y = array([Origin[1],P2[1]])
        Z = array([Origin[2],P2[2]])
                
        return [X,Y,Z]
    
    def plotme(self):
        '''Plots itself
        inactive until ShowPlots is called'''        
        import pylab as p
        import mpl_toolkits.mplot3d.axes3d as p3
        
        X = self.coordz[0]
        Y = self.coordz[1]
        Z = self.coordz[2]
        
        fig=p.figure(None)
        ax = p3.Axes3D(fig)
        
        ax.plot(X,Y,Z)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        self.axis_equal(ax)
        
        p.draw()
        
        return ax
    
    def axis_equal(self,ax):
        '''Makes axis the same size'''
        X = self.coordz[0]
        Y = self.coordz[1]
        Z = self.coordz[2]
        
        Xmax = max(X)
        Ymax = max(Y)
        Zmax = max(Z)
        
        maxx = max(Xmax,Ymax,Zmax)
        
        ax.set_xlim3d(min(X),maxx)
        ax.set_ylim3d(min(Y),maxx)
        ax.set_zlim3d(min(Z),maxx)
        return
    
    def ShowPlots(self):
        '''Triggers pylab.show()'''
        import pylab as p
        
        p.show()
        return
    