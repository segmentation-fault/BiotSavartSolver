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

The script below tests a toroidal solenoid percurred by a current
'''

from numpy import r_, zeros, pi, ones, array, linspace
import myShapes
import Biot_Savart
from math import sqrt

myShape = myShapes.Wire()

#Toroidal coil along the z axis
R1 = 10;
R2 = 1;
N = 100;
step = 0.001

myShape.Create_Toroidal_Coil(R1, R2, N, step)
myShape.Set_Current(1, 0)

#Points in the middle of the solenoid
myCircle = myShapes.Wire();
myCircle.Create_Loop([0,0],R1,5,'xy')
points = myCircle.coordz

#Solve and get the norm of the B field
Solver = Biot_Savart.BSSolver()
B=Solver.Solve(myShape, points, 0.1)
Bm = Solver.absVector(B)

#Analytical solution for a toroidal solenoid
mu0=4*pi*1e-7;
BmAnalytical = abs(mu0*myShape.I*N/(2*pi*R1)*ones(len(points[0])))

#Plot numerical vs analytical
#import pylab as p
#p.plot(linspace(0,2*pi,len(points[0])),abs(BmAnalytical-Bm))
#p.show()

RMSE = sqrt(sum((BmAnalytical-Bm)**2.0)/len(Bm))
print 'RMSE:',RMSE
print 'Analytical B:',BmAnalytical[0],'Tesla'