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

The script below tests a simple solenoid percurred by a current
'''
from numpy import r_, zeros, pi, ones, array
import myShapes
import Biot_Savart
from math import sqrt

myShape = myShapes.Wire()

#Solenoid with N turns and length l toward z axis
N = 100;
l = 1;
R = 1;
step = 0.1

myShape.Create_Solenoid(R, N, l, step)

#Points in the z axis (at the center axis of the solenoid) internal enough
pointsZ = array([0.4,0.5,0.6])
pointsY = zeros(len(pointsZ))
pointsX = zeros(len(pointsZ))
points = [pointsX,pointsY,pointsZ]

#Solve and get the norm of the B field
Solver = Biot_Savart.BSSolver()
B=Solver.Solve(myShape, points, 0.001)
Bm = Solver.absVector(B)

#Analytical solution for a solenoid
mu0=4*pi*1e-7;
BmAnalytical = abs(mu0*myShape.I*N/l*ones(len(pointsZ)))

#Plot numerical vs analytical
#import pylab as p
#p.plot(pointsZ,BmAnalytical,pointsZ,Bm,'x')
#p.show()

RMSE = sqrt(sum((BmAnalytical-Bm)**2.0)/len(Bm))
print 'RMSE:',RMSE
print 'Analytical B:',BmAnalytical[0],'Tesla'