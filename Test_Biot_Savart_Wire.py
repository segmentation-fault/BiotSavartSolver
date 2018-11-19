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

The script below tests a simple long wire percurred by a current
'''
from numpy import array, r_, zeros, pi, mean
import myShapes
import Biot_Savart
from math import sqrt

Solver = Biot_Savart.BSSolver()

#Very long wire toward the z-axis
myShape = myShapes.Wire()
myShape.coordz=[]
myShape.AugmentWire( 0, 0, 100,array([0,0,0]))
#10 points at the equator
pointsY = r_[1:11]
pointsX = zeros(len(pointsY))
pointsZ = zeros(len(pointsY))+50
points = [pointsX,pointsY,pointsZ]

#Solve and get the norm of the B field
B=Solver.Solve(myShape, points, 0.1)
Bm = Solver.absVector(B)

#Analytical solution for an infinite wire
mu0 = 4*pi*1e-7;
BmAnalytical = abs(mu0*myShape.I/(2*pi*pointsY))

#Plot numerical vs analytical
import pylab as p
p.plot(pointsY,BmAnalytical,pointsY,Bm,'x')
p.show()

RMSE = sqrt(sum((BmAnalytical-Bm)**2.0)/len(Bm))
print 'RMSE:',RMSE
print 'Analytical B (mean):',mean(BmAnalytical),'Tesla'