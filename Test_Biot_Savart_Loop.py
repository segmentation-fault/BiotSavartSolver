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

The script below tests a simple loop percurred by a current
'''

from numpy import r_, zeros, pi, mean
import myShapes
import Biot_Savart
from math import sqrt

myShape = myShapes.Wire()

#Loop of radius 1m in the yz plane
R = 1
C = [0 , 0]
NOP = 100
myShape.Create_Loop(C, R, NOP,'yz')

#Points in the x axis (at the center axis of the loop)
pointsX = r_[1:10.5:0.5]
pointsY = zeros(len(pointsX))
pointsZ = zeros(len(pointsX))
points = [pointsX,pointsY,pointsZ]

#Solve and get the norm of the B field
Solver = Biot_Savart.BSSolver()
B=Solver.Solve(myShape, points, 0.1)
Bm = Solver.absVector(B)

#Analytical solution for a loop
mu0=4*pi*1e-7;
BmAnalytical = abs(mu0*myShape.I*R**2/(2*(pointsX**2+R**2)**(3.0/2.0)))

#Plot numerical vs analytical
import pylab as p
p.plot(pointsX,BmAnalytical,pointsX,Bm,'x')
p.show()

RMSE = sqrt(sum((BmAnalytical-Bm)**2.0)/len(Bm))
print 'RMSE:',RMSE
print 'Analytical B (mean):',mean(BmAnalytical),'Tesla'