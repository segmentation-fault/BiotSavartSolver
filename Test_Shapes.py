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

The script below tests simple shapes and plots them
'''

from numpy import pi,array
import myShapes

myShape = myShapes.Wire()

#Toroidal solenoid
myShape.plotme()

#Solenoid
N = 100;
l = 1;
R = 0.1;
step = 0.1
myShape.Create_Solenoid(R, N, l, step)
myShape.plotme()

#loop
R = 2
C = [1 , 1]
NOP = 100
myShape.Create_Loop(C, R, NOP)
myShape.plotme()

#segmented wire
myShape.coordz=[]
myShape.AugmentWire( pi/2, 0, 1,array([0,0,0]))
myShape.AugmentWire( 0, 0, 1)
myShape.AugmentWire( pi/4, pi/4, 1)
myShape.plotme()

myShape.ShowPlots()