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

The script below tests discretization
'''

from numpy import array, pi
import myShapes
import Discretizer

myShape = myShapes.Wire()
Discr = Discretizer.Discretizer()

#segmented wire
myShape.coordz=[]
myShape.AugmentWire( 0, 0, 1,array([0,0,0]))
dSegment=Discr.Discretize_Uniform(myShape, 0.1, 0)
ax=myShape.plotme()
Discr.plotme(ax, dSegment)
myShape.ShowPlots()