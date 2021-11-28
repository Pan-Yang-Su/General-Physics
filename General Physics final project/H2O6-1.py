from vpython import *
import math
size1,size2, m_o, m_h, k_bond = 60E-12, 53E-12, 16.0/6E23, 1.0/6E23, 10000.0 # These numbers are all real
size=(2*size1+size2)
m=(2*m_h+m_o)/3
theta=104.45/360*2*pi
d = 4*95.84E-12 # four times the real bond
d2 = d*(sin(theta/2))*2
g=9.8

class H2O_molecule:
      def __init__(self, pos, axis1 = vector(d, 0, 0), axis2 = vector(d*math.cos(theta), d*math.sin(theta), 0), make_trail = False):
            self.O = sphere(pos = pos, radius = size, color = color.red)
            self.H1 = sphere(pos = pos+axis1, radius = size, color = color.blue)
            self.H2 = sphere(pos = pos+axis2, radius = size, color = color.blue)
            self.bond1 = cylinder(pos = pos, axis = axis1, radius = (size1)/2.0, color = color.white)
            self.bond2 = cylinder(pos = pos, axis = axis2, radius = (size2)/2.0, color = color.white)
            self.bond3 = cylinder(pos = pos+axis1, axis = axis2-axis1, radius = (size2)/2.0, color = color.white, opacity = 0.2)
            self.O.m = m
            self.H1.m,self.H2.m= m,m
           
            self.bond1.k = k_bond
            self.bond2.k = k_bond
            self.bond3.k = k_bond
            
      def bond_force_on_O(self): # return bond force acted on the O atom
            force1=self.bond1.k*(mag(self.bond1.axis)-d)*norm(self.bond1.axis)
            force2=self.bond2.k*(mag(self.bond2.axis)-d)*norm(self.bond2.axis)
            return force1+force2
      def bond_force_on_H1(self):
            force1=-self.bond1.k*(mag(self.bond1.axis)-d)*norm(self.bond1.axis)
            force2=self.bond3.k*(mag(self.bond3.axis)-d2)*norm(self.bond3.axis)
            return force1+force2
      def bond_force_on_H2(self):
            force1=-self.bond2.k*(mag(self.bond2.axis)-d)*norm(self.bond2.axis)
            force2=-self.bond3.k*(mag(self.bond3.axis)-d2)*norm(self.bond3.axis)
            return force1+force2
      
      def Total_U(self): #return elastic potential energy
            U=0.5*self.bond1.k*(mag(self.bond1.axis)-d)**2+0.5*self.bond2.k*(mag(self.bond2.axis)-d)**2+0.5*self.bond3.k*(mag(self.bond3.axis)-d2)**2
           
            return U
      def G_Force(self):
            return self.H1.m*g,self.H2.m*g,self.O.m*g
