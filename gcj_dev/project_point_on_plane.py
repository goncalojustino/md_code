# Python3 program to find  
# foot of perpendicular  
# of a point in a 3 D plane.  
  
# Plane is ax + by +cz +d = 0
# point is x1,y1,z1
# projection is x2,y2, z2

# Given a normal vector (a,b,c) and a centroid (a1,b1,c1)
# d = - a * a1 - b * b1 - c * c1

import numpy
from numpy import linalg

# Function to find foot of perpendicular  
def foot(a, b, c, a1,b1,c1, x1, y1, z1) : 
    d = - a * a1 - b * b1 - c * c1
    k = (-a * x1 - b * y1 - c * z1 - d) / (a * a + b * b + c * c);  
    x2 = a * k + x1;  
    y2 = b * k + y1;  
    z2 = c * k + z1;  
  
    print("x2 =",round(x2,1))  
    print("y2 =",round(y2,1)) 
    print("z2 =",round(z2,1)) 

    centroid1 = [a1,b1,c1]
    centroid2 = [x2,y2,z2]
    centroid1 = numpy.array(centroid1)
    centroid2 = numpy.array(centroid2)
    distance = numpy.linalg.norm(centroid1-centroid2)
    print ('distance=',distance)
  
  
# Driver Code  
if __name__ == "__main__" :  
  
    a = 1
    b = -2 
    c = 0
    a1 = 0
    b1 = 3
    c1 = 4 
    x1 = -1 
    y1 = 3
    z1 = 4 
  
    # function call  
    foot(a, b, c, a1,b1,c1, x1, y1, z1)  