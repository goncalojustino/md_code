# Python program to find the Perpendicular(shortest) 
# distance between a point and a Plane in 3 D. 

# Plane is defined from normal vector (a, b, c)
# and a point (a1, b1, c1)

# The other point is (x1, y1, z1)


import math 
  
# Function to find distance 
def shortest_distance(x1, y1, z1, a, b, c, a1, b1, c1):  
    d = - a * a1 - b * b1 - c * c1
    e = abs((a * x1 + b * y1 + c * z1 + d))  
    f = (math.sqrt(a * a + b * b + c * c)) 
    print("Perpendicular distance is", e/f )
      
  
# Driver Code  
x1 = 10
y1 = 20
z1 = 32
a = 0
b = 0
c = 100
a1 = 0
b1 = 0
c1 = 2
  
# Function call 
shortest_distance(x1, y1, z1, a, b, c, a1, b1, c1)       