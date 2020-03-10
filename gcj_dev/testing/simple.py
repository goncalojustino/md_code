import numpy

def euclidean(array1, array2):
  x1 = float(array1[0])
  y1 = float(array1[1])
  z1 = float(array1[2])
  x2 = float(array2[0])
  y2 = float(array2[1])
  z2 = float(array2[2])
  
  euclidean = numpy.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
  return euclidean

a = [1,2,3]
b = [2,3,4]
print(euclidean(a,b))