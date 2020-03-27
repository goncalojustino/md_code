# Catarina Nascimento<!-- omit in toc --> 
- [Salt bdrige script](#salt-bdrige-script)
- [Arginine dimer script](#arginine-dimer-script)
- [Iron coordination](#iron-coordination)
  - [Distances to iron through time](#distances-to-iron-through-time)
  - [Coordination angles](#coordination-angles)
  - [Screening interactions involving aromatic residues](#screening-interactions-involving-aromatic-residues)
    - [Aromatic - ion interactions](#aromatic---ion-interactions)
      - [Steps](#steps)
      - [1. Look for point charges (cations and anions)](#1-look-for-point-charges-cations-and-anions)
      - [2. Look for aromatic rings](#2-look-for-aromatic-rings)
      - [3. Check that centroid - point charge distance is within cutt-off](#3-check-that-centroid---point-charge-distance-is-within-cutt-off)
      - [4. Compute angular displacement](#4-compute-angular-displacement)
    - [Aromatic - H-X interactions](#aromatic---h-x-interactions)
      - [1. Use same aromatics](#1-use-same-aromatics)
      - [2. Define set of HX groups anf HC groups](#2-define-set-of-hx-groups-anf-hc-groups)
      - [3. Work out the distances and the angles](#3-work-out-the-distances-and-the-angles)
  - [Katie's One Drive CysMet and CysHBond CYS-HBOND](#katies-one-drive-cysmet-and-cyshbond-cys-hbond)
    - [Aromatic - Aromatic interactions](#aromatic---aromatic-interactions)
- [Geometric snippets](#geometric-snippets)
  - [Angle between planes](#angle-between-planes)
  - [Distance between parallel planes](#distance-between-parallel-planes)
  - [Distance of a point to a plane - projection of point to plane](#distance-of-a-point-to-a-plane---projection-of-point-to-plane)
# Salt bdrige script
- how are they defined for the script
- measuring distances
# Arginine dimer script
- what parameters are defined for the script and how
- distances
- angles
# Iron coordination 
## Distances to iron through time
- Identify atoms at a distance < 2.5 A than the FE atom
- Measure distances through time
## Coordination angles
- Using the previously selected neighbours (L1, L2…) compute all Li M Lj angles on a **structure**
    L1,L2... coordinates as np.array ?
    
    for each Li,Lj pair compute L1MLj angle = angle between vectors v1 = MLi and v2 = MLj
      #note: L1,L2 gives the same angle as L2,L1  
      v1 = Li - M
      v2 = Lj - M
      cosine = numpy.dot(v1, v2) / (numpy.linalg.norm(v1) * numpy.linalg.norm(v2))
      angle = numpy.arccos(cosine)
      print(numpy.degrees(angle))
    
    store all angles - Li M Lj identification and angle amplitude
- Predicted output: for 1 M with 6 L neighbours, this gives 5 + 4 + 3 + 2 + 1 = 15 angle values, that we want to compute “octahedral index” on a **structure**
- All angles must be between 0 and 180
- Take the largest value (alpha) and second largest (beta) of amplitudes [no id req.]:
  ```python
    tau_six = (360 - alpha - beta) / 180
    tau_six_abs = (360 - alpha - beta) / 180
    tau_siz_normf = [360 - DEGREES(ARCCOS(5/7)) ]  / 180
    tau_six = tau_six_abs / tau_six_normf
  ```

## Screening interactions involving aromatic residues

- two goals:
  - check if the metal arom (FE, V...) participates in any interaction with aromatic rings
  - check if there are aromatic rings involved in other interactions

- Interactions to analyse
  - aromatic - cation/anion
  - aromatic - HX
  - aromatic - aromatic

### Aromatic - ion interactions

Goal: check for pyramid-like structures where M is the apex (top vertice) and the aromatic ring is the base. 

NOTE: for aromatic - HX and aromatic - aromatic the same script can be used with small code expansions

#### Steps
#### 1. Look for point charges (cations and anions)
  * `resname ARG* and name NH1`, but only if that residue has an HH11 and/or HH12 atom (this must be possible). <small>NOTE: All Arg residues have an NH1 atom; some have an HH1 H-atom bonded to that NH1 atom, while in other the same NH1 atoms is bonded to H atoms named HH11 and HH12. We only want the NH1 atoms from the second case.</small>
  * `resname ASP and name OD*`
  * `resname GLU and name OE*`
  * `resname HIS* and name ND1` but only if BOTH atoms HD1 and HE2 are present in that residue)
  * `resname HIT and name NE2`
  * `resname LYS+ and name NZ` but only if there's an HZ3 atom in that residue
  * `resname TYO abd bane OH` 
  * all metals - not sure of names, for now this functions as placeholder
    * `resname FE3 and name FE`
    * `resname V3 and name V3`
    * `resname VO and name V4`
    * `resname VO2 and vame V5`
    * `resname VO4 and name V5`

#### 2. Look for aromatic rings
For all araomatic rings, compute centroid and normal vector **(perhaps this could be stored in a dict - aromatic ring, centroid, normal, indexed by residue number ?)**
  - compute ring centroid: a point with coordinates given by average of the atom coordinates of that ring
  - vector normal to the plan of the ring, by the usual protocol:


  ```python
  # scheme for normal vectors
  aromatic_residues = MDAnalysis.Universe.select_atoms("resname PHE TRP TYR TYO HIS HSD HSE HSP HIT").residues
  # names might change

  resname = PHE OR TYR OR TYO
    aromatic_coords = positions of CG, CD1, CD2, CE1, CE2, CZ atoms (list)
    USE FUNCTIONS FOR SIX POINTS
  
  resname = HIS HSD HSE HSP HIT
    aromatic_coords = positions of CG, ND1, CE1, NE2, CD2 atoms (list)
    USE FUNCIONST FOR FIVE POINTS
  
  resname = TRP
    aromatic_coordsA = positions of CG, CD1, NE1, CE2, CD2
    USE FUNCIONST FOR FIVE POINTS
    aromatic_coordsB = positions of CD2, CE3, CZ3, CH2, CZ2, CE2
    USE FUNCTIONS FOR SIX POINTS
    aromatic_coordsC = positions of CG, CD1, NE1, CE2, CD2, CE3, CZ3, CH2, CZ2
    USE FUNCTIONS FOR NINE POINTS

  FOR SIX POINTS:
    vector1=aromatic_coords[0]-aromatic_coords[5]
    vector2=aromatic_coords[1]-aromatic_coords[4]
    vector3=aromatic_coords[2]-aromatic_coords[3]
    prod1=numpy.cross(vector1,vector2)
    prod2=numpy.cross(vector1,vector3)
    prod3=numpy.cross(vector2,vector3)
    normalVector=(prod1+prod2+prod3)/3    

  FOR FIVE POINTS:
    vector1=aromatic_coords[1]-aromatic_coords[3]
    vector2=aromatic_coords[0]-aromatic_coords[2]
    vector3=aromatic_coords[0]-((aromatic_coords[2]+aromatic_coords[3])/2)
    prod1=numpy.cross(vector1,vector2)
    prod2=numpy.cross(vector1,vector3)
    prod3=numpy.cross(vector2,vector3)
    normalVector=(prod1+prod2+prod3)/3
  
  FOR NINE POINTS:
    vector1=piCoords1[3]-piCoords1[7]
    vector2=piCoords1[0]-piCoords1[8]
    vector3=piCoords1[1]-((piCoords1[7]+piCoords1[8])*2)
    prod1=numpy.cross(vector1,vector2)
    prod2=numpy.cross(vector1,vector3)
    prod3=numpy.cross(vector2,vector3)
    normalVector=(prod1+prod2+prod3)/3
  ```

#### 3. Check that centroid - point charge distance is within cutt-off
#### 4. Compute angular displacement


```python
for each aromatic_residue
  for each point_charge
    if numpy.linalg.norm(centroid-point_charge) < cutoff
    #for now, cutoff = 6.0 angstrom
      vector_aux = (centroid - point_charge)
      angular displacement = angle between vector normal to aromatic residue and vector_aux using the dot product formula

      output = distance, angular displacement
```

### Aromatic - H-X interactions
#### 1. Use same aromatics
  * use normal vectors
  * use centroids
#### 2. Define set of HX groups anf HC groups
```python
get residues by the following resnames IF they have the mentioned H atoms (we will need the donor atom)

HX(O/N) groups #these are polar H and explicitly present 
  O donors
    resname TYR and name HH # from OH - donor atom
    resname THR and name HG1 # from OG1
    resname SER and name HG # from OG
    resname ASP and name HD2 # from OD2
    resname GLU and name HE2 # from OE2
  N donors
    resname LYS and name HZ1 # from NZ
    resname LYS and name HZ2 # from NZ
    resname LYS and name HZ3 # from NZ
    resname HIS and name HE2 # from NE2
    resname HIS and name HD1 # from ND1
    resname GLN and name HE21 # from NE2
    resname GLN and name HE22 # from NE2
    resname ASN and name HD21 # from ND2
    resname ASN and name HD22 # from ND2 
    resname ARG and name HH11 # from NH1
    resname ARG and name HH1 # from NH1
    resname ARG and name HH12 # from NH1
    resname ARG and name HH21 # from NH2
    resname ARG and name HH22 # from NH2
HS groups:
  S donors
    resname CYS and namhe HG # from SG
HCS groups:


HC groups #these are apolar H, not present
  C donors (from CH3 groups only)
    resname MET # donor = CE; next =  SD
    resname ALA # donor = CB; next = CA
    resname LEU # donor = CG2; next = CB
    resname LEU # donor = CD; next = CG1
    resname THR # donor = CG2; next = CB
    resname VAL # donor = CG1; next = CB
    resname VAL # donor = CG2; next = CB


```
#### 3. Work out the distances and the angles
```python
for HO and HN groups:
# Stojanovic2007; Plevin2010; 
  distance: centroid to donor
  if distance < 4.0:
    vector_aux1 = centroid - H atom
    vector_aux2 = donor - H atom
    alpha = angle pi-H-X; between vector_aux1 and vector_aux2
    if alpha > 120:
      beta = angle between normal vector and centroid-donor vector
      if beta < 50:
        output

for HC groups:
# Stojanovic2007; Plevin2010; 
  ditance = centroid to donor
  if 3.0 < distance < 4.5:
    beta = angle between normal vector and centroid-donor vector
    if beta < 50:
      gamma = angle between normal vector and donor-next vector
      output
```

## Katie's One Drive CysMet and CysHBond CYS-HBOND
for HS groups:
  same as for H-O/N but 3.5 < distance < 4.9
  #arene paper - also requires H to edge   
for HCS groups:
  check paper; requires also H to edge


### Aromatic - Aromatic interactions

The set of aromatic residues is computed, together with the normal vector to each plane and the centroid of each plane. 

Classification of aromatic interaction can be performed from the geometrica if parameters computed for each i,j aromatic pair:
- $d$ - centroid-centroid distance
- $R_1$ - horizontal displacement
- $R_2$ - vertical displacement
-  $\alpha$ - dihedral angle (angle between planes = angular between normal vectors)
-  $\theta$ - angular displacement (angle between $\overrightarrow{\text{CC'}}$ and normal to plane containing centroid C)

Taxonomy:
- $$d < 5.0 \wedge \alpha \in \left[\pi,\frac{\pi}{6}\right]\cup\left[\frac{5\pi}{6},\pi\right]$$
  * parallel
  * delusional classification: $V_\%$, $H_\%$ are vertical and horizontal character, which give what the largest contribution to centroid offset 
 $$V_\%=\frac{\left|R_2\right|}{d};   H_\% = \frac{\left|R_1\right|}{d}$$

- $$\alpha \in \left[\frac{\pi}{6},\frac{5\pi}{6}\right]$$
  - $R_2<3.5$ - T-orientation, edge-to-face
  - $R_2 \geq 3.5 \wedge R_1 < 3.0$ - T-orientation, face-to-edge
  - $R_2 \geq 3.5 \wedge R_1 \geq 3.0$ - L-orientation





# Geometric snippets
## Angle between planes
Given the normal vectors of two planes, v1=(a1,b1,c1) and v2=(a2,b2,c2), the angle between the planes is the angle between the vectors, and can be computed from the dot product:
```python
cosine = dot(v1,v2) / (norm(v1) * norm(v2))
angle = arccos(cosine)
```


## Distance between parallel planes

Given a plane alpha, defined by a normal vector (a,b,c) and with a centroid (x1,y1,z1): 

- Calculate the distance to another plane that contains point (centroid) (x2,y2,z2) defined by vector (e,f,g) [paralell => are multiples] using:
  
$$D=\frac{|e \times x_1 + f \times y_1 + g \times z_1 + h|}{\sqrt(d^2+e^2+f^2) }$$

```python
#define general equation of plane TWO: Ex+Fy+GZ+H=0
h = -e  * x2 - f * y2 - g * z2 

#use this equation with point of plane ONE
dist =  ABS(e*x1 + f*y1 + g*z1 + h) / SQRT (e^2 + f^2 + g^2)
```

## Distance of a point to a plane - projection of point to plane
Given a point P(x0, y0, z0), find the (minium) distance to a plane defined by a centroid C1(x1,y1,z1) and normal vector (a,b,c). This ammounts to determine point P', the projection of P onto the plane. Line PP' goes through P, P' and is perpendicular to the plane => the vector normal to the plane is a directing vector of the line.

The parametric equations of this line are:

$$ \frac{x-x_0}{a} = \frac {y-y_0}{b}=\frac{z-z_0}{c}=t \Longleftrightarrow x = at+x_0 \wedge y=bt+y_0 \wedge z=ct+z_0$$

Substituting this on the plane equation we get the value for t, which substituing back in the line equation gives the coordinates of point P' (the orthogonal projection of P on the plane).

```python
#define general equation of plane: ax+by+cz+d=0
d = -a  * x1 - b * y1 - c * z1 

#get the parametric equation of the line
#replacing t for the value gives P'(x2,y2,z2)
x2 = a*t + x0
y2 = b*t + y0
z2 = c*t + z0
 
#compute t by substituing (x,y,z) on the plane equation
#a(at+x0)+b(bt+y0)+c(ct+z0)+d=0
t = (-d - a*x0 - b*y0 - c*z0)/(a^2 + b^2 + c^2)
```
The distance from P to P' is the "vertical displacement" of the centroids, and the distance from P to C is the "horizontal displacement" of the centroids.

```python
#compute distance of P to P' = distance of point to plane
distance = ( (x2 - x0)^2 + (y2 - y0)^2 + (z2 - z0)^2 ) ^ (1/2)
#OR
pointP0 = [x0,y0,z0]
pointP1 = [x2,y2,z2]
distance = linalg.norm(pointP0 - pointP1)
```