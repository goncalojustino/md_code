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
      cosine = numpy.dot(v1, v2) / (numpy.linalg.norm(v1) - numpy.linalg.norm(v2))
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
