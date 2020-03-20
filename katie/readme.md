# Catarina Nascimento

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
    # all angles must be between 0 and 180
    
    # take the largest value (alpha) and second largest (beta) of amplitudes [no id req.]
    tau_six = (360 - alpha - beta) / 180
    tau_six_abs = (360 - alpha - beta) / 180
    tau_siz_normf = [360 - DEGREES(ARCCOS(5/7)) ]  / 180 # ISTO PERCEBE-SE ?
    tau_six = tau_six_abs / tau_six_normf

~~~~

