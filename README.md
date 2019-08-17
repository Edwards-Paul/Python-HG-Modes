# Python-HG-Modes
Modeling HG modes propagation along z given modes, optical parameters, and XY plane of calculation.

PauLisa.py contains all code for HG mode calculation, which is imported to Jupyter Notebook(JN).

## Starting
Import PauLisa.py to JN.

User defines a *Params* object of optical parameters (e.g., beam waist, q-parameter, wavelength) or uses PauLisa.defaultParams.

User defines an XY plane object, *Plane*, for calculations (i.e., x/y min/max values, x/y # of points) or uses PauLisa.defaultPlane.

User defines HG modes using nonzero coefficients(c) and HG modes(n, m) as (n, m, c).

The *Params*, *Plane*, and *Modes* are passed to *Calculate* from PauLisa.py for in-plane phase and amplitude/intensity calculations. 

*Calculate* generates a *Result* object which holds all data necessary for plotting and analysis. Alternatively, *Amplitude* and *Phase* compute respective values from x-y-z coordinates.


## Function Usage:
### OPTICAL PARAMETERS DEFINITION     
 PARAMETERS=PauLisa.Params(wavelength, w0, z0)

### XY PLANE OF CALCULATION DEFINITION     
 PLANE=PauLisa.Plane(xmin, xmax, xpoints, ymin, ymax, ypoints) 

### MODES DEFINITION AND DISPLAY     
 MODES=PauLisa.Modes((n1, m1, c1),(n2, m2, c2))     
 PauLisa.ShowModes(MODES) 

### AMPLITUDE CALCULATIONS     
 Calculate amplitude over plane: RESULT=PauLisa.Calculate(PARAMS, PLANE, MODES, z)
 
 Simple calculation from coordinates: PauLisa.Amplitude(PARAMS, x, y, z, MODES) 

### INTENSITY PLOTTING     
 PauLisa.Contour(RESULT)     
 PauLisa.IntensitySliceX(y, *RESULT)     
 PauLisa.IntensitySliceY(x, *RESULT) 

### PHASE CALCULATION     
 PauLisa.Phase(PARAMS, x, y, z, MODES) 

### PHASE PLOTTING     
 PauLisa.PhaseContour(RESULT)     
 PauLisa.PhaseSliceX(y, *RESULT)     
 PauLisa.PhaseSliceY(x, *RESULT) 
 
 #### **VARNAME represents a variable number of args of that type.*
 
## Ex. Run
All-caps are example user-defined functions.

import PauLisa as pl

params = pl.Params(1.064e-6, 0.001, 0)


plane = pl.Plane(-0.06, 0.06, 0.0001, -0.03, 0.03, 0.0001)


modes=pl.Modes((1,0,1), (0,1,1))

pl.ShowModes(modes)


f = pl.Calculate(params, plane, modes, 100)


pl.IntensitySliceX(f, -0.01)

pl.IntensitySliceY(f, 0.01)

pl.Contour(f)
