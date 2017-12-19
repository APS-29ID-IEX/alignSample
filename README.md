# alignSample
Python code for automated sample alignment at 29IDD RSXS (resonant soft x-ray 
scattering) system.

## alignRSXS.py

alignRSXS.py contains subroutines to use the BlueSky interface for the alignment 
process

### Usage

#### Sample Alignment
Iterates between x- and theta-alignment to improve sample position

```
alignRSXS(iterations = 3, alignCriteria = None, theta_kws = None, x_kws = None, 
          SimAlign = False, verbose = False)
```

#### Testing
Both translational and rotational alignment tools have a testing mode (SimAlign):

```
alignRSXS(SimAlign = True)
```
    
### Classes

#### SynErf

```
    Evaluate a point on the error function based on the value of a motor.

    Parameters
    ----------
    name            : string
    motor           : Device
    motor_field     : string
    Imax            : float, max intensity of peak
    wid             : float, optional; rough width of error function
    x0              : float, center of error function
    noise           : {'poisson', 'uniform', None}
                      Add noise to the gaussian peak.
    noise_multiplier: float; Only relevant for 'uniform' noise. Multiply the 
                       random amount of noise by 'noise_multiplier'
```

#### SynErfGauss

```
    Evaluate a point on the error function based on the value of a motor
    that also includes a small gaussian bump corresponding to the beam passing 
    under the sample holder which occurs for high x (>1 mm after sample/holder
    block beam).

    Parameters
    ----------
    name : string
    motor : Device
    motor_field : string
    Imax : number
        max intensity of peak
    wid : number, optional
        rough width of error function
    x0 : number
        center of error function
    noise : {'poisson', 'uniform', None}
        Add noise to the gaussian peak.
    noise_multiplier : float
        Only relevant for 'uniform' noise. Multiply the random amount of
        noise by 'noise_multiplier'
```

#### CurAmp
```
    ophyd device for current amplifier that requires scan state to be 
    set to passive to be triggered properly.  Sets back to original
    state after run.
    
    Parameters same as for Ophyd device
```

### Functions

#### cleanLegend
```
    Replaces messy legends produced by lmfit
    
    Parameters
    ----------
    iteration   : integer, iteration of alignment loop
    ax          : axis object whose legend is being repaced
    other       : list of dicts for other items in legend containing
                    label   :   label to use in legend
                    shape   :   shape to use in legend entry
                    color   :   color of legend entry
                    size    :   size of legend entry symbol
```
                     
#### cleanAxes
```
    Sets up multi-axis plot where each plot is a different color
    
    Parameters
    ----------
    figPlot         : Plot object
    line1, line2    : line objects for the two sets of data being displayed
    axis1, axis2    : axis objects for the two sets of data being displayed
    color1, color2  : colors of the line/axis objects
    ylabel1, ylabel2: labels for the data sets
    xlabel          : x axis label (same for both data sets)

```

#### gaussian
```
    Create gaussian for fitting and simulation
    
    Parameters
    ----------
    x       :   input upon which gaussian is evaluated
    A       :   peak of gaussian
    sigma   :   "spread" of gaussian
    x0      :   location of peak of gaussian
```

#### erfx
```
    Create error function for fitting and simulation
    
    Parameters
    ----------
    x       :   input upon which error function is evaluated
    low     :   min value of error function
    high    :   max value of error function
    sigma   :   "spread" of error function transition region
    x0      :   location of error function's "center"
```

#### color_y_axis

```
    Colors the passed axis.

    Parameters
    ----------
    ax      :   axes object whose color is to be changed
    color   :   new color of axes object
```


#### fitAndPlotBumplessData 
```
    fitAndPlotBumplessData:  
        -reduces x-scan data to eliminate bump
        -fits reduced data
        -plots (if axes supplied) fit line
        -returns center of fit's transition
        
    Parameters
    ----------
    y               : y data from scan
    x               : x data from scan
    x0              : estimated half-max from scan
    width           : estimated width of transition from scan
    fudge           : Fudge factor for x-scan data reduction to remove bump 
                      from post-scan fit (in mm).  Defaults to def_fudge
    ax              : axis object for plotting fit, 
                      Default value is None which suppresses plot
    color           : color of fit line. Defaults to red.
    linestyle       : string descriptor of linetype for fit line, 
                      Default value is solid line
```

#### xScanAndCenter
```
    align_x:  
        -scan xMotor through range
        -monitors detector looking for transition of beam/no-beam   

    Parameters
    ----------
    iteration           : overall iteration of full alignment loop, required for 
                          plot colors/legend
    detector            : detector object
    detector_name       : detector object's name
    motor               : X motor object
    motor_name          : X motor object's name
    alignRange          : list of two floats, describing range of x-alignment 
                          Default value is [-1,2] mm,
    StepRange           : list of two floats, minimum and maximum steps for 
                          adaptive scan. Default value is [0.015, 0.075] degrees 
                          Minimum thMotor step is 0.01 degrees
    targetDelta         : float, maximum jump in detector value before scan step 
                          decreased, default is set to def_xTargetDelta
                          value is 10 in detector units, 
    plotter             : axis object for x alignment data, 
                          Default value is None which later creates new object
    SimAlign            : boolean flag for simulation/testing 
                          Default value is False,
    fudge               : Fudge factor for x-scan data reduction to remove bump 
                          from post-scan fit (in mm).  Defaults to def_fudge
    verbose             : boolean for showing alignment process using LiveTable
                          callbacks. Default value is False
```

#### thetaScanAndCenter
```
    align_theta:  
        -sets tthMotor to 10 degrees
        -scan thMotor through range
        -monitors detector looking for peak
        -move thMotor to peak
        -set thMotor set point to 5 degrees
        -move both motors back to zero
        
    Parameters
    ----------
    iteration           : overall iteration of full alignment loop, required for 
                          plot colors/legend
    detector            : detector object
    detector_name       : detector object's name
    motor               : theta motor object
    motor_name          : theta motor object's name
    motor2              : 2theta motor object
    motor2_name         : 2theta motor object's name
    alignRange          : list of two floats, 
                          Default value is [4.0,6.0] degrees,
    coarsePTS           : number of measurements in coarse scan,
                          defaults to 10
    fineRadius          : float, distance from coarse peak defining region of 
                          fine scan.  Default value is 0.3 degrees,
    fStepRange          : list of two floats, minimum and maximum steps for 
                          adaptive scan. Default value is [0.015, 0.075] degrees 
                          Minimum thMotor step is 0.01 degrees
    targetDeltaFactor   : float, maximum jump in detector value before scan step 
                          decreased is determined from coarse scan's detector
                          peak value divided by the targetDeltaFactor. Default 
                          value is 10 in detector units, 
    plotter             : axis object for theta alignment data, 
                          Default value is None which later creates new object
    det4offset          : Detector's offset from 0 (for 2theta motor), defaults
                          to None leading to attempt to get it from PV
                          29idd:userCalcOut2.VAL
    SimAlign            : boolean flag for simulation/testing 
                          Default value is False,
    verbose             : boolean for showing alignment process using LiveTable
                          callbacks. Default value is False
```

#### multiScan
```
    multiScan:  
        -set tthMotor to 0 degrees
        -set thMotor to 0 degrees
        -scan x-direction for transition from no-beam to full-beam
        -scan theta direction with detector at 2theta = 10 degrees
        -plot summary (if flagged) of x(half-max), theta(peak)

    Parameters
    ----------
    detX                : X detector object (in real alignment, same as detTh)
    detX_name           : X detector object's name (in real alignment, same as 
                          detTh_name)
    detTh               : Theta detector object (in real alignment, same as 
                          detX)
    detTh_name          : Theta detector object's name (in real alignment, same 
                          as detX_name)
    motorX              : X motor object
    motorX_name         : X motor object's name
    motorTh             : Theta motor object
    mTh_name            : Theta motor object's name
    motorTTh            : 2theta motor object
    mTTh_name           : 2theta motor object's name
    iterations          : Number of iterations of the X, Theta alignment scans
                          to run through
    alignCriteria       : dict of criteria for ending alignment
                          Default value is None,
    theta_kws           : dict for arguments to be passed to theta-alignment
                          Default value is None,
    x_kws               : dict for arguments to be passed to x-alignment
                          Default value is None

    xPlot               : axes object for x-alignment plot
    thPlot              : axes object for theta-alignment plot
    sumPlot             : plot object for usmmary plot
    sumAxes             : axes for the x-/theta-alignment summaries
    SimAlign            : boolean flag for simulation/testing 
                          Default value is False,
    verbose             : boolean for showing alignment process using LiveTable
                          callbacks. Default value is False
```                       
                          
#### alignRSXS     
```
    alignRSXS:
        -sets up detector and motor objects
        -sets up summary plots
        -calls up multiScan
        
    Parameters
    ----------
    iterations      : integer number of iterations to perform [unless alignment
                      criteria is not None and met]. Default value is 3,
    alignCriteria   : dict of criteria for ending alignment
                      Default value is None,
    theta_kws       : dict for arguments to be passed to theta-alignment -- most
                      important is det4offset 
                      Default value is None,
    x_kws           : dict for arguments to be passed to x-alignment
                      Default value is None
    SimAlign        : boolean for simulated run (not implemented for full loop
                      only for individual x, theta loops via the keyword args)
    verbose         : boolean for showing alignment process using LiveTable
                      callbacks. Default value is False
    
```

## Sample Alignment stuff.ipynb

Jupyter notebook containing early tests of code snippets later incorporated in alignSample.py (and later alignRSXS.py).

## Simple Decorators on Scans.ipynb

Jupyter notebook containing tests of decorator use for callbacks in alignRSXS.py
