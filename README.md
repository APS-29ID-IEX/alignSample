# alignSample
Python code for automated sample alignment at 29IDD RSXS (resonant soft x-ray scattering) system.

## alignSample.py

alignSample.py contains subroutines to use the BlueSky interface for the alignment process

### Usage

#### Translational alignment
Move x-stage while taking measurements from detector. As the sample/sample holder pass through beam, detector signal should go to zero.  If it goes further, the beam will pass "under" the sample holder and signal could be observed again. 

'''
align_x(iterations, theta_offset = 0, alignRange = [-1.00,2.00], stepRange = [0.0015, 0.1], targetDelta = 0.0025, SimAlign = False, ax = None, bump = False, fudge = 0.75, verbose = False)
'''

#### Rotational alignment
Rotate theta state about 5 degrees to maximize signal at detector positioned a 10 degrees. The minimum theta motor step is 0.01 degrees.

'''
align_theta(iterations, x0 = 0, alignRange = [4.0,6.0], coarseStep = 0.1, fineRadius = 0.2, 	stepRange = [0.015, 0.1], targetDelta = 0.1, SimAlign = False, ax = None, verbose = False):
'''

#### Sample Alignment
Iterates between x- and theta-alignment to improve sample position

'''
align_sample(iterations = 3, alignCriteria = None, theta_kws = None, x_kws = None, SimAlign = False, plotReport = True)
'''

#### Testing
Both translational and rotational alignment tools have a testing mode (SimAlign):

'''
align_sample(theta_kws = {'SimAlign':True}, x_kws = {'SimAlign':True, 'bump':True,'fudge':0.75}, iterations
'''
    
### Classes

#### SynErf

'''
    Evaluate a point on the error function based on the value of a motor.

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
'''

#### SynErfGauss

'''
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
'''

### Functions

#### align_legend

'''
    Replaces messy legends produced by lmfit

    Parameters
    ----------
    iteration	: integer, iteration of alignment loop
    ax		: axis object whose legend is being repaced
    other	: list of dicts for other items in legend containing
			label	:   label to use in legend
			shape	:   shape to use in legend entry
			color	:   color of legend entry
			size	:   size of legend entry symbol
'''
					 
#### align_sample

'''
	Loops over align_x and align_theta until criteria met or fixed number 
	 of iterations met
		
	Parameters
        ----------
        iterations 		: integer number of iterations to perform [unless alignment
					  criteria is not None and met]. Default value is 3,
	alignCriteria	        : dict of criteria for ending alignment
					  Default value is None,
	theta_kws 		: dict for arguments to be passed to theta-alignment
					  Default value is None,
	x_kws 			: dict for arguments to be passed to x-alignment
					  Default value is None
	SimAlign		: boolean for simulated run (not implemented for full loop
					  only for individual x, theta loops via the keyword args)
	plotReport 		: boolean for plotting summary results 
					  Default value is True,
'''

#### align_theta

''' 
	sets tthMotor to 10 degrees
	set xMotor to x0 (half-max)
	scanthMotor through range
	monitors detector looking for peak

	Parameters
        ----------
	iterations		: overal iteration of full alignment loop, required for 
					  plot colors/legend
        x0			: float,  
					  Default value is 0 in mm,
	alignRange		: list of two floats, 
					  Default value is [4.0,6.0] degrees,
	coarseStep 		: float,  
					  Default value is 0.1 degrees, 
	fineRadius 		: float,  
					  Default value is 0.2 degrees,
	stepRange 		: list of two floats, 
					  Default value is [0.015, 0.1] degrees. 
					  Minimum thMotor step is 0.01 degrees
	targetDelta		: float,  
					  Default value is 0.1 in detector units, 
	SimAlign  		: boolean, 
					  Default value is False,
	ax  			: axis object for theta alignment data, 
					  Default value is None,
	verbose 		: boolean for showing alignment process
					  Default value is False):
'''

#### align_x


'''
	set tthMotor to 0 degrees
	setthMotor to 0 degrees
	scan xMotor through range
	monitors detector looking for transition of beam/no-beam	

	Parameters
        ----------
	iterations		: overal iteration of full alignment loop, required for 
					  plot colors/legend
	theta_offset		: float,  
					  Default value is 0 in degrees.
	alignRange 		: list of two floats, 
					  Default value is [-1.00,2.00] in mm.
	stepRange 		: list of two floats, 
					  Default value is [0.0015, 0.1]. 
					  Miniumum xMotor step size is 1.5 um
	targetDelta 		: float, 
					  Default value is 0.0025.
	SimAlign 		: boolean, 
					  Default value is False.
	ax 				: axis object for plotting x-alignment data, 
					  Default value is None.
	bump 			: boolean, 
					  Default value is False.
	fudge 			: float, rough estimate of sample/sample plate thickness in mm
					  Default value is 0.75.
	verbose 		: boolean, 
					  Default value is False
'''

#### color_y_axis

'''
	Colors the passed axis.

	Parameters
	----------
	ax		:	axes object whose color is to be changed
	color	:	new color of axes object
'''

#### erfx

'''
	Create error function for fitting and simulation
	
	Parameters
	----------
	x		:	input upon which error function is evaluated
	low		: 	min value of error function
	high		: 	max value of error function
	sigma		:	"spread" of error function transition region
	x0		:	location of error function's "center"
'''

#### gaussian

'''
	Create gaussian for fitting and simulation
	
	Parameters
	----------
	x		:	input upon which gaussian is evaluated
	A		: 	peak of gaussian
	sigma	:	"spread" of gaussian
	x0		:	location of peak of gaussian
'''

## Sample Alignment stuff.ipynb

Jupyter notebook containing early tests of code snippets later incorporated in alignSample.py.
