import matplotlib.pyplot as plt
plt.ion()

# Make plots live-update while scans run.
from bluesky.utils import install_qt_kicker
install_qt_kicker()
 
#Start BlueSky
from bluesky.plans import adaptive_scan, scan
from bluesky.callbacks import LiveFit, LivePlot, LiveFitPlot, LiveTable
from bluesky import RunEngine
from ophyd.signal import EpicsSignalRO, EpicsSignal
from ophyd.epics_motor import EpicsMotor
#Need for testing/troubleshooting
from ophyd.sim import SynAxis, SynGauss, SynSignal

#For fitting
import lmfit
import numpy as np
import scipy
import math 

# Change color of each axis
def color_y_axis(ax, color):
    """Color your axes."""
    for t in ax.get_yticklabels():
        t.set_color(color)
    return None

def gaussian(x, A, sigma, x0):
    return A*np.exp(-(x - x0)**2/(2 * sigma**2))
    
def erfx(x, low, high, wid, x0):
    return (high - low) * 0.5 * (1-scipy.special.erf((x-x0)/wid)) + low    

#SynErf class required to simulate translational/x-alignment
class SynErf(SynSignal):
    """
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

    Example
    -------
    motor = SynAxis(name='motor')
    det = SynGauss('det', motor, 'motor', center=0, Imax=1, sigma=1)
    """
    
    def __init__(self, name, motor, motor_field, Imax, wid, x0,
                 noise=None, noise_multiplier=1, **kwargs):
        if noise not in ('poisson', 'uniform', None):
            raise ValueError("noise must be one of 'poisson', 'uniform', None")
        self._motor = motor

        def func():
            m = motor.read()[motor_field]['value']
#            v = Imax * 0.5 * (1-math.erf(wid*(m-x0)))
            v = erfx(m, 0, Imax, wid, x0)
            if noise == 'poisson':
                v = int(np.random.poisson(np.round(v), 1))
            elif noise == 'uniform':
                v += np.random.uniform(-1, 1) * noise_multiplier
            return v

        super().__init__(func=func, name=name, **kwargs)

#SynErfGauss class required to simulate translational/x-alignment
class SynErfGauss(SynSignal):
    """
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

    Example
    -------
    motor = SynAxis(name='motor')
    det = SynGauss('det', motor, 'motor', center=0, Imax=1, sigma=1)
    """
    
    def __init__(self, name, motor, motor_field, Imax, wid, x0,
                 noise=None, noise_multiplier=1, **kwargs):
        if noise not in ('poisson', 'uniform', None):
            raise ValueError("noise must be one of 'poisson', 'uniform', None")
        self._motor = motor

        def func():
            m = motor.read()[motor_field]['value']
#            v = Imax * 0.5 * (1-math.erf(wid*(m-x0)))
            v = erfx(m, 0, Imax, wid, x0)
            v += gaussian(m, Imax*0.20, 0.1, x0 + 0.75 + 0.25*np.random.random())
            if noise == 'poisson':
                v = int(np.random.poisson(np.round(v), 1))
            elif noise == 'uniform':
                v += np.random.uniform(-1, 1) * noise_multiplier
            return v

        super().__init__(func=func, name=name, **kwargs)

def align_x(theta_offset = 0,
			alignRange = [-1.00,2.00],
			stepRange = [0.0015, 0.1],
			targetDelta = 0.0025,
			SimAlign = False,
			ax = None,
			bump = False,
			fudge = 0.5,
			verbose = False):
	"""
	align_x:  
		-set tthMotor to 0 degrees
		-setthMotor to 0 degrees
		-scan xMotor through range
		-monitors detector looking for transition of beam/no-beam	

	Parameters
    ----------
    theta_offset	: float,  
					  Default value is 0 in degrees.
	alignRange 		: list of two floats, 
					  Default value is [-1.00,2.00] in mm.
	stepRange 		: list of two floats, 
					  Default value is [0.0015, 0.1]. 
					  Miniumum xMotor step size is 1.5 um
	targetDelta 	: float, 
					  Default value is 0.0025.
	SimAlign 		: boolean, 
					  Default value is False.
	ax 				: axis object, 
					  Default value is None.
	bump 			: boolean, 
					  Default value is False.
	fudge 			: float, 
					  Default value is 0.5.
	verbose 		: boolean, 
					  Default value is False

	"""
	
	if not SimAlign:
		#set signals
		xMotor = EpicsMotor('29idd:m1', name = 'xMotor')
		thMotor = EpicsMotor('29idd:m7', name = 'thMotor')
		tthMotor = EpicsMotor('29idHydra:m1', name = 'tthMotor')
		detX =  EpicsSignalRO('29idd:ca4:read', name = 'detX')
	else:
		#set simulated signals
		xMotor = SynAxis(name='xMotor')
		simX0 = np.random.random()-0.5
		if verbose:
			print('-------------------------------')
			print('Simulated 1/2 max at ', simX0)
		if bump:
			detX = SynErfGauss('detX', xMotor, 'xMotor', 
						  Imax=0.024, wid=0.20, x0 = simX0,
						  noise='uniform', noise_multiplier=0.0005)
		else:
			detX = SynErf('detX', xMotor, 'xMotor', 
						  Imax=0.024, wid=0.20, x0 = simX0,
						  noise='uniform', noise_multiplier=0.0005)
	
	if ax is None:
		fig, ax = plt.subplots()  # explicitly create figure, axes to use below

	#Start BlueSky run engine
	RE = RunEngine({})

#	trans_model = lmfit.Model(erfx)
	comp_model = lmfit.Model(erfx, prefix="erf_")+lmfit.Model(gaussian, prefix="gau_")
	
	#set up fits

#	init_guess = {'low': 0,
#				  'high': 0.03,
#				  'wid': lmfit.Parameter('wid', value = 0.4, min=0),
#				  'x0': -0.1}
#	trans_lf = LiveFit(trans_model, 'detX', {'x': 'xMotor'}, 
#						init_guess, update_every=5)
	comp_guess = {'erf_low': 0,
				  'erf_high': 0.03,
				  'erf_wid': lmfit.Parameter('erf_wid', value = 0.4, min=0),
				  'erf_x0': -0.1,
				  'gau_A': 0.01,
				  'gau_sigma': lmfit.Parameter('gau_sigma', value = .1, min=0),
				  'gau_x0': 1.0}
	comp_lf = LiveFit(comp_model, 'detX', {'x': 'xMotor'}, 
						comp_guess, update_every=5)

#	trans_lfp = LiveFitPlot(trans_lf, ax=ax, color='r', label = 'fit')
	comp_lfp = LiveFitPlot(comp_lf, ax=ax, color='g', label = 'composite fit')

	trans_lp = LivePlot('detX',x='xMotor', ax=ax, marker='o', 
						linestyle='none', label = 'data')

	#Move motors 
	if not SimAlign:
		#movethMotor to theta_offset or 0 degrees
		thMotor.move(theta_offset)
		#move tthMotor to 0 degrees
		tthMotor.move(0.00)
		pass
		
	#run adaptive scan
	RE(adaptive_scan([detX], 'detX', xMotor, 
					start = alignRange[0],
					stop = alignRange[1], 
					min_step = stepRange[0],  
					max_step = stepRange[1], 
					target_delta = targetDelta,
					backstep = True), [trans_lp, comp_lfp])
					
	#Use composite fit data and re-fit using interval of alignRange[0] to 
	#compFit_x0+(1+fudge)*compFit_wid

	#Reduce data
	compFit_x0 = comp_lf.result.params['erf_x0'].value
	compFit_width = comp_lf.result.params['erf_wid'].value
	
	x_reduced_max = compFit_x0 + (1.0+fudge)*compFit_width
	
	xdata = np.asarray(trans_lp.x_data)
	ydata = np.asarray(trans_lp.y_data)
	ix = np.where(xdata < x_reduced_max)
	
	x_reduced = xdata[ix]
	y_reduced = ydata[ix]
		
	#fit reduced data, get center and width
	red_model = lmfit.Model(erfx, missing = 'drop')
	red_guess = {'low': min(y_reduced),
				 'high': max(y_reduced),
				 'wid': lmfit.Parameter('wid', value = compFit_width, min=0),
				 'x0': compFit_x0}
		
	params = red_model.make_params(low = min(y_reduced), high = max(y_reduced), 
			                       wid = compFit_width, x0 = compFit_x0)
	
	redFit = red_model.fit(y_reduced, params, x = x_reduced)
	redFit_x0 = redFit.result.params['x0'].value
	redFit_width = redFit.result.params['wid'].value

	#plot new fit
	redFit.plot_fit(ax = ax, data_kws = {'marker': 'x'}, fit_kws={'color':'b'},
					xlabel = 'xMotor', ylabel = 'detX')
	ax.set_title(' ')
			
	if verbose:
		print('-------------------------------')
		print('Number of measurements: ',len(trans_lp.x_data))	
		print('-------------------------------')
		#fit_x0 = trans_lf.result.params['x0'].value
		#fit_width = trans_lf.result.params['wid'].value
		#print('Erf(x) Model results')	
		#print('Half-maximum at x0 = {0:.5f}'.format(fit_x0))
		#print('Width of eft = {0:.5f}'.format(fit_width))
		#if SimAlign:
			#fitDelta = abs(fit_x0 - simX0)
			#print('Delta of {0:.1f} um'.format(1000.0*fitDelta))
			#print('Delta of {0:.1f} min resolution steps'.format(fitDelta/0.0025))
		#print('-------------------------------')
		print('Composite Model results:')
		print('Half-maximum at x0 = ',"{0:.5f}".format(compFit_x0))
		print('Width of eft = {0:.5f}'.format(compFit_width))
		if SimAlign:
			compDelta = abs(compFit_x0 - simX0)		
			print('Delta of {0:.1f} um'.format(1000.0*compDelta))
			print('Delta of {0:.1f} min resolution steps'.format(compDelta/0.0025))
		print('-------------------------------')
		print('Re-fit results:')
		print('Half-maximum at x0 = ',"{0:.5f}".format(redFit_x0))
		print('Width of eft = {0:.5f}'.format(redFit_width))
		if SimAlign:
			redDelta = abs(redFit_x0 - simX0)		
			print('Delta of {0:.1f} um'.format(1000.0*redDelta))
			print('Delta of {0:.1f} min resolution steps'.format(redDelta/0.0025))
	
	return redFit_x0
	
def align_theta(x0 = 0,
				alignRange = [4.0,6.0],
				coarseStep = 0.1, 
				fineRadius = 0.2,
				stepRange = [0.015, 0.1],	#min thMotor step is 0.01 degrees
				targetDelta = 0.1, 
				SimAlign = False,
				ax = None,
				verbose = False):
	"""
	align_theta:  
		-sets tthMotor to 10 degrees
		-set xMotor to x0 (half-max)
		-scanthMotor through range
		-monitors detector looking for peak

	Parameters
    ----------
    x0 				: float,  
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
	ax  			: plot object, 
					  Default value is None,
	verbose 		: boolean, 
					  Default value is False):
        
 	"""
	
	if not SimAlign:
		#set signals
		xMotor = EpicsMotor('29idd:m1', name = 'xMotor')
		thMotor = EpicsMotor('29idd:m7', name = 'thMotor')
		tthMotor = EpicsMotor('29idHydra:m1', name = 'tthMotor')
		detTh =  EpicsSignalRO('29idd:ca4:read', name = 'detTh')

	else:
		#set simulated signals
		thMotor = SynAxis(name='thMotor')
		simCenter = (alignRange[1]-alignRange[0])*np.random.random() + \
					alignRange[0]
		if verbose:
			print('-------------------------------')
			print('Simulated center at ',simCenter)
		detTh = SynGauss('detTh',thMotor, 'thMotor', 
						   center=simCenter, Imax=1, sigma=.03, 
						   noise='uniform', noise_multiplier=0.05)

	#Move motors 
	if not SimAlign:
		#move xMotor to x0
		
		xMotor.move(x0)
		#move tthMotor to 10 degrees
		tthMotor.move(10.0)
	
	
	#set up plot
	if ax is None:
		fig, ax = plt.subplots()  # explitly create figure, axes to use below
	
	#Start BlueSky run engine
	RE = RunEngine({})
	
	#Coarse scan to find beam
	coarsePlot = LivePlot('detTh',x='thMotor', ax=ax, marker='x', linestyle='none', label='coarse data')
	coarseSteps = math.floor((alignRange[1]-alignRange[0])/coarseStep)
	RE(scan([detTh],thMotor, alignRange[0], alignRange[1], coarseSteps), coarsePlot)
	
	coarse_th0 = coarsePlot.x_data[coarsePlot.y_data.index(max(coarsePlot.y_data))]
	if verbose:
		print('-------------------------------')
		print('Number of coarse measurements: ', coarseSteps)
		print('Coarse Peak at :', coarse_th0)
	init_guess = {'A': 1,
				  'sigma': lmfit.Parameter('sigma', .03, min=0),
				  'x0': coarse_th0}

	#set alignRange to use th0 +/- 0.1 degrees
	fineRange = [coarse_th0 - fineRadius, coarse_th0 + fineRadius]
	
	#set up fit model
	rot_model = lmfit.Model(gaussian)
		
	#set up BlueSky fit
	rot_lf = LiveFit(rot_model, 'detTh', {'x': 'thMotor'}, 
					 init_guess, update_every=5)
	rot_lfp = LiveFitPlot(rot_lf, ax=ax, color='r', label='fit')

	rot_lp = LivePlot('detTh',x='thMotor', ax=ax, marker='o', linestyle='none', label='fine data')
		
	#run adaptive scan
	RE(adaptive_scan([detTh], 'detTh',thMotor, 
					start = fineRange[0],
					stop = fineRange[1], 
					min_step = stepRange[0],  
					max_step = stepRange[1], 
					target_delta = targetDelta,
					backstep = False), [rot_lp, rot_lfp])

	if verbose:
		print('-------------------------------')
		print('Number of fine measurements: ',len(rot_lp.x_data))
		print('Maximum at theta = ',"{0:.3f}".format(rot_lf.result.params['x0'].value))
		print('-------------------------------')
		if SimAlign:
			simDelta = abs(rot_lf.result.params['x0'].value-simCenter)
			simError = simDelta/simCenter
			print('Delta of ', simDelta) 
			print('Error % of ', simError*100.0)
			print('-------------------------------')
	return rot_lf.result.params['x0'].value
	
def align_sample(iterations = 3,
				alignCriteria = None,
				theta_kws = None,
				x_kws = None,
				SimAlign = False,
				plotReport = True,
				tableReport = True,
				verbose = False):
	"""
	align_sample:
		-loops over align_x and align_theta until criteria met or fixed number 
		 of iterations met
		-other..?
		
	Parameters
    ----------
    iterations 		: integer, 
					  Default value is 3,
	alignCriteria	: dict, 
					  Default value is None,
	theta_kws 		: dict, 
					  Default value is None,
	x_kws 			: dict, 
					  Default value is None,
	SimAlign 		: boolean, 
					  Default value is False,
	plotReport 		: boolean, 
					  Default value is True,
	tableReport 	: boolean, 
					  Default value is True,
	verbose 		: boolean, 
					  Default value is False
    
	"""
				
	iter_count = 0	
	notAligned = True
	sequence = []
	
	x0 = []	
	theta0 = []
	
	fig1, ax = plt.subplots()
	fig1.suptitle('x alignment')
	
	fig2, bx = plt.subplots()
	fig2.suptitle('theta alignment')
	
	if plotReport:
		fig3, cx1 = plt.subplots()
		fig3.suptitle('Alignment summary') 
		cx2 = cx1.twinx()
	
	if x_kws is None:
		x_kws = {}
	if theta_kws is None:
		theta_kws = {}
					
	while notAligned:
		iter_count += 1
		sequence.append(iter_count)
		print('*****************************************')
		print('Starting iteration #{}'.format(iter_count))

		if iter_count == 1:
			x0.append(align_x(theta_offset = 0, ax = ax, **x_kws))
		else:
			x0.append(align_x(theta_offset = theta0[-1], ax = ax, **x_kws))

		theta0.append(align_theta(x0 = x0[-1], ax = bx, **theta_kws))
		
		if plotReport:
			cx1.plot(sequence, x0, color = 'r', marker = 'o', 
					linestyle='none', label = 'x_offset')			
			cx1.set_ylabel = 'x offset'
			cx1.set_xlabel = 'Iteration'
			cx2.plot(sequence, theta0, color = 'b', marker = 'x', 
					linestyle='none', label = 'theta_offset')
			cx2.set_ylabel = 'theta @ 2theta = 10'
			color_y_axis(cx1, 'r')
			color_y_axis(cx2, 'b')
			plt.show()
		if tableReport:
			pass

		print('Iteration #{} summary'.format(iter_count)) ## add delta_x, delta_theta
		if iter_count > 1:
			x_delta = x0[iter_count-1] - x0[iter_count-2]
			if alignCriteria is not None:
				#if (delta_x < alignCriteria['x']) AND (delta_theta < alignCritera['theta']):
					#notAligned = False
				print('Last x0 = ',x0[-1])
			else:	
				if iter_count >= iterations:
					notAligned = False
					print('Last x0 = ',x0[-1])
					print('Final step in x0, x_delta = ', x_delta)
	return
	
				
						   
						   
						   
						   
						   
