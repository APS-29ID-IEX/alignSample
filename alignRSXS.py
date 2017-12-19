"""
python modules for aligning samples at APS's 29ID RSXS system

author: Max Wyman <mwyman@anl.gov

"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
plt.ion()

# Make plots live-update while scans run.
from bluesky.utils import install_qt_kicker
install_qt_kicker()
 
#Import BlueSky components
from bluesky.plans import adaptive_scan, scan
from bluesky.callbacks import LiveFit, LivePlot, LiveFitPlot, LiveTable
from bluesky.callbacks.fitting import PeakStats
from bluesky.preprocessors import subs_decorator
from bluesky import RunEngine

#Import Ophyd components
from ophyd.signal import EpicsSignalRO, EpicsSignal
from ophyd.epics_motor import EpicsMotor
from ophyd import Device, Component as Cpt
from ophyd.sim import SynAxis, SynGauss, SynSignal

#For fitting
import lmfit
import numpy as np
import scipy
import math 

#Temp data storage
import alignSummaryData as asd

fit_colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
linestyles = ['-', '--', '-.', ':']

sequence = []
x0 = []
theta5 = []

#Default values
def_xAlignRange = [-1.00, 2.00]
def_xStepRange = [0.0015, 0.1]
def_xTargetDelta = 0.0025
def_fudge = 0.5
def_thAlignRange = [4.0,6.0]
def_thCoarseStep = 0.1 
def_thFineRadius = 0.2 
def_thStepRange = [0.015, 0.1]	#min thMotor step is 0.01 degrees

def color_y_axis(ax, color):
	"""
	Colors the passed axis.

	Parameters
	----------
	ax       :  axes object whose color is to be changed
	colors   :  new color of axes object

	"""
	for t in ax.get_yticklabels():
		t.set_color(color)
	return None

def cleanLegend(iteration, ax, other = None):
	"""
	Replaces messy legends produced by lmfit
	
	Parameters
	----------
	iteration	: integer, iteration of alignment loop
	ax			: axis object whose legend is being repaced
	other		: list of dicts for other items in legend containing
					label	:	label to use in legend
					shape	:	shape to use in legend entry
					color	:   color of legend entry
					size	:   size of legend entry symbol
					 
	"""
	
	leg = ax.legend()
	leg.remove()
	
	iterHandles = []
	for i in range(iteration):
		iterHandles.append(mpatches.Patch(color = fit_colors[i+1 % 7], 
									label = 'Iteration {} - Fit'.format(i+1)))
		
	if other is not None:
		for j in other:
			legline = mlines.Line2D([], [], markeredgecolor=j['color'],
							markerfacecolor = 'none',
							marker=j['shape'], markersize=j['size'], 
							label=j['label'], linestyle = 'none')
			iterHandles.append(legline)
			
	ax.legend(handles=iterHandles, loc = 0)

	return

def cleanAxes(figPlot, line1, line2, axis1, axis2, 
			  color1, color2, ylabel1, ylabel2, xlabel):
	"""
	Sets up multi-axis plot where each plot is a different color
	
	Parameters
	----------
	figPlot			: Plot object
	line1, line2	: line objects for the two sets of data being displayed
	axis1, axis2	: axis objects for the two sets of data being displayed
	color1, color2	: colors of the line/axis objects
	ylabel1, ylabel2: labels for the data sets
	xlabel			: x axis label (same for both data sets)
	"""			  
				  
	axis1.set_ylabel(ylabel1,color = color1)
	axis1.set_xlabel(xlabel)
	axis2.set_ylabel(ylabel2,color = color2) 
	
	axis1.relim()
	axis2.relim()
	axis1.autoscale_view()
	axis2.autoscale_view()
	
	axis1.xaxis.set_major_locator(MaxNLocator(integer=True))
	
	color_y_axis(axis1, color1)
	color_y_axis(axis2, color2)

	lns = [line1]+[line2]
	labs = [l.get_label() for l in lns]

	figPlot.legend(lns, labs)

	return

def gaussian(x, A, sigma, x0):
	"""
	Create gaussian for fitting and simulation
	
	Parameters
	----------
	x		:	input upon which gaussian is evaluated
	A		: 	peak of gaussian
	sigma	:	"spread" of gaussian
	x0		:	location of peak of gaussian
	"""
	return A*np.exp(-(x - x0)**2/(2 * sigma**2))
	
def erfx(x, low, high, wid, x0):
	"""
	Create error function for fitting and simulation
	
	Parameters
	----------
	x		:	input upon which error function is evaluated
	low		: 	min value of error function
	high	: 	max value of error function
	sigma	:	"spread" of error function transition region
	x0		:	location of error function's "center"
	"""
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
	det = SynErfGauss('det', motor, 'motor', center=0, Imax=1, sigma=1)
	"""
	
	def __init__(self, name, motor, motor_field, Imax, wid, x0,
				 noise=None, noise_multiplier=1, **kwargs):
		if noise not in ('poisson', 'uniform', None):
			raise ValueError("noise must be one of 'poisson', 'uniform', None")
		self._motor = motor

		def func():
			m = motor.read()[motor_field]['value']
			v = erfx(m, 0, Imax, wid, x0)
			v += gaussian(m, Imax*0.20, 0.1, x0 + 1 + wid/2.0 + 0.25*np.random.random())
			if noise == 'poisson':
				v = int(np.random.poisson(np.round(v), 1))
			elif noise == 'uniform':
				v += np.random.uniform(-1, 1) * noise_multiplier
			return v

		super().__init__(func=func, name=name, **kwargs)

class CurAmp(Device):
	"""
	CurAmp: ophyd device for current amplifier that requires scan state to be 
			set to passive to be triggered properly.  Sets back to original
			state after run.
	"""
	ca_value = Cpt(EpicsSignalRO, 'read')
	# trigger for detector when in it's passive state
	trig = Cpt(EpicsSignal, 'read.PROC', trigger_value=1) 
	det_state = Cpt(EpicsSignal, 'read.SCAN')
	initial_state = 0
	
	def stage(self):
		# before run get current state
		self.initial_state = self.det_state.get() 
		# put detector into passive state
		self.det_state.put(0, wait=True) 
		super().stage()

	def unstage(self):
		# return detector to pre-run state
		self.det_state.put(self.initial_state, wait=True)
		super().unstage()


def fitAndPlotBumplessData(y, x, x0, width, fudge = 0.5, ax = [], color = 'r',
						   linestyle = '-'): 
	"""
	fitAndPlotBumplessData:  
		-reduces x-scan data to eliminate bump
		-fits reduced data
		-plots (if axes supplied) fit line
		-returns center of fit's transition
		
	Parameters
    ----------
	y				: y data from scan
	x				: x data from scan
	x0				: estimated half-max from scan
	width			: estimated width of transition from scan
	fudge			: Fudge factor for x-scan data reduction to remove bump 
					  from post-scan fit (in mm).  Defaults to def_fudge
	ax 				: axis object for plotting fit, 
					  Default value is None which suppresses plot
	color			: color of fit line. Defaults to red.
	linestyle 		: string descriptor of linetype for fit line, 
					  Default value is solid line
	"""

	x_reduced_max = x0 + fudge + 0.5*width
	
	xdata = np.asarray(x)
	ydata = np.asarray(y)
	ix = np.where(xdata < x_reduced_max)
	
	x_reduced = xdata[ix]
	y_reduced = ydata[ix]
	
	red_model = lmfit.Model(erfx, missing = 'drop')
	red_guess = {'low': min(y_reduced),
				 'high': max(y_reduced),
				 'wid': lmfit.Parameter('wid', value = width, min=0),
				 'x0': x0}

	params = red_model.make_params(low = min(y_reduced), high = max(y_reduced), 
								   wid = width, x0 = x0)

	redFit = red_model.fit(y_reduced, params, x = x_reduced)
	redFit_x0 = redFit.result.params['x0'].value
	redFit_width = redFit.result.params['wid'].value
	
	#plot new fit
	if ax is not None:
		redFit.plot_fit(ax = ax, data_kws = {'visible':False}, 
						fit_kws={'color': color, 
								 'linestyle' : linestyle},
						init_kws={'visible':False},
						xlabel = 'motorX', ylabel = 'detX')
		ax.set_title(' ')

	return redFit_x0


def xScanAndCenter(iteration, detector, det_name, motor, motor_name, 
				   alignRange = [-1, 2], stepRange = [0.01, 0.10], 
				   targetDelta = def_xTargetDelta, plotter = [],
				   SimAlign = False, fudge = def_fudge, verbose = False):
	"""
	align_x:  
		-scan xMotor through range
		-monitors detector looking for transition of beam/no-beam	

	Parameters
    ----------
	iteration			: overall iteration of full alignment loop, required for 
						  plot colors/legend
	detector			: detector object
	detector_name		: detector object's name
	motor				: X motor object
	motor_name			: X motor object's name
	alignRange			: list of two floats, describing range of x-alignment 
						  Default value is [-1,2] mm,
	StepRange 			: list of two floats, minimum and maximum steps for 
						  adaptive scan. Default value is [0.015, 0.075] degrees 
						  Minimum thMotor step is 0.01 degrees
	targetDelta			: float, maximum jump in detector value before scan step 
						  decreased, default is set to def_xTargetDelta
						  value is 10 in detector units, 
	plotter  			: axis object for x alignment data, 
						  Default value is None which later creates new object
	SimAlign  			: boolean flag for simulation/testing 
						  Default value is False,
	fudge				: Fudge factor for x-scan data reduction to remove bump 
						  from post-scan fit (in mm).  Defaults to def_fudge
	verbose 			: boolean for showing alignment process using LiveTable
						  callbacks. Default value is False
	"""


	if plotter is None:
		figX, xAx = plt.subplots()
	else:
		xAx = plotter	
		
	cur_color = fit_colors[iteration % len(fit_colors)]
	cur_linestyle = linestyles[iteration // len(fit_colors)]

	scanPlot = LivePlot(det_name,x=motor_name, markeredgecolor = cur_color, 
						markerfacecolor = 'none', ax = xAx, linestyle = 'none', 
						marker = '^', label = '{} - data'.format(iteration))
		
	comp_model = lmfit.Model(erfx, prefix="erf_") + lmfit.Model(gaussian, 
							 prefix = "gau_")
	comp_guess = {'erf_low': 0,
				  'erf_high': 0.03,
				  'erf_wid': lmfit.Parameter('erf_wid', value = 0.4, min=0),
				  'erf_x0': -0.1,
				  'gau_A': 0.01,
				  'gau_sigma': lmfit.Parameter('gau_sigma', value = .1, min=0),
				  'gau_x0': 1.0}
	xLiveFit = LiveFit(comp_model, det_name, {'x': motor_name}, comp_guess, 
					   update_every=5)
	
	lt = LiveTable([detector, motor])
	
	if verbose:
		cbs = [scanPlot, xLiveFit, lt]
	else:
		cbs = [scanPlot, xLiveFit]
				
	@subs_decorator(cbs)
	def preFitScan(detectors, motor, alignRange, stepRange, targetDelta):
		yield from adaptive_scan([detectors], det_name, motor, 
								 start = alignRange[0], stop = alignRange[1], 
								 min_step = stepRange[0], 
								 max_step = stepRange[1], 
								 target_delta = targetDelta,
								 backstep = True)

	yield from preFitScan(detector, motor, alignRange, stepRange, targetDelta)

	x0_rough = xLiveFit.result.params['erf_x0'].value
	width = xLiveFit.result.params['erf_wid'].value
		
	new_x0 = fitAndPlotBumplessData(np.asarray(scanPlot.y_data), 
									np.asarray(scanPlot.x_data), 
									x0_rough, width, ax = xAx, 
									color = cur_color, 
									linestyle = cur_linestyle,
									fudge = fudge)

	other_data = []
	leg_raw	= {'label'  : 'Data',
			   'shape'  : '^',
			   'color'  : 'k',
			   'size'   : '10'}
	
	other_data = [leg_raw]

	cleanLegend(iteration, xAx, other_data)
	
	asd.x0.append(new_x0)
	
	if not SimAlign:
		if (asd.x0[-1] <= min(alignRange) or asd.x0[-1] >= max(alignRange)):
			print("Peak estimated to be outside of alignment range")
		else:
			motor.move(asd.x0[-1])


def thetaScanAndCenter(iteration, detector, detector_name, 
					   motor, motor_name, motor2, motor2_name,
					   alignRange = [4, 6], coarsePTS=10, 
					   fineRadius=0.3, fStepRange = [0.015, 0.075], 
					   targetDeltaFactor = 10, plotter = [], det4offset = None
					   SimAlign = False, verbose = False):
	"""
	align_theta:  
		-sets tthMotor to 10 degrees
		-scan thMotor through range
		-monitors detector looking for peak
		-move thMotor to peak
		-set thMotor set point to 5 degrees
		-move both motors back to zero
		
	Parameters
    ----------
	iteration			: overall iteration of full alignment loop, required for 
						  plot colors/legend
	detector			: detector object
	detector_name		: detector object's name
	motor				: theta motor object
	motor_name			: theta motor object's name
	motor2				: 2theta motor object
	motor2_name			: 2theta motor object's name
	alignRange			: list of two floats, 
						  Default value is [4.0,6.0] degrees,
	coarsePTS			: number of measurements in coarse scan,
						  defaults to 10
	fineRadius 			: float, distance from coarse peak defining region of 
						  fine scan.  Default value is 0.3 degrees,
	fStepRange 			: list of two floats, minimum and maximum steps for 
						  adaptive scan. Default value is [0.015, 0.075] degrees 
						  Minimum thMotor step is 0.01 degrees
	targetDeltaFactor	: float, maximum jump in detector value before scan step 
						  decreased is determined from coarse scan's detector
						  peak value divided by the targetDeltaFactor. Default 
						  value is 10 in detector units, 
	plotter  			: axis object for theta alignment data, 
						  Default value is None which later creates new object
	det4offset			: Detector's offset from 0 (for 2theta motor), defaults
					      to None leading to attempt to get it from PV
					      29idd:userCalcOut2.VAL
	SimAlign  			: boolean flag for simulation/testing 
						  Default value is False,
	verbose 			: boolean for showing alignment process using LiveTable
						  callbacks. Default value is False
	"""

	
	# detector offset angle PV: 29idd:userCalcOut2.VAL
	# motor2 move to offset + 10 degrees
	if det4offset is None:
		d4offset =  EpicsSignalRO('29idd:userCalcOut2.VAL', name = 'd4offset')
		det4offset = d4offset.get()
		
	if not SimAlign:
		motor2.move(det4offset+10.0)
		motor.move(alignRange[0])
	    
	if plotter is None:
		figTh, thAx = plt.subplots()
	else:
		thAx = plotter	
		
	cur_color = fit_colors[iteration % len(fit_colors)]
	   
	coarsePlot = LivePlot(detector_name,x=motor_name, linestyle = 'none', 
						  marker = '^', markerfacecolor = 'none',
						  markeredgecolor = cur_color, ax = thAx, 
						  label = '{} - coarse'.format(iteration))
	coarsePeak = PeakStats(motor_name,detector_name)

	lt = LiveTable([detector, motor])
	
	if verbose: 
		coarse_cbs = [coarsePlot,coarsePeak, lt]
	else:
		coarse_cbs = [coarsePlot,coarsePeak]
		
	@subs_decorator(coarse_cbs)
	def coarseScan(detector, motor, cRange, pts = 10):
		#Coarse scan to find region of peak
		yield from scan([detector], motor, cRange[0], cRange[1], num = pts)
		
	yield from coarseScan(detector, motor, alignRange, pts = coarsePTS)

	finePlot = LivePlot(detector_name,x=motor_name, linestyle = 'none', 
						marker = 'o', markerfacecolor = 'none', 
						markeredgecolor = cur_color, ax = thAx, 
						label = '{} - fine'.format(iteration))
	fineThetaModel = lmfit.Model(gaussian)
	fineThetaInitGuess = {'A': coarsePeak.max[1], 
						  'sigma': lmfit.Parameter('sigma', .03, min=0),
						  'x0': coarsePeak.max[0]}
	fineThetaLiveFit = LiveFit(fineThetaModel, detector_name, {'x': motor_name}, 
							   fineThetaInitGuess, update_every=1)
	fineThetaLiveFitPlot = LiveFitPlot(fineThetaLiveFit, ax = thAx, label='fit', 
									   color = fit_colors[iteration % 7],
									   linestyle = linestyles[iteration // 7])
		
	fRange = [coarsePeak.max[0]-fineRadius, coarsePeak.max[0]+fineRadius]
	   
	thTargetDelta = coarsePeak.max[1]/targetDeltaFactor
	   
	if not SimAlign:
		if coarsePeak.max[0] > alignRange[1]-1:
			motor.move(4.00)
	
	if verbose:
		fine_cbs = [finePlot, fineThetaLiveFitPlot, lt]
	else:
		fine_cbs = [finePlot, fineThetaLiveFitPlot]
		
	@subs_decorator(fine_cbs)
	def fineScan(detectors, motor, fRange, pts = 50):
		#Adaptive scan on region localized near peak
		yield from adaptive_scan([detector], detector_name, motor, 
								 start = fRange[0], stop = fRange[1], 
								 min_step = fStepRange[0], 
								 max_step = fStepRange[1], 
								 target_delta = thTargetDelta, 
								 backstep = False)
	
	yield from fineScan(detector, motor, fRange)  
	
	asd.theta5.append(fineThetaLiveFit.result.params['x0'].value)
	
	if not SimAlign:
		if (asd.theta5[-1] <= min(alignRange) or asd.theta5[-1] >= max(alignRange)):
			print("Iteration #{}: Peak estimated to be outside of alignment range".format(iteration))
		else:
			motor.move(asd.theta5[-1])
			motor.set_current_position(5.00)
		motor.move(0)
		motor2.move(det4offset)
	
	return

def multiScan(detX, detX_name, detTh, detTh_name, 
			  motorX, mX_name, motorTh, mTh_name, motorTTh, mTTh_name, 
			  iterations = 1, alignCriteria = {}, theta_kws = None, 
			  x_kws = None, SimAlign = False, xPlot = [], thPlot = [],
			  sumPlot = None, sumAxes = None, verbose = False):
	"""
	multiScan:  
		-set tthMotor to 0 degrees
		-set thMotor to 0 degrees
		-scan x-direction for transition from no-beam to full-beam
		-scan theta direction with detector at 2theta = 10 degrees
		-plot summary (if flagged) of x(half-max), theta(peak)

	Parameters
    ----------
	detX				: X detector object (in real alignment, same as detTh)
	detX_name			: X detector object's name (in real alignment, same as 
						  detTh_name)
	detTh				: Theta detector object (in real alignment, same as 
						  detX)
	detTh_name			: Theta detector object's name (in real alignment, same 
						  as detX_name)
	motorX				: X motor object
	motorX_name			: X motor object's name
	motorTh				: Theta motor object
	mTh_name			: Theta motor object's name
	motorTTh			: 2theta motor object
	mTTh_name			: 2theta motor object's name
	iterations			: Number of iterations of the X, Theta alignment scans
						  to run through
	alignCriteria		: dict of criteria for ending alignment
						  Default value is None,
	theta_kws 			: dict for arguments to be passed to theta-alignment
					      Default value is None,
	x_kws 				: dict for arguments to be passed to x-alignment
						  Default value is None

	xPlot				: axes object for x-alignment plot
	thPlot				: axes object for theta-alignment plot
	sumPlot				: plot object for usmmary plot
	sumAxes				: axes for the x-/theta-alignment summaries
	SimAlign  			: boolean flag for simulation/testing 
						  Default value is False,
	verbose 			: boolean for showing alignment process using LiveTable
						  callbacks. Default value is False
	"""
	if not SimAlign:
		motorTh.move(0)
		try:
			motorTTh.move(theta_kws['det40ffset'])
		except:
			d4offset =  EpicsSignalRO('29idd:userCalcOut2.VAL', 
									  name = 'd4offset')
			det4offset = d4offset.get()
			motorTTh.move(det4offset)

	for i in range(iterations):
		asd.sequence = np.arange(i+1)+1
		yield from xScanAndCenter(i+1, detX, detX_name, 
								  motorX, mX_name, plotter = xPlot,
								  SimAlign = SimAlign, verbose = verbose,
								  **x_kws)
		yield from thetaScanAndCenter(i+1, detTh, detTh_name, 
									  motorTh, mTh_name, 
									  motorTTh, mTTh_name,
									  plotter = thPlot,
									  SimAlign = SimAlign, verbose = verbose,
									   **theta_kws)
				
		if sumPlot is not None:
			if i == 0:
				xLine, = sumAxes[0].plot(asd.sequence, asd.x0, 
										markeredgecolor = 'r', 
										markerfacecolor = 'none', 
										marker = '^', linestyle='none', 
										label = 'x_offset')
				thetaLine, = sumAxes[1].plot(asd.sequence, asd.theta5, 
											markeredgecolor = 'b', 
											markerfacecolor = 'none', 
											marker = 'o', linestyle='none', 
											label = 'theta_offset')
			else:
				xLine.set_data(asd.sequence, asd.x0)
				thetaLine.set_data(asd.sequence, asd.theta5)
			cleanAxes(sumPlot, xLine, thetaLine, sumAxes[0], sumAxes[1], 
					  'r', 'b', 'x offset (mm)', 'theta @ 2theta = 10', 
					  'Iteration')
					
		
def alignRSXS(iterations = 3, alignCriteria = None, theta_kws = None, 
			  x_kws = None, SimAlign = False, verbose = False):     
	"""
	alignRSXS:
		-sets up detector and motor objects
		-sets up summary plots
		-calls up multiScan
		
	Parameters
    ----------
    iterations 		: integer number of iterations to perform [unless alignment
					  criteria is not None and met]. Default value is 3,
	alignCriteria	: dict of criteria for ending alignment
					  Default value is None,
	theta_kws 		: dict for arguments to be passed to theta-alignment -- most
					  important is det4offset 
					  Default value is None,
	x_kws 			: dict for arguments to be passed to x-alignment
					  Default value is None
	SimAlign		: boolean for simulated run (not implemented for full loop
					  only for individual x, theta loops via the keyword args)
	verbose 		: boolean for showing alignment process using LiveTable
					  callbacks. Default value is False
    
	"""
	
	asd.sequence = []
	asd.x0 = []
	asd.theta5 = []

	if alignCriteria is None:
		alignCriteria = {}  #TODO meant as way of introducing variable number of 
							#steps for alignment process, where the number of 
							#iterations would based off of the change between 
							#iterations
	if x_kws is None:
		x_kws = {}
	if theta_kws is None:
		theta_kws = {}

	detX_name = 'ca4det_ca_value'
	detTh_name = detX_name 
	mX_name = 'motorX'
	mTh_name = 'motorTh'
	mTTh_name = 'motorTTh'

	#setup ophyd objects
	if not SimAlign:
		motorX = EpicsMotor('29idd:m1', name = mX_name)
		motorTh = EpicsMotor('29idd:m7', name = mTh_name)
		motorTTh = EpicsMotor('29idHydra:m1', name = mTTh_name)
		
		detX = CurAmp('29idd:ca4:', name = 'ca4det', read_attrs = ['ca_value'])
		detTh = detX		
	else:
		motorX = SynAxis(name = mX_name)
		motorTh = SynAxis(name = mTh_name)
		motorTTh = SynAxis(name = mTTh_name)

		try:
			thAlignRange = theta_kws['alignRange'] 
		except:
			thAlignRange = def_thAlignRange 

		simCenter = (thAlignRange[1]-thAlignRange[0])*np.random.random() + \
					thAlignRange[0]
		try:
			xAlignRange = x_kws['alignRange'] 
		except:
			xAlignRange = def_xAlignRange 

		simX0 = (xAlignRange[1]-xAlignRange[0])/3*(1+np.random.random()) + \
					xAlignRange[0]

		detX = SynErfGauss(detX_name, motorX, mX_name, 
						  Imax=0.024, wid=0.20, x0 = simX0,
						  noise='uniform', noise_multiplier=0.0005)
						  
		detTh = SynGauss(detTh_name, motorTh, mTh_name, 
						   center=simCenter, Imax=1, sigma=.03, 
						   noise='uniform', noise_multiplier=0.05)
		
	#setup plots

	fig1, ax = plt.subplots() 	# x-alignment plot
	fig2, bx = plt.subplots() 	# theta-alignment plot
	fig3, cx1 = plt.subplots() 	# summary plot
	fig3.suptitle('Alignment summary')
	cx2 = cx1.twinx()

	#start run engine

	RE = RunEngine({})
	RE(multiScan(detX, detX_name, detTh, detTh_name, 
				 motorX, mX_name, motorTh, mTh_name, motorTTh, mTTh_name,
				 iterations = iterations, alignCriteria = alignCriteria, 
				 theta_kws = theta_kws, x_kws = x_kws, 
				 sumPlot = fig3, sumAxes = [cx1, cx2], xPlot = ax, thPlot = bx,
				 SimAlign = SimAlign, verbose = verbose))

					  

