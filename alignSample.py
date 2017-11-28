import matplotlib.pyplot as plt
plt.ion()

# Make plots live-update while scans run.
from bluesky.utils import install_qt_kicker
install_qt_kicker()
 
#Start BlueSky
from bluesky.plans import adaptive_scan
from bluesky.callbacks import LiveFit, LivePlot, LiveFitPlot, LiveTable
from bluesky import RunEngine
from ophyd.signal import EpicsSignalRO, EpicsSignal

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
    return (high - low) * 0.5 * (1-scipy.special.erf((x-x0)*wid)) + low    

#SynErf class required to simulate translational/x-alignment
class SynErf(SynSignal):
    """
    Evaluate a point on the error function based on the value of a motor.

    Parameters
    ----------
    name : string
    motor : Device
    motor_field : string
    x0 : number
        center of error function
    Imax : number
        max intensity of peak
    wid : number, optional
        rough width of error function
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

def align_x(theta_offset = 0,
			alignRange = [-1.00,1.00], 
			stepRange = [0.05, 0.5],
			targetDelta = 0.0025, 
			SimAlign = False):
	"""
	align_x:  
		-set two_theta_motor to 0 degrees
		-set theta_motor to 0 degrees
		-scan x_motor through range
		-monitors detector looking for transition of beam/no-beam	
	"""
	
	if not SimAlign:
		#set signals
		x_motor = EpicsSignal('29idd:m1', name = 'x_motor')
		theta_motor = EpicsSignal('29idd:m7', name = 'theta_motor')
		two_theta_motor = EpicsSignal('PV_name', name = 'two_theta_motor')
		det_tran =  EpicsSignalRO('29idd:ca4', name = 'name')
	else:
		#set simulated signals
		x_motor = SynAxis(name='x_motor')
		det_tran = SynErf('det_tran', x_motor, 'x_motor', 
						  Imax=0.024, wid=5.0, x0 = 0.2,
						  noise='uniform', noise_multiplier=0.0005)
	
	trans_model = lmfit.Model(erfx)
	init_guess = {'low': 0,
              'high': 0.03,
              'wid': lmfit.Parameter('wid', 10, min=0),
              'x0': -0.1}
	trans_lf = LiveFit(trans_model, 'det_tran', {'x': 'x_motor'}, 
						init_guess, update_every=5)
	fig, ax = plt.subplots()  # explitly create figure, axes to use below
	trans_lfp = LiveFitPlot(trans_lf, ax=ax, color='r', label = 'fit')
	trans_lp = LivePlot('det_tran','x_motor', ax=ax, marker='o', 
						linestyle='none', label = 'data')

	#Start BlueSky run engine
	RE = RunEngine({})

	#Move motors 
	if not SimAlign:
		#move theta_motor to theta_offset or 0 degrees
		#theta_motor.put(theta_offset)
		#move two_theta_motor to 0 degrees
		#two_theta_motor.put(0.00)
		pass
		
	#run adaptive scan
	RE(adaptive_scan([det_tran], 'det_tran', x_motor, 
					start = alignRange[0],
					stop = alignRange[1], 
					min_step = stepRange[0],  
					max_step = stepRange[1], 
					target_delta = targetDelta,
					backstep = True), [trans_lp, trans_lfp])

	print('Half-maximum at x0 = ',"{0:.5f}".format(trans_lf.result.params['x0'].value))
	print('Number of measurements: ',len(trans_lp.x_data))	
	
	return trans_lf.result.params['x0'].value
	
	
def align_theta(x0 = 0,
				alignRange = [3.50,6.50], 
				stepRange = [0.01, 0.4],
				targetDelta = 0.25, 
				SimAlign = False):
	"""
	align_theta:  
		-sets two_theta_motor to 10 degrees
		-set x_motor to x0 (half-max)
		-scan theta_motor through range
		-monitors detector looking for peak
	"""
	
	if not SimAlign:
		#set signals
		x_motor = EpicsSignal('29idd:m1', name = 'x_motor')
		theta_motor = EpicsSignal('29idd:m7', name = 'theta_motor')
		two_theta_motor = EpicsSignal('PV_name', name = 'two_theta_motor')
		det_rot =  EpicsSignalRO('29idd:ca4', name = 'name')

	else:
		#set simulated signals
		theta_motor = SynAxis(name='theta_motor')
		det_rot = SynGauss('det_rot', theta_motor, 'theta_motor', 
						   center=4.6, Imax=1, sigma=.1, 
						   noise='uniform', noise_multiplier=0.005)
	
	#set up fit model
	rot_model = lmfit.Model(gaussian)
	init_guess = {'A': 1,
				'sigma': lmfit.Parameter('sigma', .2, min=0),
				'x0': 5}					   
	#set up BlueSky fit
	rot_lf = LiveFit(rot_model, 'det_rot', {'x': 'theta_motor'}, 
					 init_guess, update_every=5)

	#set up plots
	fig, bx = plt.subplots()  # explitly create figure, axes to use below
	rot_lfp = LiveFitPlot(rot_lf, ax=bx, color='r', label='fit')
	rot_lp = LivePlot('det_rot','theta_motor', ax=bx, marker='o', linestyle='none', label='data')

	#Start BlueSky run engine
	RE = RunEngine({})

	#Move motors 
	if not SimAlign:
		#move x_motor to x0
		#x_motor.set(x0)
		#move two_theta_motor to 10 degrees
		#theta_motor.set(10.0)
		pass
		
	#run adaptive scan
	RE(adaptive_scan([det_rot], 'det_rot', theta_motor, 
					start = alignRange[0],
					stop = alignRange[1], 
					min_step = stepRange[0],  
					max_step = stepRange[1], 
					target_delta = targetDelta,
					backstep = True), [rot_lp, rot_lfp])

	print('Maximum at theta = ',"{0:.3f}".format(rot_lf.result.params['x0'].value))
	print('Number of measurements: ',len(rot_lp.x_data))
	
	return rot_lf.result.params['x0'].value
	
def align_sample(iterations = 3,
				alignCriteria = None,
				xAlignRange = [-1.00,1.00], 
				xStepRange = [0.05, 0.5],
				xTargetDelta = 0.0025,
				thetaAlignRange = [3.50,6.50], 
				thetaStepRange = [0.01, 0.4],
				thetaTargetDelta = 0.25, 
				SimAlign = False,
				plotReport = True,
				tableReport = True):
	"""
	align_sample:
		-loops over align_x and align_theta until criteria met or fixed number 
		 of iterations met
		-other..?	
	"""
				
	iter_count = 0	
	notAligned = True
	sequence = []
	
	x0 = [0]	
	theta0 = [0]
	
	if plotReport:
		fix, cx1 = plt.subplots()
		cx2 = cx1.twinx()
					
	while notAligned:
		x0.append(align_x(theta_offset = theta0[-1]))
		theta0.append(align_theta(x0 = x0[-1]))
		
		iter_count += 1
		sequence.append(iter_count)
			
		if plotReport:
			if iter_count == 1:
				cx1.plot(sequence, x0, color = 'r', marker = 'o', 
						linestyle='none', label = 'x_offset'))
				cx2.plot(sequence, theta0, color = 'b', marker = 'x', 
						linestyle='none', label = 'theta_offset'))
				color_y_axis(cx1, 'r')
				color_y_axis(cx2, 'b')
				plt.show()
			else:
				cx1.plot(sequence, x0, color = 'r', marker = 'o', 
						linestyle='none', label = 'x_offset'))
				cx2.plot(sequence, theta0, color = 'b', marker = 'x', 
						linestyle='none', label = 'theta_offset'))
				color_y_axis(cx1, 'r')
				color_y_axis(cx2, 'b')
#				plt.show()
		if tableReport:
			pass

		if alignCriteria not None:
			# print delta_x (and criteria)
			# print delta_theta (and criteria)
			if (delta_x < alignCriteria.x AND delta_theta < alignCritera.theta):
				notAligned = False
		else:	
			# print iters
			# print delta_x
			# print delta_theta
			if iter_count >= iterations:
				notAligned = False
				
						   
						   
						   
						   
						   
