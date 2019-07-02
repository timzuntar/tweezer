import random
import math

import numpy as np
import scipy.constants

kB = scipy.constants.Boltzmann

def draw_positions(num_points, dt, trap_k, phi, trap_frequencies, trap_amplitudes, bead_radius, eta, temp=293, motion_type=1):
    """Simulates the Brownian motion of a colloidal bead trapped in an optical trap oscillating in x and y directions.
    The trap axes may be rotated. 

    Parameters
    ----------
    num_points : int
        # of data points to generate 
    dt : float
        time interval between two consecutive points [s]
    trap_k : tuple of floats
        trap stiffnesses in x- and y-directions [N/m]
    phi : float
        angle of rotation
    trap_frequencies : tuple of float
        trap oscillation frequencies in x- and y-directions [Hz]
    trap_amplitudes : tuple of float
        amplitudes of oscillation in x- and y-directions [m]
    bead_radius : float
        radius of trapped particle [m]
    eta : float
        viscosity of medium [Pa s]
    temp : float
        system temperature [K]
    motion_type : int
        governs trap motions and can take values 1 or 2.
        1: sinusoidal motion in x,y (default)
        2: linear motion in x,y; in this case, frequency parameters are ignored and the amplitudes become velocities in [m/s]

    Note
    ----
    Position values in output file are in micrometers!
    
    Returns
    -------
    times : array_like
        times of recorded points
    trap : ndarray_like
        x- and y-coordinates of trap
    trajectory : ndarray_like
        x- and y-coordinates of trapped bead

    Raises
    ------
    ValueError
        if times between consecutive output points are less than timestep of simulation
    """

    dt_internal=0.0001   #   internal time step used for simulation [s]

    if (dt <= dt_internal):
        raise ValueError("dt must be longer than time step of simulation")

    print("\nGenerating...")

    if (motion_type == 1):
        x = np.array([trap_amplitudes[0],0.])
        trap_coords=np.array([trap_amplitudes[0],0.])
    else:
        x = np.array([0.,0.])
        trap_coords=np.array([0.,0.])

    a = bead_radius    
    noise = np.array([0.,0.])
    dx = np.array([0.,0.])
    
    trajectory = np.zeros([num_points,2])   #   initial position of bead [x y] in microns
    trap = np.zeros([num_points,2])   #   initial position of trap
    times = np.zeros([num_points])

    i=0
    t=0
    last_sample_interval=dt+1e-10

    while i<num_points:
        if last_sample_interval > dt:
            last_sample_interval -= dt
            trajectory[i,0] = x[0]*1e6
            trajectory[i,1] = x[1]*1e6
            trap[i,0] = trap_coords[0]*1e6
            trap[i,1] = trap_coords[1]*1e6
            times[i] = t
            i += 1

        t += dt_internal
        last_sample_interval += dt_internal
        noise[0] = (2*random.random()-1)*math.sqrt(3)
        noise[1] = (2*random.random()-1)*math.sqrt(3)

        if (motion_type == 1):
            trap_coords[0] = trap_amplitudes[0]*math.cos(2*math.pi*trap_frequencies[0]*t)
            trap_coords[1] = trap_amplitudes[1]*math.sin(2*math.pi*trap_frequencies[1]*t)
        elif (motion_type == 2):
            trap_coords[0] = trap_amplitudes[0]*t
            trap_coords[1] = trap_amplitudes[1]*t

        coeffs = np.array(trap_k)*(np.array(trap_coords)-np.array(x))
        coeffs[0],coeffs[1] = rotate(coeffs[0],coeffs[1],phi)

        dx = (dt_internal/(6*math.pi*eta*a))*coeffs + (math.sqrt(2*kB*temp/(6*math.pi*eta*a))*math.sqrt(dt_internal))*np.array(noise)
        x += dx

    return times, trap, trajectory

def decenter(xdata, ydata, center):
    """Adds a constant offset to data."""

    return xdata + center[0], ydata + center[1]

def rotate(x, y, phi):
    """Rotates positions by angle phi in anticlockwise direction.
    
    Parameters
    ----------
    x : float
        x-coordinate
    y : float
        y-coordinate
    phi : float
        angle of rotation

    Returns
    -------
    x_rotated : float
        rotated x-coordinate
    y_rotated : float
        rotated y-coordinate
    """

    x_rotated = np.cos(phi)*x - np.sin(phi)*y
    y_rotated = np.sin(phi)*x + np.cos(phi)*y
    return x_rotated, y_rotated

def drift(xdata, ydata):
    """Adds a random sine function to the positions.

    Parameters
    ----------
    xdata : array_like
        x-coordinates
    ydata : array_like
        y-coordinates

    Returns
    -------
    x_drift : array_like
        x-coordinates with added drift
    x_drift : array_like
        y-coordinates with added drift
    """
    
    np.random.seed(seed=0)
    x_drift_parameters = np.random.rand(3)
    y_drift_parameters = np.random.rand(3)
    
    def sine(x, parameters):
        return parameters[0]*np.sin(parameters[1]*x/len(x)*5 + parameters[2])
    
    t = np.arange(len(xdata))

    x_drift = xdata + sine(t, x_drift_parameters)
    y_drift  = ydata + sine(t, y_drift_parameters)

    return x_drift, y_drift


def generate(num_points, dt, trap_k, trap_frequencies, trap_amplitudes, bead_radius, eta, phi=0., center=(0., 0.), temp=293, motion_type=1, drifting = True):
    """Simulates the Brownian motion of a colloidal bead trapped in an optical trap oscillating (or moving) in x and y directions.
    The trap axes may be rotated and an offset specified. Additionally, random drift can be added to the data. 

    Parameters
    ----------
    num_points : int
        # of data points to generate 
    dt : float
        time interval between two consecutive points [s]
    trap_k : tuple of floats
        trap stiffnesses in x- and y-directions [N/m]
    trap_frequencies : tuple of float
        trap oscillation frequencies in x- and y-directions [Hz]
    trap_amplitudes : tuple of float
        amplitudes of oscillation in x- and y-directions [m]
    bead_radius : float
        radius of trapped particle [m]
    eta : float
        viscosity of medium [Pa s]
    phi : float
        angle of rotation
    center : tuple of float
        x- and y-coordinate of offset
    temp : float
        system temperature [K]
    motion_type : int
        governs trap motions and can take values 1 or 2.
        1: sinusoidal motion in x,y (default)
        2: linear motion in x,y; in this case, frequency parameters are ignored and the amplitudes become velocities in [m/s]
    drifting : bool
        whether or not random drift is included

    Note
    ----
    Position values in output file are in micrometers!
    
    Returns
    -------
    times : array_like
        times of recorded points
    trap : ndarray_like
        x- and y-coordinates of trap
    trajectory : ndarray_like
        x- and y-coordinates of trapped bead

    Raises
    ------
    ValueError
        if times between consecutive output points are less than timestep of simulation
    """

    times,trap,trajectory = draw_positions(num_points, dt, trap_k, phi, trap_frequencies, trap_amplitudes, bead_radius, eta, temp, motion_type)

    trajectory[:,0],trajectory[:,1] = decenter(trajectory[:,0], trajectory[:,1], center)
    trap[:,0],trap[:,1] = decenter(trap[:,0], trap[:,1], center)

    if drifting == True:
        trajectory[:,0],trajectory[:,1] = drift(trajectory[:,0],trajectory[:,1])

    return times,trap,trajectory