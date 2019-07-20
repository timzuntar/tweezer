#Script for extracting force values from dataset of optical tweezer measurements
import warnings
import math

import numpy as np
import scipy.constants as constants

def rot_matrix(ks, phi):
    """
    Transforms trap coefficients into x,y-coordinates

    Parameters
    ----------
    ks : tuple of float
        trap coefficients along major and minor axis
    phi : float
        rotation angle between the two coordinate systems

    Returns
    -------
    P : ndarray_like
        2x2 tensor of coefficients in x,y-coordinates 
    """

    P = np.empty((2,2))
    
    P[0,0] = ks[0]*(math.cos(phi))**2 + ks[1]*(math.sin(phi))**2
    P[0,1] = math.sin(phi)*math.cos(phi)*(ks[0]-ks[1])
    P[1,0] = P[0,1]
    P[1,1] = ks[0]*(math.sin(phi))**2 + ks[1]*(math.cos(phi))**2
    
    return P

def rot_matrix_axial(theta=0):
    """
    Creates 2D rotation matrix.

    Parameters
    ----------
    theta : float
        angle of rotation

    Returns
    -------
    RT : ndarray_like
        2x2 rotation matrix
    """

    RT = np.empty((2,2))
    
    RT[0,0] = math.cos(theta)
    RT[0,1] = math.sin(theta)
    RT[1,0] = -RT[0,1]
    RT[1,1] = RT[0,0]
    
    return RT

def evaluate_force(trajectory_point, trap_pos_point, P, RT):
    """
    Evaluates force in a single time point.

    Parameters
    ----------
    trajectory_point : array of float
        x-,y-positions of particle
    trap_pos_point : array of float
        x-,y-positions of optical trap
    P : ndarray_like
        rotation matrix
    RT : ndarray_like
        coefficient tensor

    Returns
    -------
    RT : ndarray_like
        2x2 rotation matrix
    """

    delta = np.subtract(trajectory_point,trap_pos_point)
    force_point = np.matmul(np.matmul(RT,P),delta)
    
    return force_point*1e6

def calculate(time, trajectory, trap_position, ks, phi=0.):
    """Updated version, accounts for trap rotation.
    Provided arrays of points in time and spatial coordinates of both the optical trap
    and the trapped particle as well as trap stiffnesses, calculates forces which the trap beam exerts on the particle and their mean values.
    Assumes that trap coefficients are rotated, but trap positions aren't. 
    This is the case when calibration.calibrate is run, but x,y-data aren't overwritten. Make sure that offset has already been subtracted from trap positions.

    Parameters
    ----------
    time : array_like
        time coordinates
    trajectory : ndarray_like
        x-coordinates and y-coordinates of trapped bead
    trap_position : ndarray_like
        x-coordinates and y-coordinates of trap
    ks : tuple of floats
        trap stiffnesses in x- and y-directions [N/m]
    phi : float
        Rotation angle of long axis of trap with respect to x-axis [rad]

    Returns
    -------
    forces : ndarray_like
        n-by-2 array of forces on bead at each time point
    means : ndarray_like
        mean values of forces along each axis

    Raises
    ------
    IndexError
        if array sizes do not match
    Warning
        if any trap coefficient is less than 0 (no bound state)
    """
    
    n = len(time)
    trajectory = np.array(trajectory)
    trap_position = np.array(trap_position)

    if ( n!=len(trajectory[:, 0]) or n!=len(trap_position[:, 0])):
        raise IndexError("Array dimensions need to be identical")
    if (ks[0] < 0 or ks[1] < 0):
        warnings.warn("Value of one or more trap coefficients is negative")
        
    P = rot_matrix(ks,phi)
    RT = rot_matrix_axial()
    forces = np.empty((n, 2))

    for point in range(n):
        forces[point] = evaluate_force(trajectory[point,:],trap_position[point,:],P,RT)
        
    means = np.mean(forces, axis=0)
    print("Mean force values in pN: ", means)
    
    return forces, means

def calculate_axial(time, trajectories, trap_positions, ks_1, ks_2, phi_1=0., phi_2=0.,inter_particle = False):
    """Updated version, accounts for trap rotation.
    Provided arrays of points in time and spatial coordinates of traps, particles, as well as trap stiffnesses,
    calculates forces acting on a pair of trapped particles parallel and perpendicular to the trap-trap axis as a function of the traps' distance.

    Parameters
    ----------
    time : array_like
        array of time values
    trajectories : ndarray_like
        x- and y-coordinates of both trapped particles
    trap_positions : ndarray_like
        x- and y-coordinates of both traps
    ks_1 : tuple of floats
        stiffness of trap #1 in rotated coordinates [N/m]
    ks_2 : tuple of floats
        stiffness of trap #2 in rotated coordinates [N/m]
    phi_1 : float
        Rotation angle of long axis of 1st trap with respect to x-axis [rad]
    phi_2 : float
        Rotation angle of long axis of 2nd trap with respect to x-axis [rad]
    inter_particle : bool
        Determines whether or not the distance and angle is determined by particles' (True) or traps' (False) positions 

    Returns
    -------
    forces : ndarray_like
        parallel and perpendicular forces on both particles at each point in time
    means : array_like
        means of those forces
    distances : ndarray_like
        distances between traps -not particles- at each time point
    mean_distance : float
        mean of the above values
    """

    n = len(time)
    
    if ( n!=len(trajectories[:, 0]) or n!=len(trap_positions[:, 0])):
        raise IndexError("Array dimensions need to be identical")
    if (any(k<0 for k in ks_1) or any(k<0 for k in ks_2)):
        warnings.warn("Value of one or more trap coefficients is negative")
        
    P_1 = rot_matrix(ks_1,phi_1)
    P_2 = rot_matrix(ks_2,phi_2)
    forces_1 = np.empty((n, 2))
    forces_2 = np.empty((n, 2))
    distances = np.empty(n)
    
    for point in range(n):

        if (inter_particle == False):
            if (trap_positions[point,2] == trap_positions[point,0]):
                if (trap_positions[point,3] - trap_positions[point,1] > 0):
                    theta = math.pi/2.
                elif (trap_positions[point,3] - trap_positions[point,1] < 0):
                    theta = -math.pi/2.
            else:
                theta = np.arctan((trap_positions[point,3] - trap_positions[point,1])/(trap_positions[point,2] - trap_positions[point,0]))
            distances[point] = math.sqrt((trap_positions[point,0]-trap_positions[point,2])**2 + (trap_positions[point,1]-trap_positions[point,3])**2)
        elif (inter_particle == True):
            if (trajectories[point,2] == trajectories[point,0]):
                if (trajectories[point,3] - trajectories[point,1] > 0):
                    theta = math.pi/2.
                elif (trajectories[point,3] - trajectories[point,1] < 0):
                    theta = -math.pi/2.      
            else:
                theta = np.arctan((trajectories[point,3] - trajectories[point,1])/(trajectories[point,2] - trajectories[point,0]))      
            distances[point] = math.sqrt((trajectories[point,0]-trajectories[point,2])**2 + (trajectories[point,1]-trajectories[point,3])**2)   

        RT = rot_matrix_axial(theta)
        forces_1[point] = evaluate_force(trajectories[point,0:2],trap_positions[point,0:2],P_1,RT)
        forces_2[point] = evaluate_force(trajectories[point,2:4],trap_positions[point,2:4],P_2,RT)
        
    
    forces = np.hstack((forces_1,forces_2))
    means = np.mean(forces, axis=0)

    mean_distance = np.mean(np.fabs(distances), axis=0)

    print("Mean force values in pN: ", means)
    
    return forces,means,distances,mean_distance

def displacement_calculation(trajectory, trap_position):
    """Provided arrays of optical trap positions and particle trajectories, calculates the mean displacements and their uncertainties.

    Parameters
    ----------
    trajectory : ndarray_like
        x- and y- positions of both particles
    trap_position : ndarray_like
        x- and y- positions of both optical traps

    Returns
    -------
    means : array_like
        mean displacements of particles from trap centre in x- and y-directions
    sigmas : array_like
        uncertainties of displacements in x- and y-directions
    """

    n,m = trajectory.shape
    if (n!=len(trap_position[:, 0])):
        raise IndexError("Number of array elements needs to be identical")
    elif (trap_position.shape[1] < m):
        raise IndexError("Trap coordinates undefined for %d columns" % (m-trap_position.shape[1]))

    displacements = np.zeros((n,m))
    squares = np.zeros((n,m))
    variances = np.zeros(m)

    for point in range(n):
        for i in range(m):
            displacements[point,i] = trajectory[point,i]-trap_position[point,i]
            squares[point,i] = displacements[point,i]**2
    
    means = np.mean(displacements, axis=0)
    means_sq = np.mean(np.fabs(squares), axis=0)

    for point in range(m):
        variances[point] = math.sqrt((abs(means_sq[point] - np.square(means[point])))/n)

    return displacements, means, variances
