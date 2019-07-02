import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import calibration as cal


def read_file(path, no_of_particles):
    """Unpacks a dat file.
    Parameters
    ----------
    path : string
        location and the name of dat file
    no_of_particles : int
        number of particles in interest
    Returns
    -------
    time: array_like
        timestamps
    traps: ndarray_like
        location and strength of 4 traps
    trajectories: 2D array
        x and y trajectories of no_of_particles
    Examples
    --------
    TODO
    """
    raw_data = open(path, "r")
    columns = int(14+2*no_of_particles)
    data = np.zeros((100000, columns))
    rows = 0
    # Read file by lines
    for line in raw_data.readlines():
        # If data is missing (double tab), replace it with nan
        line = line.replace('\t\t', '\tnan')
        data[rows, :] = (line.split('\t'))[:columns]
        rows += 1
    data = data[:rows, :].astype(np.float)
    raw_data.close()
    print('Shape of initial data: ', data.shape)
    # Check data how many nans it contains in trajectories
    for i in range(rows-1, 0, -1):
        check_row = np.isnan(data[i, 14:])
        # Delete row, if it contains nan
        if(np.sum(check_row) > 0):
            data = np.delete(data, i, 0)
    print('Shape of cropped data: ', data.shape)
    return data[:, 0], data[:, 2:14], data[:, 14:columns]


def trajectory_plot(time, data, averaging_time=1.):
    """
    Creates a pair of plots showing the particle's trajectory components (x(t),y(t)). Time-averaged data is overlaid onto the raw data points.
    
    Parameters
    ----------
    time : list of floats
        times of recorded points
    data : ndarray of floats
        x- and y-positions of particle at each point in time
    averaging_time : float
        time interval over which data is averaged
    """
    data = np.array(data)
    x, x_average, _ = cal.subtract_moving_average(
        time, data[:, 0], averaging_time)
    y, y_average, time = cal.subtract_moving_average(
        time, data[:, 1], averaging_time)

    trajectory = np.stack((x, y), axis=1)
    trajectory_averaged = np.stack((x_average, y_average), axis=1)

    fig = plt.figure()
    titles = ['x', 'y']
    for i in range(2):
        ax = fig.add_subplot(1, 2, i+1)
        ax.set_title('Trajectory of a trapped particle \n in {} direction'.format(titles[i]))
        ax.grid(True)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Direction {} '.format(titles[i]) + r'[$\mu m$]')
        ax.scatter(
            time, trajectory[:, i] + trajectory_averaged[:, i],
            s=4, label=r'original'
            )
        ax.scatter(time, trajectory_averaged[:, i], s=4, label=r'averaged')
        ax.legend(loc='best')
    fig.tight_layout()
    plt.show()
    return None

def displacement_plot(time,displacements,means,sigmas,plot_means=False):
    """
    Creates plots showing components of displacements from trap center (dx(t),dy(t)) for one or multiple particles.
    Red-shaded area indicates uncertainty in mean values.
    
    Parameters
    ----------
    time : list of floats
        times of recorded points
    displacements : ndarray of floats
        x- and y-displacements of one or two particles at each time-averaging point
    means : list of floats
        mean displacements for each coordinate
    sigmas : list of floats
        uncertainties of displacements
    plot_means : bool
        whether or not mean displacements and their estimated uncertainties are to be plotted
    """

    n = len(means)//2
    fig, (ax1, ax2) = plt.subplots(1, 2)
    plt.suptitle('Displacements of trapped particles')

    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('Direction x [$\mu m$]')
    ax1.grid(True)
    ax2.set_title('')
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Direction y [$\mu m$]')
    ax2.grid(True)    

    for i in range(n):
        ax1.plot(time, displacements[:, 2*i], label=r"$\Delta x_%d$" % (i+1))
        ax2.plot(time, displacements[:, 2*i+1], label=r"$\Delta y_%d$" % (i+1))

        if (plot_means == True):
            ax1.axhline(y=means[2*i], label=r"$\Delta x_{means %d}$" % (i+1), linestyle = "--", color="black", linewidth = 1)
            ax1.axhspan(means[2*i]-sigmas[2*i], means[2*i]+sigmas[2*i], facecolor="r", alpha=0.3)

            ax2.axhline(y=means[2*i+1], label=r"$\Delta y_{means %d}$" % (i+1), linestyle = "--", color="black", linewidth = 1)
            ax2.axhspan(means[2*i+1]-sigmas[2*i+1], means[2*i+1]+sigmas[2*i+1], facecolor="r", alpha=0.3)

    ax1.legend(loc='best')
    ax2.legend(loc='best')
    fig.tight_layout()
    plt.subplots_adjust(top=0.88)
    plt.show()
    return None


def calibration_plots(time, data, averaging_time=1., temp=293.):
    """
    Creates sets of scatter (y(x)) and histogram (#points(r)) plots for the particle's position. The first plot of each set shows raw data; in the second one, trap drift is averaged out and the trapping potential's axes aligned with the coordinate system.
    
    Parameters
    ----------
    time : list of floats
        times of recorded points
    data : ndarray of floats
        x- and y-positions of particle at each point in time
    averaging_time : float
        time interval over which data is averaged
    temp : float
        temperature of system in Kelvin
    """
    data = np.array(data)
    x = cal.subtract_moving_average(time, data[:, 0], averaging_time)[0]
    y = cal.subtract_moving_average(time, data[:, 1], averaging_time)[0]
    trajectory, phi, var = cal.center_and_rotate(x, y)
    k = cal.KB*temp/var*1e12

    def scatter_plot(data, trajectory, phi):
        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.set_title('Original data')
        ax1.grid(True)
        ax1.set_xlabel('Direction x ' + r'[$\mu m$]')
        ax1.set_ylabel('Direction y ' + r'[$\mu m$]')
        ax1.scatter(data[:, 0], data[:, 1], s=4)
        ax1.set_aspect('equal')
        ax2.set_title('Centered data, phi = {:.2f} rad'.format(phi,))
        ax2.grid(True)
        ax2.set_xlabel('Direction x ' + r'[$\mu m$]')
        ax2.set_ylabel('Direction y ' + r'[$\mu m$]')
        ax2.scatter(trajectory[:, 0], trajectory[:, 1], s=4)
        ax2.set_aspect('equal')
        fig.tight_layout()
        plt.show()
        return None

    def histogram_plot(trajectory, var):
        fig = plt.figure()
        titles = ['x', 'y']
        for i in range(2):
            ax = fig.add_subplot(1, 2, i+1)
            ax.set_xlabel(('Direction {} ' + r'[$\mu m$]').format(titles[i]))
            ax.set_ylabel('Bin height')
            hist, bin_edges = np.histogram(trajectory[:, i], bins=int(
                np.sqrt(len(trajectory[:, i]))), density=True)
            ax.set_title('k_{} = {:.2e}J/m^2'.format(titles[i], k[i]))
            bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.
            ax.scatter(bin_centres, hist, s=4)
            x_model = np.linspace(min(bin_centres), max(bin_centres), 100)
            prefactor = 1./np.sqrt(2.*np.pi*var[i])
            ax.plot(x_model, prefactor*np.exp(-x_model **
                                              2./(2.*var[i])), label=r'fit')
            ax.legend(loc='best')
        fig.tight_layout()
        plt.show()
        return None

    scatter_plot(data, trajectory, phi)
    histogram_plot(trajectory, var)
    return None


def potential_plot(time, data, averaging_time=1., temp=293.):
    """
    Creates a set of plots showing the trap's potential in units of kBT along each axis (V(x),V(y)).
    
    Parameters
    ----------
    time : list of floats
        times of recorded points
    data : ndarray of floats
        x- and y-positions of particle at each point in time
    averaging_time : float
        time interval over which data is averaged
    temp : float
        temperature of system in Kelvin
    """
    positions, potential_values, _ = cal.potential(
        time, data, averaging_time=averaging_time, temp=temp)

    fig = plt.figure()
    titles = ['x', 'y']
    for i in range(2):
        ax = fig.add_subplot(1, 2, i+1)
        ax.set_title('Shape of a potential in {} direction'.format(titles[i]))
        ax.set_xlabel(('Direction {} ' + r'[$\mu m$]').format(titles[i]))
        ax.set_ylabel('Potential [kT]')
        ax.scatter(positions[i], potential_values[i], s=4)
    fig.tight_layout()
    plt.show()
    return None


def potential_polynomial_fit_plot(
    time, data, averaging_time=1., temp=293., order=2
):
    """
    Creates a set of plots showing the trap's potential in units of kBT along each axis (V(x),V(y)) and fits a polynomial curve solution.
    
    Parameters
    ----------
    time : list of floats
        times of recorded points
    data : ndarray of floats
        x- and y-positions of particle at each point in time
    averaging_time : float
        time interval over which data is averaged
    temp : float
        temperature of system in Kelvin
    order : int
        order of the polynomial to be fitted
    """
    positions, potential_values, popt, _ = cal.calibrate_by_fitting_polynomial(
        time, data, averaging_time=1, order=4)

    def f(x, *coefs):
        polynomial = 0.
        i = 0
        for coef in coefs:
            polynomial += coef * x ** (2 * i)
            i += 1
        return polynomial

    fig = plt.figure()
    titles = ['x', 'y']
    for i in range(2):
        ax = fig.add_subplot(1, 2, i+1)
        ax.set_title(
            'Shape of the potential in {} direction'.format(titles[i]))
        ax.set_xlabel(('{} ' + r'[$\mu m$]').format(titles[i]))
        ax.set_ylabel('Potential [kT]')
        ax.scatter(positions[i], potential_values[i], s=4)

        x = np.linspace(positions[i][0], positions[i][-1], 100)
        ax.plot(x, f(x, *popt[i]))
    fig.tight_layout()
    plt.show()
    return None

def force_plot(time,forces,means,plot_means=False):
    """
    Creates an F(t) plot for each axis, for one or multiple particles. Input force values should be in pN.
    
    Parameters
    ----------
    time : list of floats
        times of recorded points
    forces : ndarray of floats
        2*n-column array of forces on n trapped beads in x- and y-directions
    means : list of floats
        mean forces for each coordinate
    plot_means : bool
        whether or not mean forces and their estimated uncertainties are to be plotted
    """

    n = forces.shape[1]//2
    sq = np.mean(np.square(forces),axis=0)

    length = forces.shape[0]
    fig = plt.figure()
    plt.suptitle('Components of optical gradient force on trapped particles')

    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)

    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('$F_x$ [pN]')
    ax1.grid(True)
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('$F_y$ [pN]')
    ax2.grid(True)

    for i in range(n):
        var_x = np.sqrt((abs(sq[2*i] - np.square(means[2*i])))/length)
        var_y = np.sqrt((abs(sq[2*i+1] - np.square(means[2*i+1])))/length)

        ax1.plot(time, forces[:, 2*i], label=r"$F_{x%d}$" % (i+1))
        ax2.plot(time, forces[:, 2*i+1], label=r"$F_{y%d}$" % (i+1))

        if (plot_means == True):
            ax1.axhline(y=means[2*i], label=r"$F_{x means %d}$" % (i+1), linestyle = "--", color="black", linewidth = 1)
            ax1.axhspan(means[2*i]-var_x, means[2*i]+var_x, facecolor="r", alpha=0.3)

            ax2.axhline(y=means[2*i+1], label=r"$F_{y means %d}$" % (i+1), linestyle = "--", color="black", linewidth = 1)
            ax2.axhspan(means[2*i+1]-var_y, means[2*i+1]+var_y, facecolor="r", alpha=0.3)

    ax1.legend(loc='best')
    ax2.legend(loc='best')
    fig.tight_layout()
    plt.subplots_adjust(top=0.88)
    plt.show()
    return None

def force_plot_radial(time,forces,means,plot_means=False):
    """
    Creates an F(t) plot of effective radial forces, for one or multiple particles. Input force values should be in pN.
    
    Parameters
    ----------
    time : list of floats
        times of recorded points
    forces : ndarray of floats
        2*n-column array of forces on n trapped beads in x- and y-directions
    means : list of floats
        mean forces for each coordinate
    plot_means : bool
        whether or not mean forces and their estimated uncertainties are to be plotted
    """
    n = forces.shape[1]//2
    length = forces.shape[0]
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title('Optical gradient forces on trapped particles')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Effective radial force [pN]')
    ax.grid(True)
    for i in range(n):
        radial_forces = np.sqrt(forces[:, 2*i]**2 + forces[:, 2*i+1]**2)
        means = np.mean(radial_forces,axis=0)
        sq = np.mean(np.square(radial_forces),axis=0)
        var = np.sqrt(abs(sq - np.square(means))/length)

        ax.plot(time, radial_forces,label=r"$F_{rad%d}$" % (i+1))
        if (plot_means == True):
            ax.axhline(y=means, label=r"$F_{r_{means %d}}$" % (i+1), linestyle = "--", color="black", linewidth = 1)
            ax.axhspan(means-var, means+var, facecolor="r", alpha=0.3)

    ax.legend(loc='best')
    fig.tight_layout()
    plt.show()
    return None

def force_distance_plot(forces,means,distances,mean_distance,plot_means=False,axial=True):
    """
    Plots forces on a pair of interacting particles along each axis.

    Parameters
    ----------
    forces : ndarray of floats
        4-column array of forces on pair of trapped particles in x- and y-directions
    means : list of floats
        mean forces for each coordinate
    distances : array of floats
        inter-particle distances at each point in time
    mean_distance : float
        mean inter-particle distance
    plot_means : bool
        whether or not mean forces and distance as well as their estimated uncertainties are to be plotted
    """

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title('Sum of gradient forces on trapped particle')
    ax.set_xlabel('Inter-particle distance [$\mu m$]')
    ax.set_ylabel('$\Sigma$F [pN]')
    ax.grid(True)

    if (plot_means == True):
        sq_1 = np.mean(np.square(forces[:,0]+forces[:,2]),axis=0)
        sq_2 = np.mean(np.square(forces[:,1]+forces[:,3]),axis=0)
        length = forces.shape[0]

        var_x = np.sqrt((abs( sq_1 - np.square(means[0]+means[2]) ))/length)
        var_y = np.sqrt((abs( sq_2 - np.square(means[1]+means[3]) ))/length)
        
        ax.axvline(x=mean_distance, label=r"$r_{means}$", linestyle = "--", color="black", linewidth = 1, alpha=0.4)

        ax.axhline(y=means[0]+means[2], label=r"$\Sigma F_{\parallel avg}$", color="C0", linestyle = "--", linewidth = 1)
        ax.axhspan(means[0]+means[2]-var_x, means[0]+means[2]+var_x, facecolor="r", alpha=0.3)
        ax.axhline(y=means[1]+means[3], label=r"$\Sigma F_{\bot avg}$", color="C1", linestyle = "--", linewidth = 1)
        ax.axhspan(means[1]+means[3]-var_x, means[1]+means[3]+var_x, facecolor="r", alpha=0.3)

    if (axial == True):
        ax.plot(distances, forces[:, 0] + forces[:,2], label=r'$\Sigma F_{\parallel}$')
        ax.plot(distances, forces[:, 1] + forces[:,3], label=r'$\Sigma F_{\bot}$')
    else:
        ax.plot(distances, forces[:, 0] + forces[:,2], label=r'$\Sigma F_x$')
        ax.plot(distances, forces[:, 1] + forces[:,3], label=r'$\Sigma F_y$')     

    ax.legend(loc='best')
    fig.tight_layout()
    plt.show()
    return None