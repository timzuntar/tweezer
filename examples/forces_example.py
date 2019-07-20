import numpy as np
import tweezer.force_calc as forcecalc
import tweezer.data_generation as datagen
import tweezer.plotting as plt

# Parameters for simulation
k_values_1 = (2.5e-6,2.5e-6)
phi_1 = 0.
center_1 = (10.0,-4.)

k_values_2 = (1.e-6,1.e-6)
phi_2 = 0.
viscosity = 9.7e-4
temperature = 300
num_of_points = 1000
dt = 0.005

# Example of generating simulated motion of a bead moving in an optical trap

time, trap, trajectory = datagen.generate(num_of_points,dt,k_values_1,(0.25,0.25),(1e-5,1e-5),1.6e-6,viscosity,phi_1,center_1,293,1,False)

# Calculating forces on a single particle
f,m = forcecalc.calculate(time,trajectory[:,0:2],trap[:,0:2],k_values_1,phi_1)

displacements,means,variances = forcecalc.displacement_calculation(trajectory,trap)

plt.displacement_plot(time,displacements,means,variances,True)

# Plotting the forces
plt.force_plot(time,f,m,True)
plt.force_plot_radial(time,f,m,True)

# Example of calculating force of interaction between a pair of beads and their displacement
# here we include a test file from one of the experiments
times,traps,trajectories = plt.read_file("data/example_two_particles.dat",2)

traps_slice = np.hstack((traps[:,0:2],traps[:,3:5]))

forces,means,distances,mean_distance = forcecalc.calculate_axial(times,trajectories[:,0:4],traps_slice,k_values_1,k_values_2,phi_1,phi_2)

plt.force_plot(times,forces,means)
plt.force_distance_plot(forces,means,distances,mean_distance,True)