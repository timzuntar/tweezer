import os
from shutil import rmtree
import numpy as np
import pims
import trackpy

import tweezer.calibration as calibration
import tweezer.calibration_generate_data as calibration_generate_data
import tweezer.force_calc as force_calc
import tweezer.generate_gaussian_particles as generate_gaussian_particles
import tweezer.offset as offset
import tweezer.plotting as plotting
import tweezer.tracking as tracking


#1: Generating the synthetic data

coeffs = (0.5e-6,1.5e-6)
phi = 0.52
trap_data = np.array([25.,30.,1.])

num_points = 20
trajectory = calibration_generate_data.generate(coeffs,293,phi,(trap_data[0],trap_data[1]),num_points,False)
times = calibration_generate_data.generate_time(num_points)

#2: Using said data to generate video frames

output_folder = "generated_test_frames"
if os.path.exists(output_folder):
    rmtree(output_folder)
os.makedirs(output_folder,exist_ok=True)

for i in range(num_points):
    img = generate_gaussian_particles.make_frame((50, 50), [[trajectory[i,0], trajectory[i,1], 200, 7]])
    img.save(output_folder + "/%s.png" % str(i))

#3: Track the particle in generated frames, create .dat file
laser_powers = np.full((num_points),1.)

for i in range(num_points):
    images = pims.ImageSequence(output_folder + "/*.png")
    
features = trackpy.batch(images, 7, minmass=50, invert=False)
tracks = trackpy.link_df(features, 15, memory=10)

all_traps = [trap_data,[0,0,-1],[0,0,-1],[0,0,-1]]

traps = [[all_traps[i][:] for j in range(num_points)]for i in range(4)]

tracking.save_tracked_data_pandas(output_folder + "/tracked_data.dat",images, tracks, times, laser_powers, traps)

#4: Read .dat file back into the program

times,traps,trajectory = plotting.read_file(output_folder + "/tracked_data.dat",1)

#5: Calibration, removing offset etc.

trap_offsets = offset.single_particle_offset(traps[:,0:2],trajectory[:,::-1]) #should be close to 0 since trap isn't shifted
coeffs_estimate,rotation_angle,_ = calibration.calibrate(times,trajectory,0.005)

#6: Force analysis

print("Trap offsets:", trap_offsets)
print("Estimate of trap coefficients:", coeffs_estimate)
print("Trap rotation angle", rotation_angle)

forces,mean_forces = force_calc.calculate(times,trajectory[:,1::-1],traps[:,0:2],coeffs_estimate,rotation_angle)

#7: Create graphs
plotting.trajectory_plot(times,trajectory[:,1::-1],0.005)
plotting.force_plot(times,forces,mean_forces)