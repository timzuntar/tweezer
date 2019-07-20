Example of workflow
===================

This is a walkthrough of the data extraction and analysis process on a concrete example, namely *examples/pipeline_example.py*. It should clear things up where :ref:`ref-quickstart` didn't.

First the dependencies are imported. Then, since we do not have suitable data on hand, let us first simulate some:

.. testsetup::

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

.. testcode::

    coeffs = (0.5e-6,1.5e-6)
    phi = 0.52
    trap_data = np.array([25.,30.,1.])
    
We specify the trap parameters: stiffnesses in x- and y-directions, position, and relative strength (in this case 1).
    
.. testcode::

    num_points = 20
    trajectory = calibration_generate_data.generate(coeffs,293,phi,(trap_data[0],trap_data[1]),num_points,False)
    times = calibration_generate_data.generate_time(num_points)
    
We use these parameters to generate a particle trajectory.
    
.. testcode::

    output_folder = "generated_test_frames"
    if os.path.exists(output_folder):
        rmtree(output_folder)
    os.makedirs(output_folder,exist_ok=True)

    for i in range(num_points):
        img = generate_gaussian_particles.make_frame((50, 50), [[trajectory[i,0], trajectory[i,1], 200, 7]])
        img.save(output_folder + "/%s.png" % str(i))

Image frames are generated one by one, each corresponding to a data point. They are of course a highly idealised approximation of reality, with particles represented with "blobs" of brightness, in this case 7 pixels across.
    
.. testcode::

    for i in range(num_points):
        images = pims.ImageSequence(output_folder + "/*.png")
    
    features = trackpy.batch(images, 7, minmass=50, invert=False)
    tracks = trackpy.link_df(features, 15, memory=10)
    
We open the frames and track the particles' positions. You can see that we search for features around 7 px across and larger than 50 px in total. When running *batch* and *link_df* functions, the console should show the following output:

.. code-block:: bash

    trackpy.feature.batch:  Frame 0: 1 features
    trackpy.feature.batch:  Frame 1: 1 features
    trackpy.feature.batch:  Frame 2: 1 features
    ...
    trackpy.linking.linking.link_iter:  Frame 0: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 1: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 2: 1 trajectories present.
    ...
    
    
.. testoutput::
    :hide:

    trackpy.feature.batch:  Frame 0: 1 features
    trackpy.feature.batch:  Frame 1: 1 features
    trackpy.feature.batch:  Frame 2: 1 features
    trackpy.feature.batch:  Frame 3: 1 features
    trackpy.feature.batch:  Frame 4: 1 features
    trackpy.feature.batch:  Frame 5: 1 features
    trackpy.feature.batch:  Frame 6: 1 features
    trackpy.feature.batch:  Frame 7: 1 features
    trackpy.feature.batch:  Frame 8: 1 features
    trackpy.feature.batch:  Frame 9: 1 features
    trackpy.feature.batch:  Frame 10: 1 features
    trackpy.feature.batch:  Frame 11: 1 features
    trackpy.feature.batch:  Frame 12: 1 features
    trackpy.feature.batch:  Frame 13: 1 features
    trackpy.feature.batch:  Frame 14: 1 features
    trackpy.feature.batch:  Frame 15: 1 features
    trackpy.feature.batch:  Frame 16: 1 features
    trackpy.feature.batch:  Frame 17: 1 features
    trackpy.feature.batch:  Frame 18: 1 features
    trackpy.feature.batch:  Frame 19: 1 features
    trackpy.linking.linking.link_iter:  Frame 1: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 2: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 3: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 4: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 5: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 6: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 7: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 8: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 9: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 10: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 11: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 12: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 13: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 14: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 15: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 16: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 17: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 18: 1 trajectories present.
    trackpy.linking.linking.link_iter:  Frame 19: 1 trajectories present.

.. testcode::

    laser_powers = np.full((num_points),1.)

    all_traps = [trap_data,[0,0,-1],[0,0,-1],[0,0,-1]]
    traps = [[all_traps[i][:] for j in range(num_points)]for i in range(4)]
    
    tracking.save_tracked_data_pandas(output_folder + "/tracked_data.dat",images, tracks, times, laser_powers, traps)  
    
Laser powers and trap locations are usually read from the frame metadata. None is available in this case, so we specify it manually, and save all the data to a file.
    
.. testcode::

    times,traps,trajectory = plotting.read_file(output_folder + "/tracked_data.dat",1)
    
We read the file back. The shape of parsed data is printed out; if laser powers are not specified for a line, it is discarded.
    
.. testoutput::

    Shape of initial data:  (20, 16)
    Shape of cropped data:  (20, 16)
    
Next, the estimated trap parameters are calculated. You can see that the estimates differ very much from the originally specified values, as 20 points are not enough to compute a meaningful statistic. Trap coefficient estimates are specified in micronewtons per cm, not per m, leading to a 10^-2 scaling factor.
    
.. testcode::

    trap_offsets = offset.single_particle_offset(traps[:,0:2],trajectory[:,::-1])
    coeffs_estimate,rotation_angle,_ = calibration.calibrate(times,trajectory,0.005)
    
    print("Trap offsets in micrometers:", trap_offsets)
    print("Estimate of trap coefficients:", coeffs_estimate)
    print("Trap rotation angle", rotation_angle)
    
.. testoutput::

    Trap offsets in micrometers: [-0.23071602  0.19115445]
    Estimate of trap coefficients: (3.823906204447568e-09, 5.440606237164117e-09)
    Trap rotation angle 0.46742406927807534
    
As a last step, we calculate the estimated forces.
    
.. testcode::
    
    forces,mean_forces = force_calc.calculate(times,trajectory[:,1::-1],traps[:,0:2],coeffs_estimate,rotation_angle)
    
.. testoutput::

    Mean force values in pN:  [ 0.00108228 -0.00112729]

The values are shown in the following plots.

.. plot:: ../examples/pipeline_example.py
