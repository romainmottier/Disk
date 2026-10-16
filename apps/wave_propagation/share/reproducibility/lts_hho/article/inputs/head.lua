-- 6.7 transcranial ultrasound in a seven-tissue head model: axial T1 MRI slice (subject 1, slice 80 of
-- ~/Downloads/mri_dataset), tissues from the FSL segmentation and the T1 signal, quadtree meshes of
-- mesh_generation/brain_mri/generate_head_model.py (2 mm in the fluids, 1 mm at the tissue interfaces and in the bone;
-- level 2: titanium fixation plate (0.6 mm) and two screws on the frontal bone, meshed at 0.125 mm).
-- Bone and titanium elastic, soft tissues and water acoustic; lengths in mm, times in us, densities in g/cm^3,
-- velocities in mm/us. Initial pressure pulse (theta exp(-pi^2 r^2 / lambda^2), lambda = 7.5 mm, i.e. 200 kHz in
-- water) in the water facing the left temporal bone; homogeneous Dirichlet data on the boundary of the water box.
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "erk", final_time = 120.0, n_steps = 4000, level_measure = "h/c" }
mesh    = { type = "polygons", level = 1, files = { "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/brain_mri/meshes/head_s80_intact.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/brain_mri/meshes/head_s80_plate.txt" } }
problem = { interface = "tags", solution = "ricker",
            zones = { [1] = { medium = "elastic",  rho = 1.85,  vp = 2.8,    vs = 1.4 },   -- cortical bone
                      [2] = { medium = "elastic",  rho = 1.70,  vp = 2.3,    vs = 1.1 },   -- cancellous bone (diploe)
                      [3] = { medium = "acoustic", rho = 1.00,  vp = 1.5 },                -- water
                      [4] = { medium = "acoustic", rho = 1.09,  vp = 1.61 },               -- scalp
                      [5] = { medium = "acoustic", rho = 1.007, vp = 1.5045 },             -- cerebrospinal fluid
                      [6] = { medium = "acoustic", rho = 1.045, vp = 1.5 },                -- gray matter
                      [7] = { medium = "acoustic", rho = 1.041, vp = 1.5525 },             -- white matter
                      [8] = { medium = "elastic",  rho = 4.43,  vp = 6.1,    vs = 3.1 } }, -- titanium (Ti-6Al-4V)
            ricker = { field = "acoustic_pressure", center = {8.0, 118.5}, frequency = 0.2, velocity = 1.5, amplitude = 1.0 } }
study   = { type = "run" }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat", cpu = true, energy = true, sensor_every = 1,
            sensors = { { name = "water",        x = 14.5,  y = 118.5, fields = { "acoustic_pressure" } },
                        { name = "skull_left",   x = 29.5,  y = 118.5, fields = { "elastic_velocity" } },
                        { name = "brain_center", x = 98.5,  y = 118.5, fields = { "acoustic_pressure" } },
                        { name = "under_plate",  x = 128.5, y = 178.5, fields = { "acoustic_pressure" } },
                        { name = "skull_right",  x = 161.5, y = 118.5, fields = { "elastic_velocity" } },
                        { name = "water_right",  x = 185.5, y = 118.5, fields = { "acoustic_pressure" } } },
            snapshots = { 20.0, 40.0, 60.0, 80.0, 100.0, 120.0 } }
