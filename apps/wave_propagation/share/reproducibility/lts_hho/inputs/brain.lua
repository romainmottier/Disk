-- Elasto-acoustic application: axial head slice (MRI segmentation of ~/Downloads/mri_dataset, subject 1, slice 80),
-- quadtree meshes of mesh_generation/brain_mri (level 1: layered skull at the segmentation resolution, 8 mm background;
-- level 2: homogeneous skull at 2 mm with a temporal window refined to 0.5 mm; level 3: homogeneous skull at 2 mm and
-- a titanium implant of radius 0.65 mm at (105, 140) meshed at 0.125 mm); skull = elastic, soft tissues and water =
-- acoustic; lengths in mm, times in us, densities in g/cm^3, velocities in mm/us; initial Ricker pressure pulse in the
-- water facing the temporal window; homogeneous Dirichlet data on the boundary of the water box
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "erk", final_time = 100.0, n_steps = 3000, level_measure = "h/c" }
mesh    = { type = "polygons", level = 1, files = { "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/brain_mri/meshes/head_s80_layered_b8.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/brain_mri/meshes/head_s80_window_w8.txt", "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/brain_mri/meshes/head_s80_implant_i8.txt" } }
problem = { interface = "tags", solution = "ricker",
            zones = { [1] = { medium = "elastic",  rho = 1.85, vp = 2.8,  vs = 1.4 },    -- cortical (or homogeneous) bone
                      [2] = { medium = "elastic",  rho = 1.60, vp = 2.3,  vs = 1.1 },    -- cancellous bone
                      [3] = { medium = "acoustic", rho = 1.00, vp = 1.5 },               -- water
                      [4] = { medium = "acoustic", rho = 1.04, vp = 1.55 },              -- soft tissues
                      [5] = { medium = "elastic",  rho = 4.5,  vp = 6.1,  vs = 3.1 } },  -- titanium implant (level 3)
            ricker = { field = "acoustic_pressure", center = {28.0, 135.0}, frequency = 0.1, velocity = 1.5, amplitude = 1.0 } }
study   = { type = "run" }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat", cpu = true, energy = true, sensor_every = 1,
            sensors = { { name = "water_front", x = 28.0, y = 160.0, fields = { "acoustic_pressure" } },
                        { name = "skull_window", x = 47.0, y = 135.0, fields = { "elastic_velocity" } },
                        { name = "brain_center", x = 112.0, y = 132.0, fields = { "acoustic_pressure" } },
                        { name = "water_back", x = 205.0, y = 132.0, fields = { "acoustic_pressure" } },
                        { name = "near_implant", x = 108.0, y = 140.0, fields = { "acoustic_pressure" } } },
            snapshots = { 20.0, 40.0, 60.0, 80.0, 100.0 } }
