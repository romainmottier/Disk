-- Transcranial focused ultrasound in a seven-tissue head model with a titanium fixation plate (preprint LTS_HHO, 6.7),
-- polytopal variant: Voronoi mesh of generate_head_model_voronoi.py (no square cell: Voronoi cells of the centres of the
-- quadtree cells of generate_head_model_fitted.py jittered by 0.3 cell size, cut along the smoothed skull and the exact
-- plate and screws, --delta 0.3 --min-width 0.15), axial T1 MRI slice (subject 1, slice 80). LI with the implicit region
-- chosen by the cell width (time.fine_measure = "width": w_T < 0.75 mm). Focused transducer (arc of radius 80 mm, half
-- aperture 25 deg, axis through the plate) emitting a 3-cycle tone burst at 250 kHz towards a focus 45 mm below the outer
-- surface of the skull (frontal gray matter).
-- Lengths in mm, times in us, densities in g/cm^3, velocities in mm/us (stresses in GPa).
formulation = "first_order"
space   = { k = 3, order = "equal", stabilization = "least_squares", scaling = "1" }
time    = { scheme = "li", final_time = 80.0, n_steps = 1394, fine_size = 0.75, li_system = "constraint", li_coupling = "stage",
            li_condensation = true, fine_measure = "width", level_measure = "h/c", level_tolerance = 0.1, level_size = "width",
            level_ratios = {2, 2, 16} }
mesh    = { type = "polygons", level = 1, files = { "/home/mottie0000/Github/Diskpp/Disk/apps/wave_propagation/src/mesh_generation/brain_mri/meshes/head_s80_vor_plate.txt" } }
problem = { interface = "tags", solution = "ricker",
            zones = { [1] = { medium = "elastic",  rho = 1.85,  vp = 2.8,    vs = 1.4 },   -- cortical bone
                      [2] = { medium = "elastic",  rho = 1.70,  vp = 2.3,    vs = 1.1 },   -- cancellous bone (diploe)
                      [3] = { medium = "acoustic", rho = 1.00,  vp = 1.5 },                -- water
                      [4] = { medium = "acoustic", rho = 1.09,  vp = 1.61 },               -- scalp
                      [5] = { medium = "acoustic", rho = 1.007, vp = 1.5045 },             -- cerebrospinal fluid
                      [6] = { medium = "acoustic", rho = 1.045, vp = 1.5 },                -- gray matter
                      [7] = { medium = "acoustic", rho = 1.041, vp = 1.5525 },             -- white matter
                      [8] = { medium = "elastic",  rho = 4.43,  vp = 6.1,    vs = 3.1 } }, -- titanium (Ti-6Al-4V)
            ricker = { shape = "focused", center = {120.5, 179.6}, radius = 80.0, axis = 44.3, aperture = 25.0,
                       cycles = 3.0, frequency = 0.25, velocity = 1.5, density = 1.0, amplitude = 1.0 } }
study   = { type = "run" }
solver  = { matrix_free = true }
output  = { file = "c.txt", table = "c.dat", cpu = true, energy = true, sensor_every = 1, profile = true,
            sensors = { { name = "focus",        x = 120.5, y = 179.6, fields = { "acoustic_pressure" } },
                        { name = "water",        x = 167.0, y = 225.0, fields = { "acoustic_pressure" } },
                        { name = "brain_entry",  x = 137.0, y = 195.5, fields = { "acoustic_pressure" } },
                        { name = "beyond_focus", x = 106.2, y = 165.6, fields = { "acoustic_pressure" } },
                        { name = "skull",        x = 150.6, y = 208.9, fields = { "elastic_velocity", "elastic_stress" } } },
            snapshots = { 20.0, 28.0, 36.0, 44.0, 52.0, 60.0, 68.0, 76.0 }, snapshot_samples = 4, snapshot_elastic = "stress",
            peak_samples = 4 }
