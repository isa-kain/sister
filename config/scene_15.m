%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SISTER: Starshade Imaging Simulation Toolkit for Exoplanet Reconnaissance
%
% Copyright (c) <2019>, California Institute of Technology ("Caltech"). U.S. Government sponsorship acknowledged.
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
%
% • Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
%
% • Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
%
% • Neither the name of Caltech nor its operating division, the Jet Propulsion Laboratory, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = scene_15( opt )
% Similar example as scene_3.m but for the Remote Occulter: a 99 m diameter starshade used together with a large ground telescope. Notice that the PSF response does not include effects from atmospheric turbulence, so that this simulation is not representative of a ground telescope observation without adding those effects.
% How to add sky brightness for a ground telescope.
% How to add a particular label to the output files
% How to add an exoplanet with some given flux ratio (contrast). 
% How to choose a particular star type with its main characteristics.
% How to add the solar glint coming from the scattered Sun’s light at the petals. And
% How to turn on some noise contributions.
%
% Run it as:
% clear opt;opt.run='scene_15';sister( opt ); % Green band, with noise
%

% General optional label that will be appended to the output images and data files in addition of any other labels set by the imaging software.
  if ~isfield( opt, 'add_label' )
  opt.add_label = 'RO' ;
  end

%%%%%%%%%%%%%%%%%%
% Starshade mode %
%%%%%%%%%%%%%%%%%%
% opt.starshade.mode='non-spinning' ;  % Default is 'spinning'
opt.starshade.nominal_filename = 'UH17' ;
% Getting several properties directly from the matlab file associated with the occulter
opt = new_occulter_from_matlab_file( opt ) ;

%%%%%%%%%%%%%
% Telescope %
%%%%%%%%%%%%%
% For NI2 (WFIRST), SISTER takes into account automatically the collecting area blocked by the secondary (0.32x diameter of the primary mirror)
% and the struts (7.17% loss). Total loss of 17.41%. For TV3 (HabEx), SISTER takes into account the secondary mirror (450 mm. The primary is 4000 mm).
% There's no specs about the struts in the HabEx public report of 09/2019, although the effect should really be of second order since the secondary
% is blocking 1.3% only. These settings are consistent with  WFIRST and HabEx mission specs.

% Sampling of the pupil (16 is fast enough for creating the PSF basis on a laptop)
opt.Nx_pupil_pix = 100 ;

%%%%%%%%%%%%%%%%%%%%
% Wavelength range %
%%%%%%%%%%%%%%%%%%%%
% Instrument pass band:
opt.lambda_band_nm_min = 400 ;
opt.lambda_band_nm_max = 700 ;
opt.delta_lambda_psf_nm = 10 ;
% Spatial extension of the PSF in terms of lambda/D
opt.n_lambda_over_d = 40 ; % default is 7
% The distance from which the PSF becomes stationary. For WFIRST (NI2) and HabEx (NW2) there's a default value for each band coming from some analysis. For a new occulter, the user has to provide the value based on some own analysis.
opt.r_stationary_mas = 120 ; % Roman -> 65, from Starshade Rendezvous report

% Imaging bands:
% First wavelength, in nm, to be simulated
  if ~isfield( opt, 'lambda_imaging_1_nm' )
  opt.lambda_imaging_1_nm = 400 ; % Default value depends on the telescope. 
  end
% Last wavelength, in nm, to be simulated
  if ~isfield( opt, 'lambda_imaging_2_nm' )
  opt.lambda_imaging_2_nm = 700 ; % Default value depends on the telescope. 
  end
% Wavelength steps to be simulated
  if ~isfield( opt, 'delta_lambda_imaging_nm' )
  opt.delta_lambda_imaging_nm = 10 ; % Set it here to a constant value, or provide an array of wavelengths (see sister_imaging.m/get_lambda_scene.m)
  end

%%%%%%%%%%%%
% Detector %
%%%%%%%%%%%%
% Pixel scale of the camera. Default:
opt.pix_camera_mas = 1 ;

%%%%%%%%%%%%%%%%%%
% Detector noise %
%%%%%%%%%%%%%%%%%%
opt.noise.do = 1 ; % Shot noise, read noise and dark current generation. Default is 0, not generated.
% Total exposition time in seconds
opt.noise.exp_time_total_sec = 120 ; % Default 3600 sec
% Time per frame in seconds
opt.noise.exp_time_frame_sec = 600 ; % Default 600 sec
% Detector's gain (conversion of photons to photo-electrons)
opt.noise.gain = 1 ; % Default is 1
% Read noise in photo-electron units per detector pixel
opt.noise.read_e = 3 ; % Default is 3
% Dark current in units of photo-electrons per second, and detector pixel.
opt.noise.dark_e_s = 1e-3 ; % Default is 1e-3
% QE (of the detetcor only)
opt.qe = 1 ;
% Optical throughput
opt.optical_throughput = 0.5 ;

%%%%%%%%%%%%%%%%%%%%
% SISTER Data Cube %
%%%%%%%%%%%%%%%%%%%%
opt.cube.do = 1 ; % Default is 0. SISTER does not output a FITS file with th spectral results of the simulation.
opt.cube.fits = 0 ; % Whether the output is a FITS or Matlab file. By default, a Matlab file.
opt.cube.ra_mas = 1000 ; % By default, opt.cube.ra_mas=0; centered with the scene.
opt.cube.dec_mas = -200 ; % By default, opt.cube.dec_mas=0; centered with the scene.

%%%%%%%%%%%%%%%
% Solar glint %
%%%%%%%%%%%%%%%

opt.solar_glint.do = 1 ; % Default 0. Not included.
% Solar angle relative to Starshade and telescope
opt.solar_glint.phi_deg = 60 ; % Default is 60
opt.solar_glint.alpha_deg = 30 ; % Rotation of sun around the plane. Zero is along x axis. Default is 0 deg
% Radius (in m) used to define the 'trailing edge' of the petals
opt.solar_glint.r_m = 28.6 ; % meters. Default: 7.5 m (adequate for NI2, 26 m diameter Starshade)

%%%%%%%%%%%%%%%%%%%%
% Creating a scene %
%%%%%%%%%%%%%%%%%%%%
opt.scene.do = 1 ; % Default value is 0. It would read a pre-existing scene
% Field of View of the astrophysical scene in mas
opt.scene.fov_diam_mas = 255 ; % Default is 5000 mas
% Pixel scale on the scene
opt.pix_scene_mas = 1 ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sky model for ground telescope %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.ground_telescope.do = 1 ; % Default, 0. It is not added.
% The sky model includes scattered Moonlight, local zodiacal light, scattered starlight, emisison lines of Upper Atmosphere, airglow (residual continuum), and Molecular Emission of Lower Atmosphere (the latter essentially zero below 1.7 micron though). One has to choose among three options for the moon phase: 'full', 'half', or 'new'.
opt.ground_telescope.moon_phase = 'half' ; % default, 'half'


%%%%%%%%
% Star %
%%%%%%%%
opt.star.name = 'Sun' ;
opt.star.type = 'Sun' ; % Default is 'Sun'. Non case-sensitive.
opt.star.distance_to_earth_pc = 1.3 ; % Default is 10 pc
% Besides the star type, it's necessary to set its brightness
opt.star.app_mag_v = 5.98 ; % Sun at 17 pc. Default 4.81 (Sun at 10 pc, http://mips.as.arizona.edu/~cnaw/sun.html)

%%%%%%%%%%%
% Planets %
%%%%%%%%%%%

% Example of adding 'static' planets. That is, without an orbital motion
opt.planets.add.do = 1 ; % Default is 0
% RA/DEC are the 'x/y' coordinates on the final image. Reminder about astronomical convention: RA>0 left hand side, RA<0 right hand side, DEC>0 top, DEC<0 bottom.
% Notice that the positions are given in mas. It's possible to instead set it in AU, but the distance to the star in pc has to be known. Above there's an example where the Sun is placed at some distance, or if ExoCat is used, one may read what's the distance to the system. Alternatively, one may use the Solar System planets: 'Venus', ..., 'Neptune' (another configuration file).
% opt.planets.add.pos_arc_ra_mas = 65 / sqrt( 2 ) ;
% opt.planets.add.pos_arc_dec_mas = 65 / sqrt( 2 ) ;
% opt.planets.add.flux_ratio = 4.3e-11 ; 
opt.planets.add.planet_type = {'Earth', 'Venus', 'Jupiter'};
opt.planets.add.pos_arc_ra_au = [1, 1, 1] ;
opt.planets.add.pos_arc_dec_au = [1, 0, -1] ;
opt.planets.add.phase_angle_deg = [ 90, 120, 150 ] ;

%%%%%%%%%%%%
% Exo-zodi %
%%%%%%%%%%%%

% Times solar zodiacal light
opt.exozodi.factor = 5 ; 
% Inclination as viewed from the telescope (it does not need to be the same value as the orbital inclination)
opt.exozodi.inclination_deg = 60 ; 
% Position angle of the exoplanetary system
opt.exozodi.position_angle_deg = 52 ;

%%%%%%%%%%%%%%%%
% Plot options %
%%%%%%%%%%%%%%%%
opt.plot.do = 1 ; % Default opt.plot.do = 0 (no plotting)
% Add contour about starshade petals
opt.plot.starshade_circle.do = 1 ;
% Signal+noise (see starshade_band_imaging.m, % Combinations of signal and noise, for other options. Multiple options are fine: e.g., [ 1, 3, 4 ]. 3 figures will be created )
opt.plot.combination_list = 1 ;
% Optional labels for planets (as many as planets)
opt.plot.planet_label = { 'E', 'V', 'J' } ; % Default B, C, ...
% Positional shift for the label of each planet wrt the position of the planet (times FWHM)
opt.plot.planet_label_pos_fwhm = [ 0.8, 0.8, 0.8 ] ; % Default all 1
% A maximum limit for the linear color scale (noise is the RMS of the noise contribution over all the FOV). If the image is noiseless, the plotting script
% will choose a default min/max. Alternatively, use opt.plot.min and opt.plot.max to set the range of the plot to any specified value.
opt.plot.n_sigma_noise = 10 ; % Default is 10
% Zoom on the FOV to be plotted
opt.plot.zoom_list = 1 ; % Default is 1 (no zooom in)
% Storing the (individual) imags
opt.save_plot = 1 ; % If the image is to be stored (by default, they are. Set it to 0 not to store them.)

%%%%%%%%%%%%%%%%%%%%
% Storing products %
%%%%%%%%%%%%%%%%%%%%

% Saving the simulation
  if ~isfield( opt, 'save_output' )
  opt.save_output = 0 ;
  % Input/output with FITS files instead of matlab files (used in starshade_band_imaging.m)
  opt.fitsio = 0 ; % Default is 0, and matlab files are used for I/O.
  end
