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
function opt = scene_6( opt )
% Analogous example as scene_3.m: the only change is the solar angle phi for the solar glint (50 degrees, instead of 60).
% Example that shows how to use different features:
% How to add a particular label to the output files
% How to add an exoplanet with some given flux ratio (contrast). 
% How to choose a particular star type with its main characteristics.
% How to add a simple model of exozodiacal light.
% How to add the solar glint coming from the scattered Sun’s light at the petals. And
% How to turn on some noise contributions.
%
% Parameters shown below are the subset of all possibilities used in this run.
% Hint: exit debug mode (K> in a matlab session) typing 'dbquit all'
%

% General optional label that will be appended to the output images and data files in addition of any other labels set by the imaging software.
  if ~isfield( opt, 'add_label' )
  opt.add_label = 'Formalhaut_B' ;
  end

%%%%%%%%%%%%%%%%%%
% Starshade mode %
%%%%%%%%%%%%%%%%%%
% opt.starshade.mode='non-spinning' ;  % Default is 'spinning'

%%%%%%%%%%%%%
% Telescope %
%%%%%%%%%%%%%
% For NI2 (WFIRST), SISTER takes into account automatically the collecting area blocked by the secondary (0.32x diameter of the primary mirror)
% and the struts (7.17% loss). Total loss of 17.41%. For TV3 (HabEx), SISTER takes into account the secondary mirror (450 mm. The primary is 4000 mm).
% There's no specs about the struts in the HabEx public report of 09/2019, although the effect should really be of second order since the secondary
% is blocking 1.3% only. These settings are consistent with  WFIRST and HabEx mission specs.

% Sampling of the pupil (16 is fast enough for creating the PSF basis on a laptop)
opt.Nx_pupil_pix = 16 ;

%%%%%%%%%%%%%%%%%%%%
% Wavelength range %
%%%%%%%%%%%%%%%%%%%%

% First wavelength, in nm, to be simulated
  if ~isfield( opt, 'lambda_imaging_1_nm' )
  opt.lambda_imaging_1_nm = 615 ; % Default value depends on the telescope. 
  end
% Last wavelength, in nm, to be simulated
  if ~isfield( opt, 'lambda_imaging_2_nm' )
  opt.lambda_imaging_2_nm = 800 ; % Default value depends on the telescope. 
  end
% Wavelength steps to be simulated
  if ~isfield( opt, 'delta_lambda_imaging_nm' )
  opt.delta_lambda_imaging_nm = 10 ; % Set it here to a constant value, or provide an array of wavelengths (see sister_imaging.m/get_lambda_scene.m)
  end

%%%%%%%%%%%%
% Detector %
%%%%%%%%%%%%
% Pixel scale of the camera. Default:
% opt.pix_camera_mas = 0.4 * opt.lambda_imaging_1_nm * 1e-9 / opt.diameter_telescope_m * 180 / pi * 3600e3 ;

%%%%%%%%%%%%%%%%%%
% Detector noise %
%%%%%%%%%%%%%%%%%%
  if ~isfield( opt, 'noise' )
  opt.noise.do = 0 ; % Shot noise, read noise and dark current generation. Default is 0, not generated.
  end
% Total exposition time in seconds
opt.noise.exp_time_total_sec = 24 * 3600 ; % Default 3600 sec
% Time per frame in seconds
opt.noise.exp_time_frame_sec = 100 ; % Default 600 sec
% Detector's gain (conversion of photons to photo-electrons)
opt.noise.gain = 1 ; % Default is 1
% Read noise in photo-electron units per detector pixel
opt.noise.read_e = 3 ; % Default is 3
% Dark current in units of photo-electrons per second, and detector pixel.
opt.noise.dark_e_s = 1e-3 ; % Default is 1e-3

%%%%%%%%%%%%%%%
% Solar glint %
%%%%%%%%%%%%%%%

opt.solar_glint.do = 1 ; % Default 0. Not included.
% Solar angle relative to Starshade and telescope
opt.solar_glint.phi_deg = 50 ; % Default is 60
opt.solar_glint.alpha_deg = 0 ; % Rotation of sun around the plane. Zero is along x axis. Default is 0 deg
% Radius (in m) used to define the 'trailing edge' of the petals
opt.solar_glint.r_m = 7.5 ; % meters. Default: 7.5 m (adequate for NI2, 26 m diameter Starshade)
% Whether edges of the petals are shredded
opt.solar_glint.stealth = 0 ;
  if isfield( opt, 'solar_glint' )
    if ( opt.solar_glint.do )
    undrscr = '' ;
      if numel( opt.add_label ) > 0, undrscr = '_' ; end
    opt.add_label = [ opt.add_label undrscr 'solar_glint' ] ;
    end
  end

%%%%%%%%%%%%%%%%%%%%
% Creating a scene %
%%%%%%%%%%%%%%%%%%%%
opt.scene.do = 1 ; % Default value is 0. It would read a pre-existing scene
% Field of View of the astrophysical scene in mas
opt.scene.fov_diam_mas = 2000 ; % Default is 5000 mas

%%%%%%%%%%%%%%%%%%%%%%%%
% Local zodiacal light %
%%%%%%%%%%%%%%%%%%%%%%%%
% Recall that, by default, local zodi is not added. It should be present in any astrophysical simulation. To add it, set opt.local_zodi.do=1; in the call or configuration file.
% This option will not be used unless opt.local_zodi.do=1;
opt.local_zodi.mag_v_arcsec2 = 22.5 ; % The surface brightness of the local zodiacal light, expressed in V mag per arcsec^2. Default value is 23.

%%%%%%%%
% Star %
%%%%%%%%
opt.star.name = 'Formalhaut' ;
opt.star.type = 'A5V' ; % Default is 'Sun'. Non case-sensitive.
opt.star.distance_to_earth_pc = 7.7 ; % Default is 10 pc
% Besides the star type, it's necessary to set its brightness
opt.star.app_mag_v = 1.16 ; % Default 4.81 (Sun at 10 pc, http://mips.as.arizona.edu/~cnaw/sun.html)

%%%%%%%%%%%
% Planets %
%%%%%%%%%%%

% Example of adding 'static' planets. That is, without an orbital motion
opt.planets.add.do = 1 ; % Default is 0
% RA/DEC are the 'x/y' coordinates on the final image. Reminder about astronomical convention: RA>0 left hand side, RA<0 right hand side, DEC>0 top, DEC<0 bottom.
% Notice that the positions are given in mas. It's possible to instead set it in AU, but the distance to the star in pc has to be known. Above there's an example where the Sun is placed at some distance, or if ExoCat is used, one may read what's the distance to the system. Alternatively, one may use the Solar System planets: 'Venus', ..., 'Neptune' (another configuration file).
opt.planets.add.pos_arc_ra_mas = 154.84 ;
opt.planets.add.pos_arc_dec_mas = 506.47 ;
opt.planets.add.flux_ratio = 4.3e-11 ; % Earth-like planet with radius twice as Earth, 0.2 geometric albedo and in quadrature (Lambertian)

%%%%%%%%%%%%
% Exo-zodi %
%%%%%%%%%%%%

% Times solar zodiacal light
opt.exozodi.factor = 7 ; 
% Inclination as viewed from the telescope (it does not need to be the same value as the orbital inclination)
opt.exozodi.inclination_deg = 46.5 ; 
% Position angle of the exoplanetary system
opt.exozodi.position_angle_deg = 52 ;
% radius at which the zodiacal light becomes half its peak value
opt.exozodi.radius_au = 2 ;

%%%%%%%%%%%%%%%%
% Plot options %
%%%%%%%%%%%%%%%%
opt.plot.do = 1 ; % Default opt.plot.do = 0 (no plotting)
% Signal+noise (see starshade_band_imaging.m, % Combinations of signal and noise, for other options. Multiple options are fine: e.g., [ 1, 3, 4 ]. 3 figures will be created )
opt.plot.combination_list = 1 ;
% Optional labels for planets (as many as planets)
opt.plot.planet_label = { 'B' } ; % Default B, C, ...
% Positional shift for the label of each planet wrt the position of the planet (times FWHM)
opt.plot.planet_label_pos_fwhm = [ 0.8 ] ; % Default all 1
% A maximum limit for the linear color scale (noise is the RMS of the noise contribution over all the FOV). If the image is noiseless, the plotting script
% will choose a default min/max. Alternatively, use opt.plot.min and opt.plot.max to set the range of the plot to any specified value.
opt.plot.n_sigma_noise = 10 ; % Default is 10
opt.plot.min = 0 ;
opt.plot.max = 8000 ;
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
