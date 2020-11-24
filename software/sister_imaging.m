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
function [ scn_dt_cnv_e, scn_dt_cnv_no_plnt_e, scn_ns_e, opt_img, scn_dt_cnv_lcl_zd_e ] = sister_imaging( opt )
% Script to produce a full set of images of an scene as observed by WFIRST-S
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

% Debug mode yes/no
dbstop if error

% Code
  if ~exist( 'opt', 'var' )
  opt = [] ;
  end

% Default options
tic
opt_img = get_imaging_options( opt ) ;
toc

% Keplerian orbits: output projected RA, DEC in AU
  if ( opt_img.kplr.do ) && sum( isnan( opt_img.plnt.add.flx_rt_arry( : ) ) ) % If it is two scenes, only evaluate it once at the beginning
  opt_img = get_imaging_options( get_kepler_options( opt_img ) ) ;
  % Option to retrieve the Keplerian orbit for statistical analysis, without actually generating the scene simulation
    if ( opt_img.kplr.gt_orbt )
    scn_dt_cnv_e = NaN ; scn_dt_cnv_no_plnt_e = NaN ; scn_ns_e = NaN ;
    disp( 'Getting the Keplerian orbit ... End' )
    end
  end

if ( opt_img.return ), return ; end

% Starshade mode
disp( sprintf( 'Considering a %s Starshade', opt_img.starshade.mode ) )

% Bands to consider
opt_img.n_bnd = numel( opt_img.lmbd_img_1_nm ) ;

% Looping over bands
  for i_bnd = 1 : opt_img.n_bnd
   tm_0 = tic ;
   % Adding a label for the band in use
     if ( i_bnd  == 1 )
     tg_nm_tmp = opt_img.tg_nm ;
     end
   opt_img.tg_nm = sprintf( '%s_%04i_%04i_nm', tg_nm_tmp, opt_img.lmbd_img_1_nm( i_bnd ), opt_img.lmbd_img_2_nm( i_bnd ) ) ;
   if ~( opt_img.cluster )
   [ scn_dt_cnv_e, scn_dt_cnv_no_plnt_e, scn_ns_e, scn_dt_cnv_lcl_zd_e, opt_img ] = sister_imaging_band( i_bnd, opt_img ) ;
   else
   send_sister_imaging_band( i_bnd, opt_img ) ;
   end % opt_img.cluster
  disp( sprintf( 'Time spent imaging band %i/%i was %3.1f min', i_bnd, opt_img.n_bnd, toc( tm_0 ) / 60 ) )
  end % i_bnd
disp( 'End of the observations by Starshade' )

%%%%%%%%%%%%%%%%%
% Sub-functions %
%%%%%%%%%%%%%%%%%

% Getting the necessary data of the wavelength array associated with the scene
function [ lambda_array_scene opt ] = get_lambda_scene( opt, i_bnd )

% If it is 'file', read it off. Otherwise, it should be an array directly defined in opt.lmbd_arry_scn, or in the case of a custom scene, the *same* set of wavelength as the one used for imaging
  if ( opt.scn.do )
    if ~ischar( opt.dlt_lmbd_img_nm )
    lambda_array_scene = opt.lmbd_img_1_nm( i_bnd ) : opt.dlt_lmbd_img_nm : opt.lmbd_img_2_nm( i_bnd ) ;
    else
      if ~numel( strfind( lower( opt.dlt_lmbd_img_nm ), 'r' ) )
      disp( 'Please set opt.delta_lambda_imaging_nm to ''Rxx'' for some resolving power')
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      else
      opt.rslvng_pwr = str2num( opt.dlt_lmbd_img_nm( 2 : end ) ) ;
      lambda_array_scene = wavelength_grid( opt ) ;
      % Last wavelength is the end of the band
      lambda_array_scene = lambda_array_scene( 1 : end - 1 ) ;
      end
    end
  else
  lambda_array_scene = opt.lmbd_arry_scn ;
  end
  if strcmp( opt.lmbd_arry_scn, 'file' )
  opt.fl_lmbd_scn = [ opt.scn_dr opt.scn.nm '_misc.mat' ] ; % For miscellaneous data
    if exist( opt.fl_lmbd_scn, 'file' ) == 2
    % It should bring lambda_array_scene
    load( opt.fl_lmbd_scn )
    opt.lmbd_arry_scn = lambda_array_scene ;
    else
    disp( sprintf( 'The scene file with some miscellaneous data was not found: %s, Stopped.', opt.fl_lmbd_scn ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  end

  if ~isvector( opt.lmbd_arry_scn )
  disp( 'The array of wavelength values associated with the scene could not be read correctly. Please check the file with the wavelength array data. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

% Default options
function opt = get_imaging_options( opt )

% Getting all input fields in lowecase mode
opt = lower_opt( opt ) ;

% Saving some memory. In most applications, 7 digit precision is enough
  if ~isfield( opt, 'single_precision' )
  opt.sngl_prcsn = 1 ; % By default, store products in single precision
  else
  opt.sngl_prcsn = opt.single_precision ; 
  end


% Starshade mode: non-spinning/spinning
  if ~isfield( opt, 'starshade' )
  opt.starshade.mode = 'spinning' ;
  else
    if isfield( opt.starshade, 'mode' )
    opt.starshade.mode = lower( opt.starshade.mode ) ;
    else
    opt.starshade.mode = 'spinning' ;
    end
  end

% Check
  if ~strcmp( opt.starshade.mode, 'non-spinning' ) && ~strcmp( opt.starshade.mode, 'spinning' )
  disp( 'Either choose opt.starshade.mode=''spinning'' or opt.starshade.mode=''non-spinning''' )
  return
  end

% Nominal occulter filename. If an occulter has some small perturbations with respect another one, only the convolution with the star has noticeable changes, and that new locus is set by opt.locus.perturbed_filename.
  if ~isfield( opt.starshade, 'nominal_filename' )
  opt.starshade.nominal_filename = 'NI2' ; % NI2 corresponds to WFIRST-S, and TV3 to HabEx (latest design, after 2019)
  end

% Root directory for I/O
  if ~isfield( opt, 'root_dir' )
  opt.rt_dr = sister_installation_path() ;
  else
  opt.rt_dr = opt.root_dir ;
  end

% I/O directory of the scenes
  if ~isfield( opt, 'scene_dir' )
  opt.scn_dr = sprintf( '%input_scenes/', opt.rt_dr ) ;
  else
  opt.scn_dr = opt.scene_dir ;
  end

  if ~isfield( opt, 'path_occulter' )
  opt.path_occulter = [ opt.scene_dir 'locus/in/' ] ;
  end

% Size of the pupil data in pixels (square)
  if ~isfield( opt, 'nx_pupil_pix' ) 
  opt.nx_pupil_pix = 64 ;
  end

% PS: make the following lines consistent with get_default_options.m
% Replace by a file in FITS or Matlab format with an Nx x Nx array if you want a specific pupil
  if ~isfield( opt, 'pupil_filename' )
  opt.pupil_filename = [ opt.path_occulter 'pupil_D1Kpix_256.fits' ] ;
  end

% For now, if it is not NI2, or some perturbed locus of NI2, we will assume an ideal pupil (circularly symmetric)
  if ~strfind( opt.starshade.nominal_filename, 'NI2' )
  opt.pupil_filename = 'ideal' ;
  end

% Also for spinning starshades, we will consider an ideal pupil (otherwise the PSF basis cannot be 1-dimensional, since the telescope's pupil does not rotate whereas the starshade does.
  if strcmp( lower( opt.starshade.mode ), 'spinning' )
  opt.pupil_filename = 'ideal' ;
  end

% One may redefine the pupil as an ideal circularly symmetric pupil at any time
  if isfield( opt, 'pupil_filename' )
    if strcmp( opt.pupil_filename, 'ideal' )
    opt.pupil_filename= '0' ;
    end
  end

% Size of the secondary if a circular obscuration is fine (no struts). See makeStarshadeImage.m
  if ~isfield( opt, 'secondary_size' )
  opt.secondary_size = 0 ; % (Linear obscuration of the telescope entrance pupil (diameter ratio)
  end

% Some particular cases for the secondary obscuration (it will only have an effect if opt.pupil_filename='0', ideal pupil -see makeStarshadeImage.m). The values below must be the same as in sister_imaging.m
  if numel( strfind( lower( opt.starshade.nominal_filename ), 'ni2' ) )
  opt.secondary_size = 0.417257964788079 ; ; % WFIRST value is 0.32 but accounting for the struts, which cover almost an 7.17% of the collecting area, this is the effective value of an equivalent secondary without struts. % https://wfirst.ipac.caltech.edu/sims/Param_db.html#telescope
  end
  if numel( strfind( lower( opt.starshade.nominal_filename ), 'tv3' ) )
  opt.secondary_size = 0.1125 ; % Table 6.2-2, page 6-4 of https://www.jpl.nasa.gov/habex/pdf/HabEx-Final-Report-Public-Release-LINKED-0924.pdf. No information about struts, but it should be minor compared to the effect of the secondary, which is already a small effect (1-0.1125^2)=0.987.
  end

% Number of petals of the Starshade (default is 24 -NI2, TV3- in makeStarshadeImage.m, only updating it if it is provided)
  if ~isfield( opt.starshade, 'number_of_petals' )
  opt.n_ptl = 24 ;
  else
  opt.n_ptl = opt.starshade.number_of_petals ;
  end

% I/O directory of the PSF basis
  if ~isfield( opt, 'psf_dir' )
  opt.psf_dr = sprintf( '%ssister_basis/', opt.rt_dr ) ;
  else
  opt.psf_dr = opt.psf_dir ;
  end

% Subdirectory for the PSF data
opt.psf_dr = [ opt.psf_dr opt.starshade.mode '/' opt.starshade.nominal_filename '_' num2str( opt.n_ptl ) '_' num2str( opt.nx_pupil_pix ) '/' ] ;

% Allowing the user to set the imaging band sub-folder (this should used only in exceptional cases, where a sub-band is ananalyzed that may belong to different iamging bands and instead of the longest wavelength band (default), the use rwants to impose another band.
  if isfield( opt, 'psf_band' )
  opt.sbdr_psf = opt.psf_band ;
  end

% Output directory for the results
  if ~isfield( opt, 'output_dir' )
  opt.out_dr = sprintf( '%soutput/', opt.rt_dr ) ;
  else
  opt.out_dr = opt.output_dir ;
  end

% Storing the output arrays
  if ~isfield( opt, 'save_output' )
  opt.sv_out = 0 ;
  else
  opt.sv_out = opt.save_output ;
  end

% Storing the output arrays
  if ~isfield( opt, 'save_single_wavelength' )
  opt.sv_out_2 = 0 ;
  else
  opt.sv_out_2 = opt.save_single_wavelength ;
  end

% If FITS files are to be loaded/stored with the output products
  if ~isfield( opt, 'fitsio' )
  opt.ftsio = 0 ;
  else
  opt.ftsio = opt.fitsio ;
  end

% Sentinel to stop the execution if some initial parameter that must be set is not provided
opt.return = 0 ;

% Verbose (a way to control when displaying messages)
  if ~isfield( opt, 'verbose' )
  opt.verbose = 0 ;
  end

% Whether to re-do the computation of the accumulated image in case it exists. Default: don't repeat existing results (useful for only plotting)
  if ~isfield( opt, 'redo' )
  opt.redo = 0 ;
  end

% Whether to re-do the computation of the single wavelength images in case they exist. Default: don't repeat existing results (useful for only plotting)
  if ~isfield( opt, 'redo_2' )
  opt.redo_2 = 0 ;
  end

% Diameter of the primary mirror
  if isfield( opt, 'diameter_telescope_m' )
  opt.dmtr_tlscp_m = opt.diameter_telescope_m ;
  end

% Some well-known cases
  if numel( strfind( lower( opt.starshade.nominal_filename ), 'ni2' ) )
  opt.dmtr_tlscp_m = 2.4 ;
  end
  if numel( strfind( lower( opt.starshade.nominal_filename ), 'nw2' ) ) || numel( strfind( lower( opt.starshade.nominal_filename ), 'tv3' ) )
  opt.dmtr_tlscp_m = 4.0 ;
  end

% Cluster or serial (TBD)
  if ~isfield( opt, 'cluster' )
  opt.cluster = 0 ;
  end

% If individual wavelength slices are to be sent to a cluster (TBD)
  if ~isfield( opt, 'sub_cluster' )
  opt.sub_cluster = 0 ;
  end

% Some general parameters to create a scene
  if ~isfield( opt, 'scene' )
  opt.scn.do = 0 ;
    % Type of scene to read from (or create a custom one)
    if ~isfield( opt.scene, 'name' )
    opt.scn.nm = 'modern_cube_zodi1inc0dist10' ;
    else
    opt.scn.nm = opt.scene.name ;
    end
  else
  % Whether a scene is created instead of reading it from a file
    if ~isfield( opt.scene, 'do' )
    opt.scn.do = 0 ;
      if ~isfield( opt.scene, 'name' )
      opt.scn.nm = 'modern_cube_zodi1inc0dist10' ;
      else
      opt.scn.nm = opt.scene.name ;
      end
    else
    opt.scn.do = opt.scene.do ;
      if ( opt.scn.do )
      opt.scene.name = 'custom_scene' ;
      opt.scn.nm = opt.scene.name ;
      opt.lambda_array_scene = 'same as scene' ; 
      end
    end
    % FOV of the scene in mas (if it is created instead of read from some file)
    if ~isfield( opt.scene, 'fov_diam_mas' )
      if strfind( lower( opt.scn.nm ), 'modern_cube' )
      opt.scn.fov_diam_mas = 9999 ; % The Haystacks project has scenes with 3333x3333 pixels, each of 3 mas.
      else
      opt.scn.fov_diam_mas = 5000 ; % Diameter
      end
    else
    opt.scn.fov_diam_mas = opt.scene.fov_diam_mas ;
    end

  % Global lateral displacement of the scene
  % Horizontal direction in meters
    if ~isfield( opt.scene, 'ra_shift_m' )
    opt.scn.ra_shft_m  = NaN ;
    else
    opt.scn.ra_shft_m  = opt.scene.ra_shift_m ;
    end 
  % Vertical direction in meters
    if ~isfield( opt.scene, 'dec_shift_m' )
    opt.scn.dc_shft_m  = NaN ;
    else
    opt.scn.dc_shft_m  = opt.scene.dec_shift_m ;
    end
  % Horizontal direction in mas
    if ~isfield( opt.scene, 'ra_shift_mas' )
    opt.scn.ra_shft_mas  = NaN ;
    else
    opt.scn.ra_shft_mas  = opt.scene.ra_shift_mas ;
    end
  % Vertical direction in meters
    if ~isfield( opt.scene, 'dec_shift_mas' )
    opt.scn.dc_shft_mas  = NaN ;
    else
    opt.scn.dc_shft_mas  = opt.scene.dec_shift_mas ;
    end
  % Checks
    if ~isnan( opt.scn.ra_shft_m ) && ~isnan( opt.scn.ra_shft_mas )
    disp( 'Both opt.scene.ra_shift_m and opt.scene.ra_shift_mas are set. Please set only one option. Stopped.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
    if ~isnan( opt.scn.dc_shft_m ) && ~isnan( opt.scn.dc_shft_mas )
    disp( 'Both opt.scene.dec_shift_m and opt.scene.dec_shift_mas are set. Please set only one option. Stopped.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
    
    if ( ~isnan( opt.scn.ra_shft_m ) && isnan( opt.scn.dc_shft_m ) ) 
    disp( 'Please set both opt.scene.ra_shift_m and opt.scene.dec_shift_m. The latter is not set. Stopped.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end

    if ( isnan( opt.scn.ra_shft_m ) && ~isnan( opt.scn.dc_shft_m ) )
    disp( 'Please set both opt.scene.dec_shift_m and opt.scene.ra_shift_m. The latter is not set. Stopped.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end

    if ( ~isnan( opt.scn.ra_shft_mas ) && isnan( opt.scn.dc_shft_mas ) )
    disp( 'Please set both opt.scene.ra_shift_mas and opt.scene.dec_shift_mas. The latter is not set. Stopped.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end

    if ( isnan( opt.scn.ra_shft_mas ) && ~isnan( opt.scn.dc_shft_mas ) )
    disp( 'Please set both opt.scene.dec_shift_mas and opt.scene.ra_shift_mas. The latter is not set. Stopped.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end

  end % opt.scene

% Making sure the units of the scene are provided
% Preferred use: opt.scene.units, but it's also possible to set scene_units
  if isfield( opt, 'scene_units' )
  opt.scene.units = opt.scene_units ;
  end

  if ~isfield( opt.scene, 'units' )
    if strfind( opt.scn.nm, 'modern_cube' )
    opt.scn.unts = 'Jy' ; % The data cubes of the Haystacks project is provided in units of Jy (and cubes in wavelength). 
    else
    opt.scn.unts = 'w/m2/um' ; % Default spectral irradiance in units of W/m^2/micro-meter
    end
  else
    if ~strcmp( lower( opt.scene.units ), 'jy' ) && ~strcmp( lower( opt.scene.units ), 'w/m2/um' )
    disp( 'Set opt.scene.units=''Jy'' or opt.scene.units=''W/m2/um'' and make sure the scene flux units are as expected.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  opt.scn.unts = opt.scene.units ; 
  % Remark because the data cubes of the Haystacks project is provided in units of Jy (and cubes in wavelength).
    if strcmp( opt.scn.nm, 'modern_cube' ) && strcmp( lower( opt.scene.units ), 'w/m2/um' )
    disp( 'REMARK: Units were set to be W/m^2/micro-meter in the configuration file, but we are using data cubes from the Haystack project, which are in Jy. Continuing with Jy. Please check whether the settings of the units or the use of the Haystacks project is redundant.' )
    end
  end


% Size of the pixel of the scene in mas
% The pixel of the PSF basis is assumed to be the same (PSF files for this pixel size are required)
  if ~isfield( opt, 'pix_scene_mas' )
  opt.px_scn_mas = 3 ;
  else
  opt.px_scn_mas = opt.pix_scene_mas ;
  end
% The pixel size on the camera chip (by default same as the scene. That is no degradation)
  if ~isfield( opt, 'pix_camera_mas' )
  opt.px_cmr_mas = 0.4 * opt.lambda_imaging_1_nm * 1e-9 / opt.dmtr_tlscp_m * 180 / pi * 3600e3 ; % enough sampling
  else
  opt.px_cmr_mas = opt.pix_camera_mas ;
  end

% Wavelength array of the scenes
  if ~isfield( opt, 'lambda_array_scene' )
  opt.lmbd_arry_scn = 'file' ;
  else
  opt.lmbd_arry_scn = opt.lambda_array_scene ;
  end

% Pixel pitch of the PSF basis
  if ~isfield( opt, 'px_psf_mas' )
  l_D_tmp = opt.lambda_imaging_1_nm * 1e-9 / opt.dmtr_tlscp_m * 180 / pi * 3600e3 ;
    % For cases of large telescopes
    if ( l_D_tmp < 6 )
    opt.px_psf_mas = 1 ;
    else
    opt.px_psf_mas = 3 ;
      if ( opt.px_scn_mas < opt.px_psf_mas )
      opt.px_psf_mas = 1 ;
      end
    end
  end

% Distance between two PSF (see sister_basis.m)
  if ~isfield( opt, 'psf_spacing_mas' )
  opt.psf_spacing_mas = 1 ; 
  end

% Default imaging bands of Starshade: they have to be defined for each occulter
  if strfind( opt.starshade.nominal_filename, 'NI2' )
  opt.lambda_band_nm_min = [ 425, 606, 747, 615, 450, 590, 770  ] ; % The 615-800 is for the WFIRST-S mission probe study
  opt.lambda_band_nm_max = [ 552, 787, 970, 800, 600, 770, 1000 ] ;
  end
  if numel( strfind( opt.starshade.nominal_filename, 'nw2' ) ) || numel( strfind( opt.starshade.nominal_filename, 'tv3' ) )
  opt.lambda_band_nm_min = 300 ;
  opt.lambda_band_nm_max = 999 ;
  end
  if ~isfield( opt, 'lambda_band_nm_min' ) || ~isfield( opt, 'lambda_band_nm_max' )
  disp( 'Please, define the wavelength range for the PSF by setting opt.lambda_band_nm_min and opt.lambda_band_nm_max.' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

% Some specific bands that may be set directly
  if strfind( opt.starshade.nominal_filename, 'NI2' )
    if isfield( opt, 'band_425_552_nm' )
    opt.lambda_psf_1_nm = 425 ;
    opt.lambda_psf_2_nm = 552 ;
      if ~isfield( opt, 'lambda_imaging_1_nm' )
      opt.lambda_imaging_1_nm = opt.lambda_psf_1_nm ;
      end
      if ~isfield( opt, 'lambda_imaging_2_nm' )
      opt.lambda_imaging_2_nm = opt.lambda_psf_2_nm ;
      end
    end
    if isfield( opt, 'band_615_800_nm' )
    opt.lambda_psf_1_nm = 615 ;
    opt.lambda_psf_2_nm = 800 ;
      if ~isfield( opt, 'lambda_imaging_1_nm' )
      opt.lambda_imaging_1_nm = opt.lambda_psf_1_nm ;
      end
      if ~isfield( opt, 'lambda_imaging_2_nm' )
      opt.lambda_imaging_2_nm = opt.lambda_psf_2_nm ;
      end
    end
  end % Special bands

% Parameters that control the construction of the PSF from electric fields in the non-spinning Starshade case
% Portion of the PSF stored as times l/D
  if ~isfield( opt, 'n_lambda_over_d' )
  opt.n_lambda_over_d = 7 ;
  end

% Wavelength to consider for the image of Starshade (by default, same limits as the PSF basis, although the delta lambda may well be different. Also these wavelengths are suposed to be compatible with the wavelength simulated in the scenes. For example, Asking to image some wavelength absent in the scenes will be considered an odd case, likely from a wrong setup. The code will send a warning message when getting the scene data (convolve_with_one_wavelength.m and the sub-function get_scene_data)
  if ~isfield( opt, 'lambda_imaging_1_nm' )
    if strfind( opt.starshade.nominal_filename, 'NI2' )
    opt.lmbd_img_1_nm = [ 425, 606, 747 ] ;
    end
    if strfind( opt.starshade.nominal_filename, 'nw2' ) || strfind( opt.starshade.nominal_filename, 'tv3' )
    opt.lmbd_img_1_nm = 300 ;
    end
  else
  opt.lmbd_img_1_nm = opt.lambda_imaging_1_nm ;
  end
  if ~isfield( opt, 'lambda_imaging_2_nm' )
    if strfind( opt.starshade.nominal_filename, 'NI2' )
    opt.lmbd_img_2_nm = [ 552, 787, 970 ] ;
    end
    if strfind( opt.starshade.nominal_filename, 'nw2' ) || strfind( opt.starshade.nominal_filename, 'tv3' )
    opt.lmbd_img_2_nm = 1000 ;
    end
  else
  opt.lmbd_img_2_nm = opt.lambda_imaging_2_nm ;
  end

% Used to identify the band under consideration
  if ~isfield( opt, 'lambda_1_nm' )
  opt.lambda_1_nm = opt.lmbd_img_1_nm ;
  end
  if ~isfield( opt, 'lambda_2_nm' )
  opt.lambda_2_nm = opt.lmbd_img_2_nm ;
  end

% Spacing of the imaging band (in nanometers)
  if ~isfield( opt, 'delta_lambda_imaging_nm' )
  % By default, same as in the scenes provided
  opt.dlt_lmbd_img_nm = 'scene' ;
  else
  opt.dlt_lmbd_img_nm = opt.delta_lambda_imaging_nm ;
  end

% Arrays of wavelength to be considered 
% 1) Wavelength array in the scene
opt.n_bnd = numel( opt.lmbd_img_1_nm ) ;
  for i_bnd = 1 : opt.n_bnd
  [ scn_wl_lst opt ] = get_lambda_scene( opt, i_bnd ) ;
  % Guessing the units of the wavelength and bringing them all to nanometer
  opt.units_wavelength = 'nanometer' ;
  opt.units_wavelength_si = 1e-9 ; % m
  fct_wl = - 1 ;
    if ( mean( scn_wl_lst ) > 100 ) && ( mean( scn_wl_lst ) < 4000 )
    fct_wl = 1 ; % nm
    end
    if ( mean( scn_wl_lst ) > 1e-3 ) && ( mean( scn_wl_lst ) < 10 )
    fct_wl = 1e3 ; % micro meter
    end
    if ( mean( scn_wl_lst ) > 1e-7 ) && ( mean( scn_wl_lst ) < 1e-5 )
    fct_wl = 1e9 ; % m
    end
    if ( fct_wl == -1 )
    disp( 'Units of the wavelength array unidentified. Stopped' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end

  % When getting the array of wavelengths from an external scene, the range may well be larger than the imaging one.
    if ~( opt.scn.do )
    q_1 = find( scn_wl_lst * fct_wl >= opt.lmbd_img_1_nm ) ;
      if ~numel( q_1 )
      disp( sprintf( 'The wavelength range provided in the external scene (%04.2f nm, %04.2f nm) does not cover the lower imaging wavelength %04.2f nm'), scn_wl_lst( 1 ) * fct_wl, scn_wl_lst( end ) * fct_wl, opt.lmbd_img_1_nm )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      else
      q_2 = find( scn_wl_lst( q_1 ) * fct_wl <= opt.lmbd_img_2_nm ) ;
        if ~numel( q_2 )
        disp( sprintf( 'The wavelength range provided in the external scene (%04.2f nm, %04.2f nm) does not cover the higher imaging wavelength %04.2f nm'), scn_wl_lst( 1 ) * fct_wl, scn_wl_lst( end ) * fct_wl, opt.lmbd_img_2_nm )
        disp( ' ' )
        disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
        else
        % Final cut
        scn_wl_lst = scn_wl_lst( q_1( 1 ) : q_1( q_2( end ) ) ) ;
        end
      end
    end

  opt.lmbd_arry_scn_nm( i_bnd, 1 : numel( scn_wl_lst ) ) = scn_wl_lst * fct_wl ;

  % 2) Proposed wavelength to be imaged. It is trimed to be within the band of the instrument right after this step.
    if ischar( opt.dlt_lmbd_img_nm )
      if strcmp( opt.dlt_lmbd_img_nm, 'scene' )
      opt.lmbd_arry_img_nm = opt.lmbd_arry_scn_nm ;
      end
      if numel( strfind( lower( opt.dlt_lmbd_img_nm ), 'r' ) )
      opt.rslvng_pwr = str2num( opt.dlt_lmbd_img_nm( 2 : end ) ) ;
      opt.lmbd_arry_img_nm = wavelength_grid( opt ) ;
      % Last wavelength is the end of the band
      opt.lmbd_arry_img_nm = opt.lmbd_arry_img_nm( 1 : end -1 ) ;
      end
    else
    lmbd_arry_img_nm_tmp = opt.lmbd_img_1_nm( i_bnd ) : opt.dlt_lmbd_img_nm : opt.lmbd_img_2_nm( i_bnd ) ;
    opt.lmbd_arry_img_nm( i_bnd, 1 : numel( lmbd_arry_img_nm_tmp ) ) = lmbd_arry_img_nm_tmp ;
    end
  end % i_bnd

% Wavelength range to consider for the PSF of Starshade taking into account the imaging options (imaging may be a subset of the PSF basis)
% Checking if the imaging wavelength is within the PSF wavelength range
n_bnd_tmp = numel( opt.lambda_band_nm_min ) ;
    for i_bnd = 1 : opt.n_bnd
      if numel( opt.lambda_band_nm_max ) ~= n_bnd_tmp
      disp( sprintf( 'The number of bands is inconsistent in sister_imaging.m. Fix it. Stopped.' ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end

    % Defining the imaging band as that one where the wavelength fits within one of the official bands (it should be the same, but one may try some modifications within them)
      out_of_band = 1 ;
      lbl_bnd = '' ;
      for i_tmp = 1 : n_bnd_tmp
        if ( opt.lambda_band_nm_min( i_tmp ) <= opt.lmbd_arry_img_nm( i_bnd, 1 ) ) && ( opt.lmbd_arry_img_nm( i_bnd, end ) <= opt.lambda_band_nm_max( i_tmp ) )
        out_of_band = 0 ;
        idx_bnd = i_tmp ;
        end
      lbl_bnd = sprintf( '%s  [%04.2f, %04.2f]', lbl_bnd, opt.lambda_band_nm_min( i_tmp ), opt.lambda_band_nm_max( i_tmp ) ) ;
      end
    % If unidentified, return
      if ( out_of_band )
      disp( sprintf( 'The imaging band [%04.2f, %04.2f] nm does not correspond to any of the expected instrumental bands %s nm. Check your input band. If you need to modify the instrumental band, update sister_imaging to your desired band before continuing. Stopping.', opt.lmbd_arry_img_nm( i_bnd, 1 ), opt.lmbd_arry_img_nm( i_bnd, end ), lbl_bnd ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    clear out_of_band
    opt.lmbd_psf_1_nm( i_bnd ) = opt.lambda_band_nm_min( idx_bnd ) ;
    opt.lmbd_psf_2_nm( i_bnd ) = opt.lambda_band_nm_max( idx_bnd ) ;
    end % i_bnd

% Wavelength spacing of the PSF basis (comes from get_default_options.m)
  if ~isfield( opt, 'delta_lambda_nm' )
  opt.dlt_lmbd_psf_nm = 10 ;
  else
  opt.dlt_lmbd_psf_nm = opt.delta_lambda_nm  ;
  end
% Alternative equivalent notation (consistent with get_default_options.m as well)
  if isfield( opt, 'delta_lambda_psf_nm' )
  opt.dlt_lmbd_psf_nm = opt.delta_lambda_psf_nm ;
  end

% Precision for the unit conversion
  if ~ischar( opt.dlt_lmbd_img_nm )
    if opt.dlt_lmbd_img_nm > ( opt.lmbd_img_2_nm - opt.lmbd_img_1_nm ) / 10
    disp( 'Unit conversion accuracy to photon/s requires running the scenes with steps of wavelength equal or less than 1/10 of the bandwidth' )
    end
  end

% Distance Starshade telescope
  if ~isfield( opt, 'distance_starshade_telescope_m' )
  % Some default cases that depend on some well known configurations
  [ dst_m geo_iwa_mas ] = set_starshade_distance_and_geometric_iwa( opt, 1 ) ;
  opt.dst_strshd_tlscp_m = dst_m ;
  else
  opt.dst_strshd_tlscp_m = opt.distance_starshade_telescope_m ;
  end

% Geometric IWA
  if ~isfield( opt, 'geo_iwa_mas' )
  [ dst_m geo_iwa_mas ] = set_starshade_distance_and_geometric_iwa( opt, 1 ) ;
  opt.geo_iwa_mas = geo_iwa_mas ;
  end

% Radius of the non-stationary region. That is, distance in mas from the center of the Starshade for which the PSF bacomes stationary. This parameter is et for WFIRST and HabEx in set_r_stationary_mas.m. However, for new occulters, the user must set it up in the running configuration file.
  if isfield( opt, 'r_stationary_mas' )
  opt.r_st_mas = opt.r_stationary_mas ;
  disp( sprintf( 'PSF becomes stationary at distances of %2.2f mas or greater', opt.r_st_mas ) )
  end

% Plotting the results (most properties set in get_plot_options in sister_imaging_band)
  if ~isfield( opt, 'plot' )
  opt.plt.do = 0 ;
  else
  opt.plt.do = opt.plot.do ;
  end

% Output format for the images
  if ~isfield( opt, 'format_fig' )
  opt.frmt_fg = 'png' ;
  else
  opt.frmt_fg = opt.format_fig ;
  end

% Storing the plots
  if ~isfield( opt, 'save_plot' )
  opt.sv_plt = 1 ;
  else
  opt.sv_plt = opt.save_plot ;
  end

% Video option
  if ~isfield( opt, 'video' )
  opt.vd.do = 0 ;
  else
    if ~isfield( opt.video, 'do' )
    opt.vd.do = 0 ;
    else
    opt.vd.do = opt.video.do ;
    end
    if ~isfield( opt.video, 'delay' )
    opt.vd.dly = 40 ;
    else
    opt.vd.dly = opt.video.delay ;
    end
    if ( opt.vd.do )
    opt.sv_plt = 1 ;
    end
  end

% Directory where the figures are stored
  if ~isfield( opt, 'fig_dir' )
  opt.fg_dr = [ opt.out_dr 'fig/' ] ;
  else
  opt.fg_dr = opt.fig_dir ;
  end

% Directory where the figures are stored
  if ~isfield( opt, 'video_dir' )
  opt.vd_dr = [ opt.out_dr 'video/' ] ;
  else
  opt.vd_dr = opt.video_dir ;
  end

% Delay in seconds between two images
  if ~isfield( opt, 'video_delay_s' )
  opt.vd_dly_s = 0.5 ; % default 0.5 seconds
  else
  opt.vd_dly_s = opt.video_delay_s ;
  end

% Creating the directory if it does not exist
  if ~isdir( opt.fg_dr )
  mkdir( opt.fg_dr ) ;
  end

% Star
  if ~isfield( opt, 'star' ) % || ~( opt.scn.do )
  opt.star.do = 0 ;
  opt.str.tp = '' ; 
  opt.str.flx = NaN ;
  end
  if isfield( opt, 'star' )
  % Type of star: 'sun', or 'a0v', 'a5v', 'f5v', 'g0v', 'g5v' (closer to sun), 'k0v', 'k5v', 'm0v', 'm5v'
    if ~isfield( opt.star, 'type' )
    opt.str.tp = 'sun' ;
    else
    opt.str.tp = lower( opt.star.type ) ;
    end
  % Whether using the Willmer (2018), AM0 (2000) or WMO (1985) solar irradiance (C.N.A. Willmer, Astrophysical Journal Supplements, 2018, 236, 47, or http://mips.as.arizona.edu/~cnaw/sun.html, https://www.nrel.gov/grid/solar-resource/spectra-astm-e490.html, https://www.nrel.gov/grid/solar-resource/spectra-wehrli.html)
    if ~isfield( opt.star, 'am0' )
    opt.str.am0 = 0 ; % by default, use the solar's irradiance from Willmer 2018.
    else
    opt.str.am0 = opt.star.am0 ;
    end
    if ~isfield( opt.star, 'wmo' )
    opt.str.wmo = 0 ; % by default, use the solar's irradiance from Willmer 2018.
    else
    opt.str.wmo = opt.star.wmo ;
    end
  % Check
    if ( opt.str.am0 ) && ( opt.str.wmo )
    disp( 'Select only one solar spectrum. You have selected both AM0 and WMO options. Choose on of both, and make sure you do not want to use the default case: Willmer 2018.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end

  % Star flux (rare: should usually be in the scene file or computed with the star type and its magnitude. However, it may also be set here if a specific value wants to be used)
    if ~isfield( opt.star, 'flux' ) && ~isfield( opt.str, 'flx' )
    opt.str.flx = NaN ; % not used
    else
      if isfield( opt.star, 'flux' )
      opt.str.flx = opt.star.flux ;
      end
    end  

 % Distance to Earth
    if ~isfield( opt.star, 'distance_to_earth_pc' )
    opt.str.dst_pc = 10 ;
    else
    opt.str.dst_pc = opt.star.distance_to_earth_pc ;
    end

  % Magnitude. Recall for the Sun: apparent (-26.76 V), 5 pc (3.33 V) and absolute 10 pc (4.81 V). http://mips.as.arizona.edu/~cnaw/sun.html)
    if ~isfield( opt.star, 'app_mag_v' )
    opt.str.mg = 4.81 - 2.5 * log( 10 / opt.str.dst_pc ) ; % V4.81 Sun at 10 pc
    else
    opt.str.mg = opt.star.app_mag_v ;
    end
  % Name
    if ~isfield( opt.star, 'name' )
    opt.str.nm = 'Sun-like' ;
    else
    opt.str.nm = opt.star.name ;
    end

  % Mass in solar masses (for orbits)
    if ~isfield( opt.star, 'times_sun_mass' )
    opt.str.mss_sn = 1 ; % 1 solar mass
    else
    opt.str.mss_sn = opt.star.times_sun_mass ;
    end

  % RA of the star (J2000)
    if ~isfield( opt.star, 'ra_deg' )
    opt.str.ra_dg = 0 ;
    else
    opt.str.ra_dg = opt.star.ra_deg ;
    end

  % DEC of the star (J2000)
    if ~isfield( opt.star, 'dec_deg' )
    opt.str.dc_dg = 0 ;
    else
    opt.str.dc_dg = opt.star.dec_deg ;
    end
 
  % Proper motion (Background motion at different epochs)
    if ~isfield( opt.star, 'pm_ra_mas_yr' )
    opt.str.pm_ra_mas_yr = 0 ;
    else
    opt.str.pm_ra_mas_yr = opt.star.pm_ra_mas_yr ;
    end
    if ~isfield( opt.star, 'pm_dec_mas_yr' )
    opt.str.pm_dc_mas_yr = 0 ;
    else
    opt.str.pm_dc_mas_yr = opt.star.pm_dec_mas_yr ;
    end

  % Some stored stars (it will overwrite any parameter set before if it is found in the ExoCat database)
  opt = get_star_properties( opt ) ;

  % Conversion between AU and mas
  opt.au2mas = 1000 / opt.str.dst_pc ;
  opt.mas2au = 1 / opt.au2mas ;

  % Some checks
    if ~strcmp( opt.str.tp, '' )
    % Range of wavelengths in the HST_stars_0mag_flux_... file from external sources (see prepare_scene.m).
    lmbd_str_1_nm = 250 ;
    lmbd_str_2_nm = 1050 ;
      if strcmp( lower( opt.str.tp ), 'sun' )
      lmbd_str_1_nm = 100 ; % This is the minimum of the solar spectrum in the WMO file (which is not the default, but it's low enough. The default is Willmer 2018)
      lmbd_str_2_nm = 10000 ; % Enough for any starshade simulations
      end
      
      if ( opt.lmbd_img_1_nm < lmbd_str_1_nm )
      disp( sprintf( 'The minimum wavelength to be simulated (%04.2f nm) is less than the available one (%04.2f nm) in the files with the star spectrum. Change it to be greater or equal than the one in the file. PS: if you can choose opt.star.type=''sun'', then the available range is 119-10,000 nm.', opt.lmbd_img_1_nm, lmbd_str_1_nm ) ) ;
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
      if ( opt.lmbd_img_2_nm > lmbd_str_2_nm )
      disp( sprintf( 'The maximum wavelength to be simulated (%04.2f nm) is greater than the available one (%04.2f nm) in the files with the star spectrum. Change it to be less or equal than the one in the file. PS: if you can choose opt.star.type=''sun'', then the available range is 119-10,000 nm.', opt.lmbd_img_2_nm, lmbd_str_2_nm ) ) ;
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end
  end % opt.star

% Local zodiacal light for the scenes

  if ~isfield( opt, 'local_zodi' )
  opt.lcl_zd.do = 0 ; % Default, not added.
  opt.local_zodi.do = 0 ;
  else
    if isfield( opt.local_zodi, 'do' )
    opt.lcl_zd.do = opt.local_zodi.do ;
    else
    opt.lcl_zd.do = 0 ; % * If there are other opt.local_zodi. options set, but opt.local_zodi.do is not, don't include it.
    end
  end
  
  if ( opt.lcl_zd.do )
  % Relative magnitude of the local zodi
  % 1) Set directly
  lcl_zd_ext = 0 ;
    if ~isfield( opt.local_zodi, 'mag_v_arcsec2' )
    opt.lcl_zd.mg_v_arcs2 = 23 ; % Standard value used in Exposure Time Calculators if there's no information about the Helio-Ecliptic coordinates. It corresponds to helio-ecliptic longitudes greater than 90 deg, or some high ecliptic latitudes (>60 deg). See get_local_zodi_mag_v_arcsec2.m below.
    opt.local_zodi.mag_v_arcsec2 = opt.lcl_zd.mg_v_arcs2 ; % For the record
    else
    opt.lcl_zd.mg_v_arcs2 = opt.local_zodi.mag_v_arcsec2 ;
    % Check of consistency (it should be within 20-23.5)
      if ( opt.lcl_zd.mg_v_arcs2 < 20 ) || ( opt.lcl_zd.mg_v_arcs2 > 23.5 )
      disp( sprintf( 'WARNING: the surface brightness of the local zodiacal light (m_V/arcsec^2=%3.2f) is outside the expected range of 20-23.5. Make sure this is intentional. Continuing.', opt.lcl_zd.mg_v_arcs2 ) )
      end
   lcl_zd_ext = 1 ;
   end

  % If the heliocentric coordinates are provided, the magnitude is derived from STScI compiled data for HST. NOTE: helio-ecliptic longitude is different from Earth's ecliptic longitude of the target star. See appendix A in STScI, TIR-CRDS-2015-01.pdf (R. Diaz, 04/01/2015, http://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/reference-data-for-calibration-and-tools/documentation/_documents/TIR-CRDS-2015-01.pdf)
    if ~isfield( opt.local_zodi, 'eclipitic_latitude_deg' )
    opt.lcl_zd.ecl_lt_dg = NaN ; % Making sure it is not used 
    else
    opt.lcl_zd.ecl_lt_dg = opt.local_zodi.eclipitic_latitude_deg ;
    end

    if ~isfield( opt.local_zodi, 'helio_eclipitic_longitude_deg' )
    opt.lcl_zd.hl_ecl_ln_dg = NaN ; % Making sure it is not used
    else
    opt.lcl_zd.hl_ecl_ln_dg = opt.local_zodi.helio_eclipitic_longitude_deg ;
    end

  % Check of consistency
    if ~isnan( opt.lcl_zd.ecl_lt_dg ) && isnan( opt.lcl_zd.hl_ecl_ln_dg ) 
    disp( sprintf( 'WARNING: the ecliptic latitude was set to %3.2f deg, but the helio-ecliptic longitude was not set. One must provide a value. Stopped.', opt.lcl_zd.ecl_lt_dg ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end

    if isnan( opt.lcl_zd.hl_ecl_ln_dg ) && ~isnan( opt.lcl_zd.ecl_lt_dg )
    disp( sprintf( 'WARNING: the helio-ecliptic longitude was set to %3.2f deg, but the helio-ecliptic latitude was not set. One must provide a value. Stopped.', opt.lcl_zd.ecl_lt_dg ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end

  % If the both ecliptic coordinates are set
    if ~isnan( opt.lcl_zd.hl_ecl_ln_dg ) &&  ~isnan( opt.lcl_zd.ecl_lt_dg )
      if lcl_zd_ext
      disp( sprintf( 'WARNING: the zodiacal light was set to %2.2f V mag/arcsec^2, but helio-ecliptic coordinates have been provided.', opt.lcl_zd.mg_v_arcs2 ) )
      end
    opt.lcl_zd.mg_v_arcs2 = get_local_zodi_mag_v_arcsec2( opt.lcl_zd.hl_ecl_ln_dg, opt.lcl_zd.ecl_lt_dg ) ;
    opt.local_zodi.mag_v_arcsec2 = opt.lcl_zd.mg_v_arcs2 ; % For the record

      if lcl_zd_ext
      disp( sprintf( ' The actual value will be substituted by %2.2f, that corresponds to the ecliptic coordinates choice. Verify this is what is intended. Continuing ...', opt.local_zodi.mag_v_arcsec2 ) )
      end

      if isnan( opt.lcl_zd.mg_v_arcs2 )
      disp( sprintf( 'WARNING: the corresponding value of the local zodiacal surface brightness for the ecliptic coordinates (%3.1f deg,%3.1f deg) is NaN, due to the proximity of the line of sight to the Sun. Check this was intentional. Stopped.', opt.lcl_zd.hl_ecl_ln_dg, opt.lcl_zd.ecl_lt_dg ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end
  disp( sprintf( 'Local zodiacal light. Surface brigthness is %2.2f V mag/arcsec^2', opt.local_zodi.mag_v_arcsec2 ) )
  end % End of local zodiacal light

% Number of loops (modified if kepler orbits are set)
  if ~isfield( opt, 'n_lp' )
  opt.n_lp = 1 ;
  end

% Planets
opt.au2km = 149597870.7 ; % defined once
  if ~isfield( opt, 'planets' )
  % Removing planets
  opt.plnt.rmv.do = 0 ;
  % First/Second positions of the planets with respect the star in an 2-dimensional array.
  % PS: matlab uses row/column if row means horizontal, 'y level', and column means vertical, 'x level'. 
  opt.plnt.rmv.arc_dc_mas = NaN ;
  opt.plnt.rmv.arc_ra_mas = NaN ;
  % Adding planets
  opt.plnt.add.do = 0 ;
  opt.plnt.add.arc_dc_mas = NaN ;
  opt.plnt.add.arc_ra_mas = NaN ;
  % When added, a flux ratio needs to be provided
  opt.plnt.add.flx_rt = NaN ;
  opt.plnt.add.flx_rt_arry =  NaN ;
  % The flux of the star should be given in some file (see opt.str.flx_fl or the scene data cubes should have that information. See the add_planets sub-function.
  else

    if ~isfield( opt.planets, 'remove' )
    opt.planets.remove.dummy = 0 ;
    end

    if ~isfield( opt.planets.remove, 'do' )
    opt.plnt.rmv.do = 0 ;
    else
    opt.plnt.rmv.do = opt.planets.remove.do ;
    end

    if ~isfield( opt.planets.remove, 'pos_arc_dec_mas' )
    opt.plnt.rmv.arc_dc_mas = NaN ;
    else
    opt.plnt.rmv.arc_dc_mas = opt.planets.remove.pos_arc_dec_mas ;
    end
   
    if ~isfield( opt.planets.remove, 'pos_arc_ra_mas' )
    opt.plnt.rmv.arc_ra_mas = NaN ;
    else
    opt.plnt.rmv.arc_ra_mas = opt.planets.remove.pos_arc_ra_mas ;
    end

  % Minimum checks of consistency
    if numel( opt.plnt.rmv.arc_dc_mas ) ~= numel( opt.plnt.rmv.arc_ra_mas )
    disp( 'Number of ''ra/dec'' planet positions to remove is different. Non-sense. Stopped.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end

  % Conversion into pixels
  opt.plnt.rmv.arc_dc_px = round( opt.plnt.rmv.arc_dc_mas / opt.px_scn_mas ) ;
  opt.plnt.rmv.arc_ra_px = round( opt.plnt.rmv.arc_ra_mas / opt.px_scn_mas ) ;

    if ~isfield( opt.planets, 'add' )
    opt.planets.add.dummy = 0 ;
    end

    if ~isfield( opt.planets.add, 'do' )
    opt.plnt.add.do = 0 ;
    else
    opt.plnt.add.do = opt.planets.add.do ;
    end

  % The positions of the planets are with respect the star center
  % Pos_1 will be vertical (row)
    if ~isfield( opt.planets.add, 'pos_arc_dec_mas' )
    opt.plnt.add.arc_dc_mas = NaN ;
    else
    opt.plnt.add.arc_dc_mas = opt.planets.add.pos_arc_dec_mas ;
    end

  % Pos_2 will be horizontal (column)
    if ~isfield( opt.planets.add, 'pos_arc_ra_mas' )
    opt.plnt.add.arc_ra_mas = NaN ;
    else
    opt.plnt.add.arc_ra_mas = opt.planets.add.pos_arc_ra_mas ;
    end

  % The positions can also be given in AU
    if isfield( opt.planets.add, 'pos_arc_dec_au' ) || isfield( opt.planets.add, 'pos_arc_ra_au' )
      if isfield( opt.planets.add, 'pos_arc_dec_au' )
      opt.plnt.add.arc_dc_mas = opt.planets.add.pos_arc_dec_au * opt.au2mas ;
      else
      disp( 'The dec of the planet(s) is not provided while the ra is. Fix it. Stopped.' )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    % The RA should also be given
      if isfield( opt.planets.add, 'pos_arc_ra_au' )
      opt.plnt.add.arc_ra_mas = opt.planets.add.pos_arc_ra_au * opt.au2mas ;
      else
      disp( 'The ra of the planet(s) is not provided while the dec is. Fix it. Stopped.' )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end

  % Conversion into pixels
  opt.plnt.add.arc_dc_px = round( opt.plnt.add.arc_dc_mas / opt.px_scn_mas ) ;
  opt.plnt.add.arc_ra_px = round( opt.plnt.add.arc_ra_mas / opt.px_scn_mas ) ;
  
  % Number of planets
    if ~isfield( opt, 'n_pl' )
      if ~isnan( opt.plnt.add.arc_ra_px )
      opt.n_pl = numel( opt.plnt.add.arc_ra_px ) ;
      end

      if isfield( opt.planets.add, 'r_orb_au' )
      opt.n_pl = numel( opt.planets.add.r_orb_au ) ;
      end

    % If it is not set
      if ~isfield( opt, 'n_pl' )
      disp( 'WARNING: the number of planets could not be determined. Stopped.' )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end

    if ~isfield( opt.planets.add, 'flux_ratio' )
    opt.plnt.add.flx_rt = NaN ;
    % In case that, instead of the flux ratio, a phase angle and planet type is provided
      if isfield( opt.planets.add, 'phase_angle_deg' )
      % If there is only 1 value, replicate.
        if numel( opt.planets.add.phase_angle_deg ) == 1
        phs_ang_dg_tmp = opt.planets.add.phase_angle_deg ;
          for i_pl = 1 : opt.n_pl
          opt.plnt.add.phs_ang_dg( i_pl ) = phs_ang_dg_tmp ;
          end
        else
        opt.plnt.add.phs_ang_dg = opt.planets.add.phase_angle_deg ;
        end
      % If the number of entries is inconsistent, stop
        if numel( opt.plnt.add.phs_ang_dg ) ~= opt.n_pl
        disp( sprintf( 'WARNING: the number of values set for opt.planets.add.phase_angle_deg is %i, and should be the same as the number of planets %i, or just 1 value, the same for all. Stopped', numel( opt.plnt.add.phs_ang_dg ), opt.n_pl ) )
        disp( ' ' )
        disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
        end
      % Check no value corresponds to face on
        for i_pl = 1 : opt.n_pl
        ok = 1 ;
          if ~( opt.plnt.add.phs_ang_dg( i_pl ) )
          disp( sprintf( 'WARNING: the phase angle for planet %i is zero, which corresponds to face on, and that would set the planet on top of the star. Choose a different value for the phase angle.', i_pl ) )
          ok = 0 ;
          end
          if ~ok 
          disp( ' ' )
          disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
          end
        end
      opt.plnt.add.phs_ang_rd = pi / 180 * opt.plnt.add.phs_ang_dg ;
      end
    else
    opt.plnt.add.flx_rt = opt.planets.add.flux_ratio ;
    % Checking the sizes are correct
      if ( size( opt.plnt.add.flx_rt, 2 ) ~= numel( opt.lmbd_arry_img_nm ) )  || ( size( opt.plnt.add.flx_rt, 1 ) ~= numel( opt.plnt.add.arc_dc_mas ) )
      % Transpose case
        if ( size( opt.plnt.add.flx_rt, 1 ) == numel( opt.lmbd_arry_img_nm ) )  && ( size( opt.plnt.add.flx_rt, 2 ) == numel( opt.plnt.add.arc_dc_mas ) )
        opt.plnt.add.flx_rt = opt.plnt.add.flx_rt' ;
        else
         % Replicating the dimension with the planet flux ratios
         disp( 'Replicating the flux ratio values for all wavelengths' )
           if ( size( opt.plnt.add.flx_rt, 1 ) == numel( opt.plnt.add.arc_dc_mas ) )
              if  ( size( opt.plnt.add.flx_rt, 2 ) ~= 1 )
              disp( sprintf( 'The option opt.planets.add.flux_ratio has dimensions %ix%i, which is inconsistent with the expected number of elements of the  imaging wavelength array (%i) and the number of planets (%i). Make sure opt.planets.add.flux_ratio is consistent.',  size( opt.plnt.add.flx_rt, 2 ),  size( opt.plnt.add.flx_rt, 1 ), numel( opt.lmbd_arry_img_nm ), numel( opt.plnt.add.arc_dc_mas ) ) )
              disp( ' ' )
              disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
              else
              disp( 'Replicating the values in opt.planet.flux_ratio for each imaging wavelength. That is, assuming constant albedo across wavelength.' )
              tmp = opt.plnt.add.flx_rt ;
              opt.plnt.add = rmfield( opt.plnt.add, 'flx_rt' ) ;
                for i_wl = 1 : numel( opt.lmbd_arry_img_nm )
                opt.plnt.add.flx_rt( :, i_wl ) = tmp ;
                end
               clear tmp
              end
           end
         % Similar case, but if it would be 1xn_pl for some reason (transposed)
           if ( size( opt.plnt.add.flx_rt, 2 ) == numel( opt.plnt.add.arc_dc_mas ) )
              if  ( size( opt.plnt.add.flx_rt, 1 ) ~= 1 )
              disp( sprintf( 'The option opt.planets.add.flux_ratio has dimensions %ix%i, which is inconsistent with the expected number of elements of the  imaging wavelength array (%i) and the number of planets (%i). Make sure opt.planets.add.flux_ratio is consistent.',  size( opt.plnt.add.flx_rt, 1 ),  size( opt.plnt.add.flx_rt, 2 ), numel( opt.lmbd_arry_img_nm ), numel( opt.plnt.add.arc_dc_mas ) ) )
              disp( ' ' )
              disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
              else
              disp( 'Replicating the values in opt.planet.flux_ratio for each imaging wavelength. That is, assuming constant albedo across wavelength.' )
              tmp = opt.plnt.add.flx_rt ;
              opt.plnt.add = rmfield( opt.plnt.add, 'flx_rt' ) ;
                for i_wl = 1 : numel( opt.lmbd_arry_img_nm )
                opt.plnt.add.flx_rt( :, i_wl ) = tmp ;
                end
               clear tmp
              end
           end
        end
      end
    % To be used in sister_imaging_band
    opt.planets.add.flux_ratio_array = opt.plnt.add.flx_rt ;
    end

  % For the case when new planets are added externally, without Keplerian orbits
  if isfield( opt.planets.add, 'radius' )
    for i_pl = 1 : opt.n_pl
    fct_km = NaN ;
    str_tmp = opt.planets.add.radius{ i_pl } ;
      if ~isnan( str_tmp )
      str_tmp = lower( str_tmp ) ;
        if strfind( str_tmp, 'km' )
        q_str = findstr( str_tmp, 'km' ) ;
        opt.planets.add.radius_km( i_pl ) = str2num( str_tmp( 1 : q_str - 1 ) ) ;
        end
        if numel( strfind( str_tmp, 'r_e' ) ) || numel( strfind( str_tmp, 're' ) )
        fct_km = 6371 ; % mean radius
        q_str = findstr( str_tmp, 'r' ) ;
        opt.planets.add.radius_km( i_pl ) = fct_km * str2num( str_tmp( 1 : q_str - 1 ) ) ;
        end
        if numel( strfind( str_tmp, 'r_j' ) ) || numel( strfind( str_tmp, 'rj' ) )
        fct_km = 71492 ;
        q_str = findstr( str_tmp, 'r' ) ;
        opt.planets.add.radius_km( i_pl ) = fct_km * str2num( str_tmp( 1 : q_str - 1 ) ) ;
        end
      % Check of consistency
        if isnan( fct_km )
        disp( sprintf( 'Unidentified units for planet %i in opt.planets.add.radius. Stopped.', i_pl ) )
        disp( ' ' )
        disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
        end
      end % ~isnan( str_tmp )
    end % i_pl
  end % If radius is introduced

  if isfield( opt.planets.add, 'radius_km' )
    if numel( opt.planets.add.radius_km ) == 1
    opt.plnt.add.rd_km = opt.planets.add.radius_km * ones( 1, opt.n_pl ) ;
    else
      if numel( opt.planets.add.radius_km ) < opt.n_pl
      disp( sprintf( 'The number of planets is %i, but opt.planets.add.radius_km has only %i elements. Please, fix it. Stopped.', opt.n_pl, numel( opt.planets.add.radius_km ) ) ) ;
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
      if numel( opt.planets.add.radius_km ) > opt.n_pl
      disp( sprintf( 'WARNING: the number of planets is %i, but opt.kepler.planets.add_km has more elements (%i). Only the first %i values will be used. Please, make sure this is intentional. Continuing.', opt.n_pl, numel( opt.planets_add.radius_km ), opt.n_pl ) ) ;
      end
      for i_pl = 1 : opt.n_pl
        if ~isnan( opt.planets.add.radius_km( i_pl ) )
        opt.plnt.add.rd_km( i_pl ) = opt.planets.add.radius_km( i_pl ) ;
        end
      end
    end
  opt.kplr.rd_km = opt.plnt.add.rd_km ;
  end % opt.planets.add.radius_km

  % For Kepler orbits
    if ~isfield( opt.planets.add, 'flux_ratio_array' )
    opt.plnt.add.flx_rt_arry = NaN ;
    else
    opt.plnt.add.flx_rt_arry = opt.planets.add.flux_ratio_array ;
    end

    if ~isfield( opt.planets.add, 'phase_function' )
      if ( opt.n_pl )
        for i_pl = 1 : opt.n_pl
        opt.plnt.add.phs_fnctn{ i_pl } = 'Lambertian' ;
        end
      end
    else
      if ~iscell( opt.planets.add.phase_function )
      disp( 'The option opt.planets.add.phase_function is a cell: opt.planets.add.phase_function = { option, ... }. Fix it. Stopped.' )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    n_phs_f = numel( opt.planets.add.phase_function ) ;
      if n_phs_f == 1
        for i_pl = 1 : opt.n_pl
        opt.plnt.add.phs_fnctn{ i_pl } = opt.planets.add.phase_function ;
        end
      else
        if n_phs_f ~= opt.n_pl
        disp( sprintf( 'The number of elements in the phase function (%i) is greater than 1, but inconsistent with the number of planets (%i). Fix it. Stopped.', n_phs_f, opt.n_pl ) ) ;
        else
        opt.plnt.add.phs_fnctn = opt.planets.add.phase_function ;
        end
      end
    end
  % Assigning the flux ratio when the phase angle was set
    if isfield( opt.planets.add, 'phase_angle_deg' )
      if ~isfield( opt.planets.add, 'planet_type' )
      disp( 'WARNING: if the option opt.planets.add.phase_angle_deg was set, it is necessary to set opt.planets.add.planet_type. See SISTER Handbook for available options. Stopped.' )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop 
      else
      opt.plnt.add.plnt_tp = opt.planets.add.planet_type ;
      % Replicating the value if there's only 1
        if numel( opt.plnt.add.plnt_tp ) == 1
        plnt_tp_tmp = opt.plnt.add.plnt_tp ;
          for i_pl = 1 : opt.n_pl
          opt.plnt.add.plnt_tp( i_pl ) = plnt_tp_tmp ;
          end
        end
      % For the case when the actual distance between the planet and the host star is set by the user
      if isfield( opt.planets.add, 'r_orb_au' )
      opt.plnt.add.r_orb_au = opt.planets.add.r_orb_au ;
       % Replicating the value if there's only 1
        if numel( opt.plnt.add.r_orb_au ) == 1
        r_orb_au_tmp = opt.plnt.add.r_orb_au ;
          for i_pl = 1 : opt.n_pl
          opt.plnt.add.r_orb_au( i_pl ) = r_orb_au_tmp ;
          end
        end
      end % planets.add.r_orb_au

      % Using Kepler temporarily
        if ~isfield( opt.plnt.add, 'rd_km' )
        opt.kepler.planet_type = opt.plnt.add.plnt_tp ;
        opt.kplr.dummy = 1 ; % Only necessary because we are identifyig planets with the ones coded up in get_data_from_planet_type, but not suing full Keplerian orbits.
        opt = get_data_from_planet_type( opt ) ;
        % Checking the sentinel
          for i_pl = 1 : opt.n_pl
            if isnan( opt.kplr.rd_km( i_pl ) ) 
            disp( sprintf( 'Planet %s not identified. Fix it. Stopped', opt.kepler.planet_type{ i_pl } ) )
            disp( ' ' )
            disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
            else
            % Field that needs be parsed to planets.add
            opt.plnt.add.rd_km( i_pl ) = opt.kplr.rd_km( i_pl ) ;
            end
          end
        end

      % For planets externally introduced (not in the list in get_data_from_planet_type.m)
        if ~isfield( opt.plnt.add, 'gm_albd' )
          for i_pl = 1 : opt.n_pl
          opt.plnt.add.gm_albd{ i_pl } = lower( opt.plnt.add.plnt_tp{ i_pl } ) ;
          % To be used in get_planet_flux_ratio.m
          opt.kplr.gm_albd{ i_pl } = opt.plnt.add.gm_albd{ i_pl } ;
          end
        end

      % Assigning the flux array for each planet
      opt = get_planet_flux_ratio( opt ) ;
      % To be used in sister_imaging_band
      opt.planets.add.flux_ratio_array = opt.plnt.add.flx_rt_arry ;
      % Turning off Kepler option
      opt.kepler.do = 0 ;
      opt.kplr.do = 0 ;
      end % planet_type
    end % phase_angle_deg
  end % Planets field

% Kepler orbits (most options are set in get_kepler_options, get_planet_data and get_kepler_imaging_options, in sister_imaging_band)
  if ~isfield( opt, 'kepler' )
  opt.kplr.do = 0 ;
  else
  opt.kplr.do = opt.kepler.do ;
  % Add the planets if Kepler option is set on
    if ( opt.kplr.do )
    opt.plnt.add.do = 1 ;
    % And removing the planets if Haystacks input scene is used (because it has an instantation of the Solar System planets)
      if strcmp( opt.scn.nm( 1 :11 ), 'modern_cube' ) 
      opt.plnt.rmv.do = 1 ;
      end
    end
    if ~isfield( opt.kepler, 'system' )
    opt.kplr.sys = 'custom' ; % None in particular
    else
    opt.kplr.sys = opt.kepler.system ;
    end
  % Option to retrieve the Keplerian orbit for statistical analysis, without actually generating the scene simulation
    if ~isfield( opt.kepler, 'get_orbit' )
    opt.kplr.gt_orbt = 0 ; % By default, don't do it
    else
    opt.kplr.gt_orbt = opt.kepler.get_orbit ;
    end
    % Returning after getting the orbit
      if ( opt.kplr.gt_orbt )
      opt.return = 1 ;
      end
  end

% Add extragalactic background field
  if ~isfield( opt, 'extragalactic_background' )
  opt.bckgrnd.add = 0 ;
  else
    if ~isfield( opt.extragalactic_background, 'add' )
    opt.bckgrnd.add = 0 ;
    else
    opt.bckgrnd.add = opt.extragalactic_background.add ;
    end
    % Default extragalactic_background is from Haystacks project (some peculiarities as explained in the subfunction add_extragalactic_background in convolve_with_one_wavelength.m)
    if ~isfield( opt.extragalactic_background, 'image' )
    opt.bckgrnd.img = 'GALAXIES_10lat' ;
    else
    opt.bckgrnd.img = opt.extragalactic_background.image ;
    end
    % Center of the Haystacks extragalactic_background to be selected for the scenes (Haystacks extragalactic_background is 3600x3600 pix, each 10 mas, much bigger than the default 3333x3333 pix of 3 mas each of the Haystacks solar system scenes)
    if ~isfield( opt.extragalactic_background, 'pos_arc_dec_mas' )
    opt.bckgrnd.arc_dc_mas = 0 ; % Same center as external extragalactic background data
    else
    opt.bckgrnd.arc_dc_mas = opt.extragalactic_background.pos_arc_dec_mas ;
    end
    if ~isfield( opt.extragalactic_background, 'pos_arc_ra_mas' )
    opt.bckgrnd.arc_ra_mas = 0 ; % Same center as external extragalactic background data
    else
    opt.bckgrnd.arc_ra_mas = opt.extragalactic_background.pos_arc_ra_mas ;
    end
  end % Background field

% Add background stars
  if ~isfield( opt, 'background_stars' )
  opt.bckgrnd_str.add = 0 ;
  else 
  opt.bckgrnd_str.add = opt.background_stars.add ;
  % For each star in the background, one needs to specify its: type, flux or apparent magnitude, distance, proper motion, relative location to the center of the image

  % The number of background stars may be specified beforehand
    if isfield( opt.background_stars, 'number' )
    opt.bckgrnd_str.n = opt.background_stars.number ;
    end

  % Type of star: 'a0v', 'a5v', 'f5v', 'g0v', 'g5v' (closer to sun), 'k0v', 'k5v', 'm0v', 'm5v', or 'sun'
    if ~isfield( opt.background_stars, 'type' )
    opt.bckgrnd_str.tp = 'sun' ;
    else
    opt.bckgrnd_str.tp = lower( opt.background_stars.type ) ;
    end

  % If the number of backgroud stars was not assigned, set it now
    if ~isfield( opt.bckgrnd_str, 'n' )
    opt.bckgrnd_str.n = numel( opt.bckgrnd_str.tp ) ;
    end

  % Check there are no more background stars than types if the number was set beforehand
    if isfield( opt.background_stars, 'number' )
      if ( opt.bckgrnd_str.n < numel( opt.bckgrnd_str.tp ) )
      disp( sprintf( 'There are %i types of stars set for %i stars. Please, set the number of types to be less or equal to the number of stars. Stopped.', opt.bckgrnd_str.n, numel( opt.bckgrnd_str.tp ) ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end

  % Assigning random types if the number of stars was set beforehand, *and* there are not enough star types to assign a one-to-one correspondance
  if isfield( opt.background_stars, 'number' )
    if ( opt.bckgrnd_str.n > numel( opt.bckgrnd_str.tp ) )
    % temporary copy
    bckgrnd_str_tp_in = opt.bckgrnd_str.tp ;
      for i_str = 1 : opt.bckgrnd_str.n
        for i_str = 1 : opt.bckgrnd_str.n
        opt.bckgrnd_str.tp{ i_str } = bckgrnd_str_tp_in{ ceil( numel( bckgrnd_str_tp_in ) * rand( 1, 1 ) ) } ;
        end
      end
    end
  end

% Adding background stars (TB continued)
disp( ' ' )
disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop

  % Star flux (rare: usually one will set the apparent magnitude)
    if ~isfield( opt.star, 'flux' ) && ~isfield( opt.str, 'flx' )
    opt.str.flx = NaN ; % not used
    else
      if isfield( opt.star, 'flux' )
      opt.str.flx = opt.star.flux ;
      end
    end

  % Magnitude. Recall for the Sun: apparent (-26.76 V), 5 pc (3.33 V) and absolute 10 pc (4.83 V). http://mips.as.arizona.edu/~cnaw/sun.html)
    if ~isfield( opt.star, 'app_mag_v' )
    opt.str.mg = 4.81 ; % Sun at 10 pc
    else
    opt.str.mg = opt.star.app_mag_v ;
    end
  % Name
    if ~isfield( opt.star, 'name' )
    opt.str.nm = 'Sun-like' ;
    else
    opt.str.nm = opt.star.name ;
    end

  % Mass in solar masses (for orbits)
    if ~isfield( opt.star, 'times_sun_mass' )
    opt.str.mss_sn = 1 ; % 1 solar mass
    else
    opt.str.mss_sn = opt.star.times_sun_mass ;
    end
  end % Background stars

% Create an exo-zodi emission for the scenes
  if ~isfield( opt, 'exozodi' )
  opt.zd.fct = 0 ;
  else
    if ~isfield( opt.exozodi, 'radius_au' )
    opt.zd.rd_au = 0.1 ; % This is the radius at which Bertrand Menesson's exo-zodi data fall to half its peak value. This zodi was obtained by running zodipic (see convolve_with_one_wavelength.m and search for radius_au). Therefore, it follows the solar system zodiacal light in its spatial structure. We assume it here as some sort of reference, until zodipic or an exozodiacal emission code is added to this suite.
    else
    opt.zd.rd_au = opt.exozodi.radius_au ;
    % Setting the value in mas
      if ~isfield( opt.exozodi, 'radius_mas' )
      opt.exozodi.radius_mas = opt.zd.rd_au * opt.au2mas ;
      else
        if ( opt.exozodi.radius_mas == opt.exozodi.radius_au * opt.au2mas )
        disp( sprintf( 'The half radius intensity radius of the exozodiacal emission has been set twice to an equivalent value. opt.exozodi.radius_au=%f, and opt.exozodi.radius_mas=%f. Make sure to set only one. Continuing', opt.exozodi.radius_au, opt.exozodi.radius_mas ) )
        else
        disp( sprintf( 'The half radius intensity radius of the exozodiacal emission has been set twice to different values. opt.exozodi.radius_au=%f, and opt.exozodi.radius_mas=%f. Choose one value. Stopped', opt.exozodi.radius_au, opt.exozodi.radius_mas ) )
        disp( ' ' )
        disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
        end
      end
    end
    if ~isfield( opt.exozodi, 'radius_mas' )
    opt.zd.rd_mas = opt.zd.rd_au * opt.au2mas ; % 0.1 AU, see above.
    else
    opt.zd.rd_mas = opt.exozodi.radius_mas ;
    % Setting the value in AU
    opt.exozodi.radius_au = opt.zd.rd_mas * opt.mas2au ;
    opt.zd.rd_au = opt.exozodi.radius_au ;
    end

    % Inclination plane (intrinsic, wrt the invariant plane, which should be 0 deg in general. An observational inclination is defined with opt.scene.observing_inclination)
    if ~isfield( opt.exozodi, 'inclination_deg' )
    opt.zd.inc_dg = 0 ; 
% Add that if planet inclination is set, the exozodi inclination equals the one of the first planet (let outer planets have some different inclination)
    else
    opt.zd.inc_dg = opt.exozodi.inclination_deg ;
    end
    % Same but in radians
    if isfield( opt.exozodi, 'inclination_rad' )
      if  ~( opt.zd.inc_dg )
      opt.zd.inc_dg = opt.exozodi.inclination_rad * 180 / pi ;
      else
        if ( opt.zd.inc_dg * pi / 180 == opt.exozodi.inclination_rad )
        disp( sprintf( 'WARNING: the inclination of the exozodiacal emission has been set twice to the same value. In degrees is %f, and in rad is %f. Make sure to use one only. Continuing ...', opt.zd.inc_dg, opt.exozodi.inclination_rad ) )
        else
        disp( sprintf( 'The inclination of the exozodiacal emission has been set twice. In degrees is %f, and in rad is %f, which are not strictly the same. Choose one of them. Stopped.', opt.zd.inc_dg, opt.exozodi.inclination_rad ) )
        make_a_Stop
        end
      end
    end
    % Rotation about the visual axis (position angle, https://en.wikipedia.org/wiki/Position_angle. Equivalent to the Kepler parameter line of ascending node).
    if ~isfield( opt.exozodi, 'position_angle_deg' )
    opt.zd.pa_dg = 90 ; % Rotates the exo-zodi to be horizontal
    else
    opt.zd.pa_dg = opt.exozodi.position_angle_deg ;
    end
    % Same but in radians
    if isfield( opt.exozodi, 'position_angle_rad' )
    opt.zd.pa_dg = opt.exozodi.position_angle_rad * 180 / pi ;
      if isfield( opt.exozodi, 'position_angle_deg' )
        if ( opt.exozodi.position_angle_deg * pi / 180 == opt.exozodi.position_angle_rad )
        disp( sprintf( 'WARNING: the position angle of the exozodiacal emission has been set twice to the same value. In degrees is %f, and in rad is %f. Make sure to use one only. Continuing ...', opt.zd.pa_dg, opt.exozodi.position_angle_rad ) )
        else
        disp( sprintf( 'The position angle of the exozodiacal emission has been set twice. In degrees is %f, and in rad is %f, which are not strictly the same. Choose one of them. Stopped.', opt.zd.ps_dg, opt.exozodi.position_angle_rad ) )
        make_a_Stop
       end
      end
    end

    % The factor by which re-scale the zodi
    if ~isfield( opt.exozodi, 'factor' )
    opt.zd.fct = 1 ;
    else
    opt.zd.fct = opt.exozodi.factor ;
    end
  % These options are specific to some work and may not be used in general
    % Adding different extragalactic background fields
    if ~isfield( opt.exozodi, 'bck_lbl' )
    opt.zd.bck_lbl = '' ;
    else
    opt.zd.bck_lbl = opt.exozodi.bck_lbl ;
    end
  end % Exozodi emission 

% Exo-Kuiper belt
  if ~isfield( opt, 'exokuiper' )
  opt.kpr.do = 0 ; % Not included
  else
    if ~isfield( opt.exokuiper, 'do' )
    opt.kpr.do = 0 ; % Not included
    else
    opt.kpr.do = opt.exokuiper.do ;
    end

  % Kuiper interior radius in AU
    if ~isfield( opt.exokuiper, 'radius_int_au' )
    opt.kpr.int_au = 30 ;
    else
    opt.kpr.int_au = opt.exokuiper.radius_int_au ;
    end

  % Kuiper exterior radius in AU
    if ~isfield( opt.exokuiper, 'radius_ext_au' )
    opt.kpr.ext_au = 50 ;
    else
    opt.kpr.ext_au = opt.exokuiper.radius_ext_au ;
    end
  
  % Check of consistency
    if ( opt.kpr.int_au >= opt.kpr.ext_au ) 
    disp( sprintf( 'WARNING: The inner radius of the exo-Kuiper belt must be less than the exterior radius' ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end

  % Limit of the FOV of the Haystacks Project
    if ( opt.kpr.ext_au > 50 )
    disp( sprintf( 'WARNING: the value of the outer radius of the exo-Kuiper belt is %2.2f AU, greater than 50 AU. The resulting image may be cropped. Choose a smaller radius to fix any issues.', opt.kpr.ext_au ) )
    end

  % Minimum relative factor for the forward scattering
    if ~isfield( opt.exokuiper, 'minimum_scattering' )
    opt.kpr.mn_scttrng = 0.7 ;
    else
    opt.kpr.mn_scttrng = opt.exokuiper.minimum_scattering ;
    end

  % Check of consistency
    if ( opt.kpr.mn_scttrng < 0 ) || ( opt.kpr.mn_scttrng > 1 )
    disp( sprintf( 'WARNING: the minimum relative factor for the surface brightness of the exo-Kuiper belt is %2.2f, and should be within 0 and 1', opt.kpr.mn_scttrng ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end

    % Maximum relative factor for the forward scattering
    if ~isfield( opt.exokuiper, 'maximum_scattering' )
    opt.kpr.mx_scttrng = 2.1 ;
    else
    opt.kpr.mx_scttrng = opt.exokuiper.maximum_scattering ;
    end

  % Check of consistency
    if ( opt.kpr.mx_scttrng < 1 ) 
    disp( sprintf( 'WARNING: the maximum relative factor for the surface brightness of the exo-Kuiper belt is %2.2f, and should be greater than or equal to 1', opt.kpr.mx_scttrng ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end  

  end % exo-Kuiper

% Add background stars (TBW)
  if ~isfield( opt, 'stars' )
  opt.str.do = 0 ;
  else
    if ~isfield( opt.stars, 'do' )
    opt.str.do = 0 ;  
    else
    opt.str.do = opt.stars.do ;
    end
  end % Stars field

% Add solar glint
  if ~isfield( opt, 'solar_glint' )
  opt.slr_glnt.do = 0 ;
  else
    if ~isfield( opt.solar_glint, 'do' )
    opt.slr_glnt.do = 0 ;
    else
    opt.slr_glnt.do = opt.solar_glint.do ;
    end
  % Options that go with the laboratory data
    % Sub-directory where the lab data are to be found
    if ~isfield( opt.solar_glint, 'dir' )
    opt.slr_glnt.dr = [ opt.scn_dr 'solar_glint/' ] ;
    else
    opt.slr_glnt.dr = opt.solar_glint.dir ;
    end
    % Whether there's stealth included or not
    if ~isfield( opt.solar_glint, 'stealth' )
    opt.slr_glnt.stlth = 0 ;
    else
    opt.slr_glnt.stlth = opt.solar_glint.stealth ;
    end
  % The stealth option can only be used with a non-spinning starshade
    if ( opt.slr_glnt.stlth ) && ~strcmp( opt.starshade.mode, 'non-spinning' )
    disp( 'The option to add stealth edges is meant for a non-spinning starshade only. Your configuration has a spinning starshade and the stealth option. Please, fix either one. Stopped.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
    % star-starshade-sun angle. (Stuart Shaklan 10/17/18: Phi = 90 has the sun off to the side of the starshade.  Phi = 0 has the sun behind the starshade.  Missions will generally restrict phi angles between 45 and 85 degrees 0 deg means it is directly behind the starhade. For WFIRST the minimum angle is 54 deg and the maximum angle is 83 deg)
    if ~isfield( opt.solar_glint, 'phi_deg' )
    opt.slr_glnt.phi_dg = 60 ;
    else
    opt.slr_glnt.phi_dg = opt.solar_glint.phi_deg ;
    end
    if ~isfield( opt.solar_glint, 'alpha_deg' )
    opt.slr_glnt.alph_dg = 0 ; % rotation of sun around the plane.  zero is along x axis.
    else
    opt.slr_glnt.alph_dg = opt.solar_glint.alpha_deg ;
    end
    % Radius (in m) used to define the 'trailing edge' of the petals
    if ~isfield( opt.solar_glint, 'r_m' )
    opt.slr_glnt.r_m = 7.5 ; % Default value for WFIRST (Starshade radius of 13 m)
    else
    opt.slr_glnt.r_m = opt.solar_glint.r_m ;
    end
    % Basic check (reference is 7.5 m for a Starshade radius of 13 m, NI2)
      if isfield( opt, 'geo_iwa_mas' )
      r_m_mn = 7 * opt.geo_iwa_mas * opt.dst_strshd_tlscp_m / 1e3 / 3600 / 180 * pi / 13 ;
      r_m_mx = 8 * opt.geo_iwa_mas * opt.dst_strshd_tlscp_m / 1e3 / 3600 / 180 * pi / 13 ;
      r_m_avg = 7.5 * opt.geo_iwa_mas * opt.dst_strshd_tlscp_m / 1e3 / 3600 / 180 * pi / 13 ;
        if ( opt.slr_glnt.do )
          if ( opt.slr_glnt.r_m < r_m_mn ) || ( opt.slr_glnt.r_m > r_m_mx )
          disp( sprintf( 'WARNING. The value of opt.solar_glint.r_m=%2.1f used in the solar glint is not close to the expected value of %2.1f m. Continuing, but make sure this value of %2.1f m is intentional.', opt.slr_glnt.r_m, r_m_avg, opt.slr_glnt.r_m ) ) ;
          end
        end
      end

  % Lab data files
    if ~isfield( opt.solar_glint, 'zs_filename' )
    opt.slr_glnt.zs_fl = 'AverageSASCalibrationCoupon_S_20191023_nd7p18.mat' ; % 'GemRazor_02_23_2016_MeshS_0_02_26_2018_1.mat' ;  % SBS filename change 101519
    else
    opt.slr_glnt.zs_fl = opt.solar_glint.zs_filename ;
    end
    if ~isfield( opt.solar_glint, 'zp_filename' )
    opt.slr_glnt.zp_fl = 'AverageSASCalibrationCoupon_P_20191023_nd7p18.mat' ; % 'GemRazor_02_23_2016_MeshP_0_02_26_2018_1' ;
    else
    opt.slr_glnt.zp_fl = opt.solar_glint.zp_filename ;
    end
    if ~isfield( opt.solar_glint, 'zs2_filename' )
    opt.slr_glnt.zs2_fl = 'AverageSASCalibrationCoupon_S_20191023_nd7p18.mat'; % temporary until we get stealth data% % Stealth data
    else
    opt.slr_glnt.zs2_fl = opt.solar_glint.zs2_filename ;
    end
    if ~isfield( opt.solar_glint, 'zp2_filename' )
    opt.slr_glnt.zp2_fl = 'AverageSASCalibrationCoupon_P_20191023_nd7p18.mat' ; % temporary until we get stealth data% % Stealth data
    else
    opt.slr_glnt.zp2_fl = opt.solar_glint.zp2_filename ;
    end


  % Options for the case when two solar lobes are added, instead of using laboratory data. It should only be used for some tests
    % Positions are with respect the center of the scene (Stuart Shaklan's email from 02/06/18 ""Solar lobes")
    if ~isfield( opt.solar_glint, 'pos_arc_dec_mas' )
    opt.slr_glnt.arc_dc_mas = [ 48 * cos( 75 / 180 * pi ), 48 * cos( 75 / 180 * pi ) ] ;
    else
    opt.slr_glnt.arc_dc_mas = opt.solar_glint.pos_arc_dec_mas ;
    end
    if ~isfield( opt.solar_glint, 'pos_arc_ra_mas' )
    opt.slr_glnt.arc_ra_mas = [ -48 * sin( 75 / 180 * pi ), 48 * sin( 75 / 180 * pi ) ] ;
    else
    opt.slr_glnt.arc_ra_mas = opt.solar_glint.pos_arc_ra_mas ;
    end
    if ~isfield( opt.solar_glint, 'delta_mag_apparent' )
    opt.slr_glnt.dlt_mg = 53 ; % Stuart Shaklan's email 'Lobes magnitude in wavelength' 02/13/18 
    else
    opt.slr_glnt.dlt_mg = opt.solar_glint.delta_mag_apparent ;
    end
  end % Solar glint field

% Add perturbed locus
  if ~isfield( opt, 'locus' )
  opt.lcs.do = 0 ;
  else
    if ~isfield( opt.locus, 'do' )
    opt.lcs.do = 0 ;
    else
    opt.lcs.do = opt.locus.do ;
    end
      if ( opt.lcs.do )
        if ~isfield( opt.locus, 'perturbed_filename' )
        opt.lcs.prtrbd = 'NI2_test_case_1em10' ; % Perturbation consistent with WFIRST requirements
        else
        opt.lcs.prtrbd = opt.locus.perturbed_filename ;
        % Removing the mat extension if present
          if length( opt.lcs.prtrbd ) > 3
            if strcmp( opt.lcs.prtrbd( end - 3 : end ), '.mat' )
            opt.lcs.prtrbd = opt.lcs.prtrbd( 1 : end - 4 ) ;
            end
          end   
        end
        if ~isfield( opt.locus, 'unperturbed_filename' )
        opt.lcs.unprtrbd = opt.starshade.nominal_filename ;
        else
        opt.lcs.unprtrbd = opt.locus.unperturbed_filename ;
          % Removing the mat extension if present
           if strcmp( opt.lcs.unprtrbd( end - 3 : end ), '.mat' )
           opt.lcs.unprtrbd = opt.lcs.unprtrbd( 1 : end - 4 ) ;
           end
        end
        % Check of consistency. So far, baseline occulters are identified by the first 3 charcters: NI2, NW2, TV3, UH17, etc.
        if ~strcmp( opt.lcs.prtrbd( 1 : 3 ), opt.starshade.nominal_filename( 1 : 3 ) ) || ~strcmp( opt.lcs.unprtrbd( 1 : 3 ), opt.starshade.nominal_filename( 1: 3 ) )
        disp( sprintf( 'The occulter used for the locus change (%s) does not coincide with the one used for the optical imaging (%s). Fix it.', opt.lcs.prtrbd( 1 : 3 ), opt.starshade.nominal_filename( 1 : 3 ) ) )
        return
        end
        % Pupil file. Same as set before for the on-axis occulter.
        if ~isfield( opt.locus, 'pupil_filename' )
        opt.lcs.lbl_pl = opt.pupil_filename ;
        else
        opt.lcs.lbl_pl = opt.locus.pupil_filename ;
        end   
        % Main directory with the files for the locus
        if ~isfield( opt.locus, 'install_dir' )
        opt.lcs.inst_dr = opt.scn_dr ;
        else
        opt.lcs.inst_dr = opt.locus.install_dir ;
        end
        % Subdirectories where the locus files are located
        if ~isfield( opt.locus, 'in_dir' )
        opt.lcs.in_dr = 'locus/in/' ;
        else
        opt.lcs.in_dr = opt.locus.in_dir ;
        end
        % Subdirectories where the output of the locus work are stored (if they are)
        if ~isfield( opt.locus, 'out_dir' )
        opt.lcs.out_dr = 'locus/out/' ;
        else
        opt.lcs.out_dr = opt.locus.out_dir ;
        end
        % Second lambda in case a range is to be used. Otherwise, run_new_locus will set lambda_2_nm = lambda_1_nm.
        if isfield( opt.locus, 'lambda_2_nm' )
        opt.lcs.lambda_2_nm = opt.lmbd_img_1_nm ;
        end
        % Don't plot by default
        if ~isfield( opt.locus, 'do_plot' )
        opt.lcs.do_plt = 0 ;
        else
        opt.lcs.do_plt = opt.locus.do_plot ;
        end
        if ~isfield( opt.locus, 'px_res_mas') ;
        opt.lcs.px_res_mas = opt.px_psf_mas ;
        else
        opt.lcs.px_res_mas = opt.locus.px_res_mas ;
        end
        % Size of the image plane in mas
        if ~isfield( opt.locus, 'diam_img_mas')
        opt.lcs.diam_img_mas = 1001 ;
        else
        opt.lcs.diam_img_mas = opt.locus.diam_img_mas ;
        end
        % Replicating the value (in case a perturbed locus is run, this option goes with the locus options)
        if isfield( opt.starshade, 'number_of_petals' )
        opt.lcs.n_ptl = opt.starshade.number_of_petals ;
        end

      % For the case of a spinning starshade, one may choose between two options: i) set the step angle for every rotation and perform a full circle, or ii) perform as many rotations as necessary to attain a certain goal. See run_new_locus.m for details on the algorithm.

        % Independently of the previous choice, one may choose to use imrotate.m (it belongs to the Image Processing Toolbox), instead of a true propagation of the optical diffraction code. It is not recommended for studies that require some precision, since imrotate fails to reproduce the PSF response of rotated starshades for WFIRST with errors of few percent.
        if ~isfield( opt.locus, 'imrotate' )
        opt.lcs.imrotate = 0 ; % Do not use imrotate
        else
        opt.lcs.imrotate = 1 ;
        end

        % i) set the step angle for every rotation and perform a full circle
        if ~isfield( opt.locus, 'rotation_angle_step_deg' )
        opt.lcs.alph_lst_rd = NaN ; % By default don't do all the rotations
        else
        opt.lcs.alph_lst_rd = pi / 180 * ( 0 : opt.locus.rotation_angle_step_deg : 360 ) ;
        end

        % ii) perform as many rotations as necessary to attain a certain goal: for a fraction of an Earth-like planet, or for a given contrast value.
        % In case the absolute precision before is not set, one may set the relative contribution to an Earth-like planet (see run_new_locus.m for details on the algorithm)
        if ~isfield( opt.locus, 'perturbed_psf_relative_to_earth' )
        opt.lcs.rltv_earth = 0.1 ; % 10%
        else
        opt.lcs.rltv_earth = opt.locus.perturbed_psf_relative_to_earth ;
        end

        % This parameter sets the precision required when deriving the perturbed PSF for a *spinning* starshade. It's equivalent to a contrast error (see run_new_locus.m for details on the algorithm).
        if ~isfield( opt.locus, 'perturbed_psf_precision' )
        opt.lcs.epsln = NaN ;
        else
        opt.lcs.epsln = opt.locus.perturbed_psf_precision ; % For example, 1e-11
        end


        % Size of the image plane in mas. It usually happens that the non-ideal starshade has a level of starlight about 100x higher than the nominal case. The following parameter controls the spatial extent of the non-ideal on-axis PSF that will be generated and used in the imaging simulations.
        if ~isfield( opt.locus, 'times_nominal' )
        opt.lcs.tm_nmnl = 2 ; % A factor of 2 gives results for WFIRST that are 0.5% different than a higher sampling choice. For other occulters, it may need to be larger if the extension of the non-ideal PSF on the image needs to be larger thna with the default option.
        else
        opt.lcs.tm_nmnl = opt.locus.times_nominal ;
        end

        % Check of consistency
        if ( opt.lcs.tm_nmnl < 1)
        disp( sprintf( 'Please, choose a value of opt.locus.times_nominal greater than 1 to derive the SPF response of a non-ideal starshade. Right now, opt.locus.times_nominal=%1.2. Stopped.', opt.lcs.tm_nmnl ) )
        disp( ' ' )
        disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
        end

        % Make sure that when imrotate is used, there's enough room for cropping the rotated result
        if ( opt.lcs.imrotate ) && ( opt.lcs.tm_nmnl < sqrt( 2  ) ) 
        disp( sprintf( 'The value of the extension for the non-ideal PSF (%2.1f) is not large enough to use imrotate and cropping. Increasing it to sqrt(2)', op.lcs.tm_nmnl ) )
        opt.lcs.tm_nmnl = sqrt( 2 ) ;
        end

        % Additional, technical options
        % Re-doing the work? If 1, it will re-compute the electric fields at the pupil plane and at the image plane, and the PSF intensity from scratch.
        % If 0, it will check if the corresponding results exist. If they exist, they are loaded. If they don't, they will be created.
        if ~isfield( opt.locus, 'redo_reference' )
        opt.lcs.redo_reference = 0 ; % Default 0
        else
        opt.lcs.redo_reference = opt.locus.redo_reference ;
        end
        if ~isfield( opt.locus, 'redo_perturbed' )
        opt.lcs.redo_perturbed = 0 ; % Default 0
        else
        opt.lcs.redo_perturbed = opt.locus.redo_perturbed ;
        end
        % Storing the results (electric fields, pupil, image plane intensity and options)
        if ~isfield( opt.locus, 'save_work' ) && strcmp( lower( opt.starshade.mode ), 'spinning' )
        opt.lcs.save_work = 1 ;
        end
        if ~isfield( opt.locus, 'save_work' ) && strcmp( lower( opt.starshade.mode ), 'non-spinning' )
        opt.lcs.save_work = 0 ;
        end
        if isfield( opt.locus, 'save_work' )
        opt.lcs.save_work = opt.locus.save_work ;
        end
      end
  end % Locus field

% Extracting the data into a data cube
  if ~isfield( opt, 'cube' )
  opt.cb.do = 0 ;
  opt.cube.do = 0 ;
  opt.cb.plnt.do = 0 ;
  else
    if isfield( opt.cube, 'do' )
    opt.cb.do = opt.cube.do ;
      if ( opt.cube.do )
      opt.cb.do = 1 ;
      % Subdirectory where the FITS files will be stored
        if ~isfield( opt.cube, 'dir' )
        opt.cb.dr = 'cubes' ;
        else
        opt.cb.dr = opt.cube.dir ;
        end
      % FITS format or Matlab
        if ~isfield( opt.cube, 'fits' )
        opt.cb.fits = 0 ; % By default, store the output in Matlab files
        else
        opt.cb.fits = opt.cube.fits ;
        end
      % Additional label for the file name
        if ~isfield( opt.cube, 'label' )
        opt.cb.lbl = '' ;
        else
        opt.cb.lbl = opt.cube.label ;
        end
      % Cube FOV (see later the option that extracts data around eaxh planet)
        if ~isfield( opt.cube, 'fov_diam_mas' )
        opt.cb.fov_diam_mas =  1500 ; 
        end
      % Center of the FOV that will be extracted wrt center of the scene
        if ~isfield( opt.cube, 'ra_mas' )
        opt.cb.ra_mas = 0 ;
        else
        opt.cb.ra_mas = opt.cube.ra_mas ;
        end
        if ~isfield( opt.cube, 'dec_mas' )
        opt.cb.dc_mas = 0 ;
        else
        opt.cb.dc_mas = opt.cube.dec_mas ;
        end
      % Cubes around planets
        if ~isfield( opt.cube, 'planets' )
        opt.cb.plnt.do = 0 ; % By default, do not create cubes about the planets.
        opt.cube.planets.do = 0 ;
        else
        opt.cb.plnt.do = opt.cube.planets.do ;
        end
      % Area around each planet that gets stored in terms of fwhm=1.028*l/D (fwhm gets redefined at each wavelength)
        if ~isfield( opt.cube.planets, 'n_fwhm' )
        opt.cb.plnt.n_fwhm = 3 ; % Default 3 FWHM (+/- 1.5 around the planet position)
        else
        opt.cb.plnt.n_fwhm = opt.cube.planets.n_fwhm ;
        end
      % Removing the background
        if ~isfield( opt.cube, 'background' )
        opt.cb.bckgrnd = 1 ; % Default with background
        else
        opt.cb.bckgrnd = opt.cube.background ;
          if ~( opt.cb.bckgrnd )
          % Removing and adding the same planets in the case of not having removed the planets already
            if ~( opt.plnt.add.do )
            opt.plnt.rmv.do = 1 ;
            opt.plnt.add.do = 1 ;
            end
          end
        end
      end
    end
  end % SISTER data cubes

% Optical throughput (Obscuration, QE, and PSF core are accounted for by the simulation itself or other parameters of the simulations)

  % Cover both spellings
  if isfield( opt, 'optical_thruput' )
  opt.optical_throughput = opt.optical_thruput ; 
  end

  if ~isfield( opt, 'optical_throughput' )
    if strcmp( lower( opt.starshade.nominal_filename( 1: 3 ) ), 'ni2' )
    % Andrew Romero-Wolf compiled som evalues from different sources (05/29/18 Starshade Rendezvous Mission Study, Excel sheet). In particular, for WFISRT he shows a table constructed by Hong Tang (JPL) with the following values:
    % Wavelength (nm)  | Direct Imaging  |    IFS
    %      575         | CBE@EOL: 0.56   |   N/A   
    %                  | REQ@EOL: 0.46   |   N/A
    %      825         | CBE@EOL: 0.52   |   N/A
    %                  | REQ@EOL: 0.42   |   N/A
    %      660         | CBE@EOL:  N/A   |  0.43
    %                  | REQ@EOL:  N/A   |  0.29
    %      760         | CBE@EOL:  N/A   |  0.39
    %                  | REQ@EOL:  N/A   |  0.26
    
    % Use Current Best Estimates by default
      if ~isfield( opt, 'optical_cbe' )
      opt.optical_cbe = 1 ;
      end
      if ~isfield( opt, 'optical_req' )
      opt.optical_req = 0 ;
      end
      % If 'requirements' is set, then turn off 'current best estimate'
      if ( opt.optical_req ~= 0 )
      opt.optical_cbe = 0 ;
      end
      % The other way, too
      if ( opt.optical_cbe ~= 0 )
      opt.optical_req = 0 ;
      end
      
      % Direct Imaging
      opt.optical_throughput = NaN ;
      if ~( opt.cb.do )
        if ( opt.lambda_1_nm >= 425 ) && ( opt.lambda_2_nm <= 575 )
          if ( opt.optical_cbe )
          opt.optical_throughput = 0.56 ;
          end
          if ( opt.optical_req )
          opt.optical_throughput = 0.46 ;
          end
        end
        if ( opt.lambda_1_nm >= 606 ) && ( opt.lambda_2_nm <= 787 )
          if ( opt.optical_cbe )
          opt.optical_throughput = 0.54 ; % Mean of 575 and 825 nm data
          end
          if ( opt.optical_req )
          opt.optical_throughput = 0.44 ;
          end
        end
        % 615-800 nm: WFIRST-Starshade probe study
        if ( opt.lambda_1_nm >= 615 ) && ( opt.lambda_2_nm <= 800 )
          if ( opt.optical_cbe )
          opt.optical_throughput = 0.28 * 1.25 ; 
          end
          if ( opt.optical_req )
          opt.optical_throughput = 0.28 ;
          end
        end
        if ( opt.lambda_1_nm >= 747 ) && ( opt.lambda_2_nm <= 970 )
          if ( opt.optical_cbe )
          opt.optical_throughput = 0.52 ;
          end
          if ( opt.optical_req )
          opt.optical_throughput = 0.42 ;
          end
        end
      % IFS (Wavelength slices)
      else
        if ( opt.lambda_1_nm >= 615 ) && ( opt.lambda_2_nm <= 800 )
          if ( opt.optical_cbe )
          opt.optical_throughput = 0.43 ;
          end
          if ( opt.optical_req )
          opt.optical_throughput = 0.29 ;
          end
        end
        if ( opt.lambda_1_nm >= 656 ) && ( opt.lambda_2_nm <= 800 )
          if ( opt.optical_cbe )
          opt.optical_throughput = 0.43 ;
          end
          if ( opt.optical_req )
          opt.optical_throughput = 0.29 ;
          end
        end
        if ( opt.lambda_1_nm >= 820 ) && ( opt.lambda_2_nm <= 1000 )
          if ( opt.optical_cbe )
          opt.optical_throughput = 0.39 ; 
          end
          if ( opt.optical_req )
          opt.optical_throughput = 0.26 ;
          end
        end
      end
  
      % Final check
      if isnan( opt.optical_throughput )
        if ( opt.cb.do )
        lbl_inst = 'IFS' ;
        else
        lbl_inst = 'Direct Imager' ;
        end
      disp( sprintf( 'The value of the optical throughput for the wavelength range %04.2f-%04.2f nm in the %s has NOT been set. Please check the band corresponds to one of the IFS bands in WFIRST.', opt.lambda_1_nm, opt.lambda_2_nm, lbl_inst ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end % NI2

    % For any other occulter (for now assuming constant optical throughput)
    if ~( numel( strfind( lower( opt.starshade.nominal_filename ), 'ni2' ) ) )
    % Optical throughput
      if ~isfield( opt, 'optical_throughput' )
      disp( 'WARNING: please set opt.optical_throughput to some value in your configuration file. Stopped.' )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end % For any other occulter but NI2

  disp( sprintf( 'The value of the optical throughput is %1.3f', opt.optical_throughput ) ) 
  end % optical_throughput

% Quantum efficiency: do not set it here with default values. In convolve_with_one_wavelength.m it is defined for some special cases.

% In the general case, where the optical throughput and QE are set un the parameter file, we may want to include some loses (typical CBE value)
  if ~isfield( opt, 'cbe_optical_throughput_fct' )
  opt.cbe_optical_throughput_fct = 1 ;
  end
  if ~isfield( opt, 'cbe_qe_fct' )
  opt.cbe_qe_fct = 1 ;
  end

% Detector noise
  if ~isfield( opt, 'noise' )
  opt.ns.do = 0 ;
  % Detector gain (even if no noise is to be simulated, we set some gain for the signal)
  opt.ns.gn = 1 ;
  % Exposure time
  opt.ns.exp_tm_s = 3600 ;
  % Number of frame exposures
  opt.ns.n_frm = 1 ;
  else
    if ~isfield( opt.noise, 'do' )
    opt.ns.do = 0 ;
    else
    opt.ns.do = opt.noise.do  ;
    end
    if ~isfield( opt.noise, 'gain' )
    opt.ns.gn = 1 ;
    else
    opt.ns.gn = opt.noise.gain ;
    end
  % Make sure is greater or equal than 1
    if ( opt.ns.gn < 1 )
    disp( sprintf( 'WARNING: the value for the detector''s gain is %f, but it must be greater than or equal to 1. Please fix it.', opt.ns.gn ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end  
    if ~isfield( opt.noise, 'shot' )
    opt.ns.sht = 1 ;
    else
    opt.ns.sht =  opt.noise.shot ;
    end
    if ~isfield( opt.noise, 'dark_e_s' )
    opt.ns.drk_e_s = 1e-3 ; % e/pix/s
    else
    opt.ns.drk_e_s = opt.noise.dark_e_s ;
    end
    if ~isfield( opt.noise, 'cic_e' )
    opt.ns.cic_e = 1.3e-3 ; % e/pix/read. 1.3e-3 is the current best estimate for several EMCCDs with opt.noise.gain=1000 (that is, if the gain is set correspondingly, CIC becomes similar to a readout noise of 1.3 e/pix/read)
    else
    opt.ns.cic_e = opt.noise.cic_e ;
    end
    if ~isfield( opt.noise, 'read_e' )
    opt.ns.rd_e = 3 ; % e/pix/read
    else
    opt.ns.rd_e = opt.noise.read_e ;
    end
    if ~isfield( opt.noise, 'exp_time_total_sec' )
    opt.ns.exp_tm_s = 3600 ; % 1 hour
    else
    opt.ns.exp_tm_s = opt.noise.exp_time_total_sec ;
    end
    if ~isfield( opt.noise, 'exp_time_frame_sec' )
    opt.ns.n_frm = ceil( opt.ns.exp_tm_s / 600 ) ;
    else
    opt.ns.n_frm = ceil( opt.ns.exp_tm_s / opt.noise.exp_time_frame_sec ) ;
    end
  % Parameter to control when the Poisson distribution can be approximated by the Normal distribution
    if ~isfield( opt.noise, 'lambda_normal' )
    opt.ns.lmbd_nrml = 20 ;
    else
    opt.ns.lmbd_nrml = opt.noise.lambda_normal ;
    end

  end % opt.noise

% Jitter of the telescope
  if ~isfield( opt, 'jitter' )
  opt.jitter.do = 0 ;
  opt.jttr.do = 0 ;
  else
  opt.jttr = opt.jitter ;
  % RMS of the jitter in mas (treated as a circular 2-D Gaussian with sigma=RMS)
    if ~isfield( opt.jitter, 'rms_mas' )
    opt.jttr.rms_mas = 14 ; % 14 mas, mission specs for WFIRST (09/2019: https://wfirst.ipac.caltech.edu/sims/Param_db.html)
    else
    opt.jttr.rms_mas = opt.jitter.rms_mas ;
    end
  end

% Label related with imaging modes that may be recycled for adding planets, lobes, locus errors and stars.
  if ~isfield( opt, 'tag_name' )
  % First label is related with the diffuse part (long convolution time)
  opt.tg_nm_1 = sprintf( 'rp%iab%imz%ssl%ilcs%ist%i', opt.plnt.rmv.do, opt.bckgrnd.add, strrep( sprintf( '%2.2f', opt.zd.fct ), '.', 'p' ), opt.slr_glnt.do, opt.lcs.do, opt.str.do ) ;
    if ~strcmp( opt.starshade.nominal_filename, 'NI2' )
    opt.tg_nm_1 = [ opt.tg_nm_1 '_' opt.starshade.nominal_filename ] ;
    end
  % Adding nx_pupil_px and if pupil is ideal
  opt.tg_nm_1 = sprintf( '%s_Nx%i', opt.tg_nm_1, opt.nx_pupil_pix ) ;
    if isfield( opt, 'pupil_filename' )
      if strcmp( opt.pupil_filename, '0' ) == 1
      opt.tg_nm_1 = [ opt.tg_nm_1 '_ideal' ] ; 
      end
    end

  % Second label with point-like additions (short convolution time)
  opt.tg_nm_2 = sprintf( 'ap%i', opt.plnt.add.do ) ;
  opt.tg_nm = [ opt.tg_nm_1 '_' opt.tg_nm_2 ] ;
  else
  opt.tg_nm = opt.tag_name ;
  end
 
% Additional, optional label
  if isfield( opt, 'add_label' ) && length( opt.add_label )
  opt.tg_nm = [ opt.tg_nm '_' opt.add_label ] ;
  end

% Label related with a series of images (by default, only one image, no label)
  if ~isfield( opt, 'tg_lp' )
  opt.tg_lp = '' ;
  end

% This is used to keep track of the case when 2 scenes are imaged
  if ~isfield( opt, 'n_scnr' )
  opt.n_scnr = 1 ; % One scenario
  end

% Ground telescope (e.g., Remote Occulter)
% Basically, we include in SISTER three scenarios for the sky radiance and transmittance depending on the moon phase. Each of them, includes scattered Moonlight, local zodiacal light, scattered starlight, emisison lines of Upper Atmosphere, airglow (residual continuum), and Molecular Emission of Lower Atmosphere (the latter essentially zero below 1.7 micron though).
% The files were produced by Stefan Kimensberger (sic) "the ESO sky model is a product of a theoretical calculation cross calibrated with thousands of spectra. The files here are based on online version release 2.0.6, based on the papers Noll et al. (2012, A&A 543, A92) and Jones et al. (2013, A&A 560, A91). This version does not yet include the most recent minor revision of theAerosol extinction curve by Jones et al 2019.
% See also: https://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC
  if ~isfield( opt, 'ground_telescope' )
  opt.grnd_tlscp.do = 0 ; % By default, consider a space telescope

  else
    if ~isfield( opt.ground_telescope, 'do' )
    opt.grnd_tlscp.do = 0 ; % By default, consider a space telescope
    else
    opt.grnd_tlscp.do = opt.ground_telescope.do ;
    end
  % Directory where the sky model are found
    if ~isfield( opt.ground_telescope, 'dir' )
    opt.grnd_tlscp.dr = sprintf( '%sground_telescope/', opt.scn_dr ) ;
    else
    opt.grnd_tlscp.dr = opt.ground_telescope.dir ;
    end
  % Creating the directory if it does not exist
    if ~isdir( opt.grnd_tlscp.dr )
    mkdir( opt.grnd_tlscp.dr ) ;
    end
  % For now, one of these three choices: 'full', 'half', 'new'.
    if ~isfield( opt.ground_telescope, 'moon_phase' )
    opt.grnd_tlscp.mn_phs = 'half' ; % By default, half moon
    else
    opt.grnd_tlscp.mn_phs = opt.ground_telescope.moon_phase ;
    end
  % Check of consistency  
    if ~strcmp( lower( opt.grnd_tlscp.mn_phs ), 'full' ) && ~strcmp( lower( opt.grnd_tlscp.mn_phs ), 'half' ) && ~strcmp( lower( opt.grnd_tlscp.mn_phs ), 'new' )
    disp( sprintf( 'ERROR: please, choose one of the three options: full, half, or new, for the moon phase. You chose %s. Stopped', opt.grnd_tlscp.mn_phs ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  % Turning off the local zodiacal light from SISTER
    if ( opt.lcl_zd.do ) && ( opt.grnd_tlscp.do )
    disp( 'WARNING: local zodical light comes from the sky model for the ground telescope *only*.' )
    opt.lcl_zd.do = 0 ;
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-defined star systems (search by name) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = get_star_properties( opt )

% First check if it is in ExoCat
opt = list_of_targets_exocat( opt ) ;

% List of targets from Exocat
function opt = list_of_targets_exocat( opt )

% Nota Bene: Matlab cannot read binary Excel sheets unless in Windows OS. So, this is ExoCat.xlsb, then saved as xlsx.
[ num, txt, raw ] = xlsread( [ opt.scn_dr 'ExoCat1.xlsx' ] ) ; %import data from ExoCat
% IDs must be common names or Hipparcos ID #s
% Examples: 'HIP 19849', 'HIP 30920'
opt.str.nm = lower( opt.str.nm ) ;
txt = lower( txt ) ;
s1 = 'hip' ;
target = lower( opt.str.nm ) ;

%find data for HIP ID
isHIP = strncmpi( s1, target, 3 ) ;
  if isHIP == 1
  HIPID = target( findstr( s1, target ) + 3 : end ) ;
  id_num = str2double( HIPID ) ;
  ind = find( num( :,1 ) == id_num ) ;
  notfound = isempty( ind ) ;
    if notfound == 1
    disp( sprintf( 'WARNING (if using ExoCat): the star called ''%s'' was not found.', target ) )
    return
    end
  [ row, col ] = ind2sub( size( num ),ind ) ;
  opt.str.ra_dg = num( row, 12 ) ;
  opt.star.ra_deg = opt.str.ra_dg ;
  opt.str.dc_dg = num( row, 13 ) ;
  opt.star.dec_deg = opt.str.dc_dg ;
  opt.str.pm_ra_mas_yr = num( row, 14 ) ;
  opt.star.pm_ra_mas_yr = opt.str.pm_ra_mas_yr ;
  opt.str.pm_dc_mas_yr = num( row, 15 ) ;
  opt.star.pm_dec_mas_yr = opt.str.pm_dc_mas_yr ;
  opt.str.dst_pc = num( row, 18 ) ;
  opt.star.distance_to_earth_pc = opt.str.dst_pc ;
  % Star type to be identified with one of the available spectra in the code: 'a0v', 'a5v', 'f5v', 'g0v', 'g5v' (closer to sun), 'k0v', 'k5v', 'm0v', 'm5v'
  % PS: + 1, because "num" skips the first row of the Excel sheet, but "txt" does not
  str_tp = txt( row + 1, 35 ) ;
  opt.str.tp = lower( str_tp ) ;
  opt.star.type = opt.str.tp ;
  % Star Johnson V band magnitude 
  opt.str.mg = num( row, 20 ) ;
  opt.star.app_mag_v = opt.str.mg ;
  % Star mass
  opt.str.mss_sn = num( row, 46 ) ;
  opt.star.times_sun_mass = opt.str.mss_sn ;
  % Bolometric luminosity compared to the Sun
  opt.str.lmn_sn = num( row, 33 ) ;
  %find data for common names
  else
  row = find( strcmp( txt( :, 5 ), target ) ) ;
  notfound = isempty( row ) ;
    if notfound == 1
    disp( sprintf( 'Info: the star called ''%s'' was not found in ExoCat. Unless using ExoCat, disregard this message.', target ) )
    return
    end
  % The variable "num" skips the first row, whereas "txt" does not.
  row = row - 1 ;
  opt.str.ra_dg = num( row, 12 ) ;
  opt.star.ra_deg = opt.str.ra_dg ;
  opt.str.dc_dg = num( row, 13 ) ;
  opt.star.dec_deg = opt.str.dc_dg ;
  opt.str.pm_ra_mas_yr = num( row, 14 ) ;
  opt.star.pm_ra_mas_yr = opt.str.pm_ra_mas_yr ;
  opt.str.pm_dc_mas_yr = num( row, 15 ) ;
  opt.star.pm_dec_mas_yr = opt.str.pm_dc_mas_yr ;
  opt.str.dst_pc = num( row, 18 ) ;
  opt.star.distance_to_earth_pc = opt.str.dst_pc ;
  % Star type to be identified with one of the available spectra in the code: 'a0v', 'a5v', 'f5v', 'g0v', 'g5v' (closer to sun), 'k0v', 'k5v', 'm0v', 'm5v'
  % +1 to bring "row" back to its original value above
  str_tp = txt( row + 1, 35 ) ;
  opt.str.tp = lower( str_tp ) ;
  opt.star.type = opt.str.tp ;
  % Star V magnitude
  opt.str.mg = num( row, 20 ) ;
  opt.star.app_mag_v = opt.str.mg ;
  % Star mass
  opt.str.mss_sn = num( row, 46 ) ;
  opt.star.times_sun_mass = opt.str.mss_sn ;
  % Bolometric luminosity compared to the Sun
  opt.str.lmn_sn = num( row, 33 ) ;
  end

  if iscell( opt.str.tp )
  opt.str.tp = opt.str.tp{ 1 } ;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Planet flux ratio without using Keplerian orbits %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = get_planet_flux_ratio( opt )

% Finding the value of the phase function: Analytical or open a file
  for i_pl = 1 : opt.n_pl
  tmp = NaN ;
  fl_nm_tmp = lower( opt.plnt.add.phs_fnctn{ i_pl } ) ;
    if iscell( fl_nm_tmp )
    fl_nm_tmp = fl_nm_tmp{ 1 } ;
    end
    if strcmp( fl_nm_tmp( 1 : 7 ), 'lambert' )
    tmp = ( sin( opt.plnt.add.phs_ang_rd( i_pl ) ) + ( pi - opt.plnt.add.phs_ang_rd( i_pl ) ) .* cos( opt.plnt.add.phs_ang_rd( i_pl ) ) ) / pi ;
    elseif strcmp( fl_nm_tmp, 'rayleigh_w1p0' )
    % Load the data file: it brings phase_angle_deg and phase_function_value
    load( [ opt.scn_dr 'geomAlbedo/rayleigh_scalar_w1p0_phase_function.mat' ] )
    % Interpolating to the current values
    tmp = interp1( phase_angle_deg, phase_function_value, opt.plnt.add.phs_ang_dg( i_pl ) ) ;
    elseif strcmp( fl_nm_tmp, 'rayleigh_w0p3' )
    % Load the data file: it brings phase_angle_deg and phase_function_value
    load( [ opt.scn_dr 'geomAlbedo/rayleigh_scalar_w0p3_phase_function.mat' ] )
    % Interpolating to the current values
    tmp = interp1( phase_angle_deg, phase_function_value, opt.plnt.add.phs_ang_dg( i_pl ) ) ;
    end

  % Check
    if isnan( tmp )
    disp( sprintf( 'The phase function for planet %i (%s) could not be identified. Fix it. Stopped.', i_pl, opt.plnt.add.phs_fnctn{ i_pl } ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  opt.plnt.add.phs_fnctn_vl( i_pl ) = tmp ;
  
  % Geometric albedo
  % This must be fine if the geometric albedo comes from a planet type.
    if ~iscell( opt.plnt.add.gm_albd )
    disp( 'WARNING: the geometric albedo should be a cell array. Stopped.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  % It brings lambda_nm_array and geomAlbedo_array
  gm_albd_fl_nm = sprintf( '%sgeomAlbedo/%s_geomAlbedo.mat', opt.scn_dr, opt.kplr.gm_albd{ i_pl } ) ;
    if exist( gm_albd_fl_nm, 'file' ) ~= 2
    disp( sprintf( 'WARNING: the file %s could not be found. Please, verify the geometric albedo is the one expected. Stopped.', gm_albd_fl_nm ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  load( gm_albd_fl_nm )
  % Interpolating to the current wavelength under consideration (linear since it may have local, steep features)
    for i_bnd = 1 : opt.n_bnd
    % Average value over the wavelength slice
      for i_lmbd = 1 : numel( opt.lmbd_arry_img_nm )
      % Fine array of wavelengths every nm
        if ( i_lmbd < numel( opt.lmbd_arry_img_nm ) )
        lmbd_intrp_arry = opt.lmbd_arry_img_nm( i_lmbd ) : 1 : opt.lmbd_arry_img_nm( i_lmbd + 1 ) ;
        else % Next wavelength is the limit of the band
        lmbd_intrp_arry = opt.lmbd_arry_img_nm( i_lmbd ) : 1 : opt.lmbd_img_2_nm ;
        end
      % Average value
      gm_albd_arry_tmp( i_lmbd ) = sum( interp1( lambda_nm_array, geomAlbedo_array, lmbd_intrp_arry, 'linear' ) ) / numel( lmbd_intrp_arry ) ;
      end
    opt.plnt.add.gm_albd_arry( i_pl, i_bnd, 1 : numel( gm_albd_arry_tmp ) ) = gm_albd_arry_tmp ;

    % In order to use the options when a phase angle is given for the planets, instead of using Kepler, we consider two options:
    % i) the user provides the phase angle, the planet type, and the actual distance between the planet and the host star (not the SMA), or
    % ii) the user provides the phase angle, the planet type, and the two apparent positions on the image. Notice that this may yield actual distances between planet and host star that are very far away.
    % First case:
      if isfield( opt.plnt.add, 'r_orb_au' ) && ( numel( opt.plnt.add.r_orb_au ) == opt.n_pl )
      % Placing the planet in one of the axis
      % DEC = r * (c_Omega * c_omega - s_Omega * s_omega * c_inc ) = 0 --> cotan_Omega = tan_omega * c_inc
      % RA = r * ( s_Omega * c_omega + c_Omega * s_omega * c_inc ) = r * s_Omega * ( c_omega + s_omega^2 / c_omega * c_inc ) = r * s_Omega / c_omega * ( c_omega^2 + s_omega^2 * c_inc )
      % cos( phase_angle ) = s_omega * s_inc. Also, sin( phase_angle ) = sqrt( RA^2 + DEC^2 ) / r = RA / r (assume RA>0)--> RA = r * sin( phase_angle ). We do not need to find the actual values for Omega, omega and inc. One possible solution is the following: Omega=0, omega=pi/2, inc=pi/2-phase_angle_rad ( s_omega * s_inc / c_Omega = s_phase_angle, s_omega * s_inc = c_phase_angle)
      opt.plnt.add.arc_dc_mas( i_pl ) = 0 * opt.au2mas ; 
      opt.plnt.add.arc_ra_mas( i_pl ) = opt.plnt.add.r_orb_au( i_pl ) * opt.au2mas * sin( opt.plnt.add.phs_ang_rd( i_pl ) ) ;
      % Conversion into pixels
      opt.plnt.add.arc_dc_px( i_pl ) = round( opt.plnt.add.arc_dc_mas( i_pl ) / opt.px_scn_mas ) ;
      opt.plnt.add.arc_ra_px( i_pl ) = round( opt.plnt.add.arc_ra_mas( i_pl ) / opt.px_scn_mas ) ;
      else
      % Second case:
      % Setting argument of the perisatron equal to pi/2 to get the actual distance planet-star 
      % Getting the inclination
      opt.plnt.add.incl_dg( i_pl ) = 90 - opt.plnt.add.phs_ang_dg( i_pl ) ; % cos(alpha)=sin(inc)*sin(Omega=pi/2)
      opt.plnt.add.r_orb_au( i_pl ) = sqrt( ( opt.plnt.add.arc_ra_mas( i_pl )^2 + opt.plnt.add.arc_dc_mas( i_pl )^2 ) / ( cos( opt.plnt.add.incl_dg( i_pl ) * pi / 180 ) )^2 ) / opt.au2mas ;
      end
    opt.plnt.add.r_orb_rd_rt_2( i_pl ) = ( opt.plnt.add.rd_km( i_pl ) / ( opt.plnt.add.r_orb_au( i_pl ) * opt.au2km ) )^2 ;
    flx_rt_arry_tmp = opt.plnt.add.r_orb_rd_rt_2( i_pl ) * opt.plnt.add.gm_albd_arry( i_pl, : ) * opt.plnt.add.phs_fnctn_vl( i_pl ) ;
    opt.plnt.add.flx_rt_arry( i_pl, i_bnd, 1 : numel( flx_rt_arry_tmp ) ) = flx_rt_arry_tmp ;
    end % i_bnd
  end % i_pl

  if ( opt.n_pl > 1 )
  txt_plnt = '* The distances of the planets to the star are' ;
    for i_pl = 1 : opt.n_pl - 1
    txt_plnt = sprintf( '%s %0.2f, ', txt_plnt, opt.plnt.add.r_orb_au( i_pl ) ) ;
    end
  txt_plnt = sprintf( '%sand %0.2f AU.', txt_plnt, opt.plnt.add.r_orb_au( end ) ) ;
  else
  txt_plnt = sprintf( '* The distance of the planet to the star is %0.2f AU', opt.plnt.add.r_orb_au ) ;
  end
disp( txt_plnt ) 

%%%%%%%%%%%%%%%%%%%%
% Keplerian orbits %
%%%%%%%%%%%%%%%%%%%%
function opt  = get_kepler_options( opt )

% If an exo-planetary system is used
opt = get_kepler_system( opt ) ;

% Get planetary orbital and flux data
opt = get_planet_data( opt ) ;

% Timing. Array of times that will be transformed into an array of orbital phase and finally into an array of phase function values

  if ~isfield( opt.kepler, 'time_year' )
  opt.kepler.time_year.initial = 2025 ;
  opt.kepler.time_year.final = 2030 ;
  opt.kepler.time_year.interval = 0.98765 ; % Avoiding 1 year not to plot the Earth again at the same location in the enxt image
  else
    if ~isfield( opt.kepler.time_year, 'initial' )
    opt.kepler.time_year.initial = 2025 ;
    end
    if ~isfield( opt.kepler.time_year, 'final' )
    opt.kepler.time_year.final = 2030 ;
    end
    if ~isfield( opt.kepler.time_year, 'interval' )
    opt.kepler.time_year.interval = 0.98765 ; % Avoiding 1 year not to plot the Earth again at the same location in the enxt image
    end
  end

tm_yr_arry = opt.kepler.time_year.initial : opt.kepler.time_year.interval : opt.kepler.time_year.final ;
% Storing the times of observation
opt.kplr.tm_yr_obs = tm_yr_arry ;
% Referred to J2000
tm_yr_arry = tm_yr_arry - 2000 ;
% updating the number of loops
opt.n_lp = numel( tm_yr_arry ) ;

% Unless set earlier the loop should be complete
  if ~isfield( opt, 'i_lp_1' )
  opt.i_lp_1 = 1 ;
  end
  if ~isfield( opt, 'i_lp_2' )
  opt.i_lp_2 = opt.n_lp ;
  end

% Finding the planet position 
opt = get_planet_position( opt, tm_yr_arry ) ;

% Finding the phase angles and phase function values 
opt = get_phase_function( opt, tm_yr_arry ) ;

% Get the planet flux ratio for each observation time
opt = get_kepler_flux_ratio( opt ) ;

% Compute the RA/DEC change of position due to proper motion and parallax
opt = get_radec_pos_target( opt ) ;

% sv_out_2 & redo_2
% Keeping track of the intentions of the user about storing single wavelength results
opt.save_single_wavelength_user = opt.sv_out_2 ;
% Storing the single wavelength results for the Kepler case (allows one to only focus on convolving the planets and not the whole scene every time, unless the background requires to convolve it again, see next)
opt.save_single_wavelength = 1 ;

% If there's no background, no need to repeat the non-planetary convolution 
opt.redo_2_arry = zeros( 1, opt.n_lp ) ;
opt.redo_2_arry( 1 ) = 1 ;
% If there's background, modify for parallax and proper motion (ra_wpax: ra + proper motion + parallax) 
  if ( opt.bckgrnd.add )
  % Relative deviation in pixels. Notice the relative change of sign, since the star is held at the center of the image.
  opt.bckgrnd.dlt_ra_mas = -( opt.str.ra_wpax_dg - 0.5 * ( opt.str.ra_wpax_dg( end ) + opt.str.ra_wpax_dg( 1 ) ) ) .* cos( opt.str.dc_wpax_dg / 180 * pi ) * 3600 * 1e3 ;
  opt.bckgrnd.dlt_dc_mas = -( opt.str.dc_wpax_dg - 0.5 * ( opt.str.dc_wpax_dg( end ) + opt.str.dc_wpax_dg( 1 ) ) ) * 3600 * 1e3 ;
  dlt_px = sqrt( opt.bckgrnd.dlt_ra_mas.^2 + opt.bckgrnd.dlt_ra_mas.^2 ) / opt.px_scn_mas ;
  % Whether there's more than a pixel shift
  dlt_dlt_px = dlt_px( 2 : end ) - dlt_px( 1 : end - 1 ) ; 
  opt.redo_2_arry( 2 : end ) = ( abs( dlt_dlt_px ) >= 1 ) ;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some particular systems %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = get_kepler_system( opt )
  
  if numel( strfind( lower( opt.kplr.sys ), 'solar' ) )
  opt.kepler.planet_type = { 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune' } ;
    if ~isfield( opt.star, 'times_sun_mass' )
    opt.str.mss_sn = 1 ; % 1 Sun's mass
    else
    opt.str.mss_sn = opt.star.times_sun_mass ; % Allowing some variation
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Specific planet types %
%%%%%%%%%%%%%%%%%%%%%%%%%
% The following compilation is only meant to be used for generating exo-planetary systems that resemble the solar system, *not* for ephemeris.
% Radii are mean volumetric radius: http://www.smartconversion.com/otherInfo/Volumetric_mean_radius_of_planets_and_the_sun.aspx (same as the ones used in the Haystacks project)
% Orbital parameters from: http://www.met.rdg.ac.uk/~ross/Astronomy/Planets.html, extras.springer.com/2009/978-3-540-88054-7/16_vi4b_422.pdf, based upon Murray, C.D., Dermott, S.F.: Solar System Dynamics, Cambridge University Press, Cambridge 1999. See also https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
% Inclination and longitude of the ascending node defined wrt the invariable plane (conservation of L defines the axis of the 2-d solution of Kepler's orbits) https://www.aanda.org/articles/aa/pdf/2012/07/aa19011-12.pdf (DE405-JPL)
function opt = get_data_from_planet_type( opt )

  for i_pl = 1 : opt.n_pl
  % Sentinel
    if ~isfield( opt.kplr, 'rd_km' )
    opt.kplr.rd_km( i_pl ) = NaN ;
    else
      if numel( opt.kplr.rd_km ) ~= opt.n_pl
      opt.kplr.rd_km( i_pl ) = NaN ;
      end
    end

    if strcmp( lower( opt.kepler.planet_type{ i_pl } ), 'venus' )
    opt.kplr.rd_km( i_pl ) = 6051.8 ;
    opt.kplr.sma_au( i_pl ) = 0.723332 ;
    opt.kplr.ecc( i_pl ) = 0.00677323 ;
    % Prograde: same sense of rotation as the host star (viewed from the star's North Pole)
    opt.kplr.prgrd( i_pl ) = 1 ;
    % Orbital inclination (rad)
    opt.kplr.incl_rd( i_pl ) = 2.1545441 / 180 * pi;
    % Longitude of the ascending node (Omega) (rad)
    opt.kplr.lng_nd_rd( i_pl ) = 52.3081499 / 180 * pi ;
    % Argument of the periastron (omega=~omega-Omega)
    opt.kplr.arg_pr_rd( i_pl ) = 131.53298 / 180 * pi - opt.kplr.lng_nd_rd( i_pl )  ;
    % Mean anomaly J2000 (=mean longitude - longitude of the peripasis)
    opt.kplr.mn_anmly_rd_j2000( i_pl ) = ( 181.97973 - 131.53298 ) / 180 * pi ;
    % To read the file of the geometric albedo
    opt.kplr.gm_albd{ i_pl } = lower( opt.kepler.planet_type{ i_pl } ) ;
    end
    if strcmp( lower( opt.kepler.planet_type{ i_pl } ), 'earth' )
    opt.kplr.rd_km( i_pl ) = 6371 ;
    opt.kplr.sma_au( i_pl ) = 1.0 ;
    opt.kplr.ecc( i_pl ) =  0.01671022 ;
    % Prograde: same sense of rotation as the host star (viewed from the star's North Pole)
    opt.kplr.prgrd( i_pl ) = 1 ;
    % Orbital inclination (rad)
    opt.kplr.incl_rd( i_pl ) = 1.5717094 / 180 * pi ;
    % Longitude of the ascending node (Omega) (rad)
    opt.kplr.lng_nd_rd( i_pl ) = 284.5053506 / 180 * pi ;
    % Argument of the periastron (omega)
    opt.kplr.arg_pr_rd( i_pl ) = -4.28988266 ;
    % Mean anomaly J2000 (=mean longitude - longitude of the peripasis)
    opt.kplr.mn_anmly_rd_j2000( i_pl ) = ( 100.46435 - 102.94719 ) / 180 * pi ;
    % To read the file of the geometric albedo
    opt.kplr.gm_albd{ i_pl } = lower( opt.kepler.planet_type{ i_pl } ) ;
    end
    if strcmp( lower( opt.kepler.planet_type{ i_pl } ), 'mars' )
    opt.kplr.rd_km( i_pl ) = 3389.5 ;
    opt.kplr.sma_au( i_pl ) = 1.523662 ;
    opt.kplr.ecc( i_pl ) =  0.09341233 ;
    % Prograde: same sense of rotation as the host star (viewed from the star's North Pole)
    opt.kplr.prgrd( i_pl ) = 1 ;
    % Orbital inclination (rad)
    opt.kplr.incl_rd( i_pl ) = 1.6311858 / 180 * pi ;
    % Longitude of the ascending node (Omega) (rad)
    opt.kplr.lng_nd_rd( i_pl ) = 352.9528964 / 180 * pi ;
    % Argument of the periastron (omega)
    opt.kplr.arg_pr_rd( i_pl ) = 4.99971008  ;
    % Mean anomaly J2000 (=mean longitude - longitude of the peripasis)
    opt.kplr.mn_anmly_rd_j2000( i_pl ) = ( 355.45332 - 336.04084 ) / 180 * pi ;
    % To read the file of the geometric albedo
    opt.kplr.gm_albd{ i_pl } = lower( opt.kepler.planet_type{ i_pl } ) ;
    end
    if strcmp( lower( opt.kepler.planet_type{ i_pl } ), 'jupiter' )
    opt.kplr.rd_km( i_pl ) = 69911 ;
    opt.kplr.sma_au( i_pl ) = 5.203363 ;
    opt.kplr.ecc( i_pl ) = 0.04839266 ;
    % Prograde: same sense of rotation as the host star (viewed from the star's North Pole)
    opt.kplr.prgrd( i_pl ) = 1 ;
    % Orbital inclination (rad)
    opt.kplr.incl_rd( i_pl ) = 0.3219652 / 180 * pi ;
    % Longitude of the ascending node (Omega) (rad)
    opt.kplr.lng_nd_rd( i_pl ) = 306.9167004 / 180 * pi ;
    % Argument of the periastron (omega)
    opt.kplr.arg_pr_rd( i_pl ) = -1.49753261  ;
    % Mean anomaly J2000 (=mean longitude - longitude of the peripasis)
    opt.kplr.mn_anmly_rd_j2000( i_pl ) = ( 34.40438 - 14.75385 ) / 180 * pi ;
    % To read the file of the geometric albedo
    opt.kplr.gm_albd{ i_pl } = lower( opt.kepler.planet_type{ i_pl } ) ;
    end
    if strcmp( lower( opt.kepler.planet_type{ i_pl } ), 'saturn' )
    opt.kplr.rd_km( i_pl ) = 58232 ;
    opt.kplr.sma_au( i_pl ) = 9.537070 ;
    opt.kplr.ecc( i_pl ) = 0.0541506 ;
    % Prograde: same sense of rotation as the host star (viewed from the star's North Pole)
    opt.kplr.prgrd( i_pl ) = 1 ;
    % Orbital inclination (rad)
    opt.kplr.incl_rd( i_pl ) = 0.9254704 / 180 * pi ;
    % Longitude of the ascending node (Omega) (rad)
    opt.kplr.lng_nd_rd( i_pl ) = 122.2651836 / 180 * pi ;
    % Argument of the periastron (omega)
    opt.kplr.arg_pr_rd( i_pl ) = -0.3714602 ;
    % Mean anomaly J2000 (=mean longitude - longitude of the peripasis)
    opt.kplr.mn_anmly_rd_j2000( i_pl ) = ( 49.94432 - 92.43194 ) / 180 * pi ;
    % To read the file of the geometric albedo
    opt.kplr.gm_albd{ i_pl } = lower( opt.kepler.planet_type{ i_pl } ) ;
    end
    if strcmp( lower( opt.kepler.planet_type{ i_pl } ), 'uranus' )
    opt.kplr.rd_km( i_pl ) = 25362 ;
    opt.kplr.sma_au( i_pl ) = 19.191263 ;
    opt.kplr.ecc( i_pl ) = 0.0471677 ;
    % Prograde: same sense of rotation as the host star (viewed from the star's North Pole)
    opt.kplr.prgrd( i_pl ) = 1 ;
    % Orbital inclination (rad)
    opt.kplr.incl_rd( i_pl ) = 0.9946692 / 180 * pi ;
    % Longitude of the ascending node (Omega) (rad)
    opt.kplr.lng_nd_rd( i_pl ) = 308.4427476 / 180 * pi ;
    % Argument of the periastron (omega)
    opt.kplr.arg_pr_rd( i_pl ) = 1.68833303 ;
    % Mean anomaly J2000 (=mean longitude - longitude of the peripasis)
    opt.kplr.mn_anmly_rd_j2000( i_pl ) = ( 313.23218 - 170.96424 ) / 180 * pi ;
    % To read the file of the geometric albedo
    opt.kplr.gm_albd{ i_pl } = lower( opt.kepler.planet_type{ i_pl } ) ;
    end
    if strcmp( lower( opt.kepler.planet_type{ i_pl } ), 'neptune' )
    opt.kplr.rd_km( i_pl ) = 24622 ;
    opt.kplr.sma_au( i_pl ) = 30.068964 ;
    opt.kplr.ecc( i_pl ) = 0.00858587 ;
    % Prograde: same sense of rotation as the host star (viewed from the star's North Pole)
    opt.kplr.prgrd( i_pl ) = 1 ;
    % Orbital inclination (rad)
    opt.kplr.incl_rd( i_pl ) = 0.7354155 / 180 * pi ;
    % Longitude of the ascending node (Omega) (rad)
    opt.kplr.lng_nd_rd( i_pl ) = 189.2848872/ 180 * pi ;
    % Argument of the periastron (omega)
    opt.kplr.arg_pr_rd( i_pl ) = -1.51407897 ;
    % Mean anomaly J2000 (=mean longitude - longitude of the peripasis)
    opt.kplr.mn_anmly_rd_j2000( i_pl ) = ( 304.88003 - 44.97135 ) / 180 * pi ;
    % To read the file of the geometric albedo
    opt.kplr.gm_albd{ i_pl } = lower( opt.kepler.planet_type{ i_pl } ) ;
    end
    if strcmp( lower( opt.kepler.planet_type{ i_pl } ), 'pluto' )
    opt.kplr.rd_km( i_pl ) = 1195 ;
    opt.kplr.sma_au( i_pl ) = 39.48168677 ;
    opt.kplr.ecc( i_pl ) = 0.24880766 ;
    % Prograde: same sense of rotation as the host star (viewed from the star's North Pole)
    opt.kplr.prgrd( i_pl ) = 1 ;
    % Orbital inclination (rad)
    opt.kplr.incl_rd( i_pl ) = 15.5541473 / 180 * pi ;
    % Longitude of the ascending node (Omega) (rad)
    opt.kplr.lng_nd_rd( i_pl ) = 107.0600642 / 180 * pi ;
    % Argument of the periastron (omega)
    opt.kplr.arg_pr_rd( i_pl ) = 112.59714 / 180 * pi ;
    % Mean anomaly J2000 (=mean longitude - longitude of the peripasis)
    opt.kplr.mn_anmly_rd_j2000( i_pl ) = ( 238.92881 - 224.06676 ) / 180 * pi ;
    % To read the file of the geometric albedo
    opt.kplr.gm_albd{ i_pl } = lower( opt.kepler.planet_type{ i_pl } ) ;
    end

  % If this option is used, and a SMA is not specified (see next), then the SMA is rescaled according to the square root of the luminosity of the star (keeping habitable zone rescaled).
    if isfield( opt.str, 'lmn_sn' )
    opt.kplr.sma_au( i_pl ) = opt.kplr.sma_au( i_pl ) * sqrt( opt.str.lmn_sn ) ;
    end
   end % i_pl

% Function that defines each type of planet (Data defined wrt J2000)
function opt = get_planet_data( opt )

% Finding the number of planets. One of these should be set
n_pl = 0 ;
  if isfield( opt.kepler, 'sma' ) 
  n_pl_1 = numel( opt.kepler.sma ) ;
  end
  if isfield( opt.kepler, 'sma_au' )
  n_pl_2 = numel( opt.kepler.sma_au ) ;
  end  
  if isfield( opt.kepler, 'planet_type' )
  n_pl_3 = numel( opt.kepler.planet_type ) ;
  end
  if exist( 'n_pl_1', 'var' ) 
  n_pl = n_pl_1 ;
    if exist( 'n_pl_2', 'var' )
      if n_pl_1 ~= n_pl_2
      disp( sprintf( 'The number pf planets from opt.kepler.sma (%i) and opt.kepler.sma_au (%i) is inconsistent. Fix it. Stopped', n_pl_1, n_pl_2 ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end
    if exist( 'n_pl_3', 'var' )
      if n_pl_1 ~= n_pl_3
      disp( sprintf( 'The number pf planets from opt.kepler.sma (%i) and opt.kepler.planet_type (%i) is inconsistent. Fix it. Stopped', n_pl_1, n_pl_3 ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end
  end
  if exist( 'n_pl_2', 'var' )
  n_pl = n_pl_2 ;
    if exist( 'n_pl_3', 'var' )
      if n_pl_2 ~= n_pl_3
      disp( sprintf( 'The number pf planets from opt.kepler.sma_au (%i) and opt.kepler.planet_type (%i) is inconsistent. Fix it. Stopped', n_pl_2, n_pl_3 ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
    end
  end
  if exist( 'n_pl_3', 'var' )
  n_pl = n_pl_3 ;
  end

  if ~n_pl
  disp( 'Either opt.kepler.sma, sma_au or planet_type should be set. Fix it. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end
opt.n_pl = n_pl ;

  if isfield( opt.kepler, 'planet_type' )
  opt = get_data_from_planet_type( opt ) ;
  end % if planet_type is a field of opt.kepler

% Replacing the above default values by any user supplied ones (e.g., moving Venus elsewhere)
  if isfield( opt.kepler, 'sma' )
  n_pl = numel( opt.kepler.sma ) ;
    for i_pl = 1 : n_pl
    fct_au = NaN ;
    str_tmp = opt.kepler.sma{ i_pl } ;
      if ~isnan( str_tmp )
      str_tmp = lower( str_tmp ) ;
        if strfind( str_tmp, 'au' )
        fct_au = 1 ;
        q_str = findstr( str_tmp, 'au' ) ;
        opt.kepler.sma_au( i_pl ) = str2num( str_tmp( 1 : q_str - 1 ) ) ;
        end
        if strfind( str_tmp, 'km' )
        fct_au = 1 / opt.au2km ;
        q_str = findstr( str_tmp, 'km' ) ;
        opt.kepler.sma_au( i_pl ) = fct_au * str2num( str_tmp( 1 : q_str - 1 ) ) ;
        end
        if strfind( str_tmp, 'mas' )
        fct_au = opt.mas2au ;
        q_str = findstr( str_tmp, 'mas' ) ;
        opt.kepler.sma_au( i_pl ) = fct_au * str2num( str_tmp( 1 : q_str - 1 ) ) ;
        end
      % Check of consistency (if opt.kepler.sma_au( i_pl ) is NaN is fine, it means there's a default value in get_kepler_system to be used)
        if isnan( fct_au ) && ~isnan( str_tmp )
        disp( sprintf( 'Unidentified units for planet %i in opt.kepler.sma. Stopped.', i_pl ) )
        disp( ' ' )
        disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
        end
      else % ~isnan( str_tmp )
        if isfield( opt.kplr, 'sma_au' )
          if ~isnan( opt.kplr.sma_au( i_pl ) )
          opt.kepler.sma_au( i_pl ) = opt.kplr.sma_au( i_pl ) ;
          end
        end
      end
    end % i_pl
  end % If SMA is introduced

  % Planet radii
  if isfield( opt.kepler, 'radius' )
    for i_pl = 1 : n_pl
    fct_km = NaN ;
    str_tmp = opt.kepler.radius{ i_pl } ;

      if ~isnan( str_tmp )
      str_tmp = lower( str_tmp ) ;
        if strfind( str_tmp, 'km' )
        q_str = findstr( str_tmp, 'km' ) ;
        opt.kepler.radius_km( i_pl ) = str2num( str_tmp( 1 : q_str - 1 ) ) ;
        end
        if numel( strfind( str_tmp, 'r_e' ) ) || numel( strfind( str_tmp, 're' ) )
        fct_km = 6371 ; % mean radius
        q_str = findstr( str_tmp, 'r' ) ;
        opt.kepler.radius_km( i_pl ) = fct_km * str2num( str_tmp( 1 : q_str - 1 ) ) ;
        end
        if numel( strfind( str_tmp, 'r_j' ) ) || numel( strfind( str_tmp, 'rj' ) )
        fct_km = 71492 ;
        q_str = findstr( str_tmp, 'r' ) ;
        opt.kepler.radius_km( i_pl ) = fct_km * str2num( str_tmp( 1 : q_str - 1 ) ) ;
        end
      % Check of consistency
        if isnan( fct_km )
        disp( sprintf( 'Unidentified units for planet %i in opt.kepler.radius. Stopped.', i_pl ) )
        disp( ' ' )
        disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
        end
      end % ~isnan( str_tmp )
    end % i_pl
  end % If radius is introduced

  if ~isfield( opt.kepler, 'radius_km' )
    if ~isfield( opt.kplr, 'rd_km' )
    opt.kplr.rd_km = 6371 * ones( 1, opt.n_pl ) ; % Default 1 Earth mean radius
    else
    opt.kepler.radius_km = opt.kplr.rd_km ;
    end
  else
    if numel( opt.kepler.radius_km ) == 1
    opt.kplr.rd_km = opt.kepler.radius_km * ones( 1, n_pl ) ;
    else
      if numel( opt.kepler.radius_km ) < n_pl
      disp( sprintf( 'The number of planets is %i, but opt.kepler.radius_km has only %i elements. Please, fix it. Stopped.', n_pl, numel( opt.kepler.radius_km ) ) ) ;
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      end
      if numel( opt.kepler.radius_km ) > n_pl
      disp( sprintf( 'WARNING: the number of planets is %i, but opt.kepler.radius_km has more elements (%i). Only the first %i values will be used. Please, make sure this is intentional. Continuing.', n_pl, numel( opt.kepler.radius_km ), n_pl ) ) ;
      end
      for i_pl = 1 : n_pl
        if ~isnan( opt.kepler.radius_km( i_pl ) )
        opt.kplr.rd_km( i_pl ) = opt.kepler.radius_km( i_pl ) ;
        end
      end
    end
  end

  if isfield( opt.kepler, 'sma_au' )
    for i_pl =1 : n_pl
      if numel( opt.kepler.sma_au ) == 1
        if ~isnan( opt.kepler.sma_au )
        opt.kplr.sma_au( i_pl ) = opt.kepler.sma_au ;
        end
      elseif ~isnan( opt.kepler.sma_au( i_pl ) )
      opt.kplr.sma_au( i_pl ) = opt.kepler.sma_au( i_pl ) ;
      end
    end
  end

% Period of the planet (Earth = 1 year). T = 2pi sqrt( a^3/(GM_reduced)), M_reduced~M_star for exoplanets.
  for i_pl = 1 : n_pl
  opt.kplr.prd_yr( i_pl ) = sqrt( opt.kplr.sma_au( i_pl )^3 / opt.str.mss_sn ) ;
  end

  if isfield( opt.kepler, 'eccentricity' )
    for i_pl =1 : n_pl
      if numel( opt.kepler.eccentricity ) == 1
        if ~isnan( opt.kepler.eccentricity )
          if ( abs( opt.kepler.eccentricity ) > 1 )
          disp( sprintf( 'The eccentricity was set to %f and it should be between 0 and 1. Fix it.', opt.kepler.eccentricity ) )
          disp( ' ' )
          disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
          end
        opt.kplr.ecc( i_pl ) = opt.kepler.eccentricity ;
        end
      elseif ~isnan( opt.kepler.eccentricity( i_pl ) )
        if ( abs( opt.kepler.eccentricity( i_pl ) ) > 1 )
        disp( sprintf( 'The eccentricity for planet %i was set to %f and it should be between 0 and 1. Fix it.', i_pl, opt.kepler.eccentricity( i_pl ) ) )
        disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
        disp( ' ' )
        end
      opt.kplr.ecc( i_pl ) = opt.kepler.eccentricity( i_pl ) ;
      end
    end
  end
  if ~isfield( opt.kplr, 'ecc' )
    for i_pl =1 : n_pl
    opt.kplr.ecc( i_pl ) = 0 ;
    end
  end

% Whether the mean anomaly of the planet increases in time, or the reverse. All planets in the solar system have a prograde orbit.
  if isfield( opt.kepler, 'prograde' )
    for i_pl =1 : n_pl
      if numel( opt.kepler.prograde ) == 1
        if ~isnan( opt.kepler.prograde )
        opt.kplr.prgrd( i_pl ) = opt.kepler.prograde ;
        end
      elseif ~isnan( opt.kepler.prograde( i_pl ) )
      opt.kplr.prgrd( i_pl ) = opt.kepler.prograde( i_pl ) ;
      end
    end
  end
  if ~isfield( opt.kplr, 'prgrd' )
    for i_pl =1 : n_pl
    opt.kplr.prgrd( i_pl ) = 1 ;
    end
  end

% Inclination with respect the invariant plane (or some intrinsic reference plane of the extra-solar system)
  if isfield( opt.kepler, 'inclination_rad' )
    for i_pl = 1 : n_pl
      if numel( opt.kepler.inclination_rad ) == 1
        if ~isnan( opt.kepler.inclination_rad )
        opt.kplr.incl_rd( i_pl ) = opt.kepler.inclination_rad ;
        end
      elseif ~isnan( opt.kepler.inclination_rad( i_pl ) )
      opt.kplr.incl_rd( i_pl ) = opt.kepler.inclination_rad( i_pl ) ;
      end
    end
  end
  if ~isfield( opt.kplr, 'incl_rd' ) && ~isfield( opt.kepler, 'inclination_deg' )
    for i_pl = 1 : n_pl
    opt.kplr.incl_rd( i_pl ) = 0 ;
    end
  end
 
% Keep the inclination between [-pi/2,pi/2]
  opt.kplr.incl_rd = asin( sin( opt.kplr.incl_rd ) ) ;

% Adding the option in degrees
  if isfield( opt.kepler, 'inclination_deg' )
    for i_pl = 1 : n_pl
      if numel( opt.kepler.inclination_deg ) == 1
        if ~isnan( opt.kepler.inclination_deg )
        opt.kplr.incl_dg( i_pl ) = opt.kepler.inclination_deg ;
        opt.kplr.incl_rd( i_pl ) = opt.kepler.inclination_deg * pi / 180 ;
        end
      elseif ~isnan( opt.kepler.inclination_deg( i_pl ) )
      opt.kplr.incl_dg( i_pl ) = opt.kepler.inclination_deg( i_pl ) ;
      opt.kplr.incl_rd( i_pl ) = opt.kepler.inclination_deg( i_pl ) * pi / 180 ;
      end
    end
  end
  if ~isfield( opt.kplr, 'incl_dg' )
    for i_pl = 1 : n_pl
    opt.kplr.incl_dg( i_pl ) = opt.kplr.incl_rd( i_pl ) * 180 / pi ;
    end
  end
  % Check of consistency (to 0.001 deg)
    if isfield( opt.kplr, 'incl_rd' ) && isfield( opt.kplr, 'incl_dg' )
      for i_pl = 1 : n_pl
        if ( round( 1000 * opt.kplr.incl_rd( i_pl ) * 180 / pi ) ~= round( 1000 * opt.kplr.incl_dg( i_pl ) ) )
        disp( sprintf( 'The orbital inclination has been set both in degrees and radians, but they look different: opt.kepler.inclination_rad=%2.3f (%2.3f deg) and opt.kepler.inclination_deg=%2.3f. Set one of them only. Stopped.', opt.kplr.incl_rd( i_pl ), opt.kplr.incl_rd( i_pl ) * 180 / pi, opt.kplr.incl_dg( i_pl ) ) )
        disp( ' ' )
        disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
        end
    end
  end

  if isfield( opt.kepler, 'longitude_ascending_node_rad' )
    for i_pl =1 : n_pl
      if numel( opt.kepler.longitude_ascending_node_rad ) == 1
        if ~isnan( opt.kepler.longitude_ascending_node_rad )
        opt.kplr.lng_nd_rd( i_pl ) = opt.kepler.longitude_ascending_node_rad ;
        end
      elseif ~isnan( opt.kepler.longitude_ascending_node_rad( i_pl ) )
      opt.kplr.lng_nd_rd( i_pl ) = opt.kepler.longitude_ascending_node_rad( i_pl ) ;
      end
    end
  end
  if ~isfield( opt.kplr, 'lng_nd_rd' )
    for i_pl =1 : n_pl
    opt.kplr.lng_nd_rd( i_pl ) = pi / 2 ; % Orbit is horizontal (parallel to RA)
    end
  end

% Check: if the inclination is different, the longitude of the ascending node should be changed accordingly. 
  if ( abs( opt.kplr.incl_rd ) > 10 / 180 * pi )
  % Checking if the angular distances are too big
  lng_nd_ok = 1 ;
    for i_pl_1 = 1 : n_pl
      for i_pl_2 = i_pl_1 : n_pl
      % Somewhat arbitrary sentinel
      lng_nd_dg_sntnl = 20 ;
      dt_prdct = ( cos( opt.kplr.lng_nd_rd( i_pl_1 ) ) * cos( opt.kplr.lng_nd_rd( i_pl_2 ) ) + sin( opt.kplr.lng_nd_rd( i_pl_1 ) ) * sin( opt.kplr.lng_nd_rd( i_pl_2 ) ) ) ;
      % Some rounding-off of Matlab generates some 0 + 1e-6 i numbers 
        if ( real( acos( dt_prdct ) ) > lng_nd_dg_sntnl * pi / 180 ) 
        lng_nd_ok = 0 ;
        end
      end
    end
    if ~lng_nd_ok
    disp( sprintf( 'WARNING: some of the values of the longitude of the ascending node are different to each other by more than %3.1f deg. Although correct, the orbits may look odd. If necessary set opt.kepler.longitude_ascending_node_rad to similar values (recall the meaning of the Kepler parameters for further details).', lng_nd_dg_sntnl ) )
    end
  end 

  if isfield( opt.kepler, 'argument_periapsis_rad' )
    for i_pl =1 : n_pl
      if numel( opt.kepler.argument_periapsis_rad ) == 1
        if ~isnan( opt.kepler.argument_periapsis_rad )
        opt.kplr.arg_pr_rd( i_pl ) = opt.kepler.argument_periapsis_rad ;
        end
      elseif ~isnan( opt.kepler.argument_periapsis_rad( i_pl ) )
      opt.kplr.arg_pr_rd( i_pl ) = opt.kepler.argument_periapsis_rad( i_pl ) ;
      end
    end
  end
  if ~isfield( opt.kplr, 'arg_pr_rd' )
    for i_pl = 1 : n_pl
    opt.kplr.arg_pr_rd( i_pl ) = 0 ;
    end
  end

  % If the mean anomaly in J2000 has been set, update opt.kepler.time_periapsis_passage_year
  if isfield( opt.kplr, 'mn_anmly_rd_j2000' ) && ~isfield( opt, 'mean_anomaly_rad_j2000' )
  opt.kepler.mean_anomaly_rad_j2000 = opt.kplr.mn_anmly_rd_j2000 ;
  end
  if isfield( opt, 'mean_anomaly_rad_j2000' )
    for i_pl = 1 : n_pl
      if numel( opt.kepler.mean_anomaly_rad_j2000 ) == 1
        if ~isnan( opt.kepler.mean_anomaly_rad_j2000 )
        opt.kplr.tm_pr_pssg_yr( i_pl ) = 2000 - opt.kepler.mean_anomaly_rad_j2000( i_pl ) / 2 / pi * opt.kplr.prd_yr( i_pl ) * opt.kplr.prgrd( i_pl ) ;
        end
      elseif  ~isnan( opt.kepler.mean_anomaly_rad_j2000( i_pl ) )
      opt.kplr.tm_pr_pssg_yr( i_pl ) = 2000 - opt.kepler.mean_anomaly_rad_j2000( i_pl ) / 2 / pi * opt.kplr.prd_yr( i_pl ) * opt.kplr.prgrd( i_pl ) ;
      end
    end
  end
  if ~isfield( opt.kplr, 'mn_anmly_rd_j2000' )
    for i_pl = 1 : n_pl
    opt.kplr.mn_anmly_rd_j2000( i_pl ) = 2 * pi * ( 2000 - opt.kepler.time_year.initial ) / opt.kplr.prd_yr( i_pl ) * opt.kplr.prgrd( i_pl ) ; % Consistent with the default value for time_periapsis_passage_year: the planets start at their periapsis
    end
  end

  if isfield( opt.kepler, 'time_periapsis_passage_year' )
    for i_pl = 1 : n_pl
      if numel( opt.kepler.time_periapsis_passage_year ) == 1
        if ~isnan( opt.kepler.time_periapsis_passage_year )
        opt.kplr.tm_pr_pssg_yr( i_pl ) = opt.kepler.time_periapsis_passage_year ;
        end
      elseif ~isnan( opt.kepler.time_periapsis_passage_year( i_pl ) )
      opt.kplr.tm_pr_pssg_yr( i_pl ) = opt.kepler.time_periapsis_passage_year( i_pl ) ;
      end
    end
  end
  if ~isfield( opt.kplr, 'tm_pr_pssg_yr' )
    for i_pl = 1 : n_pl
    opt.kplr.tm_pr_pssg_yr( i_pl ) = opt.kepler.time_year.initial ; % By default, the planets start at their periapsis
    end
  end

  if isfield( opt.kepler, 'geomAlbedo' )
  % Needs to be a cell (numeric or string)
    if ( numel( opt.kepler.geomAlbedo ) > 1 ) && ~iscell( opt.kepler.geomAlbedo )
    disp( 'opt.kepler.geomAlbedo has to be a cell. Syntax is opt.kepler.geomAlbedo = { option, ... }. Fix it. Stopped.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
    for i_pl = 1 : n_pl
      if numel( opt.kepler.geomAlbedo ) == 1
        if ~isnan( opt.kepler.geomAlbedo{ 1 } )
        opt.kplr.gm_albd{ i_pl } = opt.kepler.geomAlbedo{ 1 } ;
        end
      elseif ~isnan( opt.kepler.geomAlbedo{ i_pl } )
      opt.kplr.gm_albd{ i_pl } = opt.kepler.geomAlbedo{ i_pl } ;
      end
    end
  end
  if ~isfield( opt.kplr, 'gm_albd' )
    for i_pl = 1 : n_pl
    opt.kplr.gm_albd{ i_pl } = 'Jupiter' ; 
    end
  end

% Position angle of the exo-planetary system (https://en.wikipedia.org/wiki/Position_angle). It's equivalent to the Kepler parameter line of ascending node, although this parameter is applied to all the planets.
  if ~isfield( opt.kepler, 'position_angle_deg' )
  opt.kplr.pa_dg = 0 ;
  else
  opt.kplr.pa_dg = opt.kepler.position_angle_deg ;
  end
  % Same but in radians
  if isfield( opt.kepler, 'position_angle_rad' )
    if  ~( opt.kplr.pa_dg )
    opt.kplr.pa_dg = opt.kepler.position_angle_rad * 180 / pi ;
    else
      if ( opt.kplr.pa_dg * pi / 180 == opt.kepler.position_angle_rad )
      disp( sprintf( 'WARNING: the position angle of the orbital system has been set twice to the same value. In degrees is %f, and in rad is %f. Make sure to use one only. Continuing ...', opt.kplr.pa_dg, opt.kepler.position_angle_rad ) )
      else
      disp( sprintf( 'The position angle of the orbital system has been set twice. In degrees is %f, and in rad is %f, which are not strictly the same. Choose one of them. Stopped.', opt.kplr.ps_dg, opt.kepler.position_angle_rad ) )
      make_a_Stop
      end
    end
  end
  % Comparing with the position angle set for the exozodiacal light
  if isfield( opt, 'zd' )
    if isfield( opt.zd, 'pa_dg' )
      if ( opt.zd.pa_dg ~= opt.kplr.pa_dg )
      disp( sprintf( 'WARNING: the exozodiacal emission and the orbital system have different position angles: %f deg, and %f deg, respectively. Make sure this is intentional. Continuing ...', opt.zd.pa_dg, opt.kplr.pa_dg ) )
      end
    else
      if ~( opt.kplr.pa_dg ) && ( opt.kplr.pa_dg ) % If it is zero, it does not rotate the image
      disp( sprintf( 'WARNING: the orbital system has set a position angle of %f deg, whilst the exozodiacal emission has no position angle set. Make sure this is intentional. Continuing ...', opt.kplr.pa_dg ) ) 
      end
    end
  end

% Checks of consistency
  if ( opt.n_pl ~= numel( opt.kplr.rd_km ) )
  disp( 'The number of elements of the planet radii is not the same as the one for the semi-major axis. Fix it. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end
  if ( opt.n_pl ~= numel( opt.kplr.ecc ) )
  disp( 'The number of elements of the planet eccentricity is not the same as the one for the semi-major axis. Fix it. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end
  if ( opt.n_pl ~= numel( opt.kplr.prgrd ) )
  disp( 'The number of elements of the planet sense of rotation (prograde) is not the same as the one for the semi-major axis. Fix it. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end
  if ( opt.n_pl ~= numel( opt.kplr.incl_rd ) )
  disp( 'The number of elements of the planet inclination is not the same as the one for the semi-major axis. Fix it. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end
  if ( opt.n_pl ~= numel( opt.kplr.lng_nd_rd ) )
  disp( 'The number of elements of the planet longitude of the ascending node is not the same as the one for the semi-major axis. Fix it. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end
  if ( opt.n_pl ~= numel( opt.kplr.arg_pr_rd ) )
  disp( ' ' )
  disp( 'The number of elements of the planet argument of the periapsis is not the same as the one for the semi-major axis. Fix it. Stopped' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end
  if ( opt.n_pl ~= numel( opt.kplr.mn_anmly_rd_j2000 ) )
  disp( sprintf( 'The number of elements of the planets mean anomaly (%i) is not the same as the one for the semi-major axis (%i). Fix it. Stopped', numel( opt.kplr.mn_anmly_rd_j2000 ), opt.n_pl ) )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end
  if ( opt.n_pl ~= numel( opt.kplr.tm_pr_pssg_yr ) )
  disp( sprintf( 'The number of elements of the planets time of periapsis passage (%i) is not the same as the one for the semi-major axis (%i). Fix it. Stopped', numel( opt.kplr.tm_pr_pssg_yr ), opt.n_pl ) )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end
  if ( opt.n_pl ~= numel( opt.kplr.gm_albd ) )
  disp( 'The number of elements of the planet geometric albedos is not the same as the one for the semi-major axis. Fix it. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

% Filling out opt.kepler to avoid being overwritten
opt.kepler.radius_km = opt.kplr.rd_km ;
opt.kepler.sma_au = opt.kplr.sma_au ;
opt.kepler.eccentricity = opt.kplr.ecc ;
opt.kepler.prograde = opt.kplr.prgrd ;
opt.kepler.inclination_rad = opt.kplr.incl_rd ;
opt.kepler.longitude_ascending_node_rad = opt.kplr.lng_nd_rd ;
opt.kepler.argument_periapsis_rad = opt.kplr.arg_pr_rd ;
opt.kepler.mean_anomaly_rad_j2000 = opt.kplr.mn_anmly_rd_j2000 ;
opt.kepler.time_periapsis_passage_year = opt.kplr.tm_pr_pssg_yr ;
opt.kepler.geom_albedo = opt.kplr.gm_albd ;
opt.kepler.position_angle_deg = opt.kplr.pa_dg ;
% Adding plnt_tp
opt.kplr.plnt_tp = opt.kepler.planet_type ;

%%%%%%%%%%
% Getting the eccentric anomaly
% Fast: 0.03 sec for Solar system and 61 time stamps with pixel precision
% References: http://astrowww.phys.uvic.ca/~tatum/celmechs.html, Chapter 9, Eq. 9.6.4 and following
function opt = get_planet_position( opt, tm_yr_arry )
n_pl = numel( opt.kplr.sma_au ) ;
n_tm = numel( opt.kplr.tm_yr_obs ) ;
% For each planet:
  for i_pl = 1 : n_pl
  % Getting the mean anomaly for all observation times.
  % Standard definitions( see e.g., Eq. 9.6.4 of http://astrowww.phys.uvic.ca/~tatum/celmechs/celm9.pdf, or,
  % for the calculation of the time of the perihelion passage see e.g., https://www.math.ubc.ca/~cass/courses/m309-01a/orbits.pdf
  mn_anmly_rd = opt.kplr.prgrd( i_pl ) * 2 * pi * ( ( opt.kplr.tm_yr_obs - opt.kplr.tm_pr_pssg_yr( i_pl ) ) / opt.kplr.prd_yr( i_pl ) ) ;
  % Bringing angles within [ 0, 2pi )
  n_2pi = abs( floor( min( mn_anmly_rd / 2 / pi ) ) ) ;
  mn_anmly_rd = mod( mn_anmly_rd + n_2pi * 2 * pi, 2 * pi ) ;
  % Array of values for the eccentric anomaly (precision within a pixel. Not doing Newton-Raphson. Fast enough today to evaluate the array: 7 planets, 61 observations, with all pixel precision including Neptune -30 AU- takes 0.01 sec)
  dlt_ecc_anmly_arry = atan( opt.px_scn_mas / opt.kplr.sma_au( i_pl ) / ( 1 + opt.kplr.ecc( i_pl ) ) * opt.mas2au / 2 ) ; % Divided by 2 to always have better than a pixel precision
  ecc_anmly_rd_arry_tmp = 0 : dlt_ecc_anmly_arry : ( 2 * pi - dlt_ecc_anmly_arry ) ; 
  mn_anmly_arry_tmp = ecc_anmly_rd_arry_tmp - opt.kplr.ecc( i_pl ) * sin( ecc_anmly_rd_arry_tmp ) ;
  % For each observation time
    for i_tm = 1 : n_tm
    d_tmp = abs( mn_anmly_arry_tmp - mn_anmly_rd( i_tm ) ) ;
    opt.kplr.ecc_anmly_rd_arry( i_pl, i_tm ) = ecc_anmly_rd_arry_tmp( find( d_tmp == min( d_tmp ) ) ) ;
    % For the record
    opt.kplr.mn_anmly_rd_arry( i_pl, i_tm ) = opt.kplr.ecc_anmly_rd_arry( i_pl, i_tm ) - opt.kplr.ecc( i_pl ) * sin( opt.kplr.ecc_anmly_rd_arry( i_pl, i_tm ) ) ;
    end % i_tm
  % Getting the true anomaly from the eccentric anomaly preserving quadrant information.
  % Standard definitions between eccentric anomaly and true anomaly.
  % atan2(Y,X)
  opt.kplr.tr_anmly_rd_arry( i_pl, : ) = atan2( sqrt( 1 - opt.kplr.ecc( i_pl )^2 ) * sin( opt.kplr.ecc_anmly_rd_arry( i_pl, : ) ), cos( opt.kplr.ecc_anmly_rd_arry( i_pl, : ) ) - opt.kplr.ecc( i_pl ) ) ;
  % Keep it between [0,2pi) for simplicity
  n_2pi = abs( floor( min( opt.kplr.tr_anmly_rd_arry( i_pl, : ) / 2 / pi ) ) ) ;
  opt.kplr.tr_anmly_rd_arry( i_pl, : ) = mod( opt.kplr.tr_anmly_rd_arry( i_pl, : ) + n_2pi * 2 * pi, 2 * pi ) ;
  end % i_pl

% Finding the positions of the planets on the image plane
% Following standard conventions in Celestial mechanics. See, e.g.
% R. Fitzpatrick (2016) "Introduction to Celestial Mechanics" (http://farside.ph.utexas.edu/teaching/celestial/Celestialhtml/node33.html, http://farside.ph.utexas.edu/teaching/celestial/Celestialhtml/node34.html)
% Ben Placek and Kevin H. Knuth (2014) https://arxiv.org/pdf/1310.6764.pdf 
% Convention for orbits outside the solar system: https://en.wikipedia.org/wiki/Longitude_of_the_ascending_node
  for i_pl = 1 : n_pl
  ecc_anmly_rd_tmp = squeeze( opt.kplr.ecc_anmly_rd_arry( i_pl, : ) ) ;
  % True distance to the star
  r_true = opt.kplr.sma_au( i_pl ) * ( 1 - opt.kplr.ecc( i_pl ) * cos( ecc_anmly_rd_tmp ) ) ;
  % The convention for the tangent plane on the sky that projects the true orbit is such that (see references above):
  % i) the telescope's (observer's) line of sight is along the 'Z' axis, and pointing from the star to the observer.
  % ii) the "X" axis, which is the reference direction that defines the longitude of the ascending node, points along the vertical axis of the image. Increasing values of X point towards increasing values of declination.
  % iii) the "Y" axis is along the horizontal axis of the image plane. Increasing Y is pointing from the star towards the lhs of the image.
  % Trigonometric factors (3 rotations)
  tr_anmly_rd_tmp = squeeze( opt.kplr.tr_anmly_rd_arry( i_pl, : ) ) ;
  c_Omg = cos( opt.kplr.lng_nd_rd( i_pl ) ) ; s_Omg = sin( opt.kplr.lng_nd_rd( i_pl ) ) ;
  c_omg_tm = cos( opt.kplr.arg_pr_rd( i_pl ) + tr_anmly_rd_tmp ) ; s_omg_tm = sin(  opt.kplr.arg_pr_rd( i_pl ) + tr_anmly_rd_tmp ) ; 
  c_inc = cos( opt.kplr.incl_rd( i_pl ) ) ; s_inc = sin( opt.kplr.incl_rd( i_pl ) ) ;

  X_sky = r_true .* ( c_Omg * c_omg_tm - s_Omg * s_omg_tm * c_inc ) ;
  Y_sky = r_true .* ( s_Omg * c_omg_tm + c_Omg * s_omg_tm * c_inc ) ;
  Z_sky = r_true .* s_omg_tm * s_inc ;
  % Order to follow RA/DEC convention (when plotted with imagesc and fliplr(flipud())
  tmp_arc_ra = Y_sky * opt.au2mas ; % RA (arclength, not rue RA. That is, without cos(DEC) factor)
  tmp_arc_dc = X_sky * opt.au2mas ; % DEC
  % Position angle is redundant, but sometimes it is used to deal with, for instance, a rotated camera wrt N/E.
  % Position angle, as the longitude of the ascending node, is simply a rotation about the line of sight of the observer. It does not change the illumination or the orbital motion of the planets.
  pa_rd = opt.kplr.pa_dg * pi / 180 ;
  % Rotating DEC->RA 
  %        ^                      
  %        |                     
  %      DEC(-)                   
  %        |                      
  % <-(-)--RA--(+)->		+-----> E
  %        |			|
  %      DEC(+)			N
  %        |			|
  %        v			v
  % For the case of a single scenario, the full series of positions is computed here. If there are two scenarios, it's computed the first time only.
    if ( opt.i_lp_1 == 1 )
    opt.kplr.arc_ra_mas( i_pl, : ) = tmp_arc_ra * cos( pa_rd ) + tmp_arc_dc * sin( pa_rd ) ;
    opt.kplr.arc_dc_mas( i_pl, : ) = tmp_arc_dc * cos( pa_rd ) - tmp_arc_ra * sin( pa_rd ) ;
    end
  end

%%%%%%%%
% Getting the phase angles that will be used to derive the geometric albedo
% PS: calculation of the phase angles is fast: 0.03 sec for Solar system and 61 time stamps with pixel precision
function opt = get_phase_function( opt, tm_yr_arry )

n_pl = numel( opt.kplr.sma_au ) ;
n_tm = numel( opt.kplr.tm_yr_obs ) ;
% For each planet:
  for i_pl = 1 : n_pl
    for i_tm = 1 : n_tm
    % Finding the phase angle in degrees (common choice). The phase angle is defined between [0, 180]. No need to look for quadrants.
    opt.kplr.phs_ang_dg( i_pl, i_tm ) = 180 / pi * acos( sin( opt.kplr.tr_anmly_rd_arry( i_pl, i_tm ) + opt.kplr.arg_pr_rd( i_pl ) ) * sin( opt.kplr.incl_rd( i_pl ) ) ) ;
    % Storing it in radians, as well, for other evaluations
    opt.kplr.phs_ang_rd( i_pl, i_tm ) = opt.kplr.phs_ang_dg( i_pl, i_tm ) * pi / 180 ;
    %  Finding the ratio to convert albedo into flux
    opt.kplr.sma_rd_rt_2( i_pl, i_tm ) = ( opt.kplr.rd_km( i_pl ) / ( opt.kplr.sma_au( i_pl ) * opt.au2km * ( 1 - opt.kplr.ecc( i_pl ) * cos( opt.kplr.ecc_anmly_rd_arry( i_pl, i_tm ) ) ) ) )^2 ;
    end % i_tm
  end % i_pl

% Finding the phase function values
  if ~isfield( opt.kepler, 'phase_function' )
    for i_pl = 1 : n_pl
    opt.kplr.phs_fnctn{ i_pl } = 'Lambertian' ;
    end
  else
    if numel( opt.kepler.phase_function ) == 1 || ischar( opt.kepler.phase_function ) 
      for i_pl = 1 : 7
      opt.kplr.phs_fnctn{ i_pl } = opt.kepler.phase_function ;
      end
    else
      if numel( opt.kepler.phase_function ) ~= n_pl
      disp( sprintf( 'The number of elements in the phase function option (%i) is different to the number of planets (%i)', numel( opt.kepler.phase_function ), n_pl ) )
      disp( ' ' )
      disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
      else
      opt.kplr.phs_fnctn{ i_pl } = opt.kepler.phase_function{ i_pl } ;
      end
    end
  end

% Finding the value of the phase function: Analytical or open a file
  for i_pl = 1 : n_pl
  tmp = NaN ;
  fl_nm_tmp = lower( opt.kplr.phs_fnctn{ i_pl } ) ;
    if ~iscell( fl_nm_tmp )
    disp( 'The option opt.kepler.phase_function is a cell: opt.kepler.phase_function = { option, ... }. Fix it. Stopped.' )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end   
  fl_nm_tmp = fl_nm_tmp{ 1 } ;
    if strcmp( fl_nm_tmp( 1 : 7 ), 'lambert' ) 
    tmp = ( sin( opt.kplr.phs_ang_rd( i_pl, : ) ) + ( pi - opt.kplr.phs_ang_rd( i_pl, : ) ) .* cos( opt.kplr.phs_ang_rd( i_pl, : ) ) ) / pi ;
    elseif strcmp( fl_nm_tmp, 'rayleigh_w1p0' )
    % Load the data file: it brings phase_angle_deg and phase_function_value
    load( [ opt.scn_dr 'geomAlbedo/rayleigh_scalar_w1p0_phase_function.mat' ] ) 
    % Interpolating to the current values
    tmp = interp1( phase_angle_deg, phase_function_value, opt.kplr.phs_ang_dg( i_pl, : ) ) ;
    elseif strcmp( fl_nm_tmp, 'rayleigh_w0p3' )
    % Load the data file: it brings phase_angle_deg and phase_function_value
    load( [ opt.scn_dr 'geomAlbedo/rayleigh_scalar_w0p3_phase_function.mat' ] )
    % Interpolating to the current values
    tmp = interp1( phase_angle_deg, phase_function_value, opt.kplr.phs_ang_dg( i_pl, : ) ) ;
    end
  % Other cases to follow

  % Check
    if isnan( tmp )
    disp( sprintf( 'The phase function for planet %i (%s) could not be identified. Fix it. Stopped.', i_pl, opt.kplr.phs_fnctn{ i_pl } ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  opt.kplr.phs_fnctn_vl( i_pl, : ) = tmp ;
  end % i_pl

% Some plot
if ( 0 )
clf ;n_pl=numel(opt.kplr.sma_au); a=-500*opt.px_scn_mas*opt.mas2au;b=1500*opt.px_scn_mas*opt.mas2au;hold on ; for i = 1 : n_pl, h(i)=plot( opt.kplr.x_ps_mas( i, : )* opt.mas2au, opt.kplr.y_ps_mas( i, : ) * opt.mas2au ) ; end ; xlabel( 'AU' ) ; ylabel( 'AU' ) ; xlim( [ a, b ] ) ; ylim( [ a, b ] ) ; legend( h, opt.kepler.planet_type, 'Location', 'SouthEast' ) ; title( '' )
clf ;hold all ; for i = 1:7, h(i)=plot( sort( opt.kplr.phs_fnctn_vl( i, : ), 'descend' ), '+-' ) ; end ; grid ; grid minor ; xlabel( [ 0, 61 ] ) ; xlabel( 'Months 2025-2030' ) ; ylabel( 'Lambert Phase function' ) ; legend( h, opt.kepler.planet_type ) ; xlim( [ 0, 61 ] )
end

% Getting the flux ratio for the planets under Keplerian motion
function opt = get_kepler_flux_ratio( opt )
% Check
  if ~isfield( opt.kplr, 'gm_albd' )
  disp( 'The Kepler option must contain the field for the Geometric Albedo: opt.kepler.geomAlbedo = { option 1, ... }. Set it. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

  if ~iscell( opt.kplr.gm_albd )
  disp( 'The option opt.kepler.geomAlbedo must be a cell: opt.kepler.geomAlbedo = { option 1, ... }. Fix it . Stopped.' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

% Reading the geometric albedo or each planet (reminder: opt.plnt.add.flx_rt was initially set in the beginning)
  for i_pl = 1 : opt.n_pl
  % Sentinel
  tmp = NaN ;
  % If it is a value, just use it (constant across wavelength)
    if isfloat( opt.kplr.gm_albd{ i_pl } )
      for i_lp = 1 : opt.n_lp
       % Need to do it for each band to match the structure of a case where it depends on the wavelength (see below)
        for i_bnd = 1 : opt.n_bnd
        opt.kplr.gm_albd_arry( i_pl, i_lp, i_bnd ) = opt.kplr.gm_albd{ i_pl } ;
        opt.plnt.add.flx_rt_arry( i_pl, i_lp, i_bnd ) = opt.kplr.sma_rd_rt_2( i_pl, i_lp ) * opt.kplr.gm_albd_arry( i_pl, i_lp ) * opt.kplr.phs_fnctn_vl( i_pl, i_lp ) ;
        end
      end
    tmp = 0 ;
    end
  % If it is a character, then open the corresponding file
    if ischar( opt.kplr.gm_albd{ i_pl } )
    % It brings lambda_nm_array and geomAlbedo_array
    load( [ opt.scn_dr 'geomAlbedo/' opt.kplr.gm_albd{ i_pl } '_geomAlbedo.mat' ] )
    % Interpolating to the current wavelength under consideration (linear since it may have local, steep features)
      for i_lmbd = 1 : numel( opt.lmbd_arry_img_nm )
      % Fine array of wavelengths every nm
       if ( i_lmbd < numel( opt.lmbd_arry_img_nm ) )
        lmbd_intrp_arry = opt.lmbd_arry_img_nm( i_lmbd ) : 1 : opt.lmbd_arry_img_nm( i_lmbd + 1 ) ;
        else % Next wavelength is the limit of the band
        lmbd_intrp_arry = opt.lmbd_arry_img_nm( i_lmbd ) : 1 : opt.lmbd_img_2_nm ;
        end
      % Average value
      gm_albd_arry_tmp( i_lmbd ) = sum( interp1( lambda_nm_array, geomAlbedo_array, lmbd_intrp_arry, 'linear' ) ) / numel( lmbd_intrp_arry ) ;
      end
    % Removing unnecessary field before next assignament
      if isfield( opt.plnt.add, 'flx_rt_arry' ) && ( i_pl == 1 )
      opt.plnt.add = rmfield( opt.plnt.add, 'flx_rt_arry' ) ; 
      end
      for i_bnd = 1 : opt.n_bnd
        for i_lmbd = 1 : numel( opt.lmbd_arry_img_nm )
        flx_rt_arry_tmp = opt.kplr.sma_rd_rt_2( i_pl, : ) .* squeeze( gm_albd_arry_tmp( i_lmbd ) ) .* opt.kplr.phs_fnctn_vl( i_pl, : ) ;
        opt.plnt.add.flx_rt_arry( i_pl, :, i_bnd, i_lmbd ) = flx_rt_arry_tmp ;
        end
      end
    tmp = 0 ;
    end
    if isnan( tmp )
    disp( sprintf( 'The geometric albedo for planet %i has not been identified. Fix it. Stopped', i_pl ) )
    disp( ' ' )
    disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
    end
  end % i_pl

% Copy for iterative calls to get_imaging_options
opt.planets.add.flux_ratio_array = opt.plnt.add.flx_rt_arry ;
% For the record
opt.str.plnt_flx_arry = opt.str.flx * opt.plnt.add.flx_rt_arry ;

% Celestial coordinate changes from Ecliptic to Equatorial and viceversa. The output is 1x#coordinates.
function [ ln_dg_2, lt_dg_2 ] = celestial_coordinate_change_deg( ln_dg, lt_dg, ecl2eq, mssg )

% No messaging by default
  if ~exist( 'mssg', 'var' )
  mssg = 0 ;
  end

% Checks
  if ~exist( 'ecl2eq', 'var' )
  disp( 'Set Ecliptic->Equatorial (ecl2eq=1) or vice versa (ecl2eq=-1). Please fix it.' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

  if ( abs( ecl2eq ) ~= 1 )
  disp( 'Set Ecliptic->Equatorial (ecl2eq=1) or vice versa (ecl2eq=-1). Please fix it.' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

  if ( ecl2eq == 1 ) && ( mssg )
  disp( 'Changing celestical coordinates from Ecliptic to Equatorial.' )
  end

  if ( ecl2eq == -1 ) && ( mssg )
  disp( 'Changing celestial coordinates from Equatorial to Ecliptic.' )
  end

% Number of elements in the coordinate array
n_crd = numel( ln_dg ) ;

  if ( numel( lt_dg ) ~= n_crd )
  disp( 'The number of coordinates in longitude and latitude are different. Please fix it.' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

% Change to radians
ln_rd = ln_dg * pi / 180 ;
lt_rd = lt_dg * pi / 180 ;
C( 1, : ) = squeeze( cos( lt_rd ) .* cos( ln_rd ) ) ;
C( 2, : ) = squeeze( cos( lt_rd ) .* sin( ln_rd ) ) ;
C( 3, : ) = squeeze( sin( lt_rd ) ) ;

% Change from Equatorial to Ecliptic orthonormal basis (see Eq. 1.5.7 in Hipparcos_vol1_all_https- web.archive.org web 20160303180237 http- www.rssd.esa.int SA HIPPARCOS docs vol1_all.pdf.pdf)
epsln_rd = 23.4392911111 * pi / 180 ;
M = [ [ 1, 0,               0 ] ;
      [ 0, cos( epsln_rd ), -sin( epsln_rd ) ] ;
      [ 0, sin( epsln_rd ), cos( epsln_rd ) ] ] ;
% Actual change
  if ( ecl2eq == 1 ), C_2 = M * C ; end
  if ( ecl2eq == -1 ), C_2 = inv( M ) * C ; end

lt_dg_2 = asin( C_2( 3, : ) ) * 180 / pi ;
% In case we are exactly at the equatorial poles
q_pl = find( abs( lt_dg_2 ) == 90 ) ;
n_pl = numel( q_pl ) ;
  if ( n_pl )
  ln_dg_2( q_pl ) = atan2( 1e-13 * sign( lt_dg_2 ) .* C_2( 2, : ), 1e-13 * sign( lt_dg_2 ) .* C_2( 1, : ) ) * 180 / pi ;
  end

% All other cases
q_no_pl = find( abs( lt_dg_2 ) ~= 90 ) ;
ln_dg_2( q_no_pl ) = atan2( C_2( 2, : ), C_2( 1, : ) ) * 180 / pi ;
% making RA to be within 0 and 360 degrees
ln_dg_2 = mod( ln_dg_2 + 720, 360 ) ;

% Transposing to get a rowx1 column output (useful because ExoCat parameters are read in this order)
ln_dg_2 = ln_dg_2' ; lt_dg_2 = lt_dg_2' ;

% Deriving the parallax in RADEC (proper motion has already been added and then the parallax is calculated)
% The heliocentric position vector of a star is the same. From Starshade+Telescope's perspective, the star moves along its orbit -> parallax (pax)
% Simplifications:
% 1: Starshade-Telecope's orbit and Earth orbit are on the same plane (or ecliptic latitude Starshade+Telescope = 0. True for L2)
% 3: Starshade-Telescope's orbit locked with Earth (exact L2) in a circular (0.0167 -> ~another 1, negligible) orbit of 1 year period (this may not be good enough, since the influence of other planets on Earth's orbit might shift the position of stars by a similar amount as the parallax from J2000 to 2028, but fine for mission studies spanning a few years)
function opt = add_parallax( opt )
% Orbital phase (recall #3). NB: just remember the circular orbit assumption: it just adds a multiple of 2 pi.
phi_orb = 2 * pi * mod( opt.kplr.tm_yr_obs, 1 ) ;
% To align it with the dimensions of opt.str.ra_pm_dg and opt.str.dc_pm
phi_orb = phi_orb' ;
% Vector position of Starshade with respect the Sun amd J2000
d_L2S = (1.496e11 + 1.492e9 ) ; % m. Average distance Earth-Sun + distance Earth-L2 (1% correction only)
x_tlscp_s = d_L2S * cos( phi_orb ) ;
y_tlscp_s = d_L2S * sin( phi_orb ) ;
z_tlscp_s = 0 ; % Assuming L2 orbit on Ecliptical plane
% Vector position of the star from Starshade+Telescope
% Converting parallax to distance (opt.target.pax is in mas)
d_str = opt.str.dst_pc * 3.0857e16 ; % 1 pc = 3.0857e16 m
% Converting ra_obs, dec_obs to Ecliptic coordinates
[ ecl_lon_obs_dg ecl_lat_obs_dg ] = celestial_coordinate_change_deg( opt.str.ra_pm_dg, opt.str.dc_pm_dg, -1 ) ;
ecl_lon_obs_rd = ecl_lon_obs_dg * pi / 180 ;
ecl_lat_obs_rd = ecl_lat_obs_dg * pi / 180 ;
% Vector joining the star and the Starshade+Telescope. PS: for orbital phase equal to 0, Starshade+Telescope aligned with the reference time that defines J2000
x_str_tlscp_s = d_str * cos( ecl_lat_obs_rd ) .* cos( ecl_lon_obs_rd ) - x_tlscp_s ;
y_str_tlscp_s = d_str * cos( ecl_lat_obs_rd ) .* sin( ecl_lon_obs_rd ) - y_tlscp_s ;
z_str_tlscp_s = d_str * sin( ecl_lat_obs_rd ) ;
% Normalized position vector from WFIRST-S:
d_tlscp_s = sqrt( x_str_tlscp_s .* x_str_tlscp_s + y_str_tlscp_s .* y_str_tlscp_s + z_str_tlscp_s .* z_str_tlscp_s ) ;
x_nrm = x_str_tlscp_s ./ d_tlscp_s ;
y_nrm = y_str_tlscp_s ./ d_tlscp_s ;
z_nrm = z_str_tlscp_s ./ d_tlscp_s ;
% Deriving the new Ecliptic coordinates
ecl_lat_wpax_obs_dg = asin( z_nrm ) * 180 / pi ;
% atan2(Y,X)
ecl_lon_wpax_obs_dg = atan2( y_nrm, x_nrm ) * 180 / pi ;
% Ecliptic to Equatorial
[ ra_wpax_obs_dg dec_wpax_obs_dg ] = celestial_coordinate_change_deg( ecl_lon_wpax_obs_dg, ecl_lat_wpax_obs_dg, 1 ) ;
% Back to degrees
opt.str.ra_wpax_dg = ra_wpax_obs_dg' ;
opt.str.dc_wpax_dg = dec_wpax_obs_dg' ;
% Storing the ecliptic coordinates: longitude, latitude and helio-longitude (appendix A in http://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/reference-data-for-calibration-and-tools/documentation/_documents/TIR-CRDS-2015-01.pdf. Don't use the 180 degree shift in Eq.(7))
opt.str.ecl_ln_dg = ecl_lon_wpax_obs_dg' ;
opt.str.ecl_lt_dg = ecl_lat_wpax_obs_dg' ;
% L2 has the same period as Earth's orbit and the probe is aligned with the line going from the Sun to the Earh.
% March 20th approximated as 00:00 am 03/20: (31 + 28.25 + 19)=78.25 days out of 365.25 days/year
% Relative ecliptic longitude with respect to the Sun (elongation, or angle Sun-Earth-Target)
opt.str.hl_ecl_ln_dg = mod( ecl_lon_wpax_obs_dg' - 360 * ( mod( opt.kplr.tm_yr_obs, 1 ) - 78.25 / 365.25 ) + 720, 360 ) ;
% Deriving the corresponding surface brightness of the zodiacal light
opt.str.lcl_zd.mg_v_arcs2 = get_local_zodi_mag_v_arcsec2( opt.str.hl_ecl_ln_dg, opt.str.ecl_lt_dg ) ;
% Function to get the position of the target at some given time in Equatorial coordinates, including parallax
function opt = get_radec_pos_target( opt )
% Get proper motion and translate it into RA, DEC
% SIMBAD convention:
% pm-ra : mu-ra*cos(dec) (expressed in the ICRS system in mas/yr)
% pm-dec : mu-dec (expressed in the ICRS system in mas/yr)
% See also https://en.wikipedia.org/wiki/Proper_motion
mu_ra_mas_yr = opt.str.pm_ra_mas_yr / cos( opt.str.dc_dg * pi /180 ) ;
opt.str.ra_pm_dg = opt.str.ra_dg + mu_ra_mas_yr / 1e3 / 3600 * ( opt.kplr.tm_yr_obs - 2000 ) ;
opt.str.dc_pm_dg = opt.str.dc_dg + opt.str.pm_dc_mas_yr / 1e3 / 3600 * ( opt.kplr.tm_yr_obs - 2000 ) ;
% Parallax
opt = add_parallax( opt ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local zodiacal magnitudes                                                                 %
%                                                                                           %
% Source: STScI, TIR-CRDS-2015-01.pdf                                                       %
% (R. Diaz, 04/01/2015, http://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/reference-data-for-calibration-and-tools/documentation/_documents/TIR-CRDS-2015-01.pdf                                      %
% http://ssb.stsci.edu/cdbs_open/cdbs/work/etc/etc-cdbs/background/zodi_001.dat             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zodi_mag_v_arcsc2 = get_local_zodi_mag_v_arcsec2( helio_ecl_lon_deg, ecl_lat_deg )
% From http://ssb.stsci.edu/cdbs_open/cdbs/work/etc/etc-cdbs/background/zodi_001.dat

% Array of longitudes provided
hl_ecl_ln_dg_arry = 0 : 15 : 180 ;
% Array of latitudes provided
ecl_lt_dg_arry = 0 : 15 : 90 ;

% Table (longitude,latitude)
zodi_mag_v_arcsc2_tbl = ...
[ [   nan ,   nan ,   nan ,   nan ,   22.6,   23.0,   23.3 ] ; ...
  [   nan ,   nan ,   nan ,   nan ,   22.6,   23.1,   23.3 ] ; ...
  [   nan ,   nan ,   nan ,   22.3,   22.7,   23.1,   23.3 ] ; ...
  [   nan ,   nan ,   22.1,   22.5,   22.9,   23.1,   23.3 ] ; ...
  [   21.3,   21.9,   22.4,   22.7,   23.0,   23.2,   23.3 ] ; ...
  [   21.7,   22.2,   22.6,   22.9,   23.1,   23.2,   23.3 ] ; ...
  [   22.0,   22.3,   22.7,   23.0,   23.2,   23.3,   23.3 ] ; ...
  [   22.2,   22.5,   22.9,   23.1,   23.3,   23.3,   23.3 ] ; ...
  [   22.4,   22.6,   22.9,   23.2,   23.3,   23.3,   23.3 ] ; ...
  [   22.4,   22.6,   22.9,   23.2,   23.3,   23.4,   23.3 ] ; ...
  [   22.4,   22.6,   22.9,   23.1,   23.3,   23.4,   23.3 ] ; ...
  [   22.3,   22.5,   22.8,   23.0,   23.2,   23.4,   23.3 ] ; ...
  [   22.1,   22.4,   22.7,   23.0,   23.2,   23.4,   23.3 ] ] ;

[ X_arry Y_arry ] = meshgrid( ecl_lt_dg_arry, hl_ecl_ln_dg_arry ) ;

% Input values within limits
  if ( numel( find( helio_ecl_lon_deg < -180 ) ) ) || ( numel( find( helio_ecl_lon_deg > 360 ) ) )
  disp( 'WARNING: the helio-ecliptic longitude value(s) is/are out of the bounds (-180 deg, 180 deg). Stopped.' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end
  if ( numel( find( ecl_lat_deg < -90 ) ) ) || ( numel( find( ecl_lat_deg > 90 ) ) )
  disp( 'WARNING: the ecliptic latitude is out of the bounds (-90 deg, 90 deg). Stopped.' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end
% Respecting the symmetry between the line of sight and the Sun (towards/opposite directions)
q = find( helio_ecl_lon_deg < 0 ) ;
  if numel( q )
  helio_ecl_lon_deg( q ) = -helio_ecl_lon_deg( q ) ;
  end
q = find( helio_ecl_lon_deg > 180 ) ;
  if numel( q )
  helio_ecl_lon_deg( q ) = 360 - helio_ecl_lon_deg( q ) ;
  end
ecl_lat_deg = abs( ecl_lat_deg ) ;
% Interpolated value (checked order)
zodi_mag_v_arcsc2 = interp2( X_arry, Y_arry, zodi_mag_v_arcsc2_tbl, ecl_lat_deg, helio_ecl_lon_deg ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to send jobs to the cluster
function  send_sister_imaging_band( i_bnd, opt )
% TO BE WRITTEN
cmd = 'sister_imaging_band( i_bnd, opt )' ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wavelength grid for a given resolving power
function lmbd_arry_grd = wavelength_grid( opt )
% Multiplicative factor associated with a resolving power
fct_r = ( 2 * opt.rslvng_pwr + 1 ) / ( 2 * opt.rslvng_pwr - 1 ) ;
% Maximum wavelength
n_grd = floor( log10( opt.lmbd_img_2_nm / opt.lmbd_img_1_nm ) / log10( fct_r ) ) ;
lmbd_arry_grd( 1 ) = opt.lmbd_img_1_nm ;
  for i_grd = 1 : n_grd
  lmbd_arry_grd( i_grd + 1 ) = opt.lmbd_img_1_nm * fct_r^i_grd ;
  end
lmbd_arry_grd( end + 1 ) = opt.lmbd_img_2_nm ;
