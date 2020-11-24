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
function master_matrix( opt )
% Function to derive the master matrix for a given Starshade, band and pixel size
% July-August 2019: Sergi Hildebrandt (JPL/Caltech), and Saehui Hwang (Caltech).
%
% Examples:
% 1/ master_matrix 
% 2/ clear opt;opt.cube.planets.n_fwhm=3;opt.pix_camera_mas=21.1;master_matrix(opt)
%
% for n = [ 0.5, 1, 3 ], clear opt;opt.cube.planets.n_fwhm=n;opt.pix_camera_mas=21.1;master_matrix(opt) ; end
%

% Generic debug mode
dbstop if error 

  if ~exist( 'opt', 'var' )
  opt = [] ;
  end

opt = get_default_master_options( opt ) ;

% Master matrix: array to cover all necessary pixels within the non-stationary region (+/- 2*r_st_mas, and odd number)
x_mstr = 1 : 2 * ceil( opt.r_stationary_mas / opt.px_psf_mas ) + 1 ;
y_mstr = x_mstr ;
cntr_mstr = ( numel( x_mstr ) + 1 ) / 2 ;
ra_mas = ( cntr_mstr - x_mstr ) * opt.px_psf_mas ;
dc_mas = ( cntr_mstr - y_mstr ) * opt.px_psf_mas ;
[ ra_mstr_mas dc_mstr_mas ] = meshgrid( ra_mas, dc_mas ) ;

%%%% Running SISTER for a bunch of planets
n_ttl = numel( ra_mstr_mas ) ;
disp( sprintf( 'The total number of pixels of the Master Matrix is %i', n_ttl ) )

% Number of times that SISTER will be run 400 planets per bunch is a good balance)
n_pl = 100 * opt.px_psf_mas ; % This number may be chosen larger or smaller. The larger it is, SISTER will run less times, but then each time will need to allocate more memory for all the simultaneous planets. It depends on the computer.
n_sstr = floor( n_ttl / n_pl ) ;
  if ( n_sstr * n_pl < n_ttl )
  n_sstr = n_sstr + 1 ;
  end
disp( sprintf( 'A total of %i planets will be simulated each time', n_pl ) ) 

% Running configuration options for SISTER have been defined in get_default_master_options.m below
t_mstr_mtrx = tic ;
% Common label for all cubes
lbl_all = sprintf( 'master_%2.1f_mas_nfwhm_%2.1f', opt.px_psf_mas, opt.cube.planets.n_fwhm ) ;
  for i_str = 1 : n_sstr
  idx_1 = n_pl * ( i_str - 1 ) + 1 ; % Matlab starts with 1
  idx_2 = idx_1 + n_pl ; % each time n_pl
    if idx_2 > n_ttl
    idx_2 = n_ttl ;
    end
  % Add a label for each run
  opt.add_label = sprintf( '%s_%03i', lbl_all, i_str ) ;
  opt.planets.add.do = 1 ;
  opt.planets.add.pos_arc_ra_mas = ra_mstr_mas( idx_1 : idx_2 ) ;
  opt.planets.add.pos_arc_dec_mas = dc_mstr_mas( idx_1 : idx_2 ) ;
  opt.planets.add.flux_ratio = ones( idx_2 - idx_1 + 1, 1 ) ;
  % Timing each bunch of planets
  t_bnch = tic ;
  disp( sprintf( 'Running the bunch #%i/%i', i_str, n_sstr ) )
  [ d1 d2 d3 opt_out ] = sister_imaging( opt ) ;
  disp( sprintf( 'Time to finish the bunch #%i/%i was %2.1f minutes', i_str, n_sstr, toc( t_bnch ) / 60 ) )
  end

% Opening the files
  for i_str = 1 : n_sstr
  fl_nm = sprintf( 'output/cubes/sister_cube_%s_%03i_%inm_planets_no_background.mat', lbl_all, i_str, opt.delta_lambda_imaging_nm ) ; 
  load( fl_nm ) ;
  idx_1 = n_pl * ( i_str - 1 ) + 1 ; % Matlab starts with 1
  idx_2 = idx_1 + n_pl ; % each time n_pl
  if idx_2 > n_ttl
    idx_2 = n_ttl ;
  end
  % Pixel index, wavelength slice, x/y PSF data
  master_matrix( idx_1 : idx_2, :, :, : ) = permute( psf_cube, [ 2 1 3 4 ] ) ; % The order in the PSF cube is wavelength x #planet x "x" x "y", whereas the order in the master matrix is more convenient as #planet x wavelength x "x" x "y".
  % Remove the file
  delete( fl_nm ) ;
end

% Storing the indices of the PSF arrays in descending order (from highest to lowest)
% It takes just ~30 s.
t_idx = tic ;
% Size in pixels of the square PSF 
sz_psf_mstr_mtrx = size( master_matrix, 3 ) ;
  for i_pl = 1 : size( master_matrix, 1 )
  % The location of the maximum does not depend on the wavelength
  scn_psf = squeeze( master_matrix( i_pl, end, :, : ) ) ;
  [ dummy idx_tmp ] = sort( scn_psf( : ), 'descend' ) ;
  idx_sort_psf( i_pl, : ) = idx_tmp ;
  end

% Location and filename to store the results
dr_mstr_mtrx = sprintf( '%ssister_basis/%s/%s_%i_%i/%i_%i_nm/master_matrix/', sister_installation_path(), opt.starshade.mode, opt.starshade.nominal_filename, opt.starshade.number_of_petals, opt.Nx_pupil_pix, opt.lambda_imaging_1_nm, opt.lambda_imaging_2_nm ) ;
  if ~isdir( dr_mstr_mtrx )
  mkdir( dr_mstr_mtrx ) ;
  end
% Avoiding the string 'master' twice
  if strcmp( lbl_all( 1 : 7 ), 'master_' )
  lbl_all = lbl_all( 8 : end ) ;
  end
fl_nm_mstr_mtrx = sprintf( 'master_matrix_Nx%i_%s_geo_iwa_%03.2f_mas_%s_%04d_%04d_%04d_nm.mat', opt.Nx_pupil_pix, opt.pupil_filename, opt.geo_iwa_mas, lbl_all, opt.lambda_imaging_1_nm, opt.lambda_imaging_2_nm, opt.delta_lambda_imaging_nm ) ;
  % For delta_lambda_nm that is not an integer
  if opt.delta_lambda_imaging_nm ~= round( opt.delta_lambda_imaging_nm )
   fl_nm_mstr_mtrx = sprintf( 'master_matrix_Nx%i_%s_geo_iwa_%03.2f_mas_%s_%04d_%04d_%3.1f_nm.mat', opt.Nx_pupil_pix, opt.pupil_filename, opt.geo_iwa_mas, lbl_all, opt.lambda_imaging_1_nm, opt.lambda_imaging_2_nm, opt.delta_lambda_imaging_nm ) ;
  end

% Correspondance between RA and DEC with the index of the master matrix. Checked doing ra_mstr_mas(1:10),dc_mstr_mas(1:10) and looking that DEC changes while RA is constant. RA changes once a full column is read. Also, checked by doing the operation defined below, and looking at ra_mstr_mas and dc_mstr_mas for that master index and getting the same RA, DEC values.
ra_mas_mx = max( ra_mas ) ;
dc_mas_mx = max( dc_mas ) ; % By construction is the same value than for RA, but we label it for clarity
n_dc = numel( dc_mas ) ;
px_psf_mas = opt.px_psf_mas ;
RA_DEC_mas_to_Master_Matrix_Index = sprintf( 'round((%i-dec_mas)/%i) + 1 + %i*round((%i-ra_mas)/%i)', dc_mas_mx, px_psf_mas, n_dc, ra_mas_mx, px_psf_mas ) ;
% Saving the options used to generate the master matrix (opt_sister are store with the 'sister_cube_...' files in convolve_with_one_wavelength.m)
opt_master_matrix = opt_sister ; 
save( sprintf( '%s%s', dr_mstr_mtrx, fl_nm_mstr_mtrx ), 'master_matrix', 'opt_master_matrix', 'ra_mas_mx', 'dc_mas_mx', 'px_psf_mas', 'n_dc', 'RA_DEC_mas_to_Master_Matrix_Index', 'idx_sort_psf' ) ;
disp( sprintf( 'Master matrix generated and stored in %s%s', dr_mstr_mtrx, fl_nm_mstr_mtrx ) ) 
disp( sprintf( 'The whole process took %3.1f minutes', toc( t_mstr_mtrx ) / 60 ) )

%%%%%%%%%%%%%%%%%
% SUB-FUNCTIONS %
%%%%%%%%%%%%%%%%%
% Default options for master matrix
function opt = get_default_master_options( opt )

% Options realated with the configuraytion file for SISTER

% General optional label that will be appended to the output images and data files in addition of any other labels set by the imaging software.
  if ~isfield( opt, 'add_label' )
  opt.add_label = '' ;
  end

%%%%%%%%%%%%%%%%%%
% Starshade mode %
%%%%%%%%%%%%%%%%%%
  if ~isfield( opt, 'starshade' )
  opt.starshade.dummy = 1 ; % Dummy
      if ~isfield( opt.starshade, 'mode' )
      opt.starshade.mode='spinning' ;  % Default is 'spinning'
      end
  end

  if ~isfield( opt.starshade, 'nominal_filename' )
  opt.starshade.nominal_filename = 'NI2' ;
  end

% Getting some of the basic properties from the Matlab file with the Starshade occulter
% Basically, we need the teescope's diameter below
opt = new_occulter_from_matlab_file( opt ) ;

%%%%%%%%%%%%%
% Telescope %
%%%%%%%%%%%%%
% For NI2, use the WFIRST pupil (sister_imaging will set it up)
  if ~isfield( opt, 'pupil_filename' ) && ~strcmp( lower( opt.starshade.nominal_filename( 1 : 3) ), 'ni2' )
  opt.pupil_filename = 'ideal' ; 
  end
% Sampling of the pupil (16 is fast enough for creating the PSF basis on a laptop)
  if ~isfield( opt, 'Nx_pupil_pix' )
  opt.Nx_pupil_pix = 16 ;
  end

% Optional jitter of the telescope
  if ~isfield( opt, 'jitter' )
  opt.jitter.do = 0 ;
  else
    if ( opt.jitter.do )
      if ~isfield( opt.jitter, 'rms_mas' )
      opt.jitter.rms_mas = 15 ;
      end
    end
  end

%%%%%%%%%%%%%%%%%%%%
% Wavelength range %
%%%%%%%%%%%%%%%%%%%%
% Pixel scale of the PSF objects. It should be the same or less than the pixel scale of the scene.
  if ~isfield( opt, 'px_psf_mas' )
  opt.px_psf_mas = 3 ; % Default is 1 mas
  end

% Making the scene have the same pixel scale as the PSF
opt.pix_scene_mas = opt.px_psf_mas ;

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

opt.lambda_band_nm_min = opt.lambda_imaging_1_nm ;
opt.lambda_band_nm_max = opt.lambda_imaging_2_nm ;

opt.r_stationary_mas = set_r_stationary_mas( opt, 1 ) ;

opt.lambda_1_nm = opt.lambda_imaging_1_nm ;
opt.lambda_2_nm = opt.lambda_imaging_2_nm ;

[ opt.distance_starshade_telescope_m opt.geo_iwa_mas ] = set_starshade_distance_and_geometric_iwa( opt, 1 ) ;

%%%%%%%%%%%%
% Detector %
%%%%%%%%%%%%
% Pixel scale of the camera. Default:
  if ~isfield( opt, 'pix_camera_mas' )
  opt.pix_camera_mas = 21.1 ; % 0.42 l/D, l=575 nm, D=2.37 m, WFIRST
  end

%%%%%%%%%%%%%%%%%%%%
% Creating a scene %
%%%%%%%%%%%%%%%%%%%%
opt.scene.do = 1 ; % Default value is 0. It would read a pre-existing scene
% Scene units
opt.scene.units = 'Jy' ; % No default value is set. The user must confirm the units are correctly set to Jy
% Field of View of the astrophysical scene in mas
  if ~isfield( opt.scene, 'fov_diam_mas' )
  opt.scene.fov_diam_mas = 1000 ; % Enough to cover all planets in the stationary region
  end

%%%%%%%%%%%
% Planets %
%%%%%%%%%%%

% Example of adding 'static' planets. That is, without an orbital motion
opt.planets.add.do = 1 ; % Default is 0

%%%%%%%%%%%%%%%%%%%%
% SISTER Data Cube %
%%%%%%%%%%%%%%%%%%%%
opt.cube.do = 1 ; % Default is 0. SISTER does not output a FITS file with th spectral results of the simulation.
opt.cube.fits = 0 ; % Whether the output is a FITS or Matlab file. By default, a Matlab file.
opt.cube.background = 0 ; % Sets whether the background is included or subtracted. If 1, it is included. if 0, it is subtracted.
opt.cube.planets.do = 1 ; % Extract cubes for each planet
  if ~isfield( opt.cube.planets, 'n_fwhm' )
  opt.cube.planets.n_fwhm = 3 ; % Area to be stored with the planet at the center in terms of FWHM (array is center +/- n_fwhm/2)
  end
% Display the range of the number of pixels that will be used for building up the PSF
% This is approximate. SISTER derives th actual FWHM for each PSF in convolve_with_one_wavelength.m
n_px = 2 * ceil( opt.cube.planets.n_fwhm * 1.02 * opt.lambda_imaging_1_nm * 1e-9 / opt.diameter_telescope_m * 180 * 3600e3 / pi / opt.pix_camera_mas / 2 ) + 1 ;
disp( sprintf( 'The PSF cubes will be %i pix x %i pix for %04.2f nm', n_px, n_px, opt.lambda_imaging_1_nm ) ) ;
n_px = 2 * ceil( opt.cube.planets.n_fwhm * 1.02 * opt.lambda_imaging_2_nm * 1e-9 / opt.diameter_telescope_m * 180 * 3600e3 / pi / opt.pix_camera_mas / 2 ) + 1 ;
disp( sprintf( 'And %i pix x %i pix for %04.2f nm', n_px, n_px, opt.lambda_imaging_2_nm ) )

% Other generic options
opt.scene_dir = sprintf( '%sinput_scenes/', sister_installation_path ) ;
opt.star.type = 'sun' ;

