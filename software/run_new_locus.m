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
function int_src_prt = run_new_locus( opt_locus )

% Warning: do not output opt from this function, since it could overwrite some values

% Function to generate the on-axis response for a non-ideal starshade. Non-ideal starshades usually have much higher levels of starlight than the nominal case. Hence, the spatial extent of the on-axis PSF is larger. By default, get_locus_options below sets it to 4 times as large in both dimensions (16x area).

% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

%%%%%%%%%%%%%%
% Preliminary
% Re-doing the work? If 1, it will re-compute the electric fields at the pupil plane and at the image plane, and the PSF intensity from scratch. 
% If 0, it will check if the corresponding results exist. If they exist, they are loaded. If they don't, they will be created.
  if ~isfield( opt_locus, 'redo_perturbed' )
  redo_perturbed = 0 ;
  else
  redo_perturbed = opt_locus.redo_perturbed ;
  end

% Path where the results will be stored
opt_locus.save_path = sprintf( '%s%s', opt_locus.inst_dr, opt_locus.out_dr ) ;

% If all rotations are set to be done
all_rt = 0 ;
  if ~isnan( opt_locus.alph_lst_rd )
  opt_locus.n_rot = numel( opt_locus.alph_lst_rd ) ;
  % Sentinel that all rotations will be made
  all_rt = 1 ;
  end

% Compiling more options
opt_lcs = get_locus_options( opt_locus ) ;

% Size of the image plane in pixels
opt.nx_img = 2 * floor( opt_lcs.diam_img_mas / opt_lcs.px_res_mas / 2 ) + 1 ;
% Source location (by default, it's on-axis)
opt.x_source_mas = opt_lcs.x_source_mas ;
opt.y_source_mas = opt_lcs.y_source_mas ;
% Pupil mask
opt.pupil_filename = opt_lcs.pupil_filename ;
opt.nx_pupil_pix = opt_lcs.nx_pupil_pix ;
opt.secondary_size = opt_lcs.secondary_size ;
% Path where the locus files are to be found
opt.path_occulter = [ opt_lcs.inst_dr opt_lcs.in_dr ] ;
  if ~isdir( opt.path_occulter )
  mkdir( opt.path_occulter ) ;
  disp( [ 'Input directory created: ' opt.path_occulter ] ) ;
  end

% Other generic parameters
opt.lambda_1_nm = opt_lcs.lmbd_tmp_nm ;
% Monochromatic
opt.lambda_2_nm = opt_lcs.lmbd_tmp_nm ;
opt.delta_lambda_nm = 1 ; % No effect, since lambda_2_nm == lambda_1_nm

  if isfield( opt_lcs, 'n_ptl' )
  opt.n_ptl = opt_lcs.n_ptl ;
  end
opt.starshade = opt_lcs.starshade ;
opt.px_psf_mas = opt_lcs.px_psf_mas ;
opt.geo_iwa_mas = opt_lcs.geo_iwa_mas ;
opt.fl_nm_prt = opt_lcs.fl_nm_prt ;
opt.r_st_mas = opt_lcs.r_st_mas ;
opt.au2km = opt_lcs.au2km ;
opt.au2mas = opt_lcs.au2mas ;
opt.epsln = opt_lcs.epsln ;
opt.rltv_earth = opt_lcs.rltv_earth ;
  if strcmp( lower( opt.starshade.mode ), 'spinning' )
  opt.n_symm = opt_lcs.n_symm ;
  opt.transmission_curve_tmp = opt_lcs.transmission_curve_tmp ;
  opt.half_transmission_mas = opt_lcs.half_transmission_mas ;
  end
% Dummy
opt.diameter_telescope_m = opt_lcs.diameter_telescope_m ;
opt.starshade.nominal_filename = opt_lcs.prtrbd ;
opt.make_occulter_name = opt_lcs.prtrbd ;
opt.distance_starshade_telescope_m = opt_lcs.distance_starshade_telescope_m ;
opt.alph_lst_rd = opt_lcs.alph_lst_rd ;
opt.chck_prcsn = opt_lcs.chck_prcsn ;
opt.imrotate = opt_lcs.imrotate ;
% New file for the shape of the petals:
  if isfield( opt, 'save_filename' )
  opt = rmfield( opt, 'save_filename' ) ;
  end
opt = get_default_options( opt ) ;

% Path where the intermediate results may be stored
save_path_UTotL = [ opt_locus.save_path 'UTotL/' ] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar as in sister_basis.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For spinning, first check whether a file with the spinning PSF already exists (from a previous run)
  if strcmp( lower( opt.starshade.mode ), 'spinning' )
  % Checking if the PSF response exists
    if ( exist( opt_lcs.fl_nm_prt, 'file' ) == 2 ) && ~( redo_perturbed )
    % It will load int_src_prt, and opt_new_locus
    load( opt_lcs.fl_nm_prt ) ;
    % Whether the precision was set or estimated with an Earth-like planet
      if ~isnan( opt_lcs.epsln )
      % Distance in pc where a 10% of an Earth-like planet would have the same contrast as the precision goal with half transmission (worse case)
      dst_pc_epsln = 1000 / opt_lcs.geo_iwa_mas * ( 6371 / opt_lcs.au2km ) * sqrt( 0.2 * 0.1 / opt_new_locus.d_psf_mx * 0.5 ) ; % Multiplied by 0.5 to account for the half transmission 
      % Equivalent distance of the radius of half transmission in AU
      r_hlf_au = opt_lcs.half_transmission_mas / 1000 * dst_pc_epsln ;
      disp( sprintf( 'Loading the non-ideal response (precision %1.2e. The residual error is equivalent to 10%s of the contrast of an Earth at %2.1f pc, at an orbital distance equal to %2.1f AU): %s', opt_new_locus.d_psf_mx, '%', dst_pc_epsln, r_hlf_au, opt_lcs.fl_nm_prt ) )
      else
        if ~( all_rt ) % if not all rotations were set to be made 
        % Distance in pc where a 10% of an Earth-like planet would have the same contrast as the precision goal with half transmission (worse case)
        dst_pc_epsln = 10 * sqrt( opt_new_locus.epsln / opt_new_locus.d_psf_mx ) ; % This option found the worst radius that is greater than or equal to the half transmission IWA and used its transmission.
        % Equivalent distance of the radius of half transmission in AU
        r_hlf_au = opt_lcs.half_transmission_mas / 1000 * dst_pc_epsln ;
        disp( sprintf( 'Loading the non-ideal response (precision %1.2e. The residual error is equivalent to 10%s of the contrast of an Earth at %2.1f pc, at an orbital distance of %2.1f AU): %s', opt_new_locus.d_psf_mx, '%', dst_pc_epsln, r_hlf_au, opt_lcs.fl_nm_prt ) )
        end
      end
    % Done, so return
    return
    end
  % In case the file with the spinning PSF is not available, derive it:
  % If an angle step is not provided is there will be as many rotations as to achieve a particular precison goal
    if isnan( opt.alph_lst_rd )
    % The minimum step of the rotations is such that at the tip of the petals, the rotation is half the the pixel size set to build the PSF. 
    dlt_alph_psf = opt.px_psf_mas / opt.geo_iwa_mas / 2 ;
    % The precision on the non-ideal PSF construction controls the level of error that may be incurred. There are in fact several possible rotations that would give rise to a PSF that satisfies the criterion. In order to match the visual aspect of the non-ideal PSF with the idea of a spinning starshade, a series of alternating angles define the sequence by which the PSF precisoon goal is sought.
    % Closer power of 2 that matches dlt_aph_psf
    pwr_2_psf = ceil( log2( pi / dlt_alph_psf ) ) ;
    alph_lst_rd = [ 0 ] ; % Zero and pi rotation must be set outside the loop
      for i_alph = 0 : pwr_2_psf
      alph_tmp = unique( mod( pi / 2^i_alph + ( 0 : pi / 2^( i_alph - 1 ) : 2 * pi - pi / 2^( i_alph - 1 ) ), 2 * pi ) ) ; % Avoiding duplication
      alph_lst_rd = [ alph_lst_rd alph_tmp ] ;
      end
    opt.alph_lst_rd = alph_lst_rd ; 
    end
  else % spinning
  % For non-spinning, no rotations:
  opt.alph_lst_rd = 0 ;
  end
n_rot = numel( opt.alph_lst_rd ) ;

% Constructing the *stationary* PSF for the new configuration (no need to rotate it). Used to estimate the contrast goal when rotating the on-axis response of the non-ideal starshade.
  if ( opt.chck_prcsn )
  opt_st = opt ;
  opt_st.x_source_mas = opt.r_st_mas ; 
  opt_st.y_source_mas = 0 ;
  opt_st = get_default_options( opt_st ) ;
  fl_q = opt_st.save_filename ;
  % It is convenient to rename the file, because this code may be used different times with different specifications for the pupil
  fl_q = sprintf( '%s_pupil_%s_%i_pix_psf_pix_%i_mas', fl_q, opt_lcs.lbl_ppl, opt_st.nx_pupil_pix, opt_lcs.px_psf_mas ) ;
  fl_q_3 = sprintf( '%s%s.mat', save_path_UTotL, fl_q ) ;
    if exist( fl_q_3, 'file' ) ~= 2 || ( redo_perturbed )
      if ( exist( fl_q_3, 'file' ) ~= 2 ) && ~( opt.imrotate )
      disp( sprintf( '* The file %s... does not exist. Creating it ...', fl_q_3 ) )
      end
      if ( redo_perturbed ) && ( exist( fl_q_3, 'file' ) == 2 ) && ~( opt.imrotate )
      disp( sprintf( '* The file %s... exists, but will be created again (redo_perturbed=%i).', fl_q_3, redo_perturbed ) )
      end
      if ~( redo_perturbed )
      opt_st.save = 1 ;  % Store the longer calculations
      opt_st.save_path = save_path_UTotL ;
      end
    opt_st.save_filename = fl_q ;
    [ UTotL lambdaIn dummy pupil ] = makeStarshadeImage( opt_st ) ;
    else
      if ( opt.verbose )
      disp( sprintf( '* Loading stored results from %s (redo_perturbed=%i)', fl_q_3, redo_perturbed ) )
      end
    load( fl_q_3 )
    end
  % Removing hot pixels if required
    if ( opt_st.erase_hot_pixels ), UTotL = erase_hot_pixels( UTotL, opt_st.erase_hot_pixels ) ; end
  opt_st.lD2mas = opt_lcs.lmbd_tmp_nm / 1e9 / opt_st.dmtr_tlscp_m * 180 / pi * 3600 * 1000 ;
  % Making it half the extension of the perturbed one (this is only used to estimate the relative error between two iterations of the perturbed on-axis response, and we avoid potential issues with aliasing at the edges).
  opt_st.nx_img = 2 * floor( opt_st.nx_img / 4 ) + 1 ; % Odd number
  imagePlaneDiameterInLambdaOverD = opt_st.nx_img * opt_st.px_psf_mas / opt_st.lD2mas ;
  imagePlane = mft_shift( imagePlaneDiameterInLambdaOverD, opt_st.nx_img, UTotL .* pupil, opt_st.x_source_mas, opt_st.y_source_mas, opt_st.px_psf_mas ) ;
  % Intensity of the stationary PSF
  int_src_st = abs( imagePlane ).^2 ; % No need to normalize, since the perturbed one won't be normalized when used together, see later on.
  int_src_st = int_src_st / sqrt( 2 ) ; % Although we don't normalize the PSF, the factor 1/sqrt(2) needs to be added because the on-axis response below is derived using UTot2ImagePlane.m and there's a division by sqrt(2) added (=impeak = mft_shift( 1, Nx_img, pupil, 0, 0, 1 ) ; % Use 1 L/D square. D=diameter of the primary mirror ... peak = max(abs(impeak(:)));)
  % Storing its sum
  opt.sum_int_src_st = sum( sum( int_src_st ) ) ;  
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Constructing the spinning on-axis PSF for a non-ideal starshade. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fl_q = opt.save_filename ;
  if ( opt.imrotate )
  fl_q = sprintf( '%s_imrotate', fl_q ) ;
  end
% It is convenient to rename the file, because this code may be used different times with different specifications for the pupil
fl_q = sprintf( '%s_pupil_%s_%i_pix', fl_q, opt_lcs.lbl_ppl, opt.nx_pupil_pix ) ;
int_src_prt = 0 ;

  for i_rot = 1 : n_rot
    if ( opt.alph_lst_rd( i_rot ) )
    str_alph = strrep( sprintf( '%1.4f', opt.alph_lst_rd( i_rot ) ), '.', 'p' ) ;
    fl_q_2 = sprintf( '%s_alpha_%s', fl_q, str_alph ) ;
    else
    fl_q_2 = fl_q ;
    end
  fl_q_3 = sprintf( '%s%s.mat', save_path_UTotL, fl_q_2 ) ;
    if exist( fl_q_3, 'file' ) ~= 2 || ( redo_perturbed )
      if ( exist( fl_q_3, 'file' ) ~= 2 ) && ~( opt.imrotate )
      disp( sprintf( '* The file %s... does not exist. Creating it ...', fl_q_3 ) )
      end
      if ( redo_perturbed ) && ( exist( fl_q_3, 'file' ) == 2 ) && ~( opt.imrotate )
      disp( sprintf( '* The file %s... exists, but will be created again (redo_perturbed=%i).', fl_q_3, redo_perturbed ) )
      end
      if ~( redo_perturbed )
      opt.save = 1 ;  % Store the longer calculations
      opt.save_path = save_path_UTotL ;
      end

    % Full propagation of optical response (default)
      if ~( opt.alph_lst_rd( i_rot ) ) || ~( opt.imrotate ) % At least do full propagation for 0 degree rotation, to be used as the templae for imrotate (not default)
      opt.starshade_rotation_rad = opt.alph_lst_rd( i_rot ) ;
      opt.save_filename = fl_q_2 ;
      disp( sprintf( 'Considering the rotation %i of a maximum of %i rotations', i_rot, n_rot ) )
      [ UTotL lambdaIn dummy pupil ] = makeStarshadeImage( opt ) ;
      else
      int_src_prt_tmp = imrotate( int_src_prt_0dg, opt.alph_lst_rd( i_rot ) * 180 / pi, 'crop', 'bilinear' ) ;
        if ( ( i_rot / 20 ) == floor( i_rot / 20 ) ) % Don't print out this message too many times
        disp( '* Using imrotate instead of full propagation of the optical response (not default, set in the options).' )
        end
      end
    else
      if ( opt.verbose )
      disp( sprintf( 'Considering the rotation %i of a maximum of %i rotations', i_rot, n_rot ) )
      disp( sprintf( '* Loading stored results from %s (redo_perturbed=%i)', fl_q_3, redo_perturbed ) )
      end
    load( fl_q_3 )
    end

    % Full propagation of optical response (default)
    if ~( opt.alph_lst_rd( i_rot ) ) || ~( opt.imrotate )
    % Removing hot pixels if required
      if ( opt.erase_hot_pixels )
      UTotL = erase_hot_pixels( UTotL, opt.erase_hot_pixels ) ; 
      end
    opt.diam_img_mas = opt_lcs.diam_img_mas ;
    efDefectImg = UTot2ImagePlane( lambdaIn, opt, pupil, UTotL ) ;
    int_src_prt_tmp = abs( efDefectImg ).^2 ;
      if ( opt.imrotate ) % At least do full propagation for 0 degree rotation, to be used as the templae for imrotate (not default)
      int_src_prt_0dg = int_src_prt_tmp ;
      end
    end
  % After at least 4 rotations ('0','pi/2','pi','3pi/2'), store previous perturbed intensity
    if ( i_rot > 4 ) && ( opt.chck_prcsn )
    int_src_prt_prv = int_src_prt / ( i_rot - 1 ) ;
    end
  % Perturbed intensity
  int_src_prt = int_src_prt + int_src_prt_tmp ; % No need to normalize, since both reference and perturbed have the same number of rotations, and will be later normalized in convolve_with_one_wavelength.m
  % After at least 4 rotations ('0','pi/2','pi','3pi/2'), decide whether to stop building the PSF
  stp_psf = 0 ;
  % Store actual number of rotations
  opt.i_rot = i_rot ;
    if ( i_rot > 4 ) && ( opt.chck_prcsn )
    % Comparison between the previous and the new PSF
    [ stp_psf opt ] = compare_psf( int_src_prt_prv, int_src_prt / i_rot, int_src_st, opt ) ;
      if stp_psf, break ; end
    end

  % Generating the unperturbed on-axis response for proper calibration (1 position is enough, since the perturbed starlight is greater than the nominal one and only a first order calibration is necessary -tests documented in SISTER technical google doc)
    if ~( opt.alph_lst_rd( i_rot ) ) 
    opt_nominal = opt ;
    opt_nominal.make_occulter_name = opt_lcs.starshade.nominal_filename ;
    opt_nominal.starshade = opt_lcs.starshade ;
    opt_nominal.starshade_rotation_rad = 0 ;
    % New file for the shape of the petals:
      if isfield( opt_nominal, 'save_filename' )
      opt_nominal = rmfield( opt_nominal, 'save_filename' ) ;
      end
    opt_nominal = get_default_options( opt_nominal ) ;    
    UTotL_nominal = makeStarshadeImage( opt_nominal ) ;
      if ( opt_nominal.erase_hot_pixels )
      UTotL_nominal = erase_hot_pixels( UTotL_nominal, opt_nominal.erase_hot_pixels ) ;
      end
    efDefectImg_nominal = UTot2ImagePlane( lambdaIn, opt_nominal, pupil, UTotL_nominal ) ;
    int_src_nominal = abs( efDefectImg_nominal ).^2 ;
    end
  end % i_rot

  % In the spinning case, if the goal was not achieved, change the label to the right level of accuracy
  if strcmp( lower( opt.starshade.mode ), 'spinning' ) && ( opt.i_rot == n_rot ) && ( opt.chck_prcsn ) 
    if (  opt.d_psf_mx > opt.epsln )
    % Modifying the label in the filename to store the results
    old_lbl = strrep( strrep( sprintf( '%1.1e', opt.epsln ), '.', 'p' ), '-', 'm' ) ;
    nw_lbl = strrep( strrep( sprintf( '%1.1e', opt.d_psf_mx ), '.', 'p' ), '-', 'm' ) ;
    opt.fl_nm_prt = sprintf( '%s_%s.mat', opt.fl_nm_prt( 1 : end - 3 - numel( old_lbl ) ), nw_lbl ) ;
    disp( sprintf( 'WARNING: After all available %i rotations, the precision obtained has been %1.2e, which is not the goal set of %1.2e. Storing the results for the actual precision of %1.2e in %s', n_rot, opt.d_psf_mx, opt.epsln, opt.d_psf_mx, opt.fl_nm_prt ) ) 
    end
  end


%%%%%%%%%%%%%%%%%%%%%%%
% Normalizing the PSF %
%%%%%%%%%%%%%%%%%%%%%%%
% Number of rotations
int_src_prt = int_src_prt / opt.i_rot ;
% preserving the ratio stationary/nominal on-axis
int_src_prt = int_src_prt * opt_locus.mx_psf_nmnl / max( int_src_nominal( : ) ) ;


% Storing the PSF if required.
  if ( opt_locus.save_work ) 
  opt_new_locus = opt ;
  save( opt.fl_nm_prt, 'int_src_prt', 'opt_new_locus' ) ;
  disp( sprintf( 'File stored %s, after %i rotations', opt.fl_nm_prt, opt.i_rot ) ) ;
  % Removing the rotated individual fields (Disabled until experience shows whether it is convenient to keep them or not)
  %  delete( sprintf( '%s_alph*.*', fl_q ) ) ; 
  end

%%%%%%% SUB-FUNCTIONS
function opt_locus = get_locus_options( opt_locus )
opt_locus = lower_opt( opt_locus ) ;

  if ~isdir( opt_locus.save_path ) && ( opt_locus.save_work )
  mkdir( opt_locus.save_path ) ;
  disp( sprintf( 'Output directory created: %s', opt_locus.save_path ) ) ;
  end

% Whether the construction of the spinning, non-ideal PSF will stop at some precision goal, or consider all rotations.
  if isnan( opt_locus.alph_lst_rd )
  opt_locus.chck_prcsn = 1 ;
  else
  opt_locus.chck_prcsn = 0 ; 
  end

% Results for an on-axis source
  if ~isfield( opt_locus, 'x_source_mas' )
  opt_locus.x_source_mas = 0 ;
  end
  if ~isfield( opt_locus, 'y_source_mas' )
  opt_locus.y_source_mas = 0 ;
  end

% From the previous value we set the following two values
% Spatial extent of the PSF (the nominal PSF were rescaled to the pixel scale of the scenes)
opt_locus.diam_img_mas = 2 * floor( opt_locus.tm_nmnl * opt_locus.sz_psf_rf * opt_locus.px_scn_mas / 2 ) + 1 ; % Odd number
% Pupil mask sampling
opt_locus.nx_pupil_pix = ceil( opt_locus.tm_nmnl * opt_locus.nx_pupil_pix ) ; % integer, not necessarily power of 2

% For a spinning starshade
  if strcmp( lower( opt_locus.starshade.mode ), 'spinning' )
    if ~isfield( opt_locus, 'n_symm' )
    opt_locus.n_symm = 1 ;
    end
  end

% Label to associate with the pupil
  if strcmp( opt_locus.pupil_filename, '0' )
  opt_locus.lbl_ppl = 'ideal' ;
  else
  opt_locus.lbl_ppl = opt_locus.pupil_filename ;
    if strcmp( lower( opt_locus.filename( end - 3 : end ) ), '.mat' ) || strcmp( lower( opt_locus.filename( end - 3 : end ) ), '.fits' )
    opt_locus.lbl_ppl = opt_locus.lbl_ppl( 1 : end - 4 ) ;
    end
  end

  if ~strcmp( opt_locus.lbl_ppl, '' )
  opt_locus.lbl_ppl_2 = [ '_' opt_locus.lbl_ppl ] ;
  else
  opt_locus.lbl_ppl_2 = opt_locus.lbl_ppl ;
  end

% Filename storing the results

% I/O directory for the files with the non-ideal spinning PSF
dr_out = sprintf( '%s%s/', opt_locus.save_path, opt_locus.prtrbd ) ;
  if ( opt_locus.save_work )
    if ~isdir( dr_out )
    mkdir( dr_out )
    end
  end

lbl_non_spnng = '' ;
% If it is not spinning, add a label (usually non-spinning PSF are not stored, because the derivation is fast from the stored electric fields)
  if strcmp( opt_locus.starshade.mode, 'non-spinning' )
  lbl_non_spnng = 'non_spinning_' ;
  end
opt_locus.fl_nm_prt = sprintf( '%s%s%s_%i_mas_%i_mas_%04dnm_pupil_%ipix_psf_pix_%i_mas%s.mat', dr_out, opt_locus.prtrbd, lbl_non_spnng, opt_locus.x_source_mas, opt_locus.y_source_mas, opt_locus.lmbd_tmp_nm, opt_locus.nx_pupil_pix, opt_locus.px_psf_mas, opt_locus.lbl_ppl_2 ) ;

  % If a precision goal was set
  if ( opt_locus.chck_prcsn )
    % Adding a label in case an absolute precision for the perturbed PSF is set
    if ~isnan( opt_locus.epsln )
    opt_locus.fl_nm_prt = sprintf( '%s_%s.mat', opt_locus.fl_nm_prt( 1 : end - 4 ), strrep( strrep( sprintf( '%1.1e', opt_locus.epsln ), '.', 'p' ), '-', 'm' ) ) ;
    end
  end
 % Adding a label if the rotated PSF response was derived without the full optical propagation
  if ( opt_locus.imrotate )
  opt_locus.fl_nm_prt = sprintf( '%s_imrotate.mat', opt_locus.fl_nm_prt( 1 : end - 4 ) ) ;
  end
  if ~( opt_locus.chck_prcsn )
  opt_locus.fl_nm_prt = sprintf( '%s_nrot_%i.mat', opt_locus.fl_nm_prt( 1 : end - 4 ), opt_locus.n_rot ) ;
  end

% Function used to compare the progress of the PSF building. Once reached a target precision, the construction stops
function [ stp_psf opt ] = compare_psf( psf_0, psf_1, psf_st, opt )

% Rounding the position of the half transmission to the PSF pixel scale
r_hlf = round( opt.half_transmission_mas / opt.px_psf_mas ) ;
[ X_arry Y_arry ] = meshgrid( 1 : size( psf_0, 1 ), 1 : size( psf_0, 2 ) ) ;
X_arry = X_arry - ( size( psf_0, 1 ) + 1 ) / 2 ;
Y_arry = Y_arry - ( size( psf_0, 2 ) + 1 ) / 2 ;
R_arry_2 = X_arry .* X_arry + Y_arry .* Y_arry ;
d_psf = abs( psf_1 - psf_0 ) ;
% Zeroing all the points that are interior to the half transmission region, as a convention to compute the precision in the construction of the accumulated PSF
d_psf( R_arry_2 <= r_hlf^2 ) = 0 ;
% Maximum deviation
q_mx = find( d_psf == max( d_psf( : ) ) ) ;
q_mx = q_mx(1);  % new line added by stuart 080519
% Estimating the euivalent intensity of that feature as if it were fitted by the local PSF (local PSF ~ PSF_stationary * transmission)
sz_psf_1 = size( psf_0, 1 ) ;
sz_psf_2 = size( psf_0, 2 ) ;
cntr_psf_1 = ( sz_psf_1 - 1 ) / 2 ;
cntr_psf_2 = ( sz_psf_2 - 1 ) / 2 ;
% psf_st was derived with about half the extension of the perturbed psf (see comments above when generating int_src_st).
psf_st_tmp = psf_0 * 0 ;
cntr_psf_st = ( size( psf_st, 1 ) + 1 ) / 2 ;
hlf_psf_st = ( size( psf_st, 1 ) - 1 ) / 2 ;
psf_st_tmp( cntr_psf_st - hlf_psf_st : cntr_psf_st + hlf_psf_st, cntr_psf_st - hlf_psf_st : cntr_psf_st + hlf_psf_st) = psf_st ;
clear psf_st
% psf_1 - psf_0, without absolute value, because photometric errors may cancel eachother, and there's no need to make the calculation more restrictive (all abs). Deriving the photometry that would be at that location (outside the half transmission area, the PSF is approximated by the stationary PSF and the transmission value. Otherwise, it would be necessary to build the local PSF, but that's a second order correction to the first correction considered here, which is fine for an error estimate)
opt.d_psf_mx = abs( sum( sum( circshift( psf_1 - psf_0, [ cntr_psf_1 - mod( q_mx, sz_psf_1 ) cntr_psf_2 - ceil( q_mx / sz_psf_2 ) ]  ) .* psf_st_tmp ) ) / opt.sum_int_src_st ) ;

% Transmission at that distance
trns_r_d = 1 ;
% Radial distance where it occurs
r_d_mas = round( opt.px_psf_mas * sqrt( ( q_mx / size( psf_0, 1 ) - ( size( psf_0, 1 ) + 1 ) / 2 )^2 + ( mod( q_mx, size( psf_0, 2 ) ) - ( size( psf_0, 2 ) + 1 ) / 2 )^2 ) ) ;
  if ( r_d_mas < opt.r_st_mas )
  trnd_r_d = opt.transmission_curve_tmp( r_d_mas ) ;
  end 
  % The criterion is by default referred to a 10% of an Earth-like planet at 10 pc, or it can directly be set.
  if isnan( opt.epsln )
  % Earth-size planet contrast at this angular distance (closest approach to star assumes r_d_mas is the periastron) at 10 pc. It also assumes face-on (worst case scenario)
  r_earth_km = 6371 ;
  % The 10 pc criterion is also used when displaying the precision and loading the file above (l.116 'equivalent to an Earth')
  dst_km_d_psf = ( r_d_mas / 100 ) * opt.au2km ;
  % Assuming geometric albedo of 0.2 for Earth, which is a fair non-high average between 400-1000 nm. Relative contribution: 10%.
  opt.epsln = trns_r_d * ( r_earth_km / dst_km_d_psf )^2 * 0.2 * opt.rltv_earth ; 
  else
  opt.epsln = opt.epsln * trns_r_d ; % Getting the contrast goal into proper intensity with respect the stationary PSF
  end

% Deciding whether to stop
  if ( opt.d_psf_mx <= opt.epsln )
  stp_psf = 1 ;
  disp( sprintf( '* Non-ideal PSF precision goal achieved: %1.2e, after %i rotations', opt.d_psf_mx, opt.i_rot ) )
  else
  stp_psf = 0 ;
  disp( sprintf( 'Actual precision %1.2e (> goal of %1.2e). Continuing PSF construction ...', opt.d_psf_mx, opt.epsln ) )
  end

