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
function [ scn_dt_cnv_e_1, scn_dt_cnv_no_plnt_e_1, scn_ns_e_1, opt_1, scn_dt_cnv_e_2, scn_dt_cnv_no_plnt_e_2, scn_ns_e_2, opt_2 ] = sister( opt ) 
% Example of a user interface to run sister. sister can be run directly supplying the running  options in a structure variable. This script can handle imaging of single or double scenes.
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

% * Update the path to your local installation by editing sister_installation_path.m
  if ~isfield( opt, 'installation_path' )
  opt.installation_path = sister_installation_path() ;
  end

% Ouput variables ('1' refers to the first, or single, scene being images, and '2' to the second scene, if any):
% scn_dt_cnv_e_1,2: scene data, convolved, in count (e) units.
% scn_dt_cnv_no_plnt_e_1,2: same as before but wothout planets. That is, the starlight, exo-zodi (if any), solar glint (if any), astrophysical background (if any), etc ... It can be used to get 'planets only': scn_dt_cnv_e_1 - scn_dt_cnv_no_plnt_e_1.
% scn_ns_e_1,2: The noise associated with the convolved scene data.
% opt_1,2: (structure) all the options used to generate the output data.

% Generic debug mode (exit it with dbquit all)
dbstop if error

% Initializing output variables (only useful if you want the results in running time. Otherwise, use opt.save=1 and load the results)
scn_dt_cnv_e_1 = [] ; scn_dt_cnv_e_2 = [] ; scn_dt_cnv_no_plnt_e_1 = [] ; scn_dt_cnv_no_plnt_e_2 = [] ; scn_ns_e_1 = [] ; scn_ns_e_2 = [] ; opt_1 = [] ; opt_2 = [] ;

% Basic checks
  if ~exist( 'opt', 'var' )
  disp( 'Please provide some options. For instance, opt.run=# or ''some_case'';sister(opt);' )
  return
  end

  if ~isfield( opt, 'run' )
  disp( 'Please provide a running configuration case. For instance, opt.run=# or ''some_case'';sister(opt);' )
  return
  end

% 1 or 2 scenarios?
  if ischar( opt.run )
  n_scnr = 1 ;
  else
  n_scnr = numel( opt.run ) ;
  end
% Keeping track whether this is a double scenario observation.
opt.n_scnr = n_scnr ;
  % No more than 2 scenarios at a time
  if ( n_scnr > 2 )
    if ischar( opt.run )
    n_scnr = 1 ;
    else
    disp( 'The number of scenarios must be either 1 or 2. It can''t be more.' )
    return
    end
  end
  % Cell?
  if ( n_scnr == 1 ) && iscell( opt.run )
    disp( 'In the case of running one configuration, opt.run must be either a numerical identifier or a character vector, and not a cell.' )
  return
  end

  if ( n_scnr == 2 )
    if ischar( opt.run )
    n_scnr = 1 ;
    else
      if ~iscell( opt.run )
      disp( 'Setup opt.run as a cell: opt.run={case_1,case_2}' )
      return
      end
    end
  end

% Adding the path of the running configuration files
opt_tmp = get_running_options( opt.installation_path ) ;
addpath( opt_tmp.run_dir ) ;

%%%%%%%%%%%%%%%%%%
% First scenario %
%%%%%%%%%%%%%%%%%%
clear opt_1
  if n_scnr == 1
  run_1 = opt.run ;
  else
  run_1 = opt.run{ 1 } ;
  end
  if strcmp( run_1( end - 1 : end ), '.m' ), run_1 = run_1( 1 : end - 2 ) ; end
opt.case = run_1 ;
  if isnumeric( run_1 )
  opt_1 = id_config( opt ) ;
  end
  if ischar( run_1 )
  % Check whether it exists
    if exist( [ opt_tmp.run_dir run_1 '.m' ], 'file' ) ~= 2
    disp( 'The scenario does not correspond to any file in config/. Returning.' )
    return
    end
  fnctn_hndl = str2func( run_1 ) ; opt_1 = fnctn_hndl( opt ) ; opt_1.run = run_1 ;
  end

% Basic checks
  if ~exist( 'opt_1', 'var' )
  disp( 'The first element of opt.run must be either a numerical identifier or a character vector. It seems it was not the case.' )
  return
  end

  if numel( fieldnames( opt_1 ) ) == 2
    if isnumeric( run_1 )
    disp( 'The scenario has not been identified. Check its identifier matches one in config/id_config.m' )
    else
    disp( 'The scenario has not been identified. Check the spelling and make sure such file exists in config/' )
    end
  return
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If there's a second scenario %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ( n_scnr == 2 )
  clear opt_2
  run_2 = opt.run{ 2 } ;
  if strcmp( run_2( end - 1 : end ), '.m' ), run_2 = run_2( 1 : end - 2 ) ; end
  opt.case = run_2 ;
    if isnumeric( opt.run{ 2 } )
    opt_2 = id_config( opt ) ;
    end
    if ischar( run_2 )
      % Check whether it exists
      if exist( [ opt_tmp.run_dir run_2 '.m' ], 'file' ) ~= 2
      disp( 'The second scenario does not correspond to any file in config/. Returning.' )
      return
      end
    fnctn_hndl = str2func( run_2 ) ; opt_2 = fnctn_hndl( opt ) ; opt_2.run = run_2 ;
    end

  % Basic checks
    if ~exist( 'opt_2', 'var' )
    disp( 'Remember that the second cnfiguration in opt.run must be either a numerical identifier or a character vector. It seems it was not the case with the second scenario.' )
    return
    end  
  end % n_scnr==2

  if numel( fieldnames( opt_1 ) ) == 2
    if isnumeric( run_1 )
    disp( 'The scenario has not been identified. Check its identifier matches the one in config/id_config.m' )
    else
    disp( 'The scenario has not been identified. Check the spelling and make sure such file exists in config/' )
    end
  return
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUNNING STARSHADE_IMAGING %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % If there's only one scenario, just run it
  if ( n_scnr == 1 )
  opt_1 = get_running_options( opt.installation_path, opt_1, opt ) ;
  [ scn_dt_cnv_e_1, scn_dt_cnv_no_plnt_e_1, scn_ns_e_1, opt_1, scn_dt_cnv_lcl_zd_e_1 ] = sister_imaging( opt_1 ) ;
  end
  % Make sure the epoch is consistent between both scenarios if there's an array of dates

% If there's a second one, derive the difference (store the results if either the first or the second case were set to store the results)
  if n_scnr == 2
  opt_1 = get_running_options( opt.installation_path, opt_1, opt ) ;
  opt_2 = get_running_options( opt.installation_path, opt_2, opt ) ;
  % So far, for plotting purposes, this method of comparison can only handle one wavelength band at a time
    if ( numel( opt_1.lambda_imaging_1_nm ) > 1 ) || ( numel( opt_2.lambda_imaging_2_nm ) > 1 )
    disp( 'WARNING: using sister to analyze 2 scenes and the differences only works for *1* band. Set opt.lambda_imaging_1_nm and opt.lambda_imaging_2_nm to 1 value only in the running configuration file. Returning.' )
    return
  % Same bands
      if ( opt_1.lambda_imaging_1_nm ~= opt_2.lambda_imaging_1_nm ) || ( opt_1.lambda_imaging_2_nm ~= opt_2.lambda_imaging_2_nm )
      disp( sprintf( '* WARNING: The values of the band for the first scenario are [%04.2f,%04,2f]nm, whereas for the second scene are [%04.2f,%04.2f]nm. Usually, one would compare results for the *same* band. The code will continue though.', opt_1.lambda_imaging_1_nm, opt_1.lambda_imaging_2_nm, opt_2.lambda_imaging_1_nm, opt_2.lambda_imaging_2_nm ) )
      end
    end 
  % Checking the times are the same
    if isfield( opt_1, 'kepler' )
      if sum( abs( opt_1.kepler.time_year_obs - opt_2.kepler.time_year_obs ) )
      % In general, the time of observation should exactly be the same between both runs, so we could return the code here, without creating any images. It may happen though that sometimes, one has different targets observed at different times, but still wants to derive a comparison without writing an specific running configuration file for only changing the times of observation. We choose the latter for now. It should be an odd case, though.
      disp( sprintf( 'The epoch for the images is different between both running configurations. Choosing the second one as the reference. Make sure this is what you want.' ) )
      end
    % Copying the array of times
    opt.kepler.time_year = opt_2.kepler.time_year ;
    opt.kepler.time_year_obs = opt_2.kepler.time_year_obs ;
    opt.n_lp = opt_2.n_lp ;
    else % if kepler is not a field, it is a single epoch observation
    opt.n_lp = 1 ;
    end

  % Turning off plotting in sister, but keeping track of what's intended
  opt.plot.do = ( opt_1.plot.do ) || ( opt_2.plot.do ) ;
    if isfield( opt_1, 'video' ) && isfield( opt_2, 'video' )
    % By default choose the second options
    % turn off video of the first set
    opt_1.video.do = 0 ;
    disp( 'Using the video options from the second scene (default) for the joint imaging' )
    elseif ( isfield( opt_1, 'video' ) && ~isfield( opt_2, 'video' ) )
    disp( 'Setting the video options from the first scene for the joint imaging. Notice the second one did not have any video options.' )
    opt_2.video = opt_1.video ;
    % turn off video of the first set
    opt_1.video.do = 0 ;
    elseif ~( isfield( opt_1, 'video' ) && isfield( opt_2, 'video' ) )
    disp( 'Using the video options from the second scene (default) for the joint imaging. Notice the first one did not have any video options.' )
    end

  % Running date by date
    for i_lp = 1 : opt.n_lp
    t_strt = tic ;
    disp( sprintf( 'Considering the image %i/%i', i_lp, opt.n_lp ) )
    % Iterating data by date
    % Storing the single wavelength intermediate images while looping
    opt_1.save_single_wavelength = 1 ;
    opt_2.save_single_wavelength = 1 ;
    % Erase them at the end of the loop
      if ( i_lp == opt.n_lp )
      opt_1.save_single_wavelength = 0 ;
      opt_2.save_single_wavelength = 0 ;
      end
    
    fl_dt_out = sprintf( '%simage_data/sister_%s_and_%s_%s.mat', opt_1.output_dir, opt_1.cs_nm, opt_2.cs_nm, sprintf( '%04i', i_lp ) ) ; % the last loop tag is consistent with what's defined in sister_imaging_band.m (opt.tg_lp)
      if ( do_img_f( fl_dt_out, opt_1 ) )
      opt_1.i_lp_1 = i_lp ; opt_2.i_lp_1 = i_lp ;
      opt_1.i_lp_2 = i_lp ; opt_2.i_lp_2 = i_lp ;
        if ( i_lp > 1 )
        opt_1.planets.add.flux_ratio_array = opt_1.plnt.add.flx_rt_arry ; opt_2.planets.add.flux_ratio_array = opt_2.plnt.add.flx_rt_arry ;
        end
      % Running the first scene
      disp( sprintf( 'Considering the case %s', opt_1.cs_nm ) )
      [ scn_dt_cnv_e_1, scn_dt_cnv_no_plnt_e_1, scn_ns_e_1, opt_1, scn_dt_cnv_lcl_zd_e_1 ] = sister_imaging( opt_1 ) ;
      % Running the second scene
      disp( sprintf( 'Considering the case %s', opt_2.cs_nm ) )
      [ scn_dt_cnv_e_2, scn_dt_cnv_no_plnt_e_2, scn_ns_e_2, opt_2, scn_dt_cnv_lcl_zd_e_2 ] = sister_imaging( opt_2 ) ;
      disp( sprintf( 'Imaging %i/%i of both scenes took %3.1fs', i_lp, opt.n_lp, toc( t_strt ) ) )
      % If data want to be stored, better here with both scenarios and all the options (not the difference, to save some space)
        if ( opt_1.save_two_scenes )
        save( fl_dt_out, 'scn_dt_cnv_e_1', 'scn_dt_cnv_no_plnt_e_1', 'scn_ns_e_1', 'opt_1', 'scn_dt_cnv_e_2', 'scn_dt_cnv_no_plnt_e_2', 'scn_ns_e_2', 'opt_2' ) ;
        disp( sprintf( 'File %s saved', fl_dt_out ) )
        end
      else
      load( fl_dt_out ) 
      disp( sprintf( 'File %s read', fl_dt_out ) )
      end

    % Plotting (plotting will be one if either opt_1 or opt_1.plot.do was 1 or if it is set from the call to sister.m)
      if ( opt.plot.do )
      % Updating the real number of loops
      opt_1.n_lp = opt.n_lp ;
      opt_2.n_lp = opt.n_lp ;
      % Checking the plot options are compatible between both scenes
      flds_1 = sort( fieldnames( opt_1.plot ) ) ;
      flds_2 = sort( fieldnames( opt_2.plot ) ) ;
      n_flds_1 = numel( flds_1 ) ;
      n_flds_2 = numel( flds_2 ) ;
        if ( n_flds_1 ) ~= ( n_flds_2 )
        disp( sprintf( '* WARNING: the number of plot options in opt_1.plot, %i, is different to opt_2.plot, %i. Make sure the setup of the running configuration file is what you want. The options of opt_2.plot will be used.', n_flds_1, n_flds_2 ) )
        else
          for i_fld = 1 : n_flds_1
            if ~strcmp( flds_1{ i_fld }, flds_2{ i_fld } )
            disp( sprintf( '* WARNING: the field #%i in opt_1.plot (%s) and opt_2.plot (%s) do not coincide. Make sure the setup of the running configuration file is what you want. The options of opt_2.plot will be used.', i_fld, opt_1.plot( flds_1{ i_fld } ), opt_2.plot( flds_2{ i_fld } ) ) )
            end
          end
        end % n_flds_1 ~= n_flds_2
      % Setting both plot options to be the same
      opt_1.plot = opt_2.plot ;

      % Keeping track of whether to save the images, and turning off image storage in the individual plots
        if isfield( opt_2, 'sv_plt' )
        opt.sv_plt = opt_2.sv_plt ;
        opt_1.sv_plt =  0 ; 
        opt_2.sv_plt =  0 ; 
        end

      % If there was an offset, subtract it
        if isfield( opt_1.plt, 'bckgrnd_scn' )
          if ~isnan( opt_1.plt.bckgrnd_scn )
          opt_1.plot.min = opt_1.plt.mn - opt_1.plt.bckgrnd_scn ;
          opt_1.plt.mn = opt_1.plot.min ;
          opt_1.plot.max = opt_1.plt.mx - opt_1.plt.bckgrnd_scn ;
          opt_1.plt.mx = opt_1.plot.max ;
          end
        end

      % Getting the number of signal/noise/background & zoom combinations
      opt = get_running_options( opt.installation_path, opt_1, opt_2 ) ; % Allowing the user to overwrite the oes from the runing configuration file from the command line.
       if isfield( opt_2.plot, 'combination_list' )
       cmb_lst = opt_2.plot.combination_list ;
         if isfield( opt.plot, 'combination_list' )
         cmb_lst = opt.plot.combination_list ;
         end
       else
       cmb_lst = 1 ;
       end
       if isfield( opt_2.plot, 'zoom_list' )
       zm_lst = opt_2.plot.zoom_list ;
         if isfield( opt.plot, 'zoom_list' )
         zm_lst = opt.plot.zoom_list ;
         end
       else
       zm_lst = 1 ;
       end
       if isfield( opt_2.plot, 'add_planet_label' )
       add_plnt_lbl = opt_2.plot.add_planet_label ;
         if isfield( opt.plot, 'add_planet_label' )
         add_plnt_lbl = opt.plot.add_planet_label ;
         end
       else
       add_plnt_lbl = 1 ;
       end
       for cmb_tmp = cmb_lst
         for zm_tmp = zm_lst
           for add_plnt_lbl_tmp = add_plnt_lbl
             if isfield( opt_2.plot, 'combination_list' )
             opt_2.plot.combination_list = cmb_tmp ; opt_2.plot.combination_list = cmb_tmp ;
             end
             if isfield( opt_2.plot, 'zoom_list' )
             opt_2.plot.zoom_list = zm_tmp ; opt_2.plot.zoom_list = zm_tmp ;
             end
             if isfield( opt_2.plot, 'add_planet_label' )
             opt_2.plot.add_planet_label = add_plnt_lbl ; opt_2.plot.add_planet_label = add_plnt_lbl ;
             end
           % Keep the same min/max values for each frame
             if ( i_lp > 1 )
             opt_1.plot.min = mn_plt_i_lp_1 ; opt_2.plot.min = mn_plt_i_lp_1 ;
             opt_1.plot.max = mx_plt_i_lp_1 ; opt_2.plot.max = mx_plt_i_lp_1 ;
             opt_1.kplr.arc_ra_mas = arc_ra_mas_i_lp_1 ; opt_2.kplr.arc_ra_mas = arc_ra_mas_i_lp_1 ;
             opt_1.kplr.arc_dc_mas = arc_dc_mas_i_lp_1 ; opt_2.kplr.arc_dc_mas = arc_dc_mas_i_lp_1 ;
             end
           [ opt_1_two_scenes, opt_2_two_scenes ] = plot_two_scenarios( opt_1, opt_2, scn_dt_cnv_e_1, scn_dt_cnv_no_plnt_e_1, scn_ns_e_1, scn_dt_cnv_lcl_zd_e_1, scn_dt_cnv_e_2, scn_dt_cnv_no_plnt_e_2, scn_ns_e_2, scn_dt_cnv_lcl_zd_e_2, opt, i_lp ) ;
             if ( i_lp == 1 )
             mn_plt_i_lp_1 = opt_1_two_scenes.plot.min ;
             mx_plt_i_lp_1 = opt_1_two_scenes.plot.max ;
               if isfield( opt_1_two_scenes.kplr, 'arc_ra_mas' )
               arc_ra_mas_i_lp_1 = opt_1_two_scenes.kplr.arc_ra_mas ;
               arc_dc_mas_i_lp_1 = opt_1_two_scenes.kplr.arc_dc_mas ;
               end
             end
           % Setting the same limits for the rest of the images as the first one
           end % add_plnt_lbl_tmp
         end % zm_tmp
       end % cmb_tmp
     end % i_lp
   end % opt.plot.do
 end % n_scnr == 2

%%%%%%%%%%%%%%%%%
% SUB-FUNCTIONS %
%%%%%%%%%%%%%%%%%
% Getting some additional options that need not necessarily be in the running configuration files, 
% including the major relative paths, given the installation one. opt_2 can also be used to force some
% options, despite what is setup in the configuration file (it should only be high level options, or it
% could be confusing as to know what exactly was run)
function opt_1 = get_running_options( installation_path, opt_1, opt_2 )

  if ~exist( 'opt_1', 'var' )
  opt_1 = [] ;
  end
  
  if ~exist( 'opt_2', 'var' )
  opt_2 = [] ;
  end

% Main directories
  if ~isfield( opt_1, 'scene_dir' )
  opt_1.scene_dir = [ installation_path 'input_scenes/' ] ;
  end
  if ~isfield( opt_1, 'run_dir' )
  opt_1.run_dir = [ installation_path 'config/' ] ;
  end
  if ~isfield( opt_1, 'output_dir' )
  opt_1.output_dir = [ installation_path 'output/' ] ;
  end
  if ~isfield( opt_1, 'psf_dir' )
  opt_1.psf_dir = [ installation_path 'sister_basis/' ] ;
  end

% Related with sequence of images

% Always keeping track of the number of scenarios being imaged
  if isfield( opt_2, 'n_scnr' )
  opt_1.n_scnr = opt_2.n_scnr ;
  end

  if ~isfield( opt_1, 'case' ) && isfield( opt_1, 'run' )
  opt_1.case = opt_1.run ;
  end
  if ~isfield( opt_2, 'case' ) && isfield( opt_2, 'run' )
  opt_2.case = opt_2.run ;
  end

  if isfield( opt_1, 'case' )
    if isnumeric( opt_1.case )
    opt_1.cs_nm = sprintf( '%i', opt_1.case ) ;
    end
    if ischar( opt_1.case )
    opt_1.cs_nm = opt_1.case ;
    end
  end

  if isfield( opt_2, 'case' )
    if isnumeric( opt_2.case )
    opt_2.cs_nm = sprintf( '%i', opt_2.case ) ;
    end
    if ischar( opt_2.case )
    opt_2.cs_nm = opt_2.case ;
    end
  end

  if isfield( opt_1, 'plot' )
    if ~isfield( opt_1.plot, 'do' )
    opt_1.plot.do = 0 ;
    end
  else
  opt_1.plot.do = 0 ;
  end

  if ~isfield( opt_1, 'redo' )
  opt_1.redo = 0 ;
  end

  if isfield( opt_1, 'kepler' )
    if ( opt_1.kepler.do )
      if ~isfield( opt_1.kepler, 'time_year' )
      opt_1.kepler.time_year.initial = 2025 ;
      opt_1.kepler.time_year.final = 2030 ;
      opt_1.kepler.time_year.interval = 0.98765 ; % Avoiding 1 year not to plot the Earth again at the same location in the enxt image
      else
        if ~isfield( opt_1.kepler.time_year, 'initial' )
        opt_1.kepler.time_year.initial = 2025 ;
        end
        if ~isfield( opt_1.kepler.time_year, 'final' )
        opt_1.kepler_1.time_year.final = 2030 ;
        end
        if ~isfield( opt_1.kepler.time_year, 'interval' )
        opt_1.kepler.time_year.interval = 0.98765 ; % Avoiding 1 year not to plot the Earth again at the same location in the enxt image
        end
      end
   
    tm_yr_arry = opt_1.kepler.time_year.initial : opt_1.kepler.time_year.interval : opt_1.kepler.time_year.final ;
    % Storing the times of observation
    opt_1.kepler.time_year_obs = tm_yr_arry ;
    % updating the number of loops
    opt_1.n_lp = numel( tm_yr_arry ) ;
    end % kepler.do
  end % isfield kepler

% Options that may be changed from the call of sister
  if isfield( opt_2, 'scene_dir' )
  opt_1.scene_dir = [ installation_path 'input_scenes/' ] ;
  end
  if isfield( opt_2, 'run_dir' )
  opt_1.run_dir = [ installation_path 'config/' ] ;
  end
  if isfield( opt_2, 'output_dir' )
  opt_1.output_dir = [ installation_path 'output/' ] ;
  end
  if isfield( opt_2, 'psf_dir' )
  opt_1.psf_dir = [ installation_path 'sister_basis/' ] ;
  end

  if isfield( opt_2, 'plot' )
  fld_nm = fieldnames( opt_2.plot ) ;
    for i_nm = 1 : numel( fld_nm )
    opt_1.plot.( fld_nm{ i_nm } ) = opt_2.plot.( fld_nm{ i_nm } ) ;
    end
  end

  if isfield( opt_2, 'redo' )
  opt_1.redo = opt_2.redo ;
  end

  if isfield( opt_2, 'add_label' ) && length( opt_2.add_label )
  opt_1.add_label = opt_2.add_label ;
  else
    if isfield( opt_1, 'case' )
      if isnumeric( opt_1.case )
        if isfield( opt_1, 'add_label' ) && length( opt_1.add_label )
        opt_1.add_label = [ num2str( opt_1.case ) '_' opt_1.add_label ] ;
        else
        opt_1.add_label = num2str( opt_1.case ) ;
        end
      end
      if ischar( opt_1.case )
        if isfield( opt_1, 'add_label' ) && length( opt_1.add_label )
        opt_1.add_label = [ opt_1.case '_' opt_1.add_label ] ;
        else
        opt_1.add_label = opt_1.case ;
        end
      end
    end
  end

% Order of the difference between two scenes
  if ~isfield( opt_1, 'diff' )
  opt_1.diff = '2-1' ;
  end
  if isfield( opt_2, 'difference' )
  opt_1.diff = opt_2.difference ;
  end
% check of consistency
  if ~strcmp( opt_1.diff, '1-2' ) && ~strcmp( opt_1.diff, '2-1' )
  disp( 'Either choose ''1-2'' or ''2-1'' for the difference between the two scenes. Returning' )
make_a_stop
  return
  end

% Portion of the plot range to be used with the difference plot
  if ~isfield( opt_1, 'factor_plot_difference' )
  opt_1.fct_dff = 0.1 ; % 10%
  else
  opt_1.fct_dff = opt_1.factor_plot_difference ;
  end

% Storing the data of both scenes together
  if ~isfield( opt_1, 'save_two_scenes' )
  opt_1.save_two_scenes = 0 ;
  end

  if isfield( opt_2, 'save_two_scenes' )
  opt_1.save_two_scenes = opt_2.save_two_scenes ;
  end

% Controlling some plot combinations
  if isfield( opt_2, 'plot' )
    if isfield( opt_2.plot, 'combination_list' )
    opt_1.plot.combination_list = opt_2.plot.combination_list ;
    end
    if isfield( opt_2.plot, 'zoom_list' )
    opt_1.plot.zoom_list = opt_2.plot.zoom_list ;
    end
    if isfield( opt_2.plot, 'add_planet_label' )
    opt_1.plot.add_planet_label = opt_2.plot.add_planet_label ;
    end
    if isfield( opt_2.plot, 'starshade_circle' )
    opt_1.plot.starshade_circle = opt_2.plot.starshade_circle ;
    end
  end

function [ opt_1, opt_2 ] = plot_two_scenarios( opt_1, opt_2, scn_dt_cnv_e_1, scn_dt_cnv_no_plnt_e_1, scn_ns_e_1, scn_dt_cnv_lcl_zd_e_1, scn_dt_cnv_e_2, scn_dt_cnv_no_plnt_e_2, scn_ns_e_2, scn_dt_cnv_lcl_zd_e_2, opt, i_lp ) 

setwinsize( gcf, 1300, 400 )
%clf

% Plotting the first scenario (notice that if there were different bands, only the last one gets plotted)
subaxis( 1, 3, 1 ) ; %, 'pl', -0.01, 'pr', 0.01, 'pt', 0.0, 'pb', 0.01 ) ;
opt_1 = plot_starshade_figures( scn_dt_cnv_e_1, scn_dt_cnv_no_plnt_e_1, scn_dt_cnv_lcl_zd_e_1, scn_ns_e_1, opt_1, i_lp, numel( opt_1.lambda_imaging_1_nm ) ) ;

  % If there was an offset, subtract it
  if isfield( opt_1.plt, 'bckgrnd_scn' )
    if ~isnan( opt_1.plt.bckgrnd_scn )
    opt_1.plot.min = opt_1.plt.mn - opt_1.plt.bckgrnd_scn ;
    opt_1.plt.mn = opt_1.plot.min ;
    opt_1.plot.max = opt_1.plt.mx - opt_1.plt.bckgrnd_scn ;
    opt_1.plt.mx = opt_1.plot.max ;
    end
  end

% Plotting the difference    
% Whether 2nd minus 1st simulation is going to be plotted, or the reverse.
fct_21 = 1 ;
  if strcmp( opt_1.diff, '1-2' )
  fct_21 = -1 ;
  end
% Parsing the plot options and setting some specific ones
opt_dff = opt_2 ;
% No labels
opt_dff.plot.x_label = '' ;
opt_dff.plot.y_label = '' ;
% Simple units on colorbar
opt_dff.plot.colorbar_unit = 'e' ;
opt_dff.plot.max = max( [ -opt_1.plot.min, opt_1.plot.max ] ) * opt.fct_dff ;
opt_dff.plot.min = - opt_dff.plot.max ;
% If the comparison between the two scenes is noiseless, use the full range
  if ~isfield( opt_1.plot, 'n_sigma_noise' )
  opt_1.plot.n_sigma_noise = ( opt_dff.plot.max - opt_dff.plot.min ) ;
  end
opt_dff.plot.n_sigma_noise = opt_1.plot.n_sigma_noise * opt.fct_dff ;
% Combining the noise of both images into a common noise distribution (assuming uncorrelated noise)
ns_12 = ( scn_ns_e_2 + scn_ns_e_1 ) / sqrt( 2 ) ;
% Title
ns_12_rms = std( ns_12( : ) ) ; 
% Notice that the average noise is not the noise of the difference, but the average noise of both scenarios. That gives a better way to compare the difference with the expected instrumental noise.
opt_dff.plot.title = 'Signal difference right minus left' ;
  if ( fct_21 == - 1 )
  opt_dff.plot.title = 'Signal difference left minus right' ;
  end

subaxis( 1, 3, 2 ) ; %, 'pl', -0.01, 'pr', 0.01, 'pt', 0.0, 'pb', 0.01 ) ;
opt_dff = plot_starshade_figures( fct_21 * ( scn_dt_cnv_e_2 - scn_dt_cnv_e_1 ), fct_21 * ( scn_dt_cnv_no_plnt_e_2 - scn_dt_cnv_no_plnt_e_1 ), fct_21 * ( scn_dt_cnv_lcl_zd_e_1 - scn_dt_cnv_lcl_zd_e_2), 0 * ( scn_ns_e_1 + scn_ns_e_2 ) / sqrt( 2 ), opt_dff, i_lp, numel( opt_2.lambda_imaging_1_nm ) ) ;

subaxis( 1, 3, 3 ) ; %, 'pl', -0.01, 'pr', 0.01, 'pt', 0.0, 'pb', 0.01 ) ;
opt_2.plot = opt_1.plot ; opt_2.plt = opt_1.plt ;
opt_2 = plot_starshade_figures( scn_dt_cnv_e_2, scn_dt_cnv_no_plnt_e_2, scn_dt_cnv_lcl_zd_e_2, scn_ns_e_2, opt_2, i_lp, numel( opt_2.lambda_imaging_1_nm ) ) ;
 % If there was an offset, subtract it
  if isfield( opt_1.plt, 'bckgrnd_scn' )
    if ~isnan( opt_1.plt.bckgrnd_scn )
    opt_2.plot.min = opt_2.plt.mn - opt_2.plt.bckgrnd_scn ;
    opt_2.plt.mn = opt_2.plot.min ;
    opt_2.plot.max = opt_2.plt.mx - opt_2.plt.bckgrnd_scn ;
    opt_2.plt.mx = opt_2.plot.max ;
    end
  end
suptitle( strrep( sprintf( '%s. Run configurations: %s (left) and %s (right)', opt_2.plt.spttl, opt_1.cs_nm, opt_2.cs_nm ), '_', '\_' ) ) ;

fg_nm = sprintf( 'sister_%s_and_%s', opt_1.add_label, opt_2.add_label ) ;
img = getframe( gcf ) ;
% Creating the directory if it does not exist
  if ~isdir( opt.fg_dr )
  mkdir( opt.fg_dr ) ;
  end
  % Static image or video
  if ( opt_2.vd.do )
  % Creating the directory if it does not exist
    if ~isdir( opt.vd_dr )
    mkdir( opt.vd_dr ) ;
    end
  % Write to the GIF File
  fg_nm_vd = sprintf( '%s%s.gif', opt.vd_dr, fg_nm ) ;
  img_2 = frame2im( img );
  [ imind, cm ] = rgb2ind(img_2, 256 ) ;
    if i_lp == 1
    imwrite( imind, cm, fg_nm_vd, 'gif', 'Loopcount', inf ) ;
    else
    imwrite( imind, cm, fg_nm_vd, 'gif', 'WriteMode', 'append' ) ;
    end
  else
  imwrite( img.cdata, [ opt.fg_dr fg_nm '.' opt.frmt_fg ] ) ;
  end % opt.vd.do
 
% Checking if there's some work to do
function do_img = do_img_f( fl_dt_out, opt )
do_img = 1 ; % By default, do the imaging
  if exist( fl_dt_out, 'file' ) == 2
    if ( opt.redo )
    delete( fl_dt_out  ) ;
    else
    do_img = 0 ;
    disp( sprintf( 'File %s exists. Not imaging the scene again through Starshade', fl_dt_out ) )
    end
  end
