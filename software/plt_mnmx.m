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
function [ min_plt max_plt ] = plt_mnmx( data, opt )
% Function to help plot 2-dimensional images
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

colormap( bone )

  if ~exist( 'opt', 'var' )
  opt = [ ] ;
  end

opt_plt = get_default_options( opt ) ;

% Centering the data
sz_1 = size( data, 1 ) ;
sz_2 = size( data, 2 ) ;
cntr_1 = round( ( sz_1 + 1 ) / 2 ) ;
cntr_2 = round( ( sz_2 + 1 ) / 2 ) ;
l_1 = floor( sz_1 / 2 / opt_plt.zm ) ;
l_2 = floor( sz_2 / 2 / opt_plt.zm ) ;
% For data that have an even number of elements and zm = 1
idx_tp_1 = cntr_1 + l_1 ;
  if 2 * floor( sz_1 / 2 ) == sz_1
  idx_tp_1 = cntr_1 + l_1 - 1 ;
  end
idx_tp_2 = cntr_2 + l_2 ;
  if 2 * floor( sz_2 / 2 ) == sz_2
  idx_tp_2 = cntr_2 + l_2 - 1 ;
  end

data = data( cntr_1 - l_1 : idx_tp_1, cntr_2 - l_2 : idx_tp_2 ) ;
q = ~isnan( data ) ;
  if ~isfield( opt_plt, 'prcntl' ) && ~isfield( opt_plt, 'mnmx' )
    if ( opt_plt.log_10 )
    % Some suggestion
    q = find( data > 0 ) ;
      if ~( numel( q ) )
      disp( 'WARNING: there are no positive values to plot with log10. Returning.' )
      return
      end
    std_log = nanstd( log10( data( q ) ) ) ;
    mean_log = nanmean( log10( data( q ) ) ) ;
    min_plt = mean_log - std_log ;
    max_plt = mean_log + 8 * std_log ;
    else
      if ~( numel( q ) )
      disp( 'WARNING: all data are NaN. Nothing to plot. Returning.' )
      return
      end
    std_lnr = nanstd( data( q ) ) ;
    mean_lnr = nanmean( data( q ) ) ;
    min_plt = mean_lnr - 3 * std_lnr ;
    max_plt = mean_lnr + 3 * std_lnr ;
    end
  end
  if isfield( opt_plt, 'prcntl' )
  min_plt = prctile( data( q ), opt_plt.prcntl( 1 ) ) ;
  max_plt = prctile( data( q ), opt_plt.prcntl( 2 ) ) ;
  end
  if isfield( opt_plt, 'mnmx' )
  min_plt = opt_plt.mnmx( 1 ) ;
  max_plt = opt_plt.mnmx( 2 ) ;
  end
% If the min and max are too close
  if numel( abs( min_plt - max_plt ) / max( data( q ) ) < 1e-3 ) && ~isfield( opt_plt, 'mnmx' )
    if ~( opt_plt.log_10 )
    min_plt = min_plt / 2 ;
    max_plt = max_plt * 2 ; 
    else
    min_plt = min_plt - log10( 2 ) ;
    max_plt = max_plt + log10( 2 ) ;
    end
  end

% Removing wisely the 0s and negative data
%q_0 = find( data <= 0 ) ;
%data( q_0 ) = max(max( data ) ) / 1e10 ;
% axis
hlf_1 = ( size( data, 1 ) - 1 ) / 2 ;
hlf_2 = ( size( data, 2 ) - 1 ) / 2 ;
x_ax_arry = ( -hlf_1 : hlf_1 ) * opt_plt.px_unt ;
y_ax_arry = ( -hlf_2 : hlf_2 ) * opt_plt.px_unt ;

  if ( opt_plt.log_10 )
  % Multiplying by opt_plt.log_10 allows to plot -log10(data), that is, to revert the color scale in log plots
  imagesc( x_ax_arry, y_ax_arry, opt_plt.log_10 * log10( ( data ) ), [ min_plt max_plt ] ) 
  else
  imagesc( x_ax_arry, y_ax_arry, data, [ min_plt max_plt ] ) 
  end
% Set color of the figure box to white, instead of gray
set( gcf, 'color', 'w' ) ;
xlim( sort( [ x_ax_arry( 1 ), x_ax_arry( end ) ] ) )
ylim( sort( [ y_ax_arry( 1 ), y_ax_arry( end ) ] ) ) 
% Orientation: default N down/E right. Othersie N up and E left.
  if ~( opt_plt.noflip )
  set( gca, 'YDir', 'reverse' )
  else
  set( gca, 'XDir', 'reverse' )
  set( gca, 'YDir', 'normal' )
  end
  if ( opt_plt.clrbr )
  h_clrbr = colorbar ;
    if isfield( opt_plt, 'clrbr_lbl' )
    ylabel( h_clrbr, opt_plt.clrbr_lbl ) ;
      if isfield( opt_plt, 'clrbr_lbl_sz' )
      ylabel( h_clrbr, opt_plt.clrbr_lbl, 'FontSize', opt_plt.clrbr_lbl_sz ) ;
      end
    end
  end
  
  if ( opt_plt.pbaspect )
  pbaspect( [ 1 1 1 ] )
  end

% Sub-function
function opt = get_default_options( opt )

  if ~isfield( opt, 'zm' )
  opt.zm = 1 ;
  end

  if ~isfield( opt, 'log_10' )
  opt.log_10 = 0 ;
  end

  if ~isfield( opt, 'px_unt' )
  opt.px_unt = 1 ;
  end

  if ~isfield( opt, 'clrbr' )
  opt.clrbr = 1 ;
  end

  if ~isfield( opt, 'noflip' )
  opt.noflip = 0 ;
  end

  if ~isfield( opt, 'pbaspect' )
  opt.pbaspect = 0 ;
  end


