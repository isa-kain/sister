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

function UTotL = erase_hot_pixels( UTotL, iter, mssg )
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

% Number of times to repeat the search for hot pixels
  if ~exist( 'iter', 'var' )
  iter1 = 1 ;
  else
  iter1 = iter ;
  end
% Whether to display messages
  if ~exist( 'mssg', 'var' )
  mssg = 0 ;
  end

% Do the first iteration, but if there are no hot pixels, do not repeat the follwing ones
n_q_px = 1 ; % dummy, to start up 

  while iter1 && ( n_q_px )
  % Find hot pixels
  n_sgm = 7 ; % times the std to consider a pixel an outlier (i.e., hot pixel)
  l = 1 ; % pix around the hot pixel to compute an average
  n_lmbd = size( UTotL, 3 ) ;
  sz_ppl = size( UTotL, 1 ) ;
    for i_lmbd = 1 : n_lmbd
    UTotL_r = real( squeeze( UTotL( :, :, i_lmbd ) ) ) ;
    DUTotL_r = UTotL_r - medfilt2( UTotL_r ) ;
    % Remove edges due to median filter algorithm
    DUTotL_r( 1 : 2, : ) = NaN ; DUTotL_r( end - 2 : end, : ) = NaN ;
    DUTotL_r( :, 1 : 2 ) = NaN ; DUTotL_r( :, end - 2 : end ) = NaN ;
    q_px_r = find( abs( DUTotL_r( : ) ) > ( nanmean( DUTotL_r( : ) ) + n_sgm * nanstd( DUTotL_r( : ) ) ) ) ;
    UTotL_i = imag( squeeze( UTotL( :, :, i_lmbd ) ) ) ;
    DUTotL_i = UTotL_i - medfilt2( UTotL_i ) ;
    % Remove edges due to median filter algorithm
    DUTotL_i( 1 : 3, : ) = NaN ; DUTotL_i( end -3 : end, : ) = NaN ;
    DUTotL_i( :, 1 : 3 ) = NaN ; DUTotL_i( :, end - 3 : end ) = NaN ;
    q_px_i = find( abs( DUTotL_i( : ) ) > ( nanmean( DUTotL_i( : ) ) + n_sgm * nanstd( DUTotL_i( : ) ) ) ) ;
    % Joining both sets of hot pixels (some times they are more visible in the real or imaginary part only)
    q_px = unique( cat( 1, q_px_r, q_px_i ) ) ;
    n_q_px = numel( q_px ) ;
      if n_q_px == 0
       if ( mssg )
       disp( '(ERASE_HOT_PIXELS) No hot pixels in the electric field. Done.' )
       end
      else
        % Subsitute the hot pixel by an average of the adjacent pixels
        if ( mssg )
        disp( sprintf( '(ERASE_HOT_PIXELS) %i hot pixels erased in the electric field. Iteration %i/%i', n_q_px, iter - iter1 + 1, iter ) )
        end
        for i_px = 1 : n_q_px
        y_tmp = ceil( q_px( i_px ) / sz_ppl ) ;
        x_tmp = mod( q_px( i_px ), sz_ppl ) ;
        UTotL_r( x_tmp, y_tmp ) = NaN ;
        UTotL_r( x_tmp, y_tmp ) = mean( nanmean( UTotL_r( x_tmp - l : x_tmp + l, y_tmp - l : y_tmp + l ) ) ) ;
        UTotL_i( x_tmp, y_tmp ) = NaN ;
        UTotL_i( x_tmp, y_tmp ) = mean( nanmean( UTotL_i( x_tmp - l : x_tmp + l, y_tmp - l : y_tmp + l ) ) ) ;
        end
      end % n_q_px
    % Defining the exit result
    UTotL( :, :, i_lmbd ) = UTotL_r + i * UTotL_i ;
    end % i_lmbd
  iter1 = iter1 - 1 ;
  end % iter

