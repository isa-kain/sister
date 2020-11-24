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
function psf_out = psf_resize_odd( psf_in, fct, check )
% Script to reduce the size of a psf using imresize (with the bilinear option) for images with a well defined center (odd number of pixels in and out)
% Author: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com

% Generic debug mode
dbstop if error

% Checking that the max is centered
  if ~exist( 'check', 'var' )
  chck = 0 ;
  else
  chck = 1 ;
  end

% Checking the input image has an odd number of pixels
sz_in_1 = size( psf_in, 1 ) ;
sz_in_2 = size( psf_in, 2 ) ;
  if  2 * floor( ( sz_in_1 - 1 ) / 2 ) + 1 ~= sz_in_1
  disp( 'The input array has an even number of elements. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

  if  2 * floor( ( sz_in_2 - 1 ) / 2 ) + 1 ~= sz_in_2
  disp( 'The input array has an even number of elements. Stopped' )
  disp( ' ' )
  disp( 'PS: exit the debug mode typing dbquit all' ) ; make_a_stop
  end

cntr_in_1 = ( sz_in_1 + 1 ) / 2 ;
cntr_in_2 = ( sz_in_2 + 1 ) / 2 ;

  if ( chck )
    if ( psf_in( cntr_in_1, cntr_in_2 ) == max( psf_in( : ) ) )
    disp( 'The maximum of the input PSF is centered. Good.' )
    else
    disp( 'The maximum of the input PSF is not centered. Not good if it is in the stationary region. Continuing.' )
    end
  end

% Number of pixels of the output image (imresize uses ceil)
n_px_out_1 = ceil( sz_in_1 * fct ) ;
n_px_out_2 = ceil( sz_in_2 * fct ) ;

% If even, make it odd
  if 2 * floor( n_px_out_1 / 2 ) == n_px_out_1
  n_px_out_1 = n_px_out_1 + 1 ;
  end

  if 2 * floor( n_px_out_2 / 2 ) == n_px_out_2
  n_px_out_2 = n_px_out_2 + 1 ;
  end

% Compensating for the imsize change

psf_out = imresize( psf_in, [ n_px_out_1 n_px_out_2 ], 'bilinear' ) * ( sz_in_1 / n_px_out_1 ) * ( sz_in_2 / n_px_out_2 ) ;

  if ( chck )
    if ( psf_out( ( n_px_out_1 + 1 ) / 2, ( n_px_out_2 + 1 ) / 2 ) == max( psf_out( : ) ) )
    disp( 'The maximum of the output PSF is centered. Good.' )
    else
    disp( 'The maximum of the output PSF is not centered. Not good if it is in the stationary region. Continuing.' )
    end
  end
