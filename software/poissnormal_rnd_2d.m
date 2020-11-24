function [ rms r ] = poissnormal_rnd_2d(lambda,lambda_normal)
%
%   Based upon:
%   POISSRND Random matrices from Poisson distribution.
%   R = POISSRND(LAMBDA) returns a matrix of random numbers chosen   
%   from the Poisson distribution with parameter LAMBDA.
%
%   LAMBDA_NORMAL can be 10 or higher
%   R = POISSRND(LAMBDA,M,N,LAMBDA_NORMAL) 
%
%   POISSRND uses a waiting time method.
%   References:
%      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%      Springer-Verlag, 1986 page 504.

%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.6 $  $Date: 1997/11/29 01:46:24 $
%
%   Modified by: Sergi R. Hildebrandt, srh.jpl.caltech@gmail.com. 2-dimensional array, control when to approximate by a normal distribution (with continuity correction), and return the dispersion (first output variable), not the average+dispersion, as the original one did.

% SRH: generic debugging mode
dbstop if error

if nargin <  1, 
    error('Requires at least one input argument.'); 
end

rows = size( lambda, 1 ) ;
columns = size( lambda, 2 ) ;
% SRH: default value to 20, as Stuart Shaklan had
if ~exist( 'lambda_normal', 'var' )
    lambda_normal = 20 ;
end  

%Initialize r to zero. Anything that does not have a positive lambda will be whichever value r has in its initialization (NaN gives issues for images, chosing 0s: no electron, no noise)
r = zeros( rows, columns ) ;
% The dispersion will be stored here
rms = r ;

% SRH modification (continuity correction +0.5) 
q_normal = find( lambda >= lambda_normal ) ;
n_normal = numel( q_normal ) ;
  if n_normal
  lambda_normal_array = lambda( q_normal ) ;
  lambda_normal_array = lambda_normal_array + 0.5 ;
  r( q_normal ) =round( lambda_normal_array + sqrt( lambda_normal_array ).* randn( n_normal, 1 ) ) ; 
  rms( q_normal ) = r( q_normal ) - lambda_normal_array ;
  end % If there are data that will be approximated by a Normal distribution

q_poisson_1 = find( lambda < lambda_normal ) ;
q_poisson_2 = find( lambda( q_poisson_1 ) > 0 ) ;
q_poisson = q_poisson_1( q_poisson_2 ) ;
n_poisson = numel( q_poisson ) ;
  if n_poisson
  lambda_poisson = lambda( q_poisson ) ;
  p = zeros( n_poisson, 1 ) ;
  r_tmp = p ;
  done = ones( n_poisson, 1 ) ;

    while any( any( done ) ) ~= 0
    p = p - log( rand( n_poisson, 1 ) ) ;
    kc = [ 1 : n_poisson ]' ;
    k = find( p < lambda_poisson) ;
      if any( k )
      r_tmp( k ) = r_tmp( k ) + 1 ;
      end
    kc( k ) = [ ] ;
    done( kc ) = zeros( size( kc ) ) ;
    end
  % SRh: Putting the results back into the array
  r( q_poisson ) = r_tmp ;
  % SRH: Returning the dispersion
  rms( q_poisson ) = r_tmp - lambda_poisson ;
  end % Data that are dealt with a Poisson distribution


