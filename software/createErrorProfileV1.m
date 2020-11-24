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

function vecPetalArray = createErrorProfileV1(r, Profile, occulterDiameter, petalLength, numPetals, errorArray)
% Creates the petal perturbation corresponding to an array of static errors
% Author: Eric J. Cady. eric.j.cady@jpl.nasa.gov
% Terms carried in error budget, from Luis Marchen, 10/27/11
% (Below is direct email cut-and-paste)
% x 1.	Proportional width
% x 2.	Tip clip
% x 3.	Radial shift
% x 4.	Quadratic out-of-plane bend
% x 5.	Lateral shift
% x 6.	In plane linear bend
% x 7.	Segment rigid body (translations and rotations)
% x 8.	Tip bends radial (this you have already done)
% x 9.	Tip bends lateral (this you have already done)
% x 10.	Tip bends clocking
% x 11.	Ellipticity (elliptical truss deformation, this is one of the errors you provided to us this past summer)
% x 12.	Spatial frequency (this is what Stuart talked about yesterday) sine and cosine on segments

% Additional optional terms carried below:
% x In-plane quadratic bend
% x sine and cosine on full edges
% x petal angular and clocking shifts
% x translation and rotation of segments
% x quadratic bend of segments

%-------------------------------------------------------------

% Load up useful intermediate variables
maxR = occulterDiameter/2;
a = maxR - petalLength;

% Start with single petal
[petalX, petalY] = makePolygon(Profile, r, maxR, numPetals, 0, a, 3); % Get one petal edge
petalX = petalX(:).';
petalY = petalY(:).';
petalR = sqrt(petalX.^2 + petalY.^2);
petalT = atan2(petalY, petalX);
petalZ = zeros(size(petalX));

deltaArray = repmat([zeros(size(petalX));zeros(size(petalY))], 2*numPetals, 1);
zArray = repmat([zeros(size(petalZ))], 2*numPetals, 1);

if nargin < 6 || isempty(errorArray)
    % We're building this without errors.  Probably mostly for testing,
    % but who knows.
    
    disp('No errors, building edge from supplied profile.')
    
    vecPetalArray = cell(1, numPetals);
    for k = 1:numPetals
        tmp = [petalX  fliplr(petalX); -petalY fliplr(petalY); petalZ fliplr(petalZ)];
        tmp = [cos(2*pi/numPetals*(k-1)) -sin(2*pi/numPetals*(k-1)) 0; sin(2*pi/numPetals*(k-1)) cos(2*pi/numPetals*(k-1)) 0; 0 0 1]*tmp;
        
        vecPetalArray{k} = {tmp};
    end
    
else
    % Do errors in a specific order! Clip first, then RMS, then prop width,
    % then bending, then deployment
    
    % Future speedup possibility: slice elements out of array as used.  Not
    % clear it's worth the time to write it unless there is a lot of errors
    % applied.
    
    % Tip clip
    for j = 1:length(errorArray)
        switch errorArray{j}.name
            case 'tipclip'
                pn = errorArray{j}.number;
                if pn == 0 % global error
                    deltaArray(:, (maxR - petalR < errorArray{j}.tipclip)) = NaN;
                else
                    deltaArray((4*pn-3):4*pn, (maxR - petalR < errorArray{j}.tipclip)) = NaN;
                end
        end
    end
    
    % Bend specific segments
    for j = 1:length(errorArray)
        switch errorArray{j}.name
            case 'segbend'
                pn = errorArray{j}.number;
                pv = errorArray{j}.pv;
                dec = errorArray{j}.dec;
                
                if pn == 0 % global error
                    for k = 1:numPetals
                        if errorArray{j}.edge == 1
                            secX = petalX(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*k-3, errorArray{j}.start:1:errorArray{j}.end);
                            secY = petalY(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*k-2, errorArray{j}.start:1:errorArray{j}.end);
                            
                            neg1to1 = linspace(-1, 1, length(errorArray{j}.start:1:errorArray{j}.end));
                            n1to1dist = sqrt((secX(end) - secX(1))^2 + (secY(end) - secY(1))^2);
                            outY = pv*((neg1to1 - dec/(n1to1dist/2)).^2 - 1/2); % positive on left
                            ang = atan2(secY(end) - secY(1), secX(end) - secX(1));
                            
                            tmpx = secX*cos(-ang) - secY*sin(-ang);
                            tmpy = secX*sin(-ang) + secY*cos(-ang);
                            
                            tmpx2 = tmpx;
                            tmpy2 = tmpy + outY;
                            
                            tmpx3 = tmpx2*cos(ang) - tmpy2*sin(ang);
                            tmpy3 = tmpx2*sin(ang) + tmpy2*cos(ang);
                            
                            deltaArray(4*k-3, errorArray{j}.start:1:errorArray{j}.end) = tmpx3 - petalX(errorArray{j}.start:1:errorArray{j}.end);
                            deltaArray(4*k-2, errorArray{j}.start:1:errorArray{j}.end) = tmpy3 - petalY(errorArray{j}.start:1:errorArray{j}.end);
                        elseif errorArray{j}.edge == 2
                            secX = petalX(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*k-1, errorArray{j}.start:1:errorArray{j}.end);
                            secY = petalY(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*k-0, errorArray{j}.start:1:errorArray{j}.end);
                            
                            neg1to1 = linspace(-1, 1, length(errorArray{j}.start:1:errorArray{j}.end));
                            n1to1dist = sqrt((secX(end) - secX(1))^2 + (secY(end) - secY(1))^2);
                            outY = -pv*((neg1to1 - dec/(n1to1dist/2)).^2 - 1/2); % negative on right
                            ang = atan2(secY(end) - secY(1), secX(end) - secX(1));
                            
                            tmpx = secX*cos(-ang) - secY*sin(-ang);
                            tmpy = secX*sin(-ang) + secY*cos(-ang);
                            
                            tmpx2 = tmpx;
                            tmpy2 = tmpy + outY;
                            
                            tmpx3 = tmpx2*cos(ang) - tmpy2*sin(ang);
                            tmpy3 = tmpx2*sin(ang) + tmpy2*cos(ang);
                            
                            deltaArray(4*k-1, errorArray{j}.start:1:errorArray{j}.end) = tmpx3 - petalX(errorArray{j}.start:1:errorArray{j}.end);
                            deltaArray(4*k-0, errorArray{j}.start:1:errorArray{j}.end) = tmpy3 - petalY(errorArray{j}.start:1:errorArray{j}.end);
                        else
                            error('Edge must be 1 (left) or 2 (right).')
                        end
                    end
                else
                    if errorArray{j}.edge == 1
                        secX = petalX(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*pn-3, errorArray{j}.start:1:errorArray{j}.end);
                        secY = petalY(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*pn-2, errorArray{j}.start:1:errorArray{j}.end);
                        
                        neg1to1 = linspace(-1, 1, length(errorArray{j}.start:1:errorArray{j}.end));
                        n1to1dist = sqrt((secX(end) - secX(1))^2 + (secY(end) - secY(1))^2);
                        outY = pv*((neg1to1 - dec/(n1to1dist/2)).^2 - 1/2); % positive on left
                        ang = atan2(secY(end) - secY(1), secX(end) - secX(1));
                        
                        tmpx = secX*cos(-ang) - secY*sin(-ang);
                        tmpy = secX*sin(-ang) + secY*cos(-ang);
                        
                        tmpx2 = tmpx;
                        tmpy2 = tmpy + outY;
                        
                        tmpx3 = tmpx2*cos(ang) - tmpy2*sin(ang);
                        tmpy3 = tmpx2*sin(ang) + tmpy2*cos(ang);
                        
                        deltaArray(4*pn-3, errorArray{j}.start:1:errorArray{j}.end) = tmpx3 - petalX(errorArray{j}.start:1:errorArray{j}.end);
                        deltaArray(4*pn-2, errorArray{j}.start:1:errorArray{j}.end) = tmpy3 - petalY(errorArray{j}.start:1:errorArray{j}.end);
                    elseif errorArray{j}.edge == 2
                        secX = petalX(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*pn-1, errorArray{j}.start:1:errorArray{j}.end);
                        secY = petalY(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*pn-0, errorArray{j}.start:1:errorArray{j}.end);
                        
                        neg1to1 = linspace(-1, 1, length(errorArray{j}.start:1:errorArray{j}.end));
                        n1to1dist = sqrt((secX(end) - secX(1))^2 + (secY(end) - secY(1))^2);
                        outY = pv*((neg1to1 - dec/(n1to1dist/2)).^2 - 1/2); % negative on right
                        ang = atan2(secY(end) - secY(1), secX(end) - secX(1));
                        
                        tmpx = secX*cos(-ang) - secY*sin(-ang);
                        tmpy = secX*sin(-ang) + secY*cos(-ang);
                        
                        tmpx2 = tmpx;
                        tmpy2 = tmpy + outY;
                        
                        tmpx3 = tmpx2*cos(ang) - tmpy2*sin(ang);
                        tmpy3 = tmpx2*sin(ang) + tmpy2*cos(ang);
                        
                        deltaArray(4*pn-1, errorArray{j}.start:1:errorArray{j}.end) = tmpx3 - petalX(errorArray{j}.start:1:errorArray{j}.end);
                        deltaArray(4*pn-0, errorArray{j}.start:1:errorArray{j}.end) = tmpy3 - petalY(errorArray{j}.start:1:errorArray{j}.end);
                    else
                        error('Edge must be 1 (left) or 2 (right).')
                    end
                end
        end
    end
    
    % Random sinusoids on specific segments
    for j = 1:length(errorArray)
        switch errorArray{j}.name
            case 'segrandsin'
                if (errorArray{j}.seed < 0)
                    seed = sum(100*clock);
                else
                    seed = errorArray{j}.seed;
                end
                rng(seed);
                
                zeroToOne = linspace(0, 1, length(errorArray{j}.start:1:errorArray{j}.end));
                outY = zeros(size(errorArray{j}.start:1:errorArray{j}.end));
                fn = errorArray{j}.freq;
                A = errorArray{j}.amp;
                
                for q = 1:length(fn)
                    outY = outY + A(q)*sin(2*pi*(zeroToOne*fn(q) + rand(1,1)));
                end
                pn = errorArray{j}.number;
                if pn == 0 % global error
                    for k = 1:numPetals
                        deltaArray(4*k-2, errorArray{j}.start:1:errorArray{j}.end) = deltaArray(4*k-2, errorArray{j}.start:1:errorArray{j}.end) + outY;
                        deltaArray(4*k-0, errorArray{j}.start:1:errorArray{j}.end) = deltaArray(4*k-0, errorArray{j}.start:1:errorArray{j}.end) - outY;
                    end
                else
                    if errorArray{j}.edge == 1
                        deltaArray(4*pn-2, errorArray{j}.start:1:errorArray{j}.end) = deltaArray(4*pn-2, errorArray{j}.start:1:errorArray{j}.end) + outY;
                    elseif errorArray{j}.edge == 2
                        deltaArray(4*pn-0, errorArray{j}.start:1:errorArray{j}.end) = deltaArray(4*pn-0, errorArray{j}.start:1:errorArray{j}.end) - outY;
                    else
                        error('Edge must be 1 (left) or 2 (right).')
                    end
                end
        end
    end
    
    % Specific sinusoids on specific segments
    for j = 1:length(errorArray)
        switch errorArray{j}.name
            case 'segsin'
                zeroToOne = linspace(0, 1, length(errorArray{j}.start:1:errorArray{j}.end));
                fn = errorArray{j}.freq;
                A = errorArray{j}.amp;
                p = errorArray{j}.phase;
                outY = A*sin(2*pi*(zeroToOne*fn) + p);
                
                pn = errorArray{j}.number;
                if pn == 0 % global error
                    for k = 1:numPetals
                        deltaArray(4*k-2, errorArray{j}.start:1:errorArray{j}.end) = deltaArray(4*k-2, errorArray{j}.start:1:errorArray{j}.end) + outY;
                        deltaArray(4*k-0, errorArray{j}.start:1:errorArray{j}.end) = deltaArray(4*k-0, errorArray{j}.start:1:errorArray{j}.end) - outY;
                    end
                else
                    if errorArray{j}.edge == 1
                        deltaArray(4*pn-2, errorArray{j}.start:1:errorArray{j}.end) = deltaArray(4*pn-2, errorArray{j}.start:1:errorArray{j}.end) + outY;
                    elseif errorArray{j}.edge == 2
                        deltaArray(4*pn-0, errorArray{j}.start:1:errorArray{j}.end) = deltaArray(4*pn-0, errorArray{j}.start:1:errorArray{j}.end) - outY;
                    else
                        error('Edge must be 1 (left) or 2 (right).')
                    end
                end
        end
    end
    
    % Translate and rotate specific segments
    for j = 1:length(errorArray)
        switch errorArray{j}.name
            case 'segtr'
                pn = errorArray{j}.number;
                if pn == 0 % global error
                    for k = 1:numPetals
                        secX = petalX(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*k-3, errorArray{j}.start:1:errorArray{j}.end);
                        secY = petalY(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*k-2, errorArray{j}.start:1:errorArray{j}.end);
                        deltaArray(4*k-3, errorArray{j}.start:1:errorArray{j}.end) = ((secX-mean(secX))*cos(errorArray{j}.dtheta) ...
                            - (secY-mean(secY))*sin(errorArray{j}.dtheta) - errorArray{j}.dx + mean(secX)) - petalX(errorArray{j}.start:1:errorArray{j}.end);
                        deltaArray(4*k-2, errorArray{j}.start:1:errorArray{j}.end) = ((secX-mean(secX))*sin(errorArray{j}.dtheta) ...
                            + (secY-mean(secY))*cos(errorArray{j}.dtheta) - errorArray{j}.dy + mean(secY)) - petalY(errorArray{j}.start:1:errorArray{j}.end);
                        
                        secX = petalX(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*k-1, errorArray{j}.start:1:errorArray{j}.end);
                        secY = -petalY(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*k-0, errorArray{j}.start:1:errorArray{j}.end);
                        deltaArray(4*k-1, errorArray{j}.start:1:errorArray{j}.end) = ((secX-mean(secX))*cos(errorArray{j}.dtheta) ...
                            - (secY-mean(secY))*sin(errorArray{j}.dtheta) - errorArray{j}.dx + mean(secX)) - petalX(errorArray{j}.start:1:errorArray{j}.end);
                        deltaArray(4*k-0, errorArray{j}.start:1:errorArray{j}.end) = ((secX-mean(secX))*sin(errorArray{j}.dtheta) ...
                            + (secY-mean(secY))*cos(errorArray{j}.dtheta) - errorArray{j}.dy + mean(secY)) + petalY(errorArray{j}.start:1:errorArray{j}.end);
                    end
                else
                    if errorArray{j}.edge == 1
                        secX = petalX(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*pn-3, errorArray{j}.start:1:errorArray{j}.end);
                        secY = petalY(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*pn-2, errorArray{j}.start:1:errorArray{j}.end);
                        deltaArray(4*pn-3, errorArray{j}.start:1:errorArray{j}.end) = ((secX-mean(secX))*cos(errorArray{j}.dtheta) ...
                            - (secY-mean(secY))*sin(errorArray{j}.dtheta) - errorArray{j}.dx + mean(secX)) - petalX(errorArray{j}.start:1:errorArray{j}.end);
                        deltaArray(4*pn-2, errorArray{j}.start:1:errorArray{j}.end) = ((secX-mean(secX))*sin(errorArray{j}.dtheta) ...
                            + (secY-mean(secY))*cos(errorArray{j}.dtheta) - errorArray{j}.dy + mean(secY)) - petalY(errorArray{j}.start:1:errorArray{j}.end);
                    elseif errorArray{j}.edge == 2
                        secX = petalX(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*pn-1, errorArray{j}.start:1:errorArray{j}.end);
                        secY = -petalY(errorArray{j}.start:1:errorArray{j}.end) + deltaArray(4*pn-0, errorArray{j}.start:1:errorArray{j}.end);
                        deltaArray(4*pn-1, errorArray{j}.start:1:errorArray{j}.end) = ((secX-mean(secX))*cos(errorArray{j}.dtheta) ...
                            - (secY-mean(secY))*sin(errorArray{j}.dtheta) - errorArray{j}.dx + mean(secX)) - petalX(errorArray{j}.start:1:errorArray{j}.end);
                        deltaArray(4*pn-0, errorArray{j}.start:1:errorArray{j}.end) = ((secX-mean(secX))*sin(errorArray{j}.dtheta) ...
                            + (secY-mean(secY))*cos(errorArray{j}.dtheta) - errorArray{j}.dy + mean(secY)) + petalY(errorArray{j}.start:1:errorArray{j}.end);
                    else
                        error('Edge must be 1 (left) or 2 (right).')
                    end
                end
        end
    end
    
    % RMS
    for j = 1:length(errorArray)
        switch errorArray{j}.name
            case 'symrms'
                if (errorArray{j}.seed < 0)
                    seed = sum(100*clock);
                else
                    seed = errorArray{j}.seed;
                end
                
                outY = rms1D(length(petalR), errorArray{j}.rms, errorArray{j}.cutoff, errorArray{j}.powerlaw, seed);
                pn = errorArray{j}.number;
                if pn == 0 % global error
                    for k = 1:numPetals
                        deltaArray(4*k-2, :) = deltaArray(4*k-2, :) + outY;
                        deltaArray(4*k-0, :) = deltaArray(4*k-0, :) - outY;
                    end
                else
                    deltaArray(4*pn-2, :) = deltaArray(4*pn-2, :) + outY;
                    deltaArray(4*pn-0, :) = deltaArray(4*pn-0, :) - outY;
                end
            case 'antisymrms'
                if (errorArray{j}.seed < 0)
                    seed = sum(100*clock);
                else
                    seed = errorArray{j}.seed;
                end
                
                outY = rms1D(length(petalR), errorArray{j}.rms, errorArray{j}.cutoff, errorArray{j}.powerlaw, seed);
                pn = errorArray{j}.number;
                if pn == 0 % global error
                    for k = 1:numPetals
                        deltaArray(4*k-2, :) = deltaArray(4*k-2, :) + outY;
                        deltaArray(4*k-0, :) = deltaArray(4*k-0, :) + outY;
                    end
                else
                    deltaArray(4*pn-2, :) = deltaArray(4*pn-2, :) + outY;
                    deltaArray(4*pn-0, :) = deltaArray(4*pn-0, :) + outY;
                end
        end
    end
    
    % Prop width
    for j = 1:length(errorArray)
        switch errorArray{j}.name
            case 'propwidth'
                pn = errorArray{j}.number;
                if pn == 0 % global error
                    for k = 1:numPetals
                        deltaArray(4*k-2, (petalR <= errorArray{j}.cutoff)) = deltaArray(4*k-2, (petalR <= errorArray{j}.cutoff)) ...
                            + errorArray{j}.propwidth*(petalY + deltaArray(4*k-2, (petalR <= errorArray{j}.cutoff)));
                        deltaArray(4*k-0, (petalR <= errorArray{j}.cutoff)) = deltaArray(4*k-0, (petalR <= errorArray{j}.cutoff)) ...
                            - errorArray{j}.propwidth*(petalY + deltaArray(4*k-0, (petalR <= errorArray{j}.cutoff)));
                    end
                else
                    deltaArray(4*pn-2, (petalR <= errorArray{j}.cutoff)) = deltaArray(4*pn-2, (petalR <= errorArray{j}.cutoff)) ...
                        + errorArray{j}.propwidth*(petalY(petalR <= errorArray{j}.cutoff) + deltaArray(4*pn-2, (petalR <= errorArray{j}.cutoff)));
                    deltaArray(4*pn-0, (petalR <= errorArray{j}.cutoff)) = deltaArray(4*pn-0, (petalR <= errorArray{j}.cutoff)) ...
                        - errorArray{j}.propwidth*(petalY(petalR <= errorArray{j}.cutoff) + deltaArray(4*pn-0, (petalR <= errorArray{j}.cutoff)));
                end
        end
    end
    
    % Tip radius/lateral/angular/clocking
    for j = 1:length(errorArray)
        switch errorArray{j}.name
            case 'tipbend'
                pn = errorArray{j}.number;
                tipBeginsAt = errorArray{j}.tipbegin;
                
                if (pn == 0)
                    for k = 1:numPetals % global
                        % One side of petal...
                        tmpX1full = petalX + deltaArray(4*k-3, :);
                        tmpY1full = petalY + deltaArray(4*k-2, :);
                        rfull = sqrt(tmpX1full.*tmpX1full + tmpY1full.*tmpY1full);
                        tmpX1 = tmpX1full(rfull >= tipBeginsAt);
                        tmpY1 = tmpY1full(rfull >= tipBeginsAt);
                        
                        newX = tmpX1 + errorArray{j}.radial;
                        newY = tmpY1 + errorArray{j}.lateral;
                        tmpX1 = cos(errorArray{j}.clocking)*(newX - tipBeginsAt - errorArray{j}.radial) + ...
                            sin(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + tipBeginsAt + errorArray{j}.radial;
                        tmpY1 = -sin(errorArray{j}.clocking)*(newX - tipBeginsAt - errorArray{j}.radial) + ...
                            cos(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + errorArray{j}.lateral;
                        newX = tmpX1;
                        newY = tmpY1;
                        
                        deltaArray(4*k-3, rfull >= tipBeginsAt) = newX - petalX(rfull >= tipBeginsAt);
                        deltaArray(4*k-2, rfull >= tipBeginsAt) = newY - petalY(rfull >= tipBeginsAt);
                        
                        % ...then the other.
                        tmpX1full = petalX + deltaArray(4*k-1, :);
                        tmpY1full = -petalY + deltaArray(4*k-0, :);
                        rfull = sqrt(tmpX1full.*tmpX1full + tmpY1full.*tmpY1full);
                        tmpX1 = tmpX1full(rfull >= tipBeginsAt);
                        tmpY1 = tmpY1full(rfull >= tipBeginsAt);
                        
                        newX = tmpX1 + errorArray{j}.radial;
                        newY = tmpY1 + errorArray{j}.lateral;
                        tmpX1 = cos(errorArray{j}.clocking)*(newX - tipBeginsAt - errorArray{j}.radial) + ...
                            sin(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + tipBeginsAt + errorArray{j}.radial;
                        tmpY1 = -sin(errorArray{j}.clocking)*(newX - tipBeginsAt - errorArray{j}.radial) + ...
                            cos(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + errorArray{j}.lateral;
                        newX = tmpX1;
                        newY = tmpY1;
                        
                        deltaArray(4*k-1, rfull >= tipBeginsAt) = newX - petalX(rfull >= tipBeginsAt);
                        deltaArray(4*k-0, rfull >= tipBeginsAt) = newY + petalY(rfull >= tipBeginsAt);
                    end
                else
                    % One side of petal...
                    tmpX1full = petalX + deltaArray(4*pn-3, :);
                    tmpY1full = petalY + deltaArray(4*pn-2, :);
                    rfull = sqrt(tmpX1full.*tmpX1full + tmpY1full.*tmpY1full);
                    tmpX1 = tmpX1full(rfull >= tipBeginsAt);
                    tmpY1 = tmpY1full(rfull >= tipBeginsAt);
                    
                    newX = tmpX1 + errorArray{j}.radial;
                    newY = tmpY1 + errorArray{j}.lateral;
                    tmpX1 = cos(errorArray{j}.clocking)*(newX - tipBeginsAt - errorArray{j}.radial) + ...
                        sin(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + tipBeginsAt + errorArray{j}.radial;
                    tmpY1 = -sin(errorArray{j}.clocking)*(newX - tipBeginsAt - errorArray{j}.radial) + ...
                        cos(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + errorArray{j}.lateral;
                    newX = tmpX1;
                    newY = tmpY1;
                    
                    deltaArray(4*pn-3, rfull >= tipBeginsAt) = newX - petalX(rfull >= tipBeginsAt);
                    deltaArray(4*pn-2, rfull >= tipBeginsAt) = newY - petalY(rfull >= tipBeginsAt);
                    
                    % ...then the other.
                    tmpX1full = petalX + deltaArray(4*pn-1, :);
                    tmpY1full = -petalY + deltaArray(4*pn-0, :);
                    rfull = sqrt(tmpX1full.*tmpX1full + tmpY1full.*tmpY1full);
                    tmpX1 = tmpX1full(rfull >= tipBeginsAt);
                    tmpY1 = tmpY1full(rfull >= tipBeginsAt);
                    
                    newX = tmpX1 + errorArray{j}.radial;
                    newY = tmpY1 + errorArray{j}.lateral;
                    tmpX1 = cos(errorArray{j}.clocking)*(newX - tipBeginsAt - errorArray{j}.radial) + ...
                        sin(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + tipBeginsAt + errorArray{j}.radial;
                    tmpY1 = -sin(errorArray{j}.clocking)*(newX - tipBeginsAt - errorArray{j}.radial) + ...
                        cos(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + errorArray{j}.lateral;
                    newX = tmpX1;
                    newY = tmpY1;
                    
                    deltaArray(4*pn-1, rfull >= tipBeginsAt) = newX - petalX(rfull >= tipBeginsAt);
                    deltaArray(4*pn-0, rfull >= tipBeginsAt) = newY + petalY(rfull >= tipBeginsAt);
                end
                
        end
    end
    
    % in-plane bending
    for j = 1:length(errorArray)
        switch errorArray{j}.name
            case 'ipbend'
                tmp = (petalR - a)/(maxR - a);
                tmppoly = errorArray{j}.constant + tmp.*(errorArray{j}.linear + errorArray{j}.quadratic.*tmp);
                dTheta = tmppoly./petalR;
                
                pn = errorArray{j}.number;
                if (pn == 0)
                    for k = 1:numPetals
                        tmpX1 = petalX + deltaArray(4*k-3, :);
                        tmpY1 = petalY + deltaArray(4*k-2, :);
                        tmpR1 = sqrt(tmpX1.*tmpX1 + tmpY1.*tmpY1);
                        tmpT1 = atan2(tmpY1, tmpX1);
                        deltaArray(4*k-3, :) = tmpR1.*cos(tmpT1 + dTheta) - petalX;
                        deltaArray(4*k-2, :) = tmpR1.*sin(tmpT1 + dTheta) - petalY;
                        
                        tmpX2 = petalX + deltaArray(4*k-1, :);
                        tmpY2 = -petalY + deltaArray(4*k-0, :);
                        tmpR2 = sqrt(tmpX2.*tmpX2 + tmpY2.*tmpY2);
                        tmpT2 = atan2(tmpY2, tmpX2);
                        deltaArray(4*k-1, :) = tmpR2.*cos(tmpT2 + dTheta) - petalX;
                        deltaArray(4*k-0, :) = tmpR2.*sin(tmpT2 + dTheta) + petalY;
                    end
                else
                    tmpX1 = petalX + deltaArray(4*pn-3, :);
                    tmpY1 = petalY + deltaArray(4*pn-2, :);
                    tmpR1 = sqrt(tmpX1.*tmpX1 + tmpY1.*tmpY1);
                    tmpT1 = atan2(tmpY1, tmpX1);
                    deltaArray(4*pn-3, :) = tmpR1.*cos(tmpT1 + dTheta) - petalX;
                    deltaArray(4*pn-2, :) = tmpR1.*sin(tmpT1 + dTheta) - petalY;
                    
                    tmpX2 = petalX + deltaArray(4*pn-1, :);
                    tmpY2 = -petalY + deltaArray(4*pn-0, :);
                    tmpR2 = sqrt(tmpX2.*tmpX2 + tmpY2.*tmpY2);
                    tmpT2 = atan2(tmpY2, tmpX2);
                    deltaArray(4*pn-1, :) = tmpR2.*cos(tmpT2 + dTheta) - petalX;
                    deltaArray(4*pn-0, :) = tmpR2.*sin(tmpT2 + dTheta) + petalY;
                end
        end
    end
    
    % out-of-plane bending
    for j = 1:length(errorArray)
        switch errorArray{j}.name
            case 'oopbend'
                pn = errorArray{j}.number;
                if (pn == 0)
                    for k = 1:numPetals
                        % One side of petal...
                        tmpX1 = petalX + deltaArray(4*k-3, :);
                        tmpY1 = petalY + deltaArray(4*k-2, :);
                        
                        xInterp = [tmpX1(1) (tmpX1(1:end-1) + tmpX1(2:end))/2 tmpX1(end)];
                        yInterp = [tmpY1(1) (tmpY1(1:end-1) + tmpY1(2:end))/2 tmpY1(end)];
                        
                        tmpSc = (tmpX1 - tmpX1(1))/(tmpX1(end)-tmpX1(1));
                        zInterp = errorArray{j}.quadratic.*tmpSc.^2 + errorArray{j}.linear.*tmpSc;
                        
                        S0 = [0 cumsum(sqrt(diff(xInterp).^2 + diff(yInterp).^2))];
                        S1 = [0 cumsum(sqrt(diff(xInterp).^2 + diff(yInterp).^2 + diff(zInterp).^2))];
                        
                        naninds = isnan(xInterp);
                        tmpX1a = interp1(S0(~naninds), xInterp(~naninds), S1, 'linear', NaN);
                        tmpY1a = interp1(S0(~naninds), yInterp(~naninds), S1, 'linear', NaN);
                        
                        deltaArray(4*k-3, :) = tmpX1a - petalX;
                        deltaArray(4*k-2, :) = tmpY1a - petalY;
                        
                        % ...then the other.
                        tmpX1 = petalX + deltaArray(4*k-1, :);
                        tmpY1 = -petalY + deltaArray(4*k-0, :);
                        
                        xInterp = [tmpX1(1) (tmpX1(1:end-1) + tmpX1(2:end))/2 tmpX1(end)];
                        yInterp = [tmpY1(1) (tmpY1(1:end-1) + tmpY1(2:end))/2 tmpY1(end)];
                        
                        tmpSc = (xInterp - xInterp(1))/(xInterp(end)-xInterp(1));
                        zInterp = errorArray{j}.quadratic.*tmpSc.^2 + errorArray{j}.linear.*tmpSc;
                        
                        S0 = [0 cumsum(sqrt(diff(xInterp).^2 + diff(yInterp).^2))];
                        S1 = [0 cumsum(sqrt(diff(xInterp).^2 + diff(yInterp).^2 + diff(zInterp).^2))];
                        
                        naninds = isnan(xInterp);
                        tmpX1a = interp1(S0(~naninds), xInterp(~naninds), S1, 'linear', NaN);
                        tmpY1a = interp1(S0(~naninds), yInterp(~naninds), S1, 'linear', NaN);
                        
                        deltaArray(4*k-1, :) = tmpX1a - petalX;
                        deltaArray(4*k-0, :) = tmpY1a + petalY;
                    end
                else
                    % One side of petal...
                    tmpX1 = petalX + deltaArray(4*pn-3, :);
                    tmpY1 = petalY + deltaArray(4*pn-2, :);
                    
                    xInterp = [tmpX1(1) (tmpX1(1:end-1) + tmpX1(2:end))/2 tmpX1(end)];
                    yInterp = [tmpY1(1) (tmpY1(1:end-1) + tmpY1(2:end))/2 tmpY1(end)];
                    
                    tmpSc = (xInterp - xInterp(1))/(maxR-xInterp(1));
                    zInterp = errorArray{j}.quadratic.*tmpSc.^2 + errorArray{j}.linear.*tmpSc;
                    
                    S0 = [0 cumsum(sqrt(diff(xInterp).^2 + diff(yInterp).^2))];
                    S1 = [0 cumsum(sqrt(diff(xInterp).^2 + diff(yInterp).^2 + diff(zInterp).^2))];
                    
                    naninds = isnan(xInterp);
                    tmpX1a = interp1(S1(~naninds), xInterp(~naninds), S0, 'linear', NaN);
                    naninds2 = isnan(tmpX1a);
                    tmpY1b = interp1(tmpX1a(~naninds2), yInterp(~naninds2), tmpX1, 'linear', NaN); % ADD TO OTHER TWO
                    
                    deltaArray(4*pn-3, :) = tmpX1 - petalX;
                    deltaArray(4*pn-2, :) = tmpY1b - petalY;
                    
                    % ...then the other.
                    tmpX1 = petalX + deltaArray(4*pn-1, :);
                    tmpY1 = -petalY + deltaArray(4*pn-0, :);
                    
                    xInterp = [tmpX1(1) (tmpX1(1:end-1) + tmpX1(2:end))/2 tmpX1(end)];
                    yInterp = [tmpY1(1) (tmpY1(1:end-1) + tmpY1(2:end))/2 tmpY1(end)];
                    
                    tmpSc = (xInterp - xInterp(1))/(maxR-xInterp(1));
                    zInterp = errorArray{j}.quadratic.*tmpSc.^2 + errorArray{j}.linear.*tmpSc;
                    
                    S0 = [0 cumsum(sqrt(diff(xInterp).^2 + diff(yInterp).^2))];
                    S1 = [0 cumsum(sqrt(diff(xInterp).^2 + diff(yInterp).^2 + diff(zInterp).^2))];
                    
                    naninds = isnan(xInterp);
                    tmpX1a = interp1(S1(~naninds), xInterp(~naninds), S0, 'linear', NaN);
                    naninds2 = isnan(tmpX1a);
                    tmpY1b = interp1(tmpX1a(~naninds2), yInterp(~naninds2), tmpX1, 'linear', NaN); % ADD TO OTHER TWO
                    
                    deltaArray(4*pn-1, :) = tmpX1 - petalX;
                    deltaArray(4*pn-0, :) = tmpY1b + petalY;
                end
        end
    end
    
    % Deployment
    for j = 1:length(errorArray)
        switch errorArray{j}.name
            case 'deployment'
                pn = errorArray{j}.number;
                if (pn == 0)
                    for k = 1:numPetals
                        tmpX1 = petalX + deltaArray(4*k-3, :);
                        tmpY1 = petalY + deltaArray(4*k-2, :);
                        
                        newX = tmpX1 + errorArray{j}.radial;
                        newY = tmpY1 + errorArray{j}.lateral;
                        tmpX1 = cos(errorArray{j}.clocking)*(newX - a -  errorArray{j}.radial) ...
                            + sin(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + a + errorArray{j}.radial;
                        tmpY1 = -sin(errorArray{j}.clocking)*(newX - a -  errorArray{j}.radial) ...
                            + cos(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + errorArray{j}.lateral;
                        newX = cos(errorArray{j}.angular)*tmpX1 + sin(errorArray{j}.angular)*tmpY1;
                        newY = -sin(errorArray{j}.angular)*tmpX1 + cos(errorArray{j}.angular)*tmpY1;
                        
                        deltaArray(4*k-3, :) = newX - petalX;
                        deltaArray(4*k-2, :) = newY - petalY;
                        
                        tmpX1 = petalX + deltaArray(4*k-1, :);
                        tmpY1 = -petalY + deltaArray(4*k-0, :);
                        
                        newX = tmpX1 + errorArray{j}.radial;
                        newY = tmpY1 + errorArray{j}.lateral;
                        tmpX1 = cos(errorArray{j}.clocking)*(newX - a -  errorArray{j}.radial) ...
                            + sin(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + a + errorArray{j}.radial;
                        tmpY1 = -sin(errorArray{j}.clocking)*(newX - a -  errorArray{j}.radial) ...
                            + cos(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + errorArray{j}.lateral;
                        newX = cos(errorArray{j}.angular)*tmpX1 + sin(errorArray{j}.angular)*tmpY1;
                        newY = -sin(errorArray{j}.angular)*tmpX1 + cos(errorArray{j}.angular)*tmpY1;
                        
                        deltaArray(4*k-1, :) = newX - petalX;
                        deltaArray(4*k-0, :) = newY + petalY;
                    end
                else
                    tmpX1 = petalX + deltaArray(4*pn-3, :);
                    tmpY1 = petalY + deltaArray(4*pn-2, :);
                    
                    newX = tmpX1 + errorArray{j}.radial;
                    newY = tmpY1 + errorArray{j}.lateral;
                    tmpX1 = cos(errorArray{j}.clocking)*(newX - a -  errorArray{j}.radial) ...
                        + sin(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + a + errorArray{j}.radial;
                    tmpY1 = -sin(errorArray{j}.clocking)*(newX - a -  errorArray{j}.radial) ...
                        + cos(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + errorArray{j}.lateral;
                    newX = cos(errorArray{j}.angular)*tmpX1 + sin(errorArray{j}.angular)*tmpY1;
                    newY = -sin(errorArray{j}.angular)*tmpX1 + cos(errorArray{j}.angular)*tmpY1;
                    
                    deltaArray(4*pn-3, :) = newX - petalX;
                    deltaArray(4*pn-2, :) = newY - petalY;
                    
                    tmpX1 = petalX + deltaArray(4*pn-1, :);
                    tmpY1 = -petalY + deltaArray(4*pn-0, :);
                    
                    newX = tmpX1 + errorArray{j}.radial;
                    newY = tmpY1 + errorArray{j}.lateral;
                    tmpX1 = cos(errorArray{j}.clocking)*(newX - a -  errorArray{j}.radial) ...
                        + sin(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + a + errorArray{j}.radial;
                    tmpY1 = -sin(errorArray{j}.clocking)*(newX - a -  errorArray{j}.radial) ...
                        + cos(errorArray{j}.clocking)*(newY - errorArray{j}.lateral) + errorArray{j}.lateral;
                    newX = cos(errorArray{j}.angular)*tmpX1 + sin(errorArray{j}.angular)*tmpY1;
                    newY = -sin(errorArray{j}.angular)*tmpX1 + cos(errorArray{j}.angular)*tmpY1;
                    
                    deltaArray(4*pn-1, :) = newX - petalX;
                    deltaArray(4*pn-0, :) = newY + petalY;
                end
        end
    end
    
    % Thermal distortions
    for j = 1:length(errorArray)
        switch errorArray{j}.name
            case 'uniformthermal'
                pn = errorArray{j}.number;
                if pn == 0 % global error
                    for k = 1:numPetals
                        tmpX1 = petalX + deltaArray(4*k-3, :);
                        tmpY1 = petalY + deltaArray(4*k-2, :);
                        
                        tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
                        tmpT1 = atan2(tmpY1, tmpX1);
                        
                        R2 = (tmpR1-a)*(1+errorArray{j}.CTE*errorArray{j}.dt)+a;
                        T2 = tmpT1*(1+errorArray{j}.CTE*errorArray{j}.dt);
                        
                        deltaArray(4*k-3, :) = R2.*cos(T2) - petalX;
                        deltaArray(4*k-2, :) = R2.*sin(T2) - petalY;
                        
                        tmpX1 = petalX + deltaArray(4*k-1, :);
                        tmpY1 = -petalY + deltaArray(4*k-0, :);
                        
                        tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
                        tmpT1 = atan2(tmpY1, tmpX1);
                        
                        R2 = (tmpR1-a)*(1+errorArray{j}.CTE*errorArray{j}.dt)+a;
                        T2 = tmpT1*(1+errorArray{j}.CTE*errorArray{j}.dt);
                        
                        deltaArray(4*k-1, :) = R2.*cos(T2) - petalX;
                        deltaArray(4*k-0, :) = R2.*sin(T2) + petalY;
                    end
                else
                    tmpX1 = petalX + deltaArray(4*pn-3, :);
                    tmpY1 = petalY + deltaArray(4*pn-2, :);
                    
                    tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
                    tmpT1 = atan2(tmpY1, tmpX1);
                    
                    R2 = (tmpR1-a)*(1+errorArray{j}.CTE*errorArray{j}.dt)+a;
                    T2 = tmpT1*(1+errorArray{j}.CTE*errorArray{j}.dt);
                    
                    deltaArray(4*pn-3, :) = R2.*cos(T2) - petalX;
                    deltaArray(4*pn-2, :) = R2.*sin(T2) - petalY;
                    
                    tmpX1 = petalX + deltaArray(4*pn-1, :);
                    tmpY1 = -petalY + deltaArray(4*pn-0, :);
                    
                    tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
                    tmpT1 = atan2(tmpY1, tmpX1);
                    
                    R2 = (tmpR1-a)*(1+errorArray{j}.CTE*errorArray{j}.dt)+a;
                    T2 = tmpT1*(1+errorArray{j}.CTE*errorArray{j}.dt);
                    
                    deltaArray(4*pn-1, :) = R2.*cos(T2) - petalX;
                    deltaArray(4*pn-0, :) = R2.*sin(T2) + petalY;
                end
            case 'radgradthermal'
                pn = errorArray{j}.number;
                if pn == 0 % global error
                    for k = 1:numPetals
                        tmpX1 = petalX + deltaArray(4*k-3, :);
                        tmpY1 = petalY + deltaArray(4*k-2, :);
                        
                        tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
                        tmpT1 = atan2(tmpY1, tmpX1);
                        
                        gtemp = (tmpR1-a)/(maxR-a)*errorArray{j}.CTE*errorArray{j}.dt;
                        R2 = (tmpR1-a).*(1+gtemp)+a;
                        T2 = tmpT1.*(1+gtemp);
                        
                        deltaArray(4*k-3, :) = R2.*cos(T2) - petalX;
                        deltaArray(4*k-2, :) = R2.*sin(T2) - petalY;
                        
                        tmpX1 = petalX + deltaArray(4*k-1, :);
                        tmpY1 = -petalY + deltaArray(4*k-0, :);
                        
                        tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
                        tmpT1 = atan2(tmpY1, tmpX1);
                        
                        gtemp = (tmpR1-a)/(maxR-a)*errorArray{j}.CTE*errorArray{j}.dt;
                        R2 = (tmpR1-a).*(1+gtemp)+a;
                        T2 = tmpT1.*(1+gtemp);
                        
                        deltaArray(4*k-1, :) = R2.*cos(T2) - petalX;
                        deltaArray(4*k-0, :) = R2.*sin(T2) + petalY;
                    end
                else
                    tmpX1 = petalX + deltaArray(4*pn-3, :);
                    tmpY1 = petalY + deltaArray(4*pn-2, :);
                    
                    tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
                    tmpT1 = atan2(tmpY1, tmpX1);
                    
                    gtemp = (tmpR1-a)/(maxR-a)*errorArray{j}.CTE*errorArray{j}.dt;
                    R2 = (tmpR1-a).*(1+gtemp)+a;
                    T2 = tmpT1.*(1+gtemp);
                    
                    deltaArray(4*pn-3, :) = R2.*cos(T2) - petalX;
                    deltaArray(4*pn-2, :) = R2.*sin(T2) - petalY;
                    
                    tmpX1 = petalX + deltaArray(4*pn-1, :);
                    tmpY1 = -petalY + deltaArray(4*pn-0, :);
                    
                    tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
                    tmpT1 = atan2(tmpY1, tmpX1);
                    
                    gtemp = (tmpR1-a)/(maxR-a)*errorArray{j}.CTE*errorArray{j}.dt;
                    R2 = (tmpR1-a).*(1+gtemp)+a;
                    T2 = tmpT1.*(1+gtemp);
                    
                    deltaArray(4*pn-1, :) = R2.*cos(T2) - petalX;
                    deltaArray(4*pn-0, :) = R2.*sin(T2) + petalY;
                end
            case 'sinethermal'
                pn = errorArray{j}.number;
                if pn == 0 % global error
                    for k = 1:numPetals
                        tmpX1 = petalX + deltaArray(4*k-3, :);
                        tmpY1 = petalY + deltaArray(4*k-2, :);
                        
                        tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
                        tmpT1 = atan2(tmpY1, tmpX1);
                        
                        gtemp = sin(2*pi*errorArray{j}.cycles*(tmpR1-a)/(maxR-a))*errorArray{j}.CTE*errorArray{j}.dt;
                        R2 = tmpR1;
                        T2 = tmpT1.*(1+gtemp);
                        
                        deltaArray(4*k-3, :) = R2.*cos(T2) - petalX;
                        deltaArray(4*k-2, :) = R2.*sin(T2) - petalY;
                        
                        tmpX1 = petalX + deltaArray(4*k-1, :);
                        tmpY1 = -petalY + deltaArray(4*k-0, :);
                        
                        tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
                        tmpT1 = atan2(tmpY1, tmpX1);
                        
                        gtemp = sin(2*pi*errorArray{j}.cycles*(tmpR1-a)/(maxR-a))*errorArray{j}.CTE*errorArray{j}.dt;
                        R2 = tmpR1;
                        T2 = tmpT1.*(1+gtemp);
                        
                        deltaArray(4*k-1, :) = R2.*cos(T2) - petalX;
                        deltaArray(4*k-0, :) = R2.*sin(T2) + petalY;
                    end
                else
                    tmpX1 = petalX + deltaArray(4*pn-3, :);
                    tmpY1 = petalY + deltaArray(4*pn-2, :);
                    
                    tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
                    tmpT1 = atan2(tmpY1, tmpX1);
                    
                    gtemp = sin(2*pi*errorArray{j}.cycles*(tmpR1-a)/(maxR-a))*errorArray{j}.CTE*errorArray{j}.dt;
                    R2 = tmpR1;
                    T2 = tmpT1.*(1+gtemp);
                    
                    deltaArray(4*pn-3, :) = R2.*cos(T2) - petalX;
                    deltaArray(4*pn-2, :) = R2.*sin(T2) - petalY;
                    
                    tmpX1 = petalX + deltaArray(4*pn-1, :);
                    tmpY1 = -petalY + deltaArray(4*pn-0, :);
                    
                    tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
                    tmpT1 = atan2(tmpY1, tmpX1);
                    
                    gtemp = sin(2*pi*errorArray{j}.cycles*(tmpR1-a)/(maxR-a))*errorArray{j}.CTE*errorArray{j}.dt;
                    R2 = tmpR1;
                    T2 = tmpT1.*(1+gtemp);
                    
                    deltaArray(4*pn-1, :) = R2.*cos(T2) - petalX;
                    deltaArray(4*pn-0, :) = R2.*sin(T2) + petalY;
                end
            case 'lengththermal'
                pn = errorArray{j}.number;
                if pn == 0 % global error
                    for k = 1:numPetals
                        deltaArray(4*k-3, :) = deltaArray(4*k-3, :) ...
                            + errorArray{j}.CTE*errorArray{j}.dt*(petalX + deltaArray(4*k-3, :)-a);
                        deltaArray(4*k-1, :) = deltaArray(4*k-1, :) ...
                            + errorArray{j}.CTE*errorArray{j}.dt*(petalX + deltaArray(4*k-1, :)-a);
                    end
                else
                    deltaArray(4*pn-3, :) = deltaArray(4*pn-3, :) ...
                        + errorArray{j}.CTE*errorArray{j}.dt*(petalX + deltaArray(4*pn-3, :)-a);
                    deltaArray(4*pn-1, :) = deltaArray(4*pn-1, :) ...
                        + errorArray{j}.CTE*errorArray{j}.dt*(petalX + deltaArray(4*pn-1, :)-a);
                end
%                 pn = errorArray{j}.number;
%                 if pn == 0 % global error
%                     for k = 1:numPetals
%                         tmpX1 = petalX + deltaArray(4*k-3, :);
%                         tmpY1 = petalY + deltaArray(4*k-2, :);
%                         
%                         tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
%                         tmpT1 = atan2(tmpY1, tmpX1);
%                         
%                         R2 = (tmpR1-a)*(1+errorArray{j}.CTE*errorArray{j}.dt)+a;
%                         T2 = tmpT1;
%                         
%                         deltaArray(4*k-3, :) = R2.*cos(T2) - petalX;
%                         deltaArray(4*k-2, :) = R2.*sin(T2) - petalY;
%                         
%                         tmpX1 = petalX + deltaArray(4*k-1, :);
%                         tmpY1 = -petalY + deltaArray(4*k-0, :);
%                         
%                         tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
%                         tmpT1 = atan2(tmpY1, tmpX1);
%                         
%                         R2 = (tmpR1-a)*(1+errorArray{j}.CTE*errorArray{j}.dt)+a;
%                         T2 = tmpT1;
%                         
%                         deltaArray(4*k-1, :) = R2.*cos(T2) - petalX;
%                         deltaArray(4*k-0, :) = R2.*sin(T2) + petalY;
%                     end
%                 else
%                     tmpX1 = petalX + deltaArray(4*pn-3, :);
%                     tmpY1 = petalY + deltaArray(4*pn-2, :);
%                     
%                     tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
%                     tmpT1 = atan2(tmpY1, tmpX1);
%                     
%                     R2 = (tmpR1-a)*(1+errorArray{j}.CTE*errorArray{j}.dt)+a;
%                     T2 = tmpT1;
%                     
%                     deltaArray(4*pn-3, :) = R2.*cos(T2) - petalX;
%                     deltaArray(4*pn-2, :) = R2.*sin(T2) - petalY;
%                     
%                     tmpX1 = petalX + deltaArray(4*pn-1, :);
%                     tmpY1 = -petalY + deltaArray(4*pn-0, :);
%                     
%                     tmpR1 = sqrt(tmpX1.^2 + tmpY1.^2);
%                     tmpT1 = atan2(tmpY1, tmpX1);
%                     
%                     R2 = (tmpR1-a)*(1+errorArray{j}.CTE*errorArray{j}.dt)+a;
%                     T2 = tmpT1;
%                     
%                     deltaArray(4*pn-1, :) = R2.*cos(T2) - petalX;
%                     deltaArray(4*pn-0, :) = R2.*sin(T2) + petalY;
%                 end
        end
    end
    
    % S-shaped structural term
    for j = 1:length(errorArray)
        switch errorArray{j}.name
            case 'sshape'
                tmp = (petalR - a)/(maxR - a);
                dTheta = 2*pi/numPetals*sin(2*pi*tmp)*errorArray{j}.smag./petalR;
                
                pn = errorArray{j}.number;
                if (pn == 0)
                    for k = 1:numPetals
                        tmpX1 = petalX + deltaArray(4*k-3, :);
                        tmpY1 = petalY + deltaArray(4*k-2, :);
                        tmpR1 = sqrt(tmpX1.*tmpX1 + tmpY1.*tmpY1);
                        tmpT1 = atan2(tmpY1, tmpX1);
                        deltaArray(4*k-3, :) = tmpR1.*cos(tmpT1 + dTheta) - petalX;
                        deltaArray(4*k-2, :) = tmpR1.*sin(tmpT1 + dTheta) - petalY;
                        
                        tmpX2 = petalX + deltaArray(4*k-1, :);
                        tmpY2 = -petalY + deltaArray(4*k-0, :);
                        tmpR2 = sqrt(tmpX2.*tmpX2 + tmpY2.*tmpY2);
                        tmpT2 = atan2(tmpY2, tmpX2);
                        deltaArray(4*k-1, :) = tmpR2.*cos(tmpT2 + dTheta) - petalX;
                        deltaArray(4*k-0, :) = tmpR2.*sin(tmpT2 + dTheta) + petalY;
                    end
                else
                    tmpX1 = petalX + deltaArray(4*pn-3, :);
                    tmpY1 = petalY + deltaArray(4*pn-2, :);
                    tmpR1 = sqrt(tmpX1.*tmpX1 + tmpY1.*tmpY1);
                    tmpT1 = atan2(tmpY1, tmpX1);
                    deltaArray(4*pn-3, :) = tmpR1.*cos(tmpT1 + dTheta) - petalX;
                    deltaArray(4*pn-2, :) = tmpR1.*sin(tmpT1 + dTheta) - petalY;
                    
                    tmpX2 = petalX + deltaArray(4*pn-1, :);
                    tmpY2 = -petalY + deltaArray(4*pn-0, :);
                    tmpR2 = sqrt(tmpX2.*tmpX2 + tmpY2.*tmpY2);
                    tmpT2 = atan2(tmpY2, tmpX2);
                    deltaArray(4*pn-1, :) = tmpR2.*cos(tmpT2 + dTheta) - petalX;
                    deltaArray(4*pn-0, :) = tmpR2.*sin(tmpT2 + dTheta) + petalY;
                end
        end
    end
    
    % Truss ellipticity
    for j = 1:length(errorArray)
        switch errorArray{j}.name
            case 'ellipticity'
                clc
                e = errorArray{j}.eccentricity;
                
                roots = zeros(4, numPetals);
                theta = zeros(1, 2*numPetals);
                d0 = zeros(1, 2*numPetals);
                
                for k = 1:numPetals
                    tmpX1 = petalX + deltaArray(4*k-3, :);
                    tmpY1 = petalY + deltaArray(4*k-2, :);
                    tmpX2 = petalX + deltaArray(4*k-1, :);
                    tmpY2 = -petalY + deltaArray(4*k-0, :);
                    
                    cn = cos(2*pi/numPetals*(k-1));
                    sn = sin(2*pi/numPetals*(k-1));
                    
                    rMat = [cn -sn 0 0; sn cn 0 0; 0 0 cn -sn; 0 0 sn cn];
                    roots(:,k) = rMat*[tmpX1(1) tmpY1(1) tmpX2(1) tmpY2(1)].';
                end
                
                for k = 1:numPetals
                    d0(2*k-1) = sqrt((roots(1,k) - roots(3,k)).^2 + (roots(2,k) - roots(4,k)).^2);
                    
                    if k == numPetals
                        m = 1;
                    else
                        m = k+1;
                    end
                    d0(2*k) = sqrt((roots(3,m) - roots(1,k)).^2 + (roots(4,m) - roots(2,k)).^2);
                    
                    theta(2*k-1) = acos((roots(1,k)*roots(3,k) + roots(2,k)*roots(4,k))/sqrt((roots(1,k)^2 + roots(2,k)^2)*(roots(3,k)^2 + roots(4,k)^2)));
                    theta(2*k) = acos((roots(3,m)*roots(1,k) + roots(4,m)*roots(2,k))/sqrt((roots(3,m)^2 + roots(4,m)^2)*(roots(1,k)^2 + roots(2,k)^2)));
                end
                
                global aGlob eGlob nGlob
                eGlob = e;
                aGlob = 1;
                nGlob = d0;
                
                options = optimset('TolFun', 1e-9, 'TolX', 1e-9, 'Algorithm', 'levenberg-marquardt', 'MaxFunEvals', 5000, 'MaxIter', 5000);
                tmplsq = lsqnonlin(@distlsq, theta, [], [], options);
                
                ecclsq = cumsum([0 tmplsq(1:end-1)]);
                ecclsq = ecclsq - ecclsq(2)/2; % Keep petal #1 centered.
                
                d1 = dists(ecclsq, aGlob, eGlob);
                newSMA = mean(d0./d1);
                
                for k = 1:numPetals
                    tmpX1 = petalX + deltaArray(4*k-3, :);
                    tmpY1 = petalY + deltaArray(4*k-2, :);
                    tmpX2 = petalX + deltaArray(4*k-1, :);
                    tmpY2 = -petalY + deltaArray(4*k-0, :);
                    
                    cn = cos(2*pi/numPetals*(k-1));
                    sn = sin(2*pi/numPetals*(k-1));
                    
                    rMat = [cn -sn 0 0; sn cn 0 0; 0 0 cn -sn; 0 0 sn cn];
                    newRoots = rMat'*[newSMA*cos(ecclsq(2*k)) newSMA*sqrt(1-e^2)*sin(ecclsq(2*k)) ...
                        newSMA*cos(ecclsq(2*k-1)) newSMA*sqrt(1-e^2)*sin(ecclsq(2*k-1))].';
                    
                    tmpmidx = (tmpX1(1) + tmpX2(1))/2;
                    tmpmidy = (tmpY1(1) + tmpY2(1))/2;
                    midx = (newRoots(1) + newRoots(3))/2;
                    midy = (newRoots(2) + newRoots(4))/2;
                    
                    thetaOld = atan2(tmpY2(1)-tmpY1(1), tmpX2(1) - tmpX1(1));
                    thetaNew = atan2(newRoots(4)-newRoots(2), newRoots(3) - newRoots(1));
                    
                    newX1 = (tmpX1 - tmpmidx)*cos(thetaNew - thetaOld) - (tmpY1 - tmpmidy)*sin(thetaNew - thetaOld) + midx;
                    newY1 = (tmpX1 - tmpmidx)*sin(thetaNew - thetaOld) + (tmpY1 - tmpmidy)*cos(thetaNew - thetaOld) + midy;
                    newX2 = (tmpX2 - tmpmidx)*cos(thetaNew - thetaOld) - (tmpY2 - tmpmidy)*sin(thetaNew - thetaOld) + midx;
                    newY2 = (tmpX2 - tmpmidx)*sin(thetaNew - thetaOld) + (tmpY2 - tmpmidy)*cos(thetaNew - thetaOld) + midy;
                    
                    deltaArray(4*k-3, :) = newX1 - petalX;
                    deltaArray(4*k-2, :) = newY1 - petalY;
                    deltaArray(4*k-1, :) = newX2 - petalX;
                    deltaArray(4*k-0, :) = newY2 + petalY;
                end
            case 'circharmonic'
                for k = 1:numPetals
                    tmpX1 = petalX + deltaArray(4*k-3, :);
                    tmpY1 = petalY + deltaArray(4*k-2, :);
                    
                    tmpR1 = sqrt(tmpX1.*tmpX1 + tmpY1.*tmpY1);
                    tmpT1 = atan2(tmpY1, tmpX1);
                    
                    R1 = tmpR1 + errorArray{j}.circeps*sin((tmpT1+ 2*pi/numPetals*(k-1))*errorArray{j}.order-errorArray{j}.dang);
                    T1 = tmpT1;
                    
                    deltaArray(4*k-3, :) = R1.*cos(T1) - petalX;
                    deltaArray(4*k-2, :) = R1.*sin(T1) - petalY;
                    
                    tmpX2 = petalX + deltaArray(4*k-1, :);
                    tmpY2 = -petalY + deltaArray(4*k-0, :);
                    
                    tmpR2 = sqrt(tmpX2.^2 + tmpY2.^2);
                    tmpT2 = atan2(tmpY2, tmpX2);
                    
                    R2 = tmpR2 + errorArray{j}.circeps*sin((tmpT1+ 2*pi/numPetals*(k-1))*errorArray{j}.order-errorArray{j}.dang);
                    T2 = tmpT2;
                    
                    deltaArray(4*k-1, :) = R2.*cos(T2) - petalX;
                    deltaArray(4*k-0, :) = R2.*sin(T2) + petalY;
                end
            case 'structmode'
                roots = zeros(4, numPetals);
                theta = zeros(1, 2*numPetals);
                d0 = zeros(1, 2*numPetals);
                
                for k = 1:numPetals
                    tmpX1 = petalX + deltaArray(4*k-3, :);
                    tmpY1 = petalY + deltaArray(4*k-2, :);
                    tmpX2 = petalX + deltaArray(4*k-1, :);
                    tmpY2 = -petalY + deltaArray(4*k-0, :);
                    
                    cn = cos(2*pi/numPetals*(k-1));
                    sn = sin(2*pi/numPetals*(k-1));
                    
                    rMat = [cn -sn 0 0; sn cn 0 0; 0 0 cn -sn; 0 0 sn cn];
                    roots(:,k) = rMat*diag([1 -1 1 -1])*[tmpX1(1) tmpY1(1) tmpX2(1) tmpY2(1)].';
                end
                
                ang1 = -pi/numPetals-0.0003314;
                R1 = [cos(ang1) -sin(ang1) 0; sin(ang1) cos(ang1) 0; 0 0 1];
                ang2 = -pi/numPetals+0.0003314;
                R2 = [cos(ang2) -sin(ang2) 0; sin(ang2) cos(ang2) 0; 0 0 1];
                z0 = 2715692/3000000;
                
                R12 = [R1 zeros(3,3); zeros(3,3) R2];
                
                %                 roots12 = R12*[roots(1,:); roots(2,:); z0*ones(1, numPetals); ...
                %                                roots(3,:); roots(4,:); z0*ones(1, numPetals)]*9.75/10;
                %                 loc = [errorArray{j}.location(:,2:4).'; errorArray{j}.location(2:end,2:4).' errorArray{j}.location(1,2:4).';];%[errorArray{j}.location(:,2:4).'; errorArray{j}.location(end,2:4).' errorArray{j}.location(1:end-1,2:4).'];
                
                location = errorArray{j}.location(:,2:4).';
                displacement = errorArray{j}.displacement(:, 2:4).';
                ld = [(location + displacement); ...
                    (location(:, 2:end) + displacement(:, 2:end)) (location(:,1) + displacement(:,1))];
                
                newr = R12'*ld*10/9.75;
                
                
                %                 n = 2;
                %                 plot(loc(1,1:end), loc(2,1:end), ...
                %                     loc(1,n), loc(2,n), 'ko', loc(4,n), loc(5,n), 'go',...
                %                     roots12(1, n), roots12(2, n), 'kx', roots12(4,n),  roots12(5,n), 'gx')
                
                
                
                for k = 1:numPetals
                    tmpX1 = petalX + deltaArray(4*k-3, :);
                    tmpY1 = petalY + deltaArray(4*k-2, :);
                    tmpX2 = petalX + deltaArray(4*k-1, :);
                    tmpY2 = -petalY + deltaArray(4*k-0, :);
                    
                    cn = cos(2*pi/numPetals*(k-1));
                    sn = sin(2*pi/numPetals*(k-1));
                    
                    rMat = [cn -sn 0 0; sn cn 0 0; 0 0 cn -sn; 0 0 sn cn];
                    newRoots = diag([1 -1 1 -1])*rMat'*[newr(1,k) newr(2,k) newr(4,k) newr(5,k)].';
                    
                    tmpmidx = (tmpX1(1) + tmpX2(1))/2;
                    tmpmidy = (tmpY1(1) + tmpY2(1))/2;
                    midx = (newRoots(1) + newRoots(3))/2;
                    midy = (newRoots(2) + newRoots(4))/2;
                    
                    thetaOld = atan2(tmpY2(1)-tmpY1(1), tmpX2(1) - tmpX1(1));
                    thetaNew = atan2(newRoots(4)-newRoots(2), newRoots(3) - newRoots(1));
                    
                    newX1 = (tmpX1 - tmpmidx)*cos(thetaNew - thetaOld) - (tmpY1 - tmpmidy)*sin(thetaNew - thetaOld) + midx;
                    newY1 = (tmpX1 - tmpmidx)*sin(thetaNew - thetaOld) + (tmpY1 - tmpmidy)*cos(thetaNew - thetaOld) + midy;
                    newX2 = (tmpX2 - tmpmidx)*cos(thetaNew - thetaOld) - (tmpY2 - tmpmidy)*sin(thetaNew - thetaOld) + midx;
                    newY2 = (tmpX2 - tmpmidx)*sin(thetaNew - thetaOld) + (tmpY2 - tmpmidy)*cos(thetaNew - thetaOld) + midy;
                    
                    deltaArray(4*k-3, :) = newX1 - petalX;
                    deltaArray(4*k-2, :) = newY1 - petalY;
                    deltaArray(4*k-1, :) = newX2 - petalX;
                    deltaArray(4*k-0, :) = newY2 + petalY;
                    zArray(2*k-1,:) = newr(3,k)-z0;
                    zArray(2*k,:) = newr(6,k)-z0;
                end
        end
    end
    
    
    
    vecPetalArray = cell(1, numPetals);
    for k = 1:numPetals
        tmp = [(petalX + deltaArray(4*k-1,:)) fliplr(petalX + deltaArray(4*k-3,:)); ...
            (-petalY + deltaArray(4*k,:)) fliplr(petalY + deltaArray(4*k-2,:)); ...
            (zArray(2*k,:)) fliplr(zArray(2*k-1,:))];
        inds1 = isnan(tmp(1,:));
        inds2 = isnan(tmp(2,:));
        inds3 = isnan(tmp(3,:));
        tmp = tmp(:, ~(inds1 | inds2 | inds3));
        
        tmp = [cos(2*pi/numPetals*(k-1)) -sin(2*pi/numPetals*(k-1)) 0; sin(2*pi/numPetals*(k-1)) cos(2*pi/numPetals*(k-1)) 0; 0 0 1]*tmp;
        
        vecPetalArray{k} = {tmp};
    end
end



