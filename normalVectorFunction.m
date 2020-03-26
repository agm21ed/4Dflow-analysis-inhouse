function [requiresTrigH, requiresTrigW, originalDimensions] = normalVectorFunction(normx,normy,normz)

requiresTrigH = 0; % height
requiresTrigW = 0; % width
originalDimensions = 0; % keep original dimensions or not


if abs(normx) > 0 && abs(normy) > 0 && abs(normz) > 0 % if the normal has 3 nonzero components e.g. [3 1 2]
    requiresTrigH = 1; requiresTrigW = 1; % return;
    
elseif normx == 0 % then the normal is in the y-z plane
    if normy == 0
        if normz == 0 % then there is an error [0 0 0]
            error('normal vector needs components!');
        else % e.g. [0 0 2]
            display('slice conforms with slice selection e.g axial') % then the normal is in the z plane only, meaning the slice is (pretty much) axial
            originalDimensions = 1;
        end
        
    elseif abs(normy) > 0 && abs(normz) > 0 % e.g. [0 4 1]
        requiresTrigH = 1; % just use trig to calculate new pixel height
        
    elseif abs(normy) > 0 && normz == 0 % e.g. [0 3 0]
        display('slice conforms with phase encoding direction e.g coronal') % then the normal is in the y plane only, meaning the slice is (pretty much) coronal
        originalDimensions = 1;  
    end
    
elseif normy == 0
    if normz == 0 % e.g. [3 0 0]
        display('slice conforms with frequency encoding direction e.g sagittal') % then the normal is in the x plane only, meaning the slice is (pretty much) sagittal
        originalDimensions = 1;
    else % has x and z components only e.g. [3 0 5]
        requiresTrigH = 1;
    end
    
else % has x and y components only e.g. [1 4 0]
        requiresTrigW = 1;
        
end

requiresTrigH
requiresTrigW
originalDimensions
        