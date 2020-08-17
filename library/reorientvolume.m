% By Ashkan Pakzad on 17th August 2020
% Reorients input volume from niftiinfo/niftiread into forced LPS
% orientation by applying necessary axes permutations and flips based on
% the metadata affine.

function output = reorientvolume(img, meta)

% get affine matrix
aff_raw_RAS = meta.Transform.T;
% remove origin information
aff_raw_RAS(4,1:3) = [0,0,0];
% remove spacing information
aff_raw_RAS = aff_raw_RAS./abs(aff_raw_RAS);
aff_raw_RAS(isnan(aff_raw_RAS)) = 0;
% convert to LPS affine
aff_raw_LPS = aff_raw_RAS*diag([-1,-1,1,1]);

% take absolute values to figure out anatomical axes
aff_raw_LPS_pos = abs(aff_raw_LPS);

% identify permutation to achieve L/R,P/A,S/I
aff_raw_LPS_pos_red = aff_raw_LPS_pos(1:3,1:3);

newaxes = [1,2,3]; % default
for i = 1:3
    vec = aff_raw_LPS_pos_red(:,i);
    newaxes(i) = find(vec);
end

% identify flips
aff_raw_LPS_red = aff_raw_LPS(1:3,1:3);
aff_LPS_red = zeros(3,3);

% construct new affine with permutation
for i = 1:3
    aff_LPS_red(:,i) = aff_raw_LPS_red(:,newaxes(i));
end

% identfy axes that need flipping
flips = [0,0,0]; % default
for i = 1:3
    if sum(aff_LPS_red(:,i)) == -1
        flips(i) = 1;
    end
end

% apply permutations and flips
output = permute(img, newaxes);
for i = 1:3
    if flips(i) == 1
    output = flip(output, i);
    end
end
end
