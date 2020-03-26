% This function checks for any pixels that are only attached to the mask
% via a corner. If so, they are removed from the mask
% Alasdair Morgan

function [pruned_mask] = maskPruningFunction(mask,xhighest3,yhighest3)

[row,col] = find(mask==1); % find coordinates of 1s in mask
count = 0;

for N = 1:size(row,1) % for every pixel in mask
    
    if [row(N) col(N)] == [xhighest3 yhighest3]; continue; % ignore seedpoint
    else
        test_coords = [row(N) col(N)];
    end
    
    % if none of the pixels directly adjacent to the current pixel are in the mask, make it 0
    if mask(test_coords(1)+1,test_coords(2)) == 0 ...
            && mask(test_coords(1)-1,test_coords(2)) == 0 ...
            && mask(test_coords(1),test_coords(2)+1) == 0 ...
            && mask(test_coords(1),test_coords(2)-1) == 0
        mask(row(N),col(N)) = 0; count = count+1;
    else
    end
end

pruned_mask = mask;

disp([num2str(count) ' edge pixel(s) removed!'])
end