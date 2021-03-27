function [smooth_cellsize,cell_smooth_transposed,rows,columns] = extract_boundary(I,po,nelx,nely,nG)

windowWidth = 45;
polynomialOrder = po;

binaryData = I;
%th = graythresh(binaryData);
binaryData = imresize(binaryData,[nely*nG nelx*nG]);
[rows,columns]=size(binaryData);
labeledData = bwlabel(binaryData,8);
whitemeasurement = regionprops(labeledData,binaryData,'all');
numberofwhite = size(whitemeasurement,1);

%---------------------------------------------------------------------------
% Extract the largest area using our custom function ExtractGeomFeatures().
numberToExtract = numberofwhite;
geom_regions = ExtractGeomFeatures(binaryData, numberToExtract);
%---------------------------------------------------------------------------
%% Boundary extrcation
% Now get the boundaries.
boundaries = bwboundaries(geom_regions);
numberOfBoundaries = size(boundaries, 1);

c = cell(1,numberOfBoundaries);

for k = 1 : numberOfBoundaries
	thisBoundary = boundaries{k};
% 	plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
    c{k} = thisBoundary; % Storing the boundaries in a cell array
end

c_transposed = c';
cellsize = size(c_transposed,1);

figure(3);
subplot(1,1,1);
axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
set(gca, 'YDir','reverse'); %% Makes the figure right side up!! 
cla;
hold on
cell_smooth = cell(1,cellsize);
% simplified_curve = cell(cellsize,1);
%% Performing S-Golay filtering on X and Y co-ordinates of individual boundaries
for j = 1 : cellsize
    cellcontents = c_transposed{j};
    xvalues = cellcontents(:,2); % Extracting X-coordinates from the cell
    yvalues = cellcontents(:,1); % Extracting Y-coordinates from the cell
    
    % Smoothing Boundaries --------------------------------------------
    smoothX = sgolayfilt(xvalues, polynomialOrder, windowWidth); % Smoothen operation on X
    smoothY = sgolayfilt(yvalues, polynomialOrder, windowWidth);% Smoothen operation on Y
    
    % Keep Corners Sharp --------------------------------------------------
    for q = 1:length(xvalues)
        if xvalues(q) == min(xvalues) || xvalues(q) == max(xvalues)
            smoothX(q) = xvalues(q);
        end
        if yvalues(q) == min(yvalues) || yvalues(q) == max(yvalues)
            smoothY(q) = yvalues(q);
        end
    end
    % Reduce Point Density---------------------------------------------
    p = [smoothX,smoothY];
    
    eps = 0.75;
    [ps] = rdp(p,eps);
    
    ps(end,:) = ps(1,:);       
        
    xr = ps(:,1); yr = ps(:,2);
    
    cell_smooth{j} = [xr yr]; % Storing all smooth boundaries in a cell
    plot(xr, yr, 'k', 'LineWidth', 1);
    hold on;
    scatter(xr,yr,'r*');
    hold on;
end
hold off;
% caption = sprintf('Outlines of smooth geometry');
% title({caption,' '});
axis off;

cell_smooth_transposed = cell_smooth';
smooth_cellsize = size(cell_smooth_transposed,1);

%==============================================================================================
% Function to return the specified number of largest or smallest voids in a binary image.
% If numberToExtract > 0 it returns the numberToExtract largest voids.
% If numberToExtract < 0 it returns the numberToExtract smallest voids.
% Example: return a binary image with only the largest void:
%   binaryImage = ExtractGeomFeatures(binaryImage, 1);
% Example: return a binary image with the 3 smallest voids:
%   binaryData = ExtractGeomFeatures(binaryImage, -3);
function binaryData = ExtractGeomFeatures(binaryData, numberToExtract)
	% Get all the void properties.  Can only pass in originalImage in version R2008a and later.
    [labeledData, ~] = bwlabel(binaryData);
	voidMeasurements = regionprops(labeledData, 'area');
	% Get all the areas
	allAreas = [voidMeasurements.Area];
	if numberToExtract > length(allAreas)
		% Limit the number they can get to the number that are there/available.
		numberToExtract = length(allAreas);
	end
	if numberToExtract > 0
		% For positive numbers, sort in order of largest to smallest.
		% Sort them.
		[~, sortIndexes] = sort(allAreas, 'descend');
	elseif numberToExtract < 0
		% For negative numbers, sort in order of smallest to largest.
		% Sort them.
		[~, sortIndexes] = sort(allAreas, 'ascend');
		% Need to negate numberToExtract so we can use it in sortIndexes later.
		numberToExtract = -numberToExtract;
	else
		% numberToExtract = 0.  Shouldn't happen.  Return no voids.
		binaryData = false(size(binaryData));
		return;
	end
	% Extract the "numberToExtract" largest void(a)s using ismember().
	geom_regions = ismember(labeledData, sortIndexes(1:numberToExtract));
	% Convert from integer labeled image into binary (logical) image.
	binaryData = geom_regions > 0;
end
end
