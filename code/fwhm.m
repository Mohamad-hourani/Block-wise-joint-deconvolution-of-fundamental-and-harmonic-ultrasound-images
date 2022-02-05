function fwhmx = fwhm(x,data)
% function width = fwhm(x,y)
% Find the half max value.

halfMax = (min(data) + max(data)) / 2;

% Find where the data first drops below half the max.

index1 = find(data >= halfMax, 1, 'first');

% Find where the data last rises above half the max.

index2 = find(data >= halfMax, 1, 'last');

fwhm = index2-index1 + 1; % FWHM in indexes.

% OR, if you have an x vector

fwhmx = x(index2) - x(index1);
end