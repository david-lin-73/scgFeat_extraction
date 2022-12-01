function matrixVals = indexMatrix(matrix, rowInd, varargin)
% Given a matrix, get the corresponding value for each rowInd without a for
% loop using 2D -> 1D matrix mapping. For example if you had an I wave
% index for all beats and want to quickly get the value of the I wave for
% each beat 
% Inputs 
%       matrix      [NxM]   N signal segments of M lenght each 
%       rowInd      [Mx1]   M row indices 
% Outputs 
%       matrixVals  [Mx1]   values at the row indices 
% Usage : matrixVals = indexMatrix(matrix, rowInd)
Plot = 0;
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp('Plot', varargin{arg}); Plot = 1; end
    end
end
    % determine if the rowInd is the same size as the number of columns 
    if size(rowInd, 1) ~= size(matrix, 2)
        error('rowInd is %d x %d, dim 1 (%d), should be the same length as the number of columns (%d) in the matrix %d x %d',...
            size(rowInd, 1), size(rowInd, 2), size(rowInd, 1), size(matrix, 2), size(matrix, 1), size(matrix, 2))
    end
    
%     %rowInd should be a row vector 
%     if size(rowInd, 2) > size(rowInd, 1)
%         error('rowInd should have indices on the rows and features on the columns')
%     end
    matrixVals = zeros(size(rowInd))*nan;
    nfeats = size(rowInd, 2);
    for feat = 1:nfeats
        % 2D -> 1D indexing and find the value 
        linIdx = sub2ind(size(matrix), rowInd(:, feat)', 1:length(rowInd(:, feat)));
        matrixVals(~isnan(linIdx), feat) = matrix(linIdx(~isnan(linIdx)));

    end
    
    if Plot
        figure; hold on;

        for feat = 1:nfeats
            scatter(rowInd(:, feat), matrixVals(:, feat)); hold on;
        end
        plot(matrix);
        
        h = get(gca, 'Children');
        set(gca, 'Children', [h(end-nfeats+1:end); h(1:end-nfeats)])
    end
    
end