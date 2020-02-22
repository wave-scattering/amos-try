function [T] = exportTmatrix( stT, format, out )
%% exportTmatrix
% Reshaping to long format and exporting T-matrix entries to a text file
%
% The output format consists of 8 columns:
% s sp m mp n np Tr Ti
% whereby
% * s  1st block index (electric/magnetic)
% * sp 2nd block index (electric/magnetic)
% * m  1st m-index
% * mp 2nd m-index (identical, due to rotational symmetry)
% * n  1st n-index
% * np 2nd n-index
% * Tr real(T_sspmmpnnp)
% * Ti imag(T_sspmmpnnp)
%
% Dependency: 
% none

if(nargin < 3)
    out = [];
end
if(nargin < 2) 
    format = '%d %d %d %d %d %d %.15g %.15g\n';
end

mMax = length(stT);
nrows = 0; % initial count
for m = 1:mMax
    
    tmpoe = stT{m}.('st4MToe');
    tmpeo = stT{m}.('st4MTeo');
    
    M11 = [tmpoe.M11(:) ; tmpeo.M11(:)];
    M12 = [tmpoe.M12(:) ; tmpeo.M12(:)];
    M21 = [tmpoe.M21(:) ; tmpeo.M21(:)];
    M22 = [tmpoe.M22(:) ; tmpeo.M22(:)];
    
    % 11 block
    [n,np] = meshgrid(tmpoe.ind1,tmpoe.ind1);
    [n2,np2] = meshgrid(tmpeo.ind1,tmpeo.ind1);
    vecn = [n(:); n2(:)];
    vecnp = [np(:); np2(:)];
    nbm = numel(tmpoe.M11);
    nbm2 = numel(tmpeo.M11);
    vecm = [repmat(tmpoe.m, nbm, 1); repmat(tmpeo.m, nbm2, 1)];
    vecs = ones(length(vecm), 1);
    Tm11 = [vecs vecs vecm vecm vecn vecnp real(M11) imag(M11)];
    
    % 12 block
    [n,np] = meshgrid(tmpoe.ind1,tmpoe.ind2);
    [n2,np2] = meshgrid(tmpeo.ind1,tmpeo.ind2);
    vecn = [n(:); n2(:)];
    vecnp = [np(:); np2(:)];
    nbm = numel(tmpoe.M12);
    nbm2 = numel(tmpeo.M12);
    vecm = [repmat(tmpoe.m, nbm, 1); repmat(tmpeo.m, nbm2, 1)];
    vecs = ones(length(vecm), 1);
    Tm12 = [vecs 2*vecs vecm vecm vecn vecnp real(M12) imag(M12)];
    
    % 21 block
    [n,np] = meshgrid(tmpoe.ind2,tmpoe.ind1);
    [n2,np2] = meshgrid(tmpeo.ind2,tmpeo.ind1);
    vecn = [n(:); n2(:)];
    vecnp = [np(:); np2(:)];
    nbm = numel(tmpoe.M21);
    nbm2 = numel(tmpeo.M21);
    vecm = [repmat(tmpoe.m, nbm, 1); repmat(tmpeo.m, nbm2, 1)];
    vecs = ones(length(vecm), 1);
    Tm21 = [2*vecs vecs vecm vecm vecn vecnp real(M21) imag(M21)];
    
    % 22 block
    [n,np] = meshgrid(tmpoe.ind2,tmpoe.ind2);
    [n2,np2] = meshgrid(tmpeo.ind2,tmpeo.ind2);
    vecn = [n(:); n2(:)];
    vecnp = [np(:); np2(:)];
    nbm = numel(tmpoe.M22);
    nbm2 = numel(tmpeo.M22);
    vecm = [repmat(tmpoe.m, nbm, 1); repmat(tmpeo.m, nbm2, 1)];
    vecs = ones(length(vecm), 1);
    Tm22 = [2*vecs 2*vecs vecm vecm vecn vecnp real(M22) imag(M22)];
    
    % assign combined values to a cell
    stT{m}.('T') = [Tm11; Tm12; Tm21; Tm22];
    stT{m}.('nrows') = size(stT{m}.('T'), 1);
    nrows = nrows + stT{m}.('nrows');
end

% combine all m T-matrix blocks
T = zeros(nrows, 8);
current = 1;
for ii = 1:mMax
    nc = stT{ii}.('nrows');
    T(current:(current + nc - 1), :) = stT{ii}.('T');
    current = current + nc;
end

% write to a file
if(~isempty(out))
    fileID = fopen(out, 'w');
    fprintf(fileID, '%d elements of T-matrix\n', nrows);
    fprintf(fileID, 's sp m mp n np Tr Ti \n');
    fprintf(fileID, format, T.');
    fclose(fileID);
end



end

