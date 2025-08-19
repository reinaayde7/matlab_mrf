function [dictc,V,S] = svd_compress_dictionary(dict,thresh)

% Normalize the dictionary
% for i=1:size(dict,2)
%     for j=1:size(dict,3)
%         dict(:,i,j) = dict(:,i,j) / norm(dict(:,i,j),2);
%     end
% end


if nargin < 2
    thresh = 1e-5;
end
dict = dict.'; % N x time
[U,S,V] = svd(dict,'econ');
% if nargin < 2
nk = max(find(diag(S)/S(1)>thresh));
fprintf('SVD dictionary compression: %d --> %d\n',size(dict,2),nk);
% end
V = V(:,1:nk);

dd = diag(S);
totE = sum(power(dd,2));
svdE = sum(power(dd(1:nk),2));
energy = svdE/totE;

figure,
plot(diag(S)/S(1),'o-');
xlim([1 nk+5])
xline(nk,'--g')
title(['Dictionary SVD compression - Energy: ' num2str(energy, '%.3f')]);

dictc = (dict*V(:,1:nk)).';