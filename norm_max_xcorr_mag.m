function [op]=norm_max_xcorr_mag(x,y)

% Calculates the NMCC between two signals x and y

% Normalized maximum cross-correlation magnitude (NMCC) between the
% signals x and y can be obtained by taking the maximum absolute value of 
% cross correlation between the two signals divided by the product of 
% Euclidean norm of the two signals.

% input
% x,y - The two input signals

%output

% op - the NMCC between the two signals
% % cs=abs(xcorr(x,y));
% % cs=sort(cs,'descend');
% % cs=mean(cs(1:10));
op=(max(abs(xcorr(x,y)))/(norm(x,2)*norm(y,2))); % original
% op=cs/(norm(x,2)*norm(y,2));

end