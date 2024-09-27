% SFMAP : h1 by w1 by 3 scene flow map
% MASK :h1 by w1 strain mask + inter-frame tool mask, == ture means invalid
    % scene flow
% IMAP : h2 by w2 index map of the merged mesh

% MASK2 : h1 by w1 mask, == true means no direct scene flow
%%
function [MASK2] = directDeformationMask(SFMAP,MASK,IMAP)

W1 = size(SFMAP,2);
H1 = size(SFMAP,1);
W2 = size(IMAP,2);
H2 = size(IMAP,1);
dH = round((H2-H1)/2);
dW = round((W2-W1)/2);

IMAP_ = IMAP(dH+1:dH+H1,dW+1:dW+W1);
omask = logical(~MASK.*IMAP_); % overlap mask, == true means valid scene flow
oflag = omask(:); % reshape(omask,[],1);

sf = reshape(SFMAP,[],3);
sf = sf(oflag,:);
len_sf = sqrt(sum(sf.^2,2));
p1 = 10;
p2 = 90;
t1 = prctile(len_sf,p1);
t2 = prctile(len_sf,p2);
lenflag = len_sf<t2+1.5*(t2-t1);
% t99 = prctile(len_sf,99);
% t1 = 0;
% lenflag = logical((len_sf<t99).*(len_sf>=t1));

flag2 = false(H1*W1,1);
flag2(oflag) = lenflag;
MASK2 = reshape(~flag2,H1,W1);

end