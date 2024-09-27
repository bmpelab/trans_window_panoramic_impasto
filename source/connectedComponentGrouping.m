% DDMASK : h1 by w1, == true means no direct defor, directDeformationMask
% IMAP : h2 by w2 hypermap.vertex_index_map

% CCMASK : ccmap whose ccflag == true
% ISO_FLAG : isolated flag, == true means isolated
%%
function [CC_MAP,CC_MASK,CC_GROUP_ID,CC_CENTERS,ISO_FLAG] = connectedComponentGrouping(DDMASK,IMAP)

H1 = size(DDMASK,1);
W1 = size(DDMASK,2);
H2 = size(IMAP,1);
W2 = size(IMAP,2);
dH = round((H2-H1)/2);
dW = round((W2-W1)/2);

% generate the connected components of IMAP
mask = false(H2,W2);
mask(dH+1:dH+H1,dW+1:dW+W1) = true;
IMAP_in = IMAP;
IMAP_in(~mask) = 0;
IMAP_out = IMAP;
IMAP_out(mask) = 0;
ccmap_in = bwlabel(logical(IMAP_in),4);
ccmap_out = bwlabel(logical(IMAP_out),4);
ccmap_out(ccmap_out(:)~=0) = ccmap_out(ccmap_out(:)~=0) + max(ccmap_in(:));
ccmap = ccmap_in + ccmap_out;
CC_MAP = ccmap;
ccflag = false(max(ccmap(:)),1);
cc_centers = zeros(length(ccflag),2); % [cs,rs]
cc_counters = zeros(length(ccflag),1);
for r = 1 : H2
    for c = 1 : W2
        cclabel = ccmap(r,c);
        if cclabel == 0
            continue;
        end
        %
        if r >= dH+1 && r <= dH+H1 && c >= dW+1 && c <= dW+W1
            % if in track
            if ~DDMASK(r-dH,c-dW)
                ccflag(cclabel) = true;
            end
        end
        %
        cc_centers(cclabel,1) = cc_centers(cclabel,1) + c;
        cc_centers(cclabel,2) = cc_centers(cclabel,2) + r;
        cc_counters(cclabel) = cc_counters(cclabel) + 1;
    end
end
cc_centers = cc_centers./cc_counters;
CC_CENTERS = cc_centers;
cc_anchors = cc_centers(ccflag);
ids = 1 : length(ccflag);
cc_anchor_ids = ids(ccflag);
group_ids = (1 : length(ccflag))';
for i = 1 : length(ccflag)
    % if the cc is an anchor
    % its group id is the same as its cclabel
    if ccflag(i)
        continue;
    end
    % find which anchor is the closest one
    cc_center = cc_centers(i,:);
    % actually, here should be Hausdorff distance !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % but the calculation of Hausdorff distance is time consuming
    distances = sqrt(sum((cc_anchors-cc_center).^2,2));
    [~,min_i] = min(distances);
    group_ids(i) = cc_anchor_ids(min_i);
end

ccmask = zeros(size(ccmap));
for i = 1 : length(ccflag)
    if ccflag(i)
        ccmask = ccmask + (ccmap==i);
    end
end
CC_MASK = ~logical(ccmask);

iso_flag = true(length(ccflag),1);
for r = 1 : H2
    for c = 1 : W2
        cclabel = ccmap(r,c);
        if cclabel == 0
            continue;
        end
        % if inside the boundary
        if (r >= dH+1 && r <= dH+H1 && c >= dW+1 && c <= dW+W1) && ~CC_MASK(r,c)
            iso_flag(cclabel) = false;
        % if outside but close to the boundary, and inside is in track
        elseif (r == dH && c >= dW+1 && c <= dW+W1)
            if ~CC_MASK(r+1,c)
                group_ids(cclabel) = CC_MAP(r+1,c);
                iso_flag(cclabel) = false;
            end
        elseif (r == dH+H1+1 && c >= dW+1 && c <= dW+W1)
            if ~CC_MASK(r-1,c)
                group_ids(cclabel) = CC_MAP(r-1,c);
                iso_flag(cclabel) = false;
            end
        elseif (c == dW && r >= dH+1 && r <= dH+H1)
            if ~CC_MASK(r,c+1)
                group_ids(cclabel) = CC_MAP(r,c+1);
                iso_flag(cclabel) = false;
            end
        elseif (c == dW+W1+1 && r >= dH+1 && r <= dH+H1)
            if ~CC_MASK(r,c-1)
                group_ids(cclabel) = CC_MAP(r,c-1);
                iso_flag(cclabel) = false;
            end
        end
    end
end
ISO_FLAG = iso_flag;
CC_GROUP_ID = group_ids;

end