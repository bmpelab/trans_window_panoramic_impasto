% merge P3MAP1 and P3MAP2
% size of P3MAP2 must be equal or smaller than P3MAP1
% P3MAP2 will be centralized at P3MAP1 for merging
% MASK1, MASK2 == true if no info
% OFFSET : image offset
% SF : smooth filter flag
% SMASK : == true means apply smoothing

% MAP3 is the merged MAP

%%
function [P3MAP3,DMAP3,CMAP3,MASK3,BIMAP3] = mergeMap_ver2(P3MAP1,P3MAP2,DMAP1,DMAP2,CMAP1,CMAP2,MASK1,MASK2,TMASK,BIMAP,SF,SMASK)

H1 = size(P3MAP1,1);
W1 = size(P3MAP1,2);
H2 = size(P3MAP2,1);
W2 = size(P3MAP2,2);

dH = round((H1-H2)/2);
dW = round((W1-W2)/2);

mask3 = MASK1;
mask3(dH+1:dH+H2,dW+1:dW+W2) = logical(MASK1(dH+1:dH+H2,dW+1:dW+W2).*MASK2);
MASK3 = mask3;

p3map3 = P3MAP1;
dmap3 = DMAP1;
bimap3 = BIMAP;
for r = dH+1 : dH+H2
    for c = dW+1 : dW+W2
        % if no info in 2, keep the info in 1
        if MASK2(r-dH,c-dW)
            continue;
        end
        % if there is info in 2 but no info in 1, keep the info in 2
        if MASK1(r,c)
            p3map3(r,c,:) = P3MAP2(r-dH,c-dW,:);
            dmap3(r,c,:) = DMAP2(r-dH,c-dW,:);
            bimap3(r,c) = 0;
        % if there is info in 2 and 1
        else
            distance = norm(squeeze(P3MAP2(r-dH,c-dW,:)-p3map3(r,c,:)));
            signed_distance = distance;
            if P3MAP2(r-dH,c-dW,3)>p3map3(r,c,3)
                signed_distance = -signed_distance;
            end
            if signed_distance < 5
                % keep info in 1
                % do nothing
            else
                % keep info in 2
                p3map3(r,c,:) = P3MAP2(r-dH,c-dW,:);
                dmap3(r,c,:) = DMAP2(r-dH,c-dW,:);
                bimap3(r,c) = 0;
            end
        end
    end
end

cmap3 = CMAP1;
for r = dH+1 : dH+H2
    for c = dW+1 : dW+W2
        % if color occluded by tool, keep the info in 1
        if TMASK(r-dH,c-dW)
            continue;
        end
        % if color is non-occluded, keep the info in 2
        cmap3(r,c,:) = CMAP2(r-dH,c-dW,:);
    end
end

% smoothing
p3map3_archived = p3map3;
if SF
    for r = 2 : H1-1
        for c = 2 : W1-1
            if ~SMASK(r,c)
                continue;
            end
            if MASK3(r,c)
                continue;
            end
            mask3_ = MASK3(r-1:r+1,c-1:c+1);
            flag = reshape(~mask3_,[],1);
            if sum(flag)<2
                continue;
            end
            map_ = p3map3_archived(r-1:r+1,c-1:c+1,:);
            eles = reshape(map_,[],size(map_,3));
            mean_ele = mean(eles(flag,:));
            direction = squeeze(map_(2,2,:));
            direction = direction/norm(direction);
            scale = (mean_ele(end)-map_(2,2,end))/direction(end);
            p3map3(r,c,:) = squeeze(map_(2,2,:))+scale*direction;
        end
    end
end

P3MAP3 = p3map3;
CMAP3 = cmap3;
DMAP3 = dmap3;
BIMAP3 = bimap3;

end