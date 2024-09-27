% P3MAP : point 3d map
% MASK : vertex_index_map == true means no information
% EDGEMASK : == true means outlier edge or undefined point

% DMAP : derivative map
% SMAP : left-right layer and top-bottom layer mask, in the first layer
    % == 0 means left, == 1 means right, == -1 means undefined,
    % in the second layer == 0 means top, == 1 means bottom, == -1 means
    % undefined
%%
function [DMAP,SMAP] = derivativeMapGenerator_ver3(P3MAP,MASK,EDGEMASK)

W = size(P3MAP,2);
H = size(P3MAP,1);
dmap = zeros(H,W,6);
smap = -ones(H,W,2);
for r = 1 : H
    for c = 1 : W
        if MASK(r,c)
            continue;
        end
        bi_flag = true;
        % search left and right
        if c < W && c > 1
            if ~MASK(r,c+1) && ~MASK(r,c-1)
                % if the right is smooth or the left is rough
                if ~EDGEMASK(r,c+1) || EDGEMASK(r,c-1)
                    dmap(r,c,1:3) = squeeze(P3MAP(r,c+1,:)-P3MAP(r,c,:));
                    smap(r,c,1) = 1;
                else
                    dmap(r,c,1:3) = squeeze(P3MAP(r,c,:)-P3MAP(r,c-1,:));
                    smap(r,c,1) = 0;
                end
            else
                bi_flag = false;
            end
        end
        if ~bi_flag
            % search right
            if c < W && ~MASK(r,c+1)
                dmap(r,c,1:3) = squeeze(P3MAP(r,c+1,:)-P3MAP(r,c,:));
                smap(r,c,1) = 1;
            % search left
            elseif c > 1 && ~MASK(r,c-1)
                dmap(r,c,1:3) = squeeze(P3MAP(r,c,:)-P3MAP(r,c-1,:));
                smap(r,c,1) = 0;
            % error
            else
                disp("error in calculating derivative map");
                DMAP = [];
                SMAP = [];
                return;
            end
        end
        bi_flag = true;
        % search bottom and top
        if r < H && r > 1
            if ~MASK(r+1,c) && ~MASK(r-1,c)
                % if the bottom is smooth or the top is rough
                if ~EDGEMASK(r+1,c) || EDGEMASK(r-1,c)
                    dmap(r,c,4:6) = squeeze(P3MAP(r+1,c,:)-P3MAP(r,c,:));
                    smap(r,c,2) = 1;
                else
                    dmap(r,c,4:6) = squeeze(P3MAP(r,c,:)-P3MAP(r-1,c,:));
                    smap(r,c,2) = 0;
                end
            else
                bi_flag = false;
            end
        end
        if ~bi_flag
            % search bottom
            if r < H && ~MASK(r+1,c)
                dmap(r,c,4:6) = squeeze(P3MAP(r+1,c,:)-P3MAP(r,c,:));
                smap(r,c,2) = 1;
            % search top
            elseif r > 1 && ~MASK(r-1,c)
                dmap(r,c,4:6) = squeeze(P3MAP(r,c,:)-P3MAP(r-1,c,:));
                smap(r,c,2) = 0;
            % error
            else
                disp("error in calculating derivative map");
                DMAP = [];
                SMAP = [];
                return;
            end
        end
    end
end
DMAP = dmap;
SMAP = smap;

end