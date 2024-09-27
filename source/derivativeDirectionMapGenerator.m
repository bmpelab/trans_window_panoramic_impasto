% IMASK : ~vertex_index_map == true, undefined
% EDGEMASK : == true, outlier edge or undefined point

% SMAP : left-right layer and top-bottom layer mask, in the first layer
    % == 0 means left, == 1 means right, == -1 means undefined,
    % in the second layer == 0 means top, == 1 means bottom, == -1 means
    % undefined
%%
function [SMAP] = derivativeDirectionMapGenerator(IMASK,EDGEMASK)

W = size(IMASK,2);
H = size(IMASK,1);
smap = -ones(H,W,2);
for r = 1 : H
    for c = 1 : W
        if IMASK(r,c)
            continue;
        end
        bi_flag = true;
        % search left and right
        if c < W && c > 1
            if ~IMASK(r,c+1) && ~IMASK(r,c-1)
                % if the right is smooth or the left is rough
                if ~EDGEMASK(r,c+1) || EDGEMASK(r,c-1)
                    smap(r,c,1) = 1;
                else
                    smap(r,c,1) = 0;
                end
            else
                bi_flag = false;
            end
        end
        if ~bi_flag
            % search right
            if c < W && ~IMASK(r,c+1)
                smap(r,c,1) = 1;
            % search left
            elseif c > 1 && ~IMASK(r,c-1)
                smap(r,c,1) = 0;
            % error
            else
                disp("error in calculating derivative direction map");
                SMAP = [];
                return;
            end
        end
        bi_flag = true;
        % search bottom and top
        if r < H && r > 1
            if ~IMASK(r+1,c) && ~IMASK(r-1,c)
                % if the bottom is smooth or the top is rough
                if ~EDGEMASK(r+1,c) || EDGEMASK(r-1,c)
                    smap(r,c,2) = 1;
                else
                    smap(r,c,2) = 0;
                end
            else
                bi_flag = false;
            end
        end
        if ~bi_flag
            % search bottom
            if r < H && ~IMASK(r+1,c)
                smap(r,c,2) = 1;
            % search top
            elseif r > 1 && ~IMASK(r-1,c)
                smap(r,c,2) = 0;
            % error
            else
                disp("error in calculating derivative direction map");
                SMAP = [];
                return;
            end
        end
    end
end
SMAP = smap;

end