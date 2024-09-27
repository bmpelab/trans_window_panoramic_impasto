% keep the edge with smaller norm of derivative
% only generate the E and DE for those vertex within the field of view

% IMASK : index mask
% BMASK : boundary mask

% E : edge differential matrix
% DE : delta coordinates

%%
function [E,DE] = dM_ver7(HYPERMAP,IMASK,BMASK_LR,BMASK_TB,R,MODE)

HC = size(IMASK,1);
WC = size(IMASK,2);

IMAP = HYPERMAP.vertex_index_map;
IMAP(IMASK) = 0;
valid_flag = logical(IMAP(:)~=0);
N = sum(valid_flag);
IMAP(valid_flag) = 1 : N;
P3MAP = HYPERMAP.point_3d_map;
DMAP = HYPERMAP.derivative_map;
SMAP = HYPERMAP.derivative_direction_map;
lr_valid_flag = true(N,1);
tb_valid_flag = true(N,1);

cols_hp1 = zeros(N,1);
cols_hn1 = zeros(N,1);
cols_vp1 = zeros(N,1);
cols_vn1 = zeros(N,1);
de_h = zeros(N,3);
de_v = zeros(N,3);
for r = 1 : HC
    for c = 1 : WC
        if IMASK(r,c)
            continue;
        end
        index = IMAP(r,c);
        % if ~ccmask(r,c), index is not zero
        % if ccmask(r,c), index may be zero
%         if index == 0
%             continue;
%         end
        % left==0 or right==1
        lr_flag = SMAP(r,c,1);
        % if right
        if lr_flag == 1
            if c == WC || BMASK_LR(r,c)
                lr_valid_flag(index) = false;
            else
                index2 = IMAP(r,c+1);
                if index2 == 0
                    disp("error in calculating E and DE");
                    E = [];
                    DE = [];
                    return;
                end
                cols_hp1(index) = index2;
                cols_hn1(index) = index;
                if MODE==1
                    de_h(index,:) = R*squeeze(P3MAP(r,c+1,:)-P3MAP(r,c,:));
                elseif MODE==2
                    de_h(index,:) = squeeze(DMAP(r,c,1:3));
                end
            end
        % if left
        elseif lr_flag == 0
            if c == 1 || BMASK_LR(r,c)
                lr_valid_flag(index) = false;
            else
                index2 = IMAP(r,c-1);
                if index2 == 0
                    disp("error in calculating E and DE");
                    E = [];
                    DE = [];
                    return;
                end
                cols_hp1(index) = index;
                cols_hn1(index) = index2;
                if MODE==1
                    de_h(index,:) = R*squeeze(P3MAP(r,c,:)-P3MAP(r,c-1,:));
                elseif MODE==2
                    de_h(index,:) = squeeze(DMAP(r,c,1:3));
                end
            end
        % error
        else
            disp("error in calculating E and DE");
            E = [];
            DE = [];
            return;
        end
        % top==0 or bottom==1
        tb_flag = SMAP(r,c,2);
        % if bottom
        if tb_flag == 1
            if r == HC || BMASK_TB(r,c)
                tb_valid_flag(index) = false;
            else
                index2 = IMAP(r+1,c);
                if index2 == 0
                    disp("error in calculating E and DE");
                    E = [];
                    DE = [];
                    return;
                end
                cols_vp1(index) = index2;
                cols_vn1(index) = index;
                if MODE==1
                    de_v(index,:) = R*squeeze(P3MAP(r+1,c,:)-P3MAP(r,c,:));
                elseif MODE==2
                    de_v(index,:) = squeeze(DMAP(r,c,4:6));
                end
            end
        % if top
        elseif tb_flag == 0
            if r == 1 || BMASK_TB(r,c)
                tb_valid_flag(index) = false;
            else
                index2 = IMAP(r-1,c);
                if index2 == 0
                    disp("error in calculating E and DE");
                    E = [];
                    DE = [];
                    return;
                end
                cols_vp1(index) = index;
                cols_vn1(index) = index2;
                if MODE==1
                    de_v(index,:) = R*squeeze(P3MAP(r,c,:)-P3MAP(r-1,c,:));
                elseif MODE==2
                    de_v(index,:) = squeeze(DMAP(r,c,4:6));
                end
            end
        % error
        else
            disp("error in calculating E and DE");
            E = [];
            DE = [];
            return;
        end
    end
end

cols_hp1 = cols_hp1(lr_valid_flag); % hotizontal positive 1
cols_hn1 = cols_hn1(lr_valid_flag); % horizontal negative 1
cols_vp1 = cols_vp1(tb_valid_flag); % vertical positive 1
cols_vn1 = cols_vn1(tb_valid_flag); % vertical negative 1
de_h = de_h(lr_valid_flag,:);
de_v = de_v(tb_valid_flag,:);

n = sum(lr_valid_flag) + sum(tb_valid_flag);

E = sparse([(1:n)';(1:n)'],[cols_hp1;cols_vp1;cols_hn1;cols_vn1],[ones(n,1);-ones(n,1)],n,N);
DE = [de_h;de_v];

end