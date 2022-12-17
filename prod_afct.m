function [ ] = prod_afct( prtn_num)
%  to make an degradation interaction
global nodes reactions reaction_type
A = nodes{prtn_num};
A_idx  = prtn_num;
%   null >> A
%% to check whether this interacton exists or not
box = zeros(size(reactions,1) , size(reactions,2));
box(1,1) = [ 0 ];
box(1,2) =   [ A_idx  ];

eq=0;
for r=1:1:size(reactions,3)
    if(isequal(box,reactions(:,:,r)))
        eq=1;
    end
end
%%
if(~eq)
    reaction_num = size(reactions,3) + 1;
    reactions(1:size(box,1),1:size(box,2),reaction_num) = box;
    reaction_type(reaction_num) = 100; % 100 : constant production
end
end