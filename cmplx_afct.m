function [ ] = cmplx_afct( reactnts , cmplx_num ,rev)
%reactnts : [A1, A2, ...]  cmplx : B  
%first : A1 + A2 + ...  -->  B
% rev added in inputs ; if 1 the reversed reactions are
%                       determined and if 0 not



%  to convert an  chea edge to its interactions and save its adress in prtns_nodes_adress
global reactions reaction_type prtns_nodes_adress prtns_activity

%to find active forms in prtns_activity
for i=1:1:length(reactnts)
    box = prtns_activity(reactnts(i),:);
    box = find(box);
    box = box(1);
    reactnts(i) = prtns_nodes_adress(reactnts(i),box); %the adress of activated protein form in nodes
end


%A1 + A2 + ... --> B
%% to check whether this interacton exists or not
box = zeros(size(reactions,1) , size(reactions,2));
box(1:length(reactnts),1) = reactnts;
box(1,2) =  cmplx_num;

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
reaction_type(reaction_num) = 10; % 10 : complex formation
end

prtns_nodes_adress(reactnts,6) = cmplx_num;
end
