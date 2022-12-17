function [ ] = dephos_auto_afct( node_num)
%prtn_num1 :  B_P  
%first : to add related reaction to reactions  B_P  -->  B

global  reactions reaction_type prtns_nodes_adress 

BP_idx = node_num;
id = find(ismember(prtns_nodes_adress(:,4),BP_idx));
B_idx = prtns_nodes_adress(id,1);

   %B_P  >> B
    %% to check whether this interacton exists or not
    box = zeros(size(reactions,1) , size(reactions,2));
    box(1,1) = [ BP_idx];
    box(1,2) = [ B_idx  ];

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
    reaction_type(reaction_num) = 50; % 50 : automatic dephosphorilation
    end

end