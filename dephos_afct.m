function [ ] = dephos_afct( prtn_num1 ,prtn_num2, step_num)
%prtn_num1 : A  prtn_num2 : B_P  
%first : to add B_p in nodes if does'nt exist
%second : to add related reaction to reactions  A->B  :  A+B  --> A + B_P
%version 2 : 1) A-> B :  Aactive+BP -> Aactive|BP  ,  Aactive|BP -> Aactive + B
%Step_num =1 > add nodes step_num = 2 > add reactions   



global nodes reactions reaction_type prtns_nodes_adress prtns_activity

A = nodes{prtn_num1};
B = nodes{prtn_num2};

%to find the active form of A
box = prtns_activity(prtn_num1,:);
box = find(box);
box = box(1);

id = prtns_nodes_adress(prtn_num1,box); %the adress of activated protein form in nodes
Aa = nodes{id};

BP =  strcat(nodes{prtn_num2},'_p');
AaBP =  strcat(Aa,' |',' ',BP);

[~,Aa_idx]  = ismember(Aa,nodes);
[~,AaBP_idx]  = ismember(AaBP,nodes);
[~,BP_idx]  = ismember(BP,nodes);



if(AaBP_idx)
    ;
else
    nodes = vertcat(nodes ,AaBP );
    AaBP_idx = length(nodes);
end
    
if(BP_idx)
    ;
else
    nodes = vertcat(nodes ,BP);
    BP_idx = length(nodes);
end

%save phosphorilated proteins adress  in prtns_nodes_adress
prtns_nodes_adress(prtn_num2,4) = BP_idx;



if(step_num==2)

    %Aa + B_P  >> Aa|B_p
    %% to check whether this interacton exists or not
    box = zeros(size(reactions,1) , size(reactions,2));
    box(1:2,1) = [ Aa_idx ; BP_idx];
    box(1,2) =   [ AaBP_idx  ];

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

    %Aa|B_P  >> Aa + B
    %% to check whether this interacton exists or not
    box = zeros(size(reactions,1) , size(reactions,2));
    box(1,1) = [AaBP_idx];
    box(1:2,2) = [ Aa_idx ; prtn_num2 ];
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
    reaction_type(reaction_num) = 11; % 11 : dephosphor
    end
end

end