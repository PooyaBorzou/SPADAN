function [ ] = phos_afct( prtn_num1 , prtn_num2,step_num)
%prtn_num1 : A  prtn_num2 : B  
%first : to add B_p in nodes if does'nt exist
%second : to add related reaction to reactions  A->B  :  A+B  --> A + B_P
%version 2 : 1) A-> B :  Aactive+B -> Aactive|B  ,  Aactive|B -> Aactive + B_p
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
 
AaB =  strcat(Aa,' |',' ',nodes{prtn_num2});
BP =  strcat(nodes{prtn_num2},'_p');
 
[~,Aa_idx]  = ismember(Aa,nodes);
[~,AaB_idx]  = ismember(AaB,nodes);
[~,BP_idx]  = ismember(BP,nodes);
 
 
if(step_num==2)
    if(AaB_idx)
        ;
    else
        nodes = vertcat(nodes ,AaB );
        AaB_idx = length(nodes);
    end
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
 
    %Aa + B  >> Aa|B
    %% to check whether this interacton exists or not
    box = zeros(size(reactions,1) , size(reactions,2));
    box(1:2,1) = [ Aa_idx ; prtn_num2];
    box(1,2) =   [ AaB_idx  ];
 
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
 
    %Aa|B  >> Aa + B_P
    %% to check whether this interacton exists or not
    box = zeros(size(reactions,1) , size(reactions,2));
    box(1,1) = [  AaB_idx];
    box(1:2,2) = [ Aa_idx ; BP_idx];
 
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
    reaction_type(reaction_num) = 12; % 12 : phosphorilation
    end
end
 
end