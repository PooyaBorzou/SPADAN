function [ ] = exp_afct( reactnts , prtn_num2,act_deact)
%prtn_nums1 : A1,A2,...  prtn_num2 : B
%first : to add A1|A2|... and A1|A2|B_g in nodes if does'nt exist
%second : to add related reaction to reactions   A-> B :  Aactive+B_g ->
%Aactive|B_g  ,  Aactive|B_g -> Aactive + B_RNA
%B_rna _> B




global nodes reactions reaction_type prtns_nodes_adress prtns_activity genes_activity

%add reactnts complex formation

%to find active forms in prtns_activity
for i=1:1:length(reactnts)
    box = prtns_activity(reactnts(i),:);
    box = find(box);
    box = box(1);
    reactnts(i) = prtns_nodes_adress(reactnts(i),box); %the adress of activated protein form in nodes
end
%to make binded reactnts 
if(length(reactnts)>1)
    
    reactnts_name = nodes(reactnts);
    cmplx = reactnts_name{1};
    for i=2:1:length(reactnts_name)
        cmplx = strcat(cmplx,' |',' ',reactnts_name{i});
    end
   Aa =  cmplx;
else
   Aa =  nodes{reactnts};
    
end


AaB_g =  strcat(Aa,' |',' ',nodes{prtn_num2},'_g');
B_g = strcat(nodes{prtn_num2},'_g');
B_rna =  strcat(nodes{prtn_num2},'_rna');

[~,AaB_g_idx]  = ismember(AaB_g,nodes);
[~,B_g_idx]   = ismember(B_g,nodes);
[~,B_rna_idx] = ismember(B_rna,nodes);


if(AaB_g_idx)
    ;
else
    nodes = vertcat(nodes ,AaB_g );
    AaB_g_idx = length(nodes);
end
if(B_g_idx)
    ;
else
    nodes = vertcat(nodes ,B_g );
    B_g_idx = length(nodes);
end
if(B_rna_idx)
    ;
else
    nodes = vertcat(nodes ,B_rna );
    B_rna_idx = length(nodes);
end

%save proteins Gene and RNA adress  in prtns_nodes_adress
prtns_nodes_adress(prtn_num2,2) = B_g_idx;
prtns_nodes_adress(prtn_num2,3) = B_rna_idx;



%% A1a + A2a + B_g  >> A1|A2|B_g
box = zeros(size(reactions,1) , size(reactions,2));
box(1:length(reactnts)+1,1) = [ reactnts ; B_g_idx];
box(1,2) = [ AaB_g_idx ];

eq=0;
for r=1:1:size(reactions,3)
    if(isequal(box,reactions(:,:,r)))
        eq=1;
    end
end
if(~eq)
    reaction_num = size(reactions,3) + 1;
    reactions(1:size(box,1),1:size(box,2),reaction_num) = box;
    reaction_type(reaction_num) = 10; % 10 : complex formation
end
%%
if(act_deact)
    %% Aa|B_g  >> A1 + A2 + B_rna + B_g
    box = zeros(size(reactions,1) , size(reactions,2));
    box(1,1) = [ AaB_g_idx];
    box(1:length(reactnts)+2,2) = [ B_g_idx ; B_rna_idx ; reactnts ];
    
    eq=0;
    for r=1:1:size(reactions,3)
        if(isequal(box,reactions(:,:,r)))
            eq=1;
        end
    end
   
    if(~eq)
        reaction_num = size(reactions,3) + 1;
        reactions(1:size(box,1),1:size(box,2),reaction_num) = box;
        reaction_type(reaction_num) = 21; % 21 : gene expression
    end
else
    %% B_g  >>  B_rna + B_g
    box = zeros(size(reactions,1) , size(reactions,2));
    box(1,1) = [ B_g_idx];
    box(1:2,2) = [ B_g_idx ; B_rna_idx ];
    
    eq=0;
    for r=1:1:size(reactions,3)
        if(isequal(box,reactions(:,:,r)))
            eq=1;
        end
    end
   
    if(~eq)
        reaction_num = size(reactions,3) + 1;
        reactions(1:size(box,1),1:size(box,2),reaction_num) = box;
        reaction_type(reaction_num) = 21; % 21 : gene expression
    end
end  
%% B_rna  >> B + B_rna
    box = zeros(size(reactions,1) , size(reactions,2));
    box(1,1) = [  B_rna_idx];
    box(1:2,2) = [ prtn_num2 , B_rna_idx];
    
    eq=0;
    for r=1:1:size(reactions,3)
        if(isequal(box,reactions(:,:,r)))
            eq=1;
        end
    end
    %
    if(~eq)
        reaction_num = size(reactions,3) + 1;
        reactions(1:size(box,1),1:size(box,2),reaction_num) = box;
        reaction_type(reaction_num) = 30; % 30 : RNA translation
    end


end
