
%v21 01/4/12
%packaging SPADAN 

clear all


%% STEP 0   :  Reading related files and experimental data proteins that exists in signor
%             and to extract related edges
disp('step 0:  Reading related files and experimental data proteins that exists in signor and extracting related edges')


[~,biomart,~]=xlsread('data files/database/biomart_uniprt_ensp_ensg_symbl.xlsx',1);

%list  of mass proteins uniprot codes
[~,Prot_ID,~]=xlsread('data files/input data/proteom_ctrl.xlsx',1);

%list of signor edges and families
[~,signor_edges,~]=xlsread('data files/database/signor_edited.xlsx',1);
signor_edges = signor_edges(2:end,:);

    %filter signor edges
    id_up       = strcmp(signor_edges(:,4), 'up-regulates' ) ;
    id_up_act   = strcmp(signor_edges(:,4), 'up-regulates activity');
    id_down     = strcmp(signor_edges(:,4), 'down-regulates' );
    id_down_act = strcmp(signor_edges(:,4), 'down-regulates activity');

    id_u = id_up|id_up_act;
    id_d =id_down| id_down_act;

    id_up = strfind(signor_edges(:,4),'up-regulates');
    id_u_t = not(cellfun('isempty',id_up));
    id_down = strfind(signor_edges(:,4),'down-regulates');
    id_d_t = not(cellfun('isempty',id_down));
    id_trans    = strcmp(signor_edges(:,5), 'transcriptional regulation' );
    id_direct  = strcmp(signor_edges(:,3), 'YES');
    id_trans = id_direct&id_trans;

    id_phos   = strcmp(signor_edges(:,5), 'phosphorylation' );
    id_dephos = strcmp(signor_edges(:,5), 'dephosphorylation' );
    id_ubiq   = strcmp(signor_edges(:,5), 'ubiquitination' );
    id_cmplx = strcmp(signor_edges(:,4), 'form complex' );
    id_tru = id_trans&id_u_t;
    id_trd = id_trans&id_d_t;

    id = id_phos|id_dephos|id_ubiq|id_cmplx|id_tru|id_trd;
    signor_edges = signor_edges(id,:);
    clear id_up id_up_act id_down id_down_act id_u id_d id_u_t id_down  id_d_t id_trans id_direct id_phos id_dephos id_ubiq id_cmplx id_tru id_trd id

signor_edges_a = signor_edges(:,1);
signor_edges_b = signor_edges(:,2);

[~,signor_family,~]=xlsread('data files/database/signor_families.xlsx',1);
signor_family_members = signor_family(2:end,3:end);
clear signor_family

%find proteins in  signor signor_edges that belong to input data too
signor_mass_id_a = ismember(signor_edges_a,Prot_ID);
signor_mass_id_b = ismember(signor_edges_b,Prot_ID);



id = bitand(signor_mass_id_a, signor_mass_id_b);
edges = signor_edges(find(id),:);
mass_prtns = edges(:,1:2);
mass_prtns = mass_prtns(:);

%find proteins in  signor protein families that belong to input data too
box = signor_family_members(find(ismember(signor_family_members,Prot_ID)));
id = find(~cellfun(@isempty,box));
box = box(id);
mass_prtns = [mass_prtns;box];

%remove duplicates in mass_prtns
%remove duplicates
[~,c,~]=unique(mass_prtns);
tmp=true(size(mass_prtns));
tmp(c)=false;
[mass_prtns{tmp}] = deal('');
mass_prtns  = mass_prtns(c);


clear edges id signor_mass_id_a signor_mass_id_b  signor_edges_a signor_edges_b signor_edges box signor_family_members c


%search complexes that their elements are in mass proteins 
[~,signor_complex,~]=xlsread('data files/database/signor_complexes.xlsx',1);
signor_complex = signor_complex(2:end,:);
complex_uni = signor_complex(:,3:end);
for i=1:1:size(complex_uni,1)
    box = complex_uni(i,:); 
    id = find(~cellfun(@isempty,box));
    box = box(id);
    id = ismember(box,mass_prtns);
    complex_mass(i) = all(id);
end

id = find(complex_mass);            %list of complexes that its elements
mass_complexes =  signor_complex(id,1);  %exists in mass data
                                   
%search families that at least a member of them is an element of  mass  
%and put a member of
%family that exists in mass data instead of whole family
[~,signor_family,~]=xlsread('data files/database/signor_families.xlsx',1);
signor_family_members = signor_family(2:end,3:end);
cnt =1;
for i=1:1:size(signor_family_members,1)
    box = signor_family_members(i,:); 
    id = find(~cellfun(@isempty,box));
    box = box(id);
    id = ismember(box,mass_prtns);
    b = find(id);
   
    if(b)
        b = b(1);
        mass_families(cnt,1) = signor_family(i+1,1);
        family_prtns(cnt,1) = box(b);
        cnt = cnt+1;
    end
    
end

mass_elements = [mass_prtns;mass_complexes;mass_families];

%remove duplicates
[~,c,~]=unique(mass_elements);
tmp=true(size(mass_elements));
tmp(c)=false;
[mass_elements{tmp}] = deal('');
mass_elements  = mass_elements(c);
clear c 




%list of signor edges
[~,signor_edges,~]=xlsread('data files/database/signor_edited.xlsx',1);
signor_edges = signor_edges(2:end,:);
signor_edges_a = signor_edges(:,1);
signor_edges_b = signor_edges(:,2);

%find mass edges in signor
signor_mass_id_a = ismember(signor_edges_a,mass_elements);
signor_mass_id_b = ismember(signor_edges_b,mass_elements);
id = bitand(signor_mass_id_a, signor_mass_id_b);

edges = signor_edges(find(id),:);
elements = edges(:,1:2);
elements = elements(:);

%remove duplicates
[~,c,~]=unique(elements);
tmp=true(size(elements));
tmp(c)=false;
[elements{tmp}] = deal('');
elements  = elements(c);
clear mass_elements b cnt complex_mass complex_uni  id signor_mass_id_a signor_mass_id_b  signor_edges_a signor_edges_b signor_edges box signor_family_members


%to convert protein families with one of their proteins in 'edges'
% and 'elements'
[~,id2] = ismember(mass_families,elements);
for i=1:1:length(id2)
    id1 = find(ismember(edges,mass_families(i)));
    if(id1)
        edges(id1) = family_prtns(i);
    end
    
    if(id2(i))
       elements(id2(i)) =  family_prtns(i);
    end
end

%remove duplicates again (a protein instead of family does not duplicated)
[~,c,~]=unique(elements);
tmp=true(size(elements));
tmp(c)=false;
[elements{tmp}] = deal('');
elements  = elements(c);

%clear all of residual variables
keep_vars = {'edges','elements','mass_complexes','Prot_ID' ,'figure_show','optimizer','init_optimizer','K_id_total','p_id_total','n_id_total'};
clearvars('-except', keep_vars{:});

global prtns_data
prtns_data = elements' ;

xlswrite('data files/output data/PPI_edges.xlsx',edges);

%% STEP 1   :  Cheking activated form of protein based in different edges
%(before starting to add them into 'reactions' )
%find all of upregulation, and down regulation edge adresses
disp('step 1:  Cheking activated form of protein based in different edges')


clear nodes reactions
global prtns_nodes_adress
%column1: indice of protein  in "elements"
%column2: indice of gene in "elements"
%column3: indice of RNA in elements
%column4: indice of phophorilated proteins in elements
%column5: indice of ubiquitinated proteins in elements
%column6-end: indice of related complex in elements


prtns_nodes_adress = zeros(length(elements),6);
prtns_nodes_adress(:,1) = 1:1:length(elements);

%a matrix that shows whether different PTmod proteins are active or not
global prtns_activity
prtns_activity = zeros(size(prtns_nodes_adress,1),size(prtns_nodes_adress,2) );
prtns_activity(:,1) = 1;
%column1: activity of protein"
%column2: null"
%column3: null"
%column4: activity of phophorilated proteins in elements
%column5: activity of ubiquitinated proteins in elements
%column6: null

%a matrix that shows whether different gene complexes are active or not
global genes_activity
genes_activity = zeros(size(prtns_nodes_adress,1),1 );
genes_activity(:,1) = 1;
%column1     : activity of gene
%column2-end : activity of different gene complexes

global nodes reactions degrad_prtns
global reaction_type %this variable contains the type of the equation of each reactio
%it is valed by using afct functions
%functions

%code_phosphorilation =1;     %1: protein phosphorilation
%code_expression = 2;         %2: gene expression
%code_translation = 3;        %3 :gene translation
%code_translation = 4;        %4 :protein ubiquitination
%code_complex_formation = 10;  %10 :complex_formation



nodes = elements;
reactions = zeros(0,0,0);

%to search between edges
id_up       = strcmp(edges(:,4), 'up-regulates' ) ;
id_up_act   = strcmp(edges(:,4), 'up-regulates activity');
id_down     = strcmp(edges(:,4), 'down-regulates' );
id_down_act = strcmp(edges(:,4), 'down-regulates activity');

id_u = id_up|id_up_act;
id_d =id_down| id_down_act;

id_up = strfind(edges(:,4),'up-regulates');
id_u_t = not(cellfun('isempty',id_up));

id_down = strfind(edges(:,4),'down-regulates');
id_d_t = not(cellfun('isempty',id_down));

id_phos   = strcmp(edges(:,5), 'phosphorylation' );
id_dephos = strcmp(edges(:,5), 'dephosphorylation' );

id_ubiq   = strcmp(edges(:,5), 'ubiquitination' );
id_bind   = strcmp(edges(:,5), 'binding' );

%to set phosphorilation prtns_activity
id = id_phos|id_dephos;
elem_phos_dephos = ismember(elements,edges(id,2));

%add phosphorilated in nodes and prtns_nodes_adress
phos_edges = edges(id_phos,1:2);
[~,phos_edges] =  ismember(phos_edges,elements);
[phos_edges,~,~] = unique(phos_edges,'rows','stable');%remove duplicates
for i=1:1:size(phos_edges,1)
    phos_afct(phos_edges(i,1),phos_edges(i,2),1);
end

dephos_edges = edges(id_dephos,1:2);
[~,dephos_edges] =  ismember(dephos_edges,elements);
[dephos_edges,~,~] = unique(dephos_edges,'rows','stable');%remove duplicates
for i=1:1:size(dephos_edges,1)
    phos_afct(dephos_edges(i,1),dephos_edges(i,2),1);
end



% dephos:inactive phos:active
id = (id_d&id_dephos)|(id_u&id_phos); %searches in edges
elem_phos_act = ismember(elements,edges(id,2));
prtns_activity(find(elem_phos_act),1) = 0;
prtns_activity(find(elem_phos_act),4) = 1;

%dephos: active  phos: inactive
elem_phos_inact = elem_phos_dephos&(~elem_phos_act);
prtns_activity(find(elem_phos_inact),1) = 1;
prtns_activity(find(elem_phos_inact),4) = 0;

%to set ubiquitination prtns_activity
elem_ubiq = ismember(elements,edges(id_ubiq,2));

%add ubiquitinated in nodes and prtns_nodes_adress
ubiq_edges = edges(id_ubiq,1:2);
[~,ubiq_edges] =  ismember(ubiq_edges,elements);
[ubiq_edges,~,~] = unique(ubiq_edges,'rows','stable');%remove duplicates
for i=1:1:size(ubiq_edges,1)
    ubiq_afct(ubiq_edges(i,1),ubiq_edges(i,2),1);
end


%normal: active(not special)  ubiq: inactive
id = id_d_t&id_ubiq;
elem_ubiq_inact = ismember(elements,edges(id,2));
%prtns_activity(find(elem_ubiq_inact),1) = 1; it's not needed because
%maybe normal is not active itself
prtns_activity(find(elem_ubiq_inact),5) = 0;

%normal: inactive ubiq: active
elem_ubiq_act = elem_ubiq&(~elem_ubiq_inact);
prtns_activity(find(elem_ubiq_act),1) = 0;
prtns_activity(find(elem_ubiq_act),5) = 1;

%% STEP 2   :  to make "reactions" matrix based on 'edges
disp('step 2 :  Making "reactions" matrix based on edges')

%to find complex forming edges
id1 = strcmp(edges(:,4), 'form complex' );
for i=1:1:length(mass_complexes)
    cmplx = mass_complexes(i);
    id2   = strcmp(edges(:,2),cmplx);
    id = id1&id2;
    reactnts = edges(id,1);
    [~,reactnts_num] = ismember(reactnts,elements);
    [~,product_num] = ismember(cmplx,elements);
    cmplx_afct(reactnts_num,product_num);
end
    
%to add phosphorilation edges to reactions
for i=1:1:size(phos_edges,1)
    id_a = phos_edges(i,1);
    id_b = phos_edges(i,2);
    phos_afct( id_a,id_b,2);
end

%to add dephosphorilation edges to reactions
for i=1:1:size(dephos_edges,1)
    id_a = dephos_edges(i,1);
    id_b = dephos_edges(i,2);
    dephos_afct( id_a,id_b,2);
end

%to add ubiquitination edges to reactions
for i=1:1:size(ubiq_edges,1)
    id_a = ubiq_edges(i,1);
    id_b = ubiq_edges(i,2);
    ubiq_afct( id_a,id_b,2);
end





%transcriptional edges
id_trans    = strcmp(edges(:,5), 'transcriptional regulation' );
id_direct  = strcmp(edges(:,3), 'YES');
id_trans = id_direct&id_trans;

    %transcription up regulation
    id = id_trans&id_u_t;
    trans_u_edges = edges(id,1:2);
    [~,trans_u_edges] =  ismember(trans_u_edges,elements);
    [trans_u_edges,~,~] = unique(trans_u_edges,'rows','stable');%remove duplicates
    
    jump = zeros(0,0);
    
    for i=1:1:size(trans_u_edges,1)
        if(~ismember(i,jump))
            box = trans_u_edges(i,2);
            box = find(ismember(trans_u_edges(:,2),box));
            jump = [jump ; box];
            reactnts = trans_u_edges(box,1);
            exp_afct(reactnts,trans_u_edges(i,2),1);
        end
        
    end
    
     %transcription down regulation
    id = id_trans&id_d_t;
    trans_d_edges = edges(id,1:2);
    [~,trans_d_edges] =  ismember(trans_d_edges,elements);
    [trans_d_edges,~,~] = unique(trans_d_edges,'rows','stable');%remove duplicates
    
    jump = zeros(0,0);
    
    for i=1:1:size(trans_d_edges,1)
        if(~ismember(i,jump))
            box = trans_d_edges(i,2);
            box = find(ismember(trans_d_edges(:,2),box));
            jump = [jump ; box];
            reactnts = trans_d_edges(box,1);
            exp_afct(reactnts,trans_d_edges(i,2),0);
        end
        
    end   
    


% to add degradation rates
%proteins and RNAs that should have degradation rate
global num_prtns
num_prtns = length(elements) - length(mass_complexes); 
degrad_prtns = [prtns_nodes_adress(1:num_prtns,1) ; prtns_nodes_adress(find(prtns_nodes_adress(:,3)),3)];
for i=1:1:length(degrad_prtns)
    deg_afct(degrad_prtns(i));
end

%to add steady production rates for uncontroleesd transcripion proteins
id = id_trans&(id_u_t|id_d_t);
box = edges(id,2);
[box,~,~] = unique(box);%remove duplicates
const_prod_id = find(~ismember(elements(1:num_prtns),box));
for i=1:1:length(const_prod_id)
    prod_afct(const_prod_id(i));
end

%to add automatic dephosphorilation edges to reactions
p = prtns_nodes_adress(:,4);
p_nodes = p(find(p));

for i=1:1:length(p_nodes)
    dephos_auto_afct( p_nodes(i));
end


%clear all of residual variables
keep_vars = {'prtns_nodes_adress','figure_show','optimizer','init_optimizer','num_prtns','edges','elements','mass_complexes','nodes','prtns_activity','prtns_data','','reaction_type','reactions','Prot_ID'};
clearvars('-except', keep_vars{:});


save v21_step2.mat

%  %}
% load('v21_step2.mat'); %skip STEP2

%% step 3  :  extracting network dyamics based on its type
disp('step 3:  Extracting network dyamics based on its type')

eq_type = {''};
global eqns Vs 
eqns =  cell(length(nodes),1);
Vs = cell(size(reactions,3) ,1);



parameters = zeros(size(reactions,3),1); %every line is about a reaction and parameters are soterd in this line


%to assign types of equation to types of reactions
for i=1:1:size(reactions,3)
    switch reaction_type(1,i)
        case 10
            eq_type{i} = 'mass_r'; %complex formation
            
        case 11
            eq_type{i} = 'mass_f'; %dephosphorilation
            
        case 0
            eq_type{i} = 'degradation';%degradation
            
        case 21
            eq_type{i} = 'mass_f'; %gene expression
            
        case 30
            eq_type{i} = 'mass_f'; %RNA translation
        
        case 12
            eq_type{i} = 'mass_f'; %phosphorilation
            
        case 40
            eq_type{i} = 'mass_f'; %ubiquitination
            
      case 100
            eq_type{i} = 'production'; % constant production
            
       case 50
            eq_type{i} = 'mass_f'; % automatic dephosphorilation    
            
        otherwise
            ;
            
    end
    
    
end

global kval_cnt kval_addr  K  X_rel_param
                           %kval_addr is a variable that stores adress of
                           %each element of matrix K in vector kval

K = zeros(size(reactions,3),2);
X_rel_param = zeros(size(nodes,1),length(K(:)) );
%every row is about a variable and every column is a parameter
%number of columns = length(reshaped K)


kval_cnt = 1; %the variable is used in mass2.m to make kval_idx


for i=1:1:length(eq_type)
    reac = reactions(find(reactions(:,1,i)),1,i)';
    prod = reactions(find(reactions(:,2,i)),2,i)';
    
    switch eq_type{i}
        case 'mass_r'
            mass2(reac,prod,i,1);
         
        case 'mass_f'
            mass2(reac,prod,i,0);
            
        case 'degradation'
            deg(reac,i);
            
        case 'production'
            product(prod,i);
            
            
    end
    
    
    
end

%to save A matrix
global A
A = zeros(length(nodes), size(reactions,3));
for i=1:1:size(reactions,3)
    
    %search for positive ones
    strp = strcat('+V(',num2str(i),')');
    box = strfind(eqns,strp);
    box = not(cellfun('isempty',box));
    [id,~] = find(box);
    A(id,i) =A(id,i)+1;
    
    %search for negative ones
    strn = strcat('-V(',num2str(i),')');
    box = strfind(eqns,strn);
    box = not(cellfun('isempty',box));
    [id,~] = find(box);
    A(id,i) =A(id,i)-1;
end
A = sparse(A);

%to make p1-p4 and n1-n4 to use instead of dynamics_auto

K_id_total = zeros(size(K,1),size(K,2));
for i=1:1:length(eq_type)
    if(strcmp(eq_type{i},'mass_f'))
       K_id_total(i,1) = 1;
       
    elseif(strcmp(eq_type{i},'mass_r'))
       K_id_total(i,1) = 1;
       K_id_total(i,2) = -1;
       
    elseif(strcmp(eq_type{i},'degradation'))
       K_id_total(i,1) = -1;
       
    elseif(strcmp(eq_type{i},'production'))
       K_id_total(i,1) = 1;   
    end
end

box=reactions(:,1,:);
for i=1:1:size(box,3)
box2(1:4,i) = box(1:4,1,i);
end
p_id_total = box2';

box=reactions(:,2,:);
for i=1:1:size(box,3)
box2(1:4,i) = box(1:4,1,i);
end
n_id_total = box2';

for i=1:1:length(eq_type)
    if(~strcmp(eq_type{i},'mass_r'))
       n_id_total(i,:) = n_id_total(i,:)*0;
       
    end
end


        
    

save v21_step3_1.mat
%  %}
% load('v21_step3_1.mat'); %skip STEP3

%% step 3.2  : to make mfile dynamics_auto
disp('step 3.2:  Generating an m file containing model dynamics ')


fileID = fopen('dynamics_auto.m','w');
fprintf(fileID,'function X_dot = dynamics_auto(t,X)\n');
fprintf(fileID,'%written automatically after v6\n');
fclose(fileID);
fileID = fopen('dynamics_auto.m','at');
fprintf(fileID, 'global  K A\n');
fprintf(fileID, 'V = zeros(size(K,1),1);\n');

for i=1:1:length(Vs)
    
    if(Vs{i})
        fprintf(fileID,Vs{i});
    else
        fprintf(fileID,'0');
    end
        
    fprintf(fileID,'\n');
    
end

fprintf(fileID, 'X_dot =A*V;\n');
fprintf(fileID,'end');
fclose(fileID);

clear   fileID i kval_cnt parameters prod reac react_node_num  reaction_type box0 box1
save v21_step3_2.mat

% 
% %}
% load('v21_step3_2.mat'); %skip STEP4

%% step 4.1 :  to save  mass experimental data in Wctrl_prtns
disp('step 4.1:  Saving  mass experimental data in Wctrl_prtns')


global  exprmnt_time  mass_name
mass_name = Prot_ID;

exprmnt_time = [0 2 6 24 48];

ensp_g_find3(); %to complete  prtns_data
                 %row1:     uniprot code
                 %row2:     gene symbol
                 %rpw3:     ensg code
                 %row4-end: ensp code for protein and ts same proteins in yasin mass data
                 
                 


                 %seprate prtns_data protein codes and save in ensps
prtns_ensps = zeros(size(prtns_data,1)-3 ,size(prtns_data,2));
for p=1:1:size(prtns_data,2)
    box =  prtns_data(:,p);
    id = find(~cellfun(@isempty,box));
    box = box(id);

box2 =0;    
for j=4:1:size(box)
    S  =  box{j};
    S(isletter(S)) = [];
    box2(j-3,1) = str2double(S);
end

if(box2)
prtns_ensps(1:length(box2),p) = box2;
end
clear box box2 
end
clear p id j S box2

xlswrite('data files\output data/Genes_data.xlsx',prtns_data);

 
%read  experimental data

%step 5.2 : read mass data 
disp('step 4.2:  Reading experimental data ')

%store adress of experimental data for each protein in mass
 for i=1:1:length(elements)

     box = strfind(mass_name,elements(i));
     box = not(cellfun('isempty',box));
     [r,~] = find(box);
     if(r)
     elem_mass_addr(i,1) = r-1;
     end
  end   


 %read control data
[mass_ctrl_r1,~,~] = xlsread('data files/input data/proteom_ctrl.xlsx',2);
[mass_ctrl_r2,~,~] = xlsread('data files/input data/proteom_ctrl.xlsx',3);
[mass_ctrl_r3,~,~] = xlsread('data files/input data/proteom_ctrl.xlsx',4);

 

  global elem_mass_ctrl 
  box = mass_ctrl_r1(elem_mass_addr,:);
  elem_mass_ctrl(1:size(box,1),1:size(box,2),1)=box;
  
  box = mass_ctrl_r2(elem_mass_addr,:);
  elem_mass_ctrl(1:size(box,1),1:size(box,2),2)=box;
  
  box = mass_ctrl_r3(elem_mass_addr,:);
  elem_mass_ctrl(1:size(box,1),1:size(box,2),3)=box;
  clear box r i  mass_ctrl_r1 mass_ctrl_r2 mass_ctrl_r3 
  
  
 %read gef data
[mass_gef_r1,~,~] = xlsread('data files/input data/proteom_g1.xlsx',2);
[mass_gef_r2,~,~] = xlsread('data files/input data/proteom_g1.xlsx',3);
[mass_gef_r3,~,~] = xlsread('data files//input data/proteom_g1.xlsx',4);
   
  global elem_mass_gef 
  box = mass_gef_r1(elem_mass_addr,:);
  elem_mass_gef(1:size(box,1),1:size(box,2),1)=box;
  
  box = mass_gef_r2(elem_mass_addr,:);
  elem_mass_gef(1:size(box,1),1:size(box,2),2)=box;
  
  box = mass_gef_r3(elem_mass_addr,:);
  elem_mass_gef(1:size(box,1),1:size(box,2),3)=box;
  clear box r i  mass_gef_r1 mass_gef_r2 mass_gef_r3 
   
  
 %read plx data
[mass_plx_r1,~,~] = xlsread('data files/input data/proteom_g2.xlsx',2);
[mass_plx_r2,~,~] = xlsread('data files/input data/proteom_g2.xlsx',3);
[mass_plx_r3,~,~] = xlsread('data files/input data/proteom_g2.xlsx',4);
   
  global elem_mass_plx 
  box = mass_plx_r1(elem_mass_addr,:);
  elem_mass_plx(1:size(box,1),1:size(box,2),1)=box;
  
  box = mass_plx_r2(elem_mass_addr,:);
  elem_mass_plx(1:size(box,1),1:size(box,2),2)=box;
  
  box = mass_plx_r3(elem_mass_addr,:);
  elem_mass_plx(1:size(box,1),1:size(box,2),3)=box;
  clear box r i  mass_plx_r1 mass_plx_r2 mass_plx_r3 
   
  
 %read gef_plx data
[mass_gef_plx_r1,~,~] = xlsread('data files/input data/proteom_g3.xlsx',2);
[mass_gef_plx_r2,~,~] = xlsread('data files/input data/proteom_g3.xlsx',3);
[mass_gef_plx_r3,~,~] = xlsread('data files/input data/proteom_g3.xlsx',4);
   
  global elem_mass_gef_plx 
  box = mass_gef_plx_r1(elem_mass_addr,:);
  elem_mass_gef_plx(1:size(box,1),1:size(box,2),1)=box;
  
  box = mass_gef_plx_r2(elem_mass_addr,:);
  elem_mass_gef_plx(1:size(box,1),1:size(box,2),2)=box;
  
  box = mass_gef_plx_r3(elem_mass_addr,:);
  elem_mass_gef_plx(1:size(box,1),1:size(box,2),3)=box;
  clear box r i  mass_gef_plx_r1 mass_gef_plx_r2 mass_gef_plx_r3 
 %-------------------------------------------------------------------
 
 %step 4.3 : read phospho data 
disp('step 4.3:  Reading phospho data ')

[~,phospho_name,~]=xlsread('data files/input data/phospho_ctrl.xlsx',1);

%store adress of experimental data for each protein in phospho
box2 = strfind(phospho_name,'CON_');
box2 = not(cellfun('isempty',box2));

 for i=1:1:length(elements)
     
     box = strfind(phospho_name,elements(i));
     box = not(cellfun('isempty',box));
     [r,~] = find(box & ~box2);
     if(r)
        if(prtns_nodes_adress(i,4)~=0)
            elem_phospho_addr(i,1) = r-1;
        end
     end
  end  




 %read control data
[phospho_ctrl_r1,~,~] = xlsread('data files/input data/phospho_ctrl.xlsx',2);
[phospho_ctrl_r2,~,~] = xlsread('data files/input data/phospho_ctrl.xlsx',3);
[phospho_ctrl_r3,~,~] = xlsread('data files/input data/phospho_ctrl.xlsx',4);

  

  global elem_phospho_ctrl 
  box = phospho_ctrl_r1(nonzeros(elem_phospho_addr),:);
  elem_phospho_ctrl(find(elem_phospho_addr),1:size(phospho_ctrl_r1,2),1)=box;
  
  box = phospho_ctrl_r2(nonzeros(elem_phospho_addr),:);
  elem_phospho_ctrl(find(elem_phospho_addr),1:size(phospho_ctrl_r2,2),2)=box;
  
  box = phospho_ctrl_r3(nonzeros(elem_phospho_addr),:);
  elem_phospho_ctrl(find(elem_phospho_addr),1:size(phospho_ctrl_r3,2),3)=box;
  clear box r i  phospho_ctrl_r1 phospho_ctrl_r2 phospho_ctrl_r3 
  
  
 %read gef data
[phospho_gef_r1,~,~] = xlsread('data files/input data/phospho_g1.xlsx',2);
[phospho_gef_r2,~,~] = xlsread('data files/input data/phospho_g1.xlsx',3);
[phospho_gef_r3,~,~] = xlsread('data files/input data/phospho_g1.xlsx',4);
   
  global elem_phospho_gef 
  box = phospho_gef_r1(nonzeros(elem_phospho_addr),:);
  elem_phospho_gef(find(elem_phospho_addr),1:size(phospho_gef_r1,2),1)=box;
  
  box = phospho_gef_r2(nonzeros(elem_phospho_addr),:);
  elem_phospho_gef(find(elem_phospho_addr),1:size(phospho_gef_r2,2),2)=box;
  
  box = phospho_gef_r3(nonzeros(elem_phospho_addr),:);
  elem_phospho_gef(find(elem_phospho_addr),1:size(phospho_gef_r3,2),3)=box;
  clear box r i  phospho_gef_r1 phospho_gef_r2 phospho_gef_r3 
   
  
 %read plx data
[phospho_plx_r1,~,~] = xlsread('data files/input data/phospho_g2.xlsx',2);
[phospho_plx_r2,~,~] = xlsread('data files/input data/phospho_g2.xlsx',3);
[phospho_plx_r3,~,~] = xlsread('data files/input data/phospho_g2.xlsx',4);
   
  global elem_phospho_plx 
   box = phospho_plx_r1(nonzeros(elem_phospho_addr),:);
  elem_phospho_plx(find(elem_phospho_addr),1:size(phospho_plx_r1,2),1)=box;
  
  box = phospho_plx_r2(nonzeros(elem_phospho_addr),:);
  elem_phospho_plx(find(elem_phospho_addr),1:size(phospho_plx_r2,2),2)=box;
  
  box = phospho_plx_r3(nonzeros(elem_phospho_addr),:);
  elem_phospho_plx(find(elem_phospho_addr),1:size(phospho_plx_r3,2),3)=box;
  clear box r i  phospho_plx_r1 phospho_plx_r2 phospho_plx_r3 
   
  
 %read gef_plx data
[phospho_gef_plx_r1,~,~] = xlsread('data files/input data/phospho_g3.xlsx',2);
[phospho_gef_plx_r2,~,~] = xlsread('data files/input data/phospho_g3.xlsx',3);
[phospho_gef_plx_r3,~,~] = xlsread('data files/input data/phospho_g3.xlsx',4);
   
  global elem_phospho_gef_plx 
   box = phospho_gef_plx_r1(nonzeros(elem_phospho_addr),:);
  elem_phospho_gef_plx(find(elem_phospho_addr),1:size(phospho_gef_plx_r1,2),1)=box;
  
  box = phospho_gef_plx_r2(nonzeros(elem_phospho_addr),:);
  elem_phospho_gef_plx(find(elem_phospho_addr),1:size(phospho_gef_plx_r2,2),2)=box;
  
  box = phospho_gef_plx_r3(nonzeros(elem_phospho_addr),:);
  elem_phospho_gef_plx(find(elem_phospho_addr),1:size(phospho_gef_plx_r3,2),3)=box;
  clear box r i  phospho_gef_plx_r1 phospho_gef_plx_r2 phospho_gef_plx_r3 
   
 
     
%   %to normalize phospho data in which the average in eevery replication be one

%elem_phospho_ctrl  
for k=1:1:3
      for i=1:1:size(elem_phospho_ctrl,1)
          if(find(elem_phospho_ctrl(i,:,k)))
              avg = mean(elem_phospho_ctrl(i,:,k));
              elem_phospho_ctrl(i,:,k) = elem_phospho_ctrl(i,:,k)/avg;
          end
       end
end
  
%elem_phospho_gef  
for k=1:1:3
      for i=1:1:size(elem_phospho_gef,1)
          if(find(elem_phospho_gef(i,:,k)))
              avg = mean(elem_phospho_gef(i,:,k));
              elem_phospho_gef(i,:,k) = elem_phospho_gef(i,:,k)/avg;
          end
       end
end
  
%elem_phospho_plx  
for k=1:1:3
      for i=1:1:size(elem_phospho_plx,1)
          if(find(elem_phospho_plx(i,:,k)))
              avg = mean(elem_phospho_plx(i,:,k));
              elem_phospho_plx(i,:,k) = elem_phospho_plx(i,:,k)/avg;
          end
       end
end
  
%elem_phospho_gef_plx  
for k=1:1:3
      for i=1:1:size(elem_phospho_gef_plx,1)
          if(find(elem_phospho_gef_plx(i,:,k)))
              avg = mean(elem_phospho_gef_plx(i,:,k));
              elem_phospho_gef_plx(i,:,k) = elem_phospho_gef_plx(i,:,k)/avg;
          end
       end
  end
  
clear idx  avg  Wctrl_int Wctrl_name

disp('step 4.4  Reading RNAs experimental data')

%% step 4.4 :  to save RNAs experimental data

[ensg,~,~]=xlsread('data files/input data/RNA_ctrl.xlsx',1);

 %read control data
 
[RNA_ctrl_r1,~,~] = xlsread('data files/input data/RNA_ctrl.xlsx',2);
[RNA_ctrl_r2,~,~] = xlsread('data files/input data/RNA_ctrl.xlsx',3);
[RNA_ctrl_r3,~,~] = xlsread('data files/input data/RNA_ctrl.xlsx',4);

%find protein ensgs in RNA data
 prtns_ensg = prtns_data(3,1:size(prtns_data,2) - length(mass_complexes));
 prtns_ensg = cell2ensp(prtns_ensg)';
[~,box] = ismember(prtns_ensg,ensg);
elem_RNA_addr = box.*(prtns_nodes_adress(1:768,3)~=0);


a = size(RNA_ctrl_r1,2);
elem_RNA_ctrl = zeros(length(elem_RNA_addr),5);
 
elem_RNA_ctrl (find(elem_RNA_addr),1:a,1)=  RNA_ctrl_r1(nonzeros(elem_RNA_addr),:);
elem_RNA_ctrl(find(elem_RNA_addr),1:a,2) =  RNA_ctrl_r2(nonzeros(elem_RNA_addr),:);
elem_RNA_ctrl(find(elem_RNA_addr),1:a,3) =  RNA_ctrl_r3(nonzeros(elem_RNA_addr),:);

clear prtns_ensg RNA_ctrl_r1 RNA_ctrl_r2 RNA_ctrl_r3 biomart i j idx Vs Wctrl_data Wctrl_val Wctrl_name ensg


 %read gef data
 
[RNA_gef_r1,~,~] = xlsread('data files/input data/RNA_g1.xlsx',2);
[RNA_gef_r2,~,~] = xlsread('data files/input data/RNA_g1.xlsx',3);
[RNA_gef_r3,~,~] = xlsread('data files/input data/RNA_g1.xlsx',4);

elem_RNA_gef = zeros(length(elem_RNA_addr),5);

elem_RNA_gef(find(elem_RNA_addr),1:a,1) =  RNA_gef_r1(nonzeros(elem_RNA_addr),:);
elem_RNA_gef(find(elem_RNA_addr),1:a,2) =  RNA_gef_r2(nonzeros(elem_RNA_addr),:);
elem_RNA_gef(find(elem_RNA_addr),1:a,3) =  RNA_gef_r3(nonzeros(elem_RNA_addr),:);

clear prtns_ensg RNA_gef_r1 RNA_gef_r2 RNA_gef_r3 biomart i j idx Vs 


 %read plx data
 
[RNA_plx_r1,~,~] = xlsread('data files/input data/RNA_g2.xlsx',2);
[RNA_plx_r2,~,~] = xlsread('data files/input data/RNA_g2.xlsx',3);
[RNA_plx_r3,~,~] = xlsread('data files/input data/RNA_g2.xlsx',4);

elem_RNA_plx = zeros(length(elem_RNA_addr),5);

elem_RNA_plx(find(elem_RNA_addr),1:a,1) =  RNA_plx_r1(nonzeros(elem_RNA_addr),:);
elem_RNA_plx(find(elem_RNA_addr),1:a,2) =  RNA_plx_r2(nonzeros(elem_RNA_addr),:);
elem_RNA_plx(find(elem_RNA_addr),1:a,3) =  RNA_plx_r3(nonzeros(elem_RNA_addr),:);

clear prtns_ensg RNA_plx_r1 RNA_plx_r2 RNA_plx_r3


 %read gef_plx data
 
[RNA_gef_plx_r1,~,~] = xlsread('data files/input data/RNA_g3.xlsx',2);
[RNA_gef_plx_r2,~,~] = xlsread('data files/input data/RNA_g3.xlsx',3);
[RNA_gef_plx_r3,~,~] = xlsread('data files/input data/RNA_g3.xlsx',4);

elem_RNA_gef_plx = zeros(length(elem_RNA_addr),5);

elem_RNA_gef_plx(find(elem_RNA_addr),1:a,1) =  RNA_gef_plx_r1(nonzeros(elem_RNA_addr),:);
elem_RNA_gef_plx(find(elem_RNA_addr),1:a,2) =  RNA_gef_plx_r2(nonzeros(elem_RNA_addr),:);
elem_RNA_gef_plx(find(elem_RNA_addr),1:a,3) =  RNA_gef_plx_r3(nonzeros(elem_RNA_addr),:);

clear prtns_ensg RNA_gef_plx_r1 RNA_gef_plx_r2 RNA_gef_plx_r3 


save v21_step4_4.mat


% %}
% load('v21_step4_4.mat');%skip above


%% step 4.5 to make exp data ready to be analyzed
disp('step 4.5  Normalizing experimental data')
%X(:,:,1): proteomics 
%X(:,:,2) phospho 
%X(:,:,3) RNA

a= size(elem_mass_ctrl,1);
b= size(elem_mass_ctrl,2);

%ctrl protein data
exp_ctrl(1:a,1:b,1) = (elem_mass_ctrl(:,:,1) + elem_mass_ctrl(:,:,2) + elem_mass_ctrl(:,:,3)) /3;
exp_ctrl_min(1:a,1:b,1) = min(elem_mass_ctrl(:,:,1), min(elem_mass_ctrl(:,:,2) ,elem_mass_ctrl(:,:,3)));
exp_ctrl_max(1:a,1:b,1) = max(elem_mass_ctrl(:,:,1), max(elem_mass_ctrl(:,:,2) ,elem_mass_ctrl(:,:,3)));
%ctrl phospho data
exp_ctrl(1:a,1:b,2) = (elem_phospho_ctrl(:,:,1) + elem_phospho_ctrl(:,:,2) + elem_phospho_ctrl(:,:,3)) /3;
exp_ctrl_min(1:a,1:b,2) = min(elem_phospho_ctrl(:,:,1), min(elem_phospho_ctrl(:,:,2) ,elem_phospho_ctrl(:,:,3)));
exp_ctrl_max(1:a,1:b,2) = max(elem_phospho_ctrl(:,:,1), max(elem_phospho_ctrl(:,:,2) ,elem_phospho_ctrl(:,:,3)));
%ctrl RNA data 
%[nM]  concentration of 1 RNA in a volume of cell in [nM]:  (1 / 6.02E14) / ( 8E-7 uL)
exp_ctrl(1:a,1:b,3)     = 2E-9 * (elem_RNA_ctrl(:,:,1) + elem_RNA_ctrl(:,:,2) + elem_RNA_ctrl(:,:,3)) /3;
exp_ctrl_min(1:a,1:b,3) = 2E-9 *min(elem_RNA_ctrl(:,:,1), min(elem_RNA_ctrl(:,:,2) ,elem_RNA_ctrl(:,:,3)));
exp_ctrl_max(1:a,1:b,3) = 2E-9 *max(elem_RNA_ctrl(:,:,1), max(elem_RNA_ctrl(:,:,2) ,elem_RNA_ctrl(:,:,3))); 


%gef protein data
exp_gef(1:a,1:b,1) = (elem_mass_gef(:,:,1) + elem_mass_gef(:,:,2) + elem_mass_gef(:,:,3)) /3;
exp_gef_min(1:a,1:b,1) = min(elem_mass_gef(:,:,1), min(elem_mass_gef(:,:,2) ,elem_mass_gef(:,:,3)));
exp_gef_max(1:a,1:b,1) = max(elem_mass_gef(:,:,1), max(elem_mass_gef(:,:,2) ,elem_mass_gef(:,:,3)));
%gef phospho data
exp_gef(1:a,1:b,2) = (elem_phospho_gef(:,:,1) + elem_phospho_gef(:,:,2) + elem_phospho_gef(:,:,3)) /3;
exp_gef_min(1:a,1:b,2) = min(elem_phospho_gef(:,:,1), min(elem_phospho_gef(:,:,2) ,elem_phospho_gef(:,:,3)));
exp_gef_max(1:a,1:b,2) = max(elem_phospho_gef(:,:,1), max(elem_phospho_gef(:,:,2) ,elem_phospho_gef(:,:,3)));
%gef RNA data
exp_gef(1:a,1:b,3)     = 2E-9 *(elem_RNA_gef(:,:,1) + elem_RNA_gef(:,:,2) + elem_RNA_gef(:,:,3)) /3;
exp_gef_min(1:a,1:b,3) = 2E-9 *min(elem_RNA_gef(:,:,1), min(elem_RNA_gef(:,:,2) ,elem_RNA_gef(:,:,3)));
exp_gef_max(1:a,1:b,3) = 2E-9 *max(elem_RNA_gef(:,:,1), max(elem_RNA_gef(:,:,2) ,elem_RNA_gef(:,:,3))); 



%plx mass data
exp_plx(1:a,1:b,1) = (elem_mass_plx(:,:,1) + elem_mass_plx(:,:,2) + elem_mass_plx(:,:,3)) /3;
exp_plx_min(1:a,1:b,1) = min(elem_mass_plx(:,:,1), min(elem_mass_plx(:,:,2) ,elem_mass_plx(:,:,3)));
exp_plx_max(1:a,1:b,1) = max(elem_mass_plx(:,:,1), max(elem_mass_plx(:,:,2) ,elem_mass_plx(:,:,3)));
%plx phospho data
exp_plx(1:a,1:b,2) = (elem_phospho_plx(:,:,1) + elem_phospho_plx(:,:,2) + elem_phospho_plx(:,:,3)) /3;
exp_plx_min(1:a,1:b,2) = min(elem_phospho_plx(:,:,1), min(elem_phospho_plx(:,:,2) ,elem_phospho_plx(:,:,3)));
exp_plx_max(1:a,1:b,2) = max(elem_phospho_plx(:,:,1), max(elem_phospho_plx(:,:,2) ,elem_phospho_plx(:,:,3)));
%plx RNA data 
exp_plx(1:a,1:b,3)     = 2E-9 *(elem_RNA_plx(:,:,1) + elem_RNA_plx(:,:,2) + elem_RNA_plx(:,:,3)) /3;
exp_plx_min(1:a,1:b,3) = 2E-9 *min(elem_RNA_plx(:,:,1), min(elem_RNA_plx(:,:,2) ,elem_RNA_plx(:,:,3)));
exp_plx_max(1:a,1:b,3) = 2E-9 *max(elem_RNA_plx(:,:,1), max(elem_RNA_plx(:,:,2) ,elem_RNA_plx(:,:,3))); 



%gef_plx mass data
exp_gef_plx(1:a,1:b,1) = (elem_mass_gef_plx(:,:,1) + elem_mass_gef_plx(:,:,2) + elem_mass_gef_plx(:,:,3)) /3;
exp_gef_plx_min(1:a,1:b,1) = min(elem_mass_gef_plx(:,:,1), min(elem_mass_gef_plx(:,:,2) ,elem_mass_gef_plx(:,:,3)));
exp_gef_plx_max(1:a,1:b,1) = max(elem_mass_gef_plx(:,:,1), max(elem_mass_gef_plx(:,:,2) ,elem_mass_gef_plx(:,:,3)));
%gef_plx phospho data
exp_gef_plx(1:a,1:b,2) = (elem_phospho_gef_plx(:,:,1) + elem_phospho_gef_plx(:,:,2) + elem_phospho_gef_plx(:,:,3)) /3;
exp_gef_plx_min(1:a,1:b,2) = min(elem_phospho_gef_plx(:,:,1), min(elem_phospho_gef_plx(:,:,2) ,elem_phospho_gef_plx(:,:,3)));
exp_gef_plx_max(1:a,1:b,2) = max(elem_phospho_gef_plx(:,:,1), max(elem_phospho_gef_plx(:,:,2) ,elem_phospho_gef_plx(:,:,3)));
%gef_plx RNA data 
exp_gef_plx(1:a,1:b,3)     = 2E-9 *(elem_RNA_gef_plx(:,:,1) + elem_RNA_gef_plx(:,:,2) + elem_RNA_gef_plx(:,:,3)) /3;
exp_gef_plx_min(1:a,1:b,3) = 2E-9 *min(elem_RNA_gef_plx(:,:,1), min(elem_RNA_gef_plx(:,:,2) ,elem_RNA_gef_plx(:,:,3)));
exp_gef_plx_max(1:a,1:b,3) = 2E-9 *max(elem_RNA_gef_plx(:,:,1), max(elem_RNA_gef_plx(:,:,2) ,elem_RNA_gef_plx(:,:,3))); 

clear elem_mass_ctrl elem_mass_gef elem_mass_plx elem_mass_gef_plx
clear elem_phospho_ctrl elem_phospho_gef elem_phospho_plx elem_phospho_gef_plx
clear elem_RNA_ctrl elem_RNA_gef elem_RNA_plx elem_RNA_gef_plx

save v21_step4_5.mat 
%}
%  load('v21_step4_5.mat');%skip above

%% Step 4.6: to set a propper initial value for parameters and nodes
%in this step the values of abundaces are in [nM]
%the volume of the ctrl cell in 0h is estimated as 8E-7 uL by dividing copy
%number to concentration in asolute proteomics
disp('step 4.6:  Finding a propper initial value for parameters and nodes')


%global initial K_G 

 

K = 1e-20 * ones(size(K,1) , size(K,2));
% id =  strcmp(eq_type, 'degradation');
% K(id , 1) = 1; %initial degradation rates



%initial
initial  = zeros(size(nodes,1),1);



initial(nonzeros(prtns_nodes_adress(:,3))) = exp_ctrl(find(prtns_nodes_adress(:,3)),1,3); %initial of RNA of ctrl

%based on papers, the percentage of phosphorylated proteins in all proteome
%is about 30% . so the initial value of phospho is adjusted as 30% of
%initial value of protein abundace
id = find(exp_ctrl(:,1,2));
box = 0*exp_ctrl(:,1,2);
box(id)  = 0.3*exp_ctrl(id,1,1); %initial value of every phosphoprotein

a = prtns_nodes_adress(:,4);
a = a(1:768);
initial( nonzeros(a)) = box(find(a));

initial(1:size(exp_ctrl,1),1) = exp_ctrl(:,1,1); %initial of proteins of ctrl
initial(id)= 0.7* initial(id);




clear initial_mass

K_G = 4E-9 ; %[nM]  concentration of 2 genes in a volume of cell in [nM]:  (2 / 6.02E14) / ( 8E-7 uL)
idx = find(prtns_nodes_adress(:,2));
initial(prtns_nodes_adress(idx,2)) = K_G;

%% Step 4.7: to calculate C_prot matrix  (y = CX)
%the addition of different complex structures of a protein should be equal
%to abundance of that protein. so C_prot will be calculated based on ids of
%elem_mass_ctrl indices
disp('step 4.7:  Calculating C matrix  (y = CX)')


C_prot = sparse(zeros(size(exp_ctrl,1),length(nodes)));
for i=1:1:size(elem_mass_addr,1)
    
        %search for complexes of the protein in nodes
        box = strfind(nodes,elements(i));
        box = not(cellfun('isempty',box));
        [id,~] = find(box);
        C_prot(i,id) =1;
        
        %search of protein ID in signor complexes
        [~,signor_complex,~]=xlsread('data files/database/signor_complexes.xlsx',1);
        [r,~]  = find(strcmp(signor_complex(:,3:end),elements(i)));
        box = signor_complex(r,1);
        id = find(ismember(elements,box));
        if(id)
            C_prot(i,id)  = 1;
        end
end


C_phos = sparse(zeros(size(elem_phospho_addr,1),size(nodes,1)));
for i=1:1:size(elem_phospho_addr,1)
    
        if(prtns_nodes_adress(i,4))
            C_phos(i,prtns_nodes_adress(i,4))  = 1;
        else
            C_phos(i,1:end)  = 0;
        end
end

C_RNA = sparse(zeros(size(elem_phospho_addr,1),length(nodes)));%should be completed
for i=1:1:size(elem_RNA_addr,1)
    
        if(prtns_nodes_adress(i,3))
            C_RNA(i,prtns_nodes_adress(i,3))  = 1;
        else
            C_RNA(i,1:end)  = 0;
        end
end

% load('v21_step6.mat','C_prot','C_phos','C_RNA')

 
%% step 4.8: to evaluate jacobian matrix of model

disp('step 4.8:  Evaluating jacobian matrix of model')

box = length(nodes);
X = sym('X',[box,1]);

% box1 = size(K,1);
% box2 = size(K,2);
% K = sym('K',[box1,box2]);
K = ones( size(K,1),size(K,2) );

p1 = sym(ones(size(K,1),1));
box = find(p_id_total(:,1));
p1(box) = X(nonzeros(p_id_total(:,1)));

%make p2
p2 = sym(ones(size(K,1),1));
box = find(p_id_total(:,2));
p2(box) = X(nonzeros(p_id_total(:,2)));

%make p3
p3 = sym(ones(size(K,1),1));
box = find(p_id_total(:,3));
p3(box) = X(nonzeros(p_id_total(:,3)));

%make p4
p4 = sym(ones(size(K,1),1));
box = find(p_id_total(:,4));
p4(box) = X(nonzeros(p_id_total(:,4)));

%make n1
n1 = sym(ones(size(K,1),1));
box = find(n_id_total(:,1));
n1(box) = X(nonzeros(n_id_total(:,1)));

%make n2
n2 = sym(ones(size(K,1),1));
box = find(n_id_total(:,2));
n2(box) = X(nonzeros(n_id_total(:,2)));

%make n3
n3 = sym(ones(size(K,1),1));
box = find(n_id_total(:,3));
n3(box) = X(nonzeros(n_id_total(:,3)));

%make n4
n4 = sym(ones(size(K,1),1));
box = find(n_id_total(:,4));
n4(box) = X(nonzeros(n_id_total(:,4)));

V = K_id_total(:,1).* K(:,1) .* p1 .*p2.* p3.* p4 + K_id_total(:,2) .* n1 .* n2 .* n3 .* n4;
X_dot_sym =A*V;

jac_sym = jacobian(X_dot_sym,X);

jpat = sparse(zeros(size(jac_sym,1),size(jac_sym,2)));
box = find(jac_sym);
jpat(box) = 1;

save ('v21_step4_8.mat');
clear jac_sym


%  load('v21_step4_8.mat');
%  clear jac_sym

%% step 4.9 performing novel parameter estimation method

global infi out
infi =0;




%load last estimated parameters before interruption in program (imadiate
%shutdwn or etc)
load('iterative_result');
load('x_idx_S_values');
infi =0;

K_string = transpose(k_rec(:,i));


% v19
optimizer =2;

while(resid_tot > 10)
   
    box = X_rel_param(x_num,:);
    check_params = find(box);
    sprintf('X%d',x_num)

    
    switch optimizer
        case 1
            [k_opt,fval,~,~] = optimizer_Ifmin5(K_string(check_params));
        case 2
            [k_opt,fval,~,~] = optimizer_simplex5_2(K_string(check_params));
        case 3
            [k_opt,fval,~,~,~] = optimizer_ILMA5(K_string(check_params));
        otherwise
            ;
    end
    
    
    %{
    %to show the editing parameters
    %------------------------
    K_show = zeros(size(K,1),size(K,2));
    K_show_l = K_show(:);
    K_show_l(check_params) = 1;
    k_l = size(K,1) * size(K,2);
    v1 = K_show_l(1:k_l);
    K_show = reshape(v1,size(K,1),size(K,2));
    [r,c] = find(K_show)
    %------------------------
    %}
    K_string(check_params) = abs(k_opt);
    cnt = cnt + 1;
    
    resid_tot = fval;
    
    
    
    % to calculate the least square
    
    if(infi)
        
        resid_tot = 10000;
        
    else
        
       ;
        
        %{
if(figure_show)
        figure(3)
        hold off
        plot( exprmnt_time, out,'linewidth',1);
        
        hold on
        plot(exprmnt_time',elem_mass_ctrl(:,:,1),'.-.','linewidth',1);
        
        plot( exprmnt_time, elem_mass_ctrl(x_num,:,1),'.','markersize',15,'color' , 'r');
        plot( exprmnt_time, out(x_num,:) ,'color','g','linewidth',2);
        title(sprintf('resid = %.15f  X%d', resid_tot,x_num));
        hold off
        drawnow
    end
        %}
    end
    
     iterative_result(1,cnt) = resid_tot;
    
    % to save the result v20
    
    
    i = size(k_rec,2)+1;
    
    t=now;
    t= datetime(t,'ConvertFrom','datenum');
    time(i) = t;
    
    k_rec(1:length(K_string),i) = K_string;
    
    values(i) = resid_tot;
    
    clear ans t x
    
    

       x_cnt = x_cnt+1;
    
    if(x_cnt > length(x_idx))
        x_cnt = 1;
    end
x_num = x_idx(x_cnt);
 save('iterative_result');
 
end

clear all
load('v20_step6_5.mat');%skip above
 load('v20_step6.mat');
%% step 7.4 to test parallerl optimizer fmincon

optimizer = 4;

k_l = size(K,1) * size(K,2);
load('iterative_result_v20.mat','k_rec');
v1 = k_rec(:,size(k_rec,2));
K = reshape(v1,size(K,1),size(K,2));


%find reaction numbers in which Braf and EGFR play a role into
[~,b] = find(strcmp(prtns_data,'BRAF'));
[~,K_BRAF_id,~]=find(ismember(reactions(:,1,:),b));

[~,b] = find(strcmp(prtns_data,'EGFR'));
[~,K_EGFR_id,~]=find(ismember(reactions(:,1,:),b)); 



prtn_num = 609;
inf_out = 0;

opt = odeset('Jpattern',jpat);



%% Step1.1 calculate output of ctrl
 
    fh = @(t,X) dynamics_auto_v20(K,A,t,X,K_id_total,p_id_total,n_id_total);
    [~,y] = ode15s(fh,exprmnt_time ,initial,opt);
    
    model_prot = C_prot*transpose(y); %calculate proteins model output
    
    model_phos = C_phos*transpose(y); %calculate phospho model output
    box = diag(5./sum(model_phos,2));
    box(isnan(box))=0;
    box(isinf(box))=1;
    model_phos = box*model_phos;  %to divide each row of phospho data to its average
    
    model_RNA =  C_RNA*transpose(y); %calculate RNA  model output
    
    model_ctrl = cat( 3 , model_prot,model_phos,model_RNA);
    
    box = 2*model_ctrl - (exp_ctrl_max+exp_ctrl_min); %negative if exp_ctrl_min<model_ctrl<exp_ctrl_max
    mid_chk = sum(box<0 ,2);
    mid_chk = ~(mid_chk==5); %if all of protein timepoints are between max and min
    
    
    
    
    dif = model_ctrl -  exp_ctrl;
    
    box = diag(5./sum(exp_ctrl(:,:,1),2));
    box(isnan(box))=0;
    box(isinf(box))=1;
    dif(1:end,1:end,1) = box*dif(:,:,1);  %to divide each row of protein dif to the average of proteins experimental data to calculate relative S
    
    box = diag(5./sum(exp_ctrl(:,:,3),2));
    box(isnan(box))=0;
    box(isinf(box))=1;
    dif(1:end,1:end,3) = box*dif(:,:,3);  %to divide each row of RNA dif to the average of RNA experimental data to calculate relative S
    
    
    box = diag(mid_chk(:,:,1)) * dif(:,:,1);% set protein difference zero if it is between borders
    box(1:end,1:end,2) = diag(mid_chk(:,:,2)) * dif(:,:,2);
    box(1:end,1:end,3) = diag(mid_chk(:,:,3)) * dif(:,:,3);
    ctrl_sqr = sum(sum(sum(box.^2)));

    for ii=1:1:768
        
        x = exprmnt_time;
        y = model_prot(ii,:);
        
        yexp = exp_ctrl(ii,:,1);
        yneg = yexp - exp_ctrl_min(ii,:,1);
        ypos = exp_ctrl_max(ii,:,1) - yexp;
        
        ctrl_prot_y(ii,1:5) = y;
        ctrl_prot_yexp(ii,1:5) = yexp;
        ctrl_prot_yneg(ii,1:5) = yneg;
        ctrl_prot_ypos(ii,1:5) = ypos;
        
        %plot phospho
%         subplot(3,1,2)
        
        x = exprmnt_time;
        y = model_phos(ii,:);
        yexp = exp_ctrl(ii,:,2);
        yneg = yexp - exp_ctrl_min(ii,:,2);
        ypos = exp_ctrl_max(ii,:,2) - yexp;

        
        ctrl_phospho_y(ii,1:5) = y;
        ctrl_phospho_yexp(ii,1:5) = yexp;
        ctrl_phospho_yneg(ii,1:5) = yneg;
        ctrl_phospho_ypos(ii,1:5) = ypos;
        
        %plot RNA
%         subplot(3,1,3)
        
        x = exprmnt_time;
        y = model_RNA(ii,:);
        yexp = exp_ctrl(ii,:,3);
        yneg = yexp - exp_ctrl_min(ii,:,3);
        ypos = exp_ctrl_max(ii,:,3) - yexp;
        
        
        
        ctrl_RNA_y(ii,1:5) = y;
        ctrl_RNA_yexp(ii,1:5) = yexp;
        ctrl_RNA_yneg(ii,1:5) = yneg;
        ctrl_RNA_ypos(ii,1:5) = ypos;
        

    end

%% Step1.2 calculate output of gef (EGFRi) copied and uncomplete
% opt = odeset('Jpattern',jpat);
K1 = K; %K values which are related to EGFR reactions became zero
K1(K_EGFR_id,1:2) = K1(K_EGFR_id,1:2)*0;

    fh = @(t,X) dynamics_auto_v20(K1,A,t,X,K_id_total,p_id_total,n_id_total);
    [~,y] = ode15s(fh,exprmnt_time ,initial,opt);
    y_gef =y;
    
    if(size(exp_gef,2) ~= size(y,1)) % if an eror occured in simulation
        inf_out = 1;
        gef_sqr = 1e100;
        
    else
        
        model_prot = C_prot*transpose(y); %calculate proteins model output
        
        model_phos = C_phos*transpose(y); %calculate phospho model output
        box = diag(5./sum(model_phos,2));
        box(isnan(box))=0;
        box(isinf(box))=0;
        model_phos = box*model_phos;  %to divide each row of phospho data to its average
        
        model_RNA =  C_RNA*transpose(y); %calculate RNA  model output
        
        model_gef = cat( 3 , model_prot,model_phos,model_RNA);
        
        
        box = 2*model_gef - (exp_gef_max+exp_gef_min); %negative if exp_gef_min<model_gef<exp_gef_max
        mid_chk = sum(box<0 ,2);
        mid_chk = ~(mid_chk==5); %if all of protein timepoints are between max and min
        
        
        model_gef = full(model_gef);
        dif = model_gef -  exp_gef;
        
        box = diag(5./sum(exp_gef(:,:,1),2));
        box(isnan(box))=0;
        box(isinf(box))=1;
        dif(1:end,1:end,1) = box*dif(:,:,1);  %to divide each row of protein dif to the average of proteins experimental data to calculate relative S
    
        box = diag(5./sum(exp_gef(:,:,3),2));
        box(isnan(box))=0;
        box(isinf(box))=1;
        dif(1:end,1:end,3) = box*dif(:,:,3);  %to divide each row of RNA dif to the average of RNA experimental data to calculate relative S
    
        
        box = diag(mid_chk(:,:,1)) * dif(:,:,1);% set protein difference zero if it is between borders
        box(1:end,1:end,2) = diag(mid_chk(:,:,2)) * dif(:,:,2);
        box(1:end,1:end,3) = diag(mid_chk(:,:,3)) * dif(:,:,3);
        
        gef_sqr = sum(sum(sum(box.^2)));
        
        
    end
    for ii=1:1:768

        x = exprmnt_time;
        y = model_prot(ii,:);
        yexp = exp_gef(ii,:,1);
        yneg = yexp - exp_gef_min(ii,:,1);
        ypos = exp_gef_max(ii,:,1) - yexp;



        gef_prot_y(ii,1:5) = y;
        gef_prot_yexp(ii,1:5) = yexp;
        gef_prot_yneg(ii,1:5) = yneg;
        gef_prot_ypos(ii,1:5) = ypos;
        

        x = exprmnt_time;
        y = model_phos(ii,:);
        yexp = exp_gef(ii,:,2);
        yneg = yexp - exp_gef_min(ii,:,2);
        ypos = exp_gef_max(ii,:,2) - yexp;

        gef_phospho_y(ii,1:5) = y;
        gef_phospho_yexp(ii,1:5) = yexp;
        gef_phospho_yneg(ii,1:5) = yneg;
        gef_phospho_ypos(ii,1:5) = ypos;
        

        x = exprmnt_time;
        y = model_RNA(ii,:);
        yexp = exp_gef(ii,:,3);
        yneg = yexp - exp_gef_min(ii,:,3);
        ypos = exp_gef_max(ii,:,3) - yexp;
        
        
        gef_RNA_y(ii,1:5) = y;
        gef_RNA_yexp(ii,1:5) = yexp;
        gef_RNA_yneg(ii,1:5) = yneg;
        gef_RNA_ypos(ii,1:5) = ypos;
        

    end


%% Step1.3 calculate output of plx (BRAFi) copied and uncomplete

K2 = K; %K values which are related to BRAF reactions became zero
K2(K_BRAF_id,1:2) = K2(K_BRAF_id,1:2)*0;


    fh = @(t,X) dynamics_auto_v20(K2,A,t,X,K_id_total,p_id_total,n_id_total);
    [~,y] = ode15s(fh,exprmnt_time ,initial,opt);
    y_plx=y;
    
    if(size(exp_plx,2) ~= size(y,1)) % if an eror occured in simulation
        inf_out = 1;
        plx_sqr = 1e100;
        
    else
        model_prot = C_prot*transpose(y); %calculate proteins model output
        
        model_phos = C_phos*transpose(y); %calculate phospho model output
        box = diag(5./sum(model_phos,2));
        box(isnan(box))=0;
        box(isinf(box))=0;
        model_phos = box*model_phos;  %to divide each row of phospho data to its average
        
        model_RNA =  C_RNA*transpose(y); %calculate RNA  model output
        
        model_plx = cat( 3 , model_prot,model_phos,model_RNA);
        
        box = 2*model_plx - (exp_plx_max+exp_plx_min); %negative if exp_plx_min<model_plx<exp_plx_max
        mid_chk = sum(box<0 ,2);
        mid_chk = ~(mid_chk==5); %if all of protein timepoints are between max and min
        
        model_plx = full(model_plx);
        dif = model_plx -  exp_plx;
        
        box = diag(5./sum(exp_plx(:,:,1),2));
        box(isnan(box))=0;
        box(isinf(box))=1;
        dif(1:end,1:end,1) = box*dif(:,:,1);  %to divide each row of protein dif to the average of proteins experimental data to calculate relative S
    
        box = diag(5./sum(exp_plx(:,:,3),2));
        box(isnan(box))=0;
        box(isinf(box))=1;
        dif(1:end,1:end,3) = box*dif(:,:,3);  %to divide each row of RNA dif to the average of RNA experimental data to calculate relative S
    
        
        box = diag(mid_chk(:,:,1)) * dif(:,:,1);% set protein difference zero if it is between borders
        box(1:end,1:end,2) = diag(mid_chk(:,:,2)) * dif(:,:,2);
        box(1:end,1:end,3) = diag(mid_chk(:,:,3)) * dif(:,:,3);
        
        plx_sqr = sum(sum(sum(box.^2)));
      
    end
    for ii=1:1:768

        x = exprmnt_time;
        y = model_prot(ii,:);
        yexp = exp_ctrl(ii,:,1);
        yneg = yexp - exp_ctrl_min(ii,:,1);
        ypos = exp_ctrl_max(ii,:,1) - yexp;


        plx_prot_y(ii,1:5) = y;
        plx_prot_yexp(ii,1:5) = yexp;
        plx_prot_yneg(ii,1:5) = yneg;
        plx_prot_ypos(ii,1:5) = ypos;



        x = exprmnt_time;
        y = model_phos(ii,:);
        yexp = exp_plx(ii,:,2);
        yneg = yexp - exp_plx_min(ii,:,2);
        ypos = exp_plx_max(ii,:,2) - yexp;


        plx_phospho_y(ii,1:5) = y;
        plx_phospho_yexp(ii,1:5) = yexp;
        plx_phospho_yneg(ii,1:5) = yneg;
        plx_phospho_ypos(ii,1:5) = ypos;
        

        x = exprmnt_time;
        y = model_RNA(ii,:);
        yexp = exp_plx(ii,:,3);
        yneg = yexp - exp_plx_min(ii,:,3);
        ypos = exp_plx_max(ii,:,3) - yexp;


        plx_RNA_y(ii,1:5) = y;
        plx_RNA_yexp(ii,1:5) = yexp;
        plx_RNA_yneg(ii,1:5) = yneg;
        plx_RNA_ypos(ii,1:5) = ypos;

    end



% to calculate square


%{
%% plot all parts of a protein togrther

    for ii=1:1:768
%%
        scrsz = get(0,'ScreenSize');
        figure('Position',scrsz);
        
        %plot ctrl
        
         suptitle(strcat('Protein',num2str(ii),':  ',prtns_data(2,ii)))
%         %-----------------------
        
%% plot protein of ctrl
        subplot(3,3,1)
        
        x = exprmnt_time;
        y = ctrl_prot_y(ii,:);
        yexp = ctrl_prot_yexp(ii,:);
        yneg = ctrl_prot_yneg(ii,:);
        ypos = ctrl_prot_ypos(ii,:);
        
        
        errorbar(x,yexp,yneg,ypos,'.', 'MarkerSize',20)
        
        hold on
        xlim([0 50])
        plot(x,y)
        hold off
        
        ntitle('Protein of CTRL')

%% plot phospho of ctrl
        subplot(3,3,2)
        
        x = exprmnt_time;
        y = ctrl_phospho_y(ii,:);
        yexp = ctrl_phospho_yexp(ii,:);
        yneg = ctrl_phospho_yneg(ii,:);
        ypos = ctrl_phospho_ypos(ii,:);
        
        
        errorbar(x,yexp,yneg,ypos,'.', 'MarkerSize',20)
        
        hold on
        xlim([0 50])
        plot(x,y)
        hold off
        
        ntitle('Phospho of CTRL')   
    
%% plot RNA of ctrl
        subplot(3,3,3)
        
        x = exprmnt_time;
        y = ctrl_RNA_y(ii,:);
        yexp = ctrl_RNA_yexp(ii,:);
        yneg = ctrl_RNA_yneg(ii,:);
        ypos = ctrl_RNA_ypos(ii,:);
        
        
        errorbar(x,yexp,yneg,ypos,'.', 'MarkerSize',20)
        
        hold on
        xlim([0 50])
        plot(x,y)
        hold off
        
        ntitle('RNA of CTRL')

%% plot protein of gef
        subplot(3,3,4)
        
        x = exprmnt_time;
        y = gef_prot_y(ii,:);
        yexp = gef_prot_yexp(ii,:);
        yneg = gef_prot_yneg(ii,:);
        ypos = gef_prot_ypos(ii,:);
        
        
        errorbar(x,yexp,yneg,ypos,'.', 'MarkerSize',20)
        
        hold on
        xlim([0 50])
        plot(x,y)
        hold off
        
        ntitle('Protein of GEF')

%% plot phospho of gef
        subplot(3,3,5)
        
        x = exprmnt_time;
        y =    gef_phospho_y(ii,:);
        yexp = gef_phospho_yexp(ii,:);
        yneg = gef_phospho_yneg(ii,:);
        ypos = gef_phospho_ypos(ii,:);
        
        
        errorbar(x,yexp,yneg,ypos,'.', 'MarkerSize',20)
        
        hold on
        xlim([0 50])
        plot(x,y)
        hold off
        
        ntitle('Phospho of GEF')   
    
%% plot RNA of gef
        subplot(3,3,6)
        
        x = exprmnt_time;
        y =    gef_RNA_y(ii,:);
        yexp = gef_RNA_yexp(ii,:);
        yneg = gef_RNA_yneg(ii,:);
        ypos = gef_RNA_ypos(ii,:);
        
        
        errorbar(x,yexp,yneg,ypos,'.', 'MarkerSize',20)
        
        hold on
        xlim([0 50])
        plot(x,y)
        hold off
        
        ntitle('RNA of GEF')        
 
%% plot protein of plx
        subplot(3,3,7)
        
        x = exprmnt_time;
        y =    plx_prot_y(ii,:);
        yexp = plx_prot_yexp(ii,:);
        yneg = plx_prot_yneg(ii,:);
        ypos = plx_prot_ypos(ii,:);
        
        
        errorbar(x,yexp,yneg,ypos,'.', 'MarkerSize',20)
        
        hold on
        xlim([0 50])
        plot(x,y)
        hold off
        
        ntitle('Protein of PLX')

%% plot phospho of plx
        subplot(3,3,8)
        
        x = exprmnt_time;
        y =    plx_phospho_y(ii,:);
        yexp = plx_phospho_yexp(ii,:);
        yneg = plx_phospho_yneg(ii,:);
        ypos = plx_phospho_ypos(ii,:);
        
        
        errorbar(x,yexp,yneg,ypos,'.', 'MarkerSize',20)
        
        hold on
        xlim([0 50])
        plot(x,y)
        hold off
        
        ntitle('Phospho of PLX')   
    
%% plot RNA of plx
        subplot(3,3,9)
        
        x = exprmnt_time;
        y =    plx_RNA_y(ii,:);
        yexp = plx_RNA_yexp(ii,:);
        yneg = plx_RNA_yneg(ii,:);
        ypos = plx_RNA_ypos(ii,:);
        
        
        errorbar(x,yexp,yneg,ypos,'.', 'MarkerSize',20)
        
        hold on
        xlim([0 50])
        plot(x,y)
        hold off
        
        ntitle('RNA of PLX')        
 
%% save   
        saveas(figure(1),strcat('plots/',num2str(ii),'_',prtns_data{2,ii},'.jpg'))
        delete(get(groot,'CurrentFigure'))
    
    
    end
%}

%% PLX + EGFR prediction
K3 = K;
K3(K_EGFR_id,1:2) = K3(K_EGFR_id,1:2)*0;
K3(K_BRAF_id,1:2) = K3(K_BRAF_id,1:2)*0;

if(inf_out)
    gef_plx_sqr = 1e100;
else
    fh = @(t,X) dynamics_auto_v20(K3,A,t,X,K_id_total,p_id_total,n_id_total);
    [~,y] = ode15s(fh,exprmnt_time ,initial,opt);
    y_gef_plx=y;
    
    if(size(exp_gef_plx,2) ~= size(y,1)) % if an eror occured in simulation
        inf_out = 1;
        gef_plx_sqr = 1e100;
        
    else
        model_prot = C_prot*transpose(y); %calculate proteins model output
        
        model_phos = C_phos*transpose(y); %calculate phospho model output
        box = diag(5./sum(model_phos,2));
        box(isnan(box))=0;
        box(isinf(box))=0;
        model_phos = box*model_phos;  %to divide each row of phospho data to its average
        
        model_RNA =  C_RNA*transpose(y); %calculate RNA  model output
        
        model_gef_plx = cat( 3 , model_prot,model_phos,model_RNA);
        
        box = 2*model_gef_plx - (exp_gef_plx_max+exp_gef_plx_min); %negative if exp_gef_plx_min<model_gef_plx<exp_gef_plx_max
        mid_chk = sum(box<0 ,2);
        mid_chk = ~(mid_chk==5); %if all of protein timepoints are between max and min
        
        model_gef_plx = full(model_gef_plx);
        dif = model_gef_plx -  exp_gef_plx;
        
        box = diag(5./sum(exp_gef_plx(:,:,1),2));
        box(isnan(box))=0;
        box(isinf(box))=1;
        dif(1:end,1:end,1) = box*dif(:,:,1);  %to divide each row of protein dif to the average of proteins experimental data to calculate relative S
    
        box = diag(5./sum(exp_gef_plx(:,:,3),2));
        box(isnan(box))=0;
        box(isinf(box))=1;
        dif(1:end,1:end,3) = box*dif(:,:,3);  %to divide each row of RNA dif to the average of RNA experimental data to calculate relative S
    
        
        box = diag(mid_chk(:,:,1)) * dif(:,:,1);% set protein difference zero if it is between borders
        box(1:end,1:end,2) = diag(mid_chk(:,:,2)) * dif(:,:,2);
        box(1:end,1:end,3) = diag(mid_chk(:,:,3)) * dif(:,:,3);
        
        gef_plx_sqr = sum(sum(sum(box.^2)));
      
    end
    for ii=1:1:768
%         
%         scrsz = get(0,'ScreenSize');
%         figure('Position',scrsz);
%         
%         
%         
%          suptitle(strcat('Protein',num2str(ii),':  ',prtns_data(2,ii)))
        
        
        
        %plot ctrl
%         suptitle('gef_plx')
        %-----------------------
        
        %plot protein
%          subplot(3,1,1)
        
        x = exprmnt_time;
        y = model_prot(ii,:);
        yexp = exp_ctrl(ii,:,1);
        yneg = yexp - exp_ctrl_min(ii,:,1);
        ypos = exp_ctrl_max(ii,:,1) - yexp;
%          errorbar(x,yexp,yneg,ypos,'o')
         
%         hold on
%         plot(x,y)
%         hold off

%         ntitle('protein')
        gef_plx_prot_y(ii,1:5) = y;
        gef_plx_prot_yexp(ii,1:5) = yexp;
        gef_plx_prot_yneg(ii,1:5) = yneg;
        gef_plx_prot_ypos(ii,1:5) = ypos;


        
%         plot phospho
%         subplot(3,1,2)
        
        x = exprmnt_time;
        y = model_phos(ii,:);
        yexp = exp_gef_plx(ii,:,2);
        yneg = yexp - exp_gef_plx_min(ii,:,2);
        ypos = exp_gef_plx_max(ii,:,2) - yexp;
%         errorbar(x,yexp,yneg,ypos,'o')
        
%         hold on
%         plot(x,y)
%         hold off

%         ntitle('phospho')

        gef_plx_phospho_y(ii,1:5) = y;
        gef_plx_phospho_yexp(ii,1:5) = yexp;
        gef_plx_phospho_yneg(ii,1:5) = yneg;
        gef_plx_phospho_ypos(ii,1:5) = ypos;
        
%         plot RNA
%         subplot(3,1,3)
        
        x = exprmnt_time;
        y = model_RNA(ii,:);
        yexp = exp_gef_plx(ii,:,3);
        yneg = yexp - exp_gef_plx_min(ii,:,3);
        ypos = exp_gef_plx_max(ii,:,3) - yexp;
%         errorbar(x,yexp,yneg,ypos,'o')
        
%         hold on
%         plot(x,y)
%         hold off
%         
%         ntitle('RNA')
        
%         saveas(figure(1),strcat('plots/',num2str(ii),'_',prtns_data{2,ii},'_gef_plx.jpg'))
%         delete(get(groot,'CurrentFigure'))
        
        
        gef_plx_RNA_y(ii,1:5) = y;
        gef_plx_RNA_yexp(ii,1:5) = yexp;
        gef_plx_RNA_yneg(ii,1:5) = yneg;
        gef_plx_RNA_ypos(ii,1:5) = ypos;

    end

end

%% plot 4 groups simultaneousely
clear all
load('v20_step6_5.mat');%skip above
load('v20_step6.mat');
 
scrsz = get(0,'ScreenSize');
h =figure('Position',scrsz);

for ii=1:1:768

suptitle(strcat('Gene',num2str(ii),':  ',prtns_data(2,ii)))

y_label_names = {'Concentration(n[M])', 'Relative Conc','Concentration(n[M])','Concentration(n[M])', 'Relative Conc','Concentration(n[M])','Concentration(n[M])', 'Relative Conc','Concentration(n[M])','Concentration(n[M])', 'Relative Conc','Concentration(n[M])'};

        

%plot ctrl protein
p1=subplot(4,3,1);
errorbar(x,ctrl_prot_yexp(ii,:),ctrl_prot_yneg(ii,:),ctrl_prot_ypos(ii,:),'.', 'MarkerSize',20);

hold on
xlim([0 50])



plot(x,ctrl_prot_y(ii,:))

xlabel('Time (hour)');
ylabel(y_label_names{1});

hold off
ntitle('CTRL Protein')

%plot ctrl phospho
subplot(4,3,2)
errorbar(x,ctrl_phospho_yexp(ii,:),ctrl_phospho_yneg(ii,:),ctrl_phospho_ypos(ii,:),'.', 'MarkerSize',20);

hold on
xlim([0 50])  

plot(x,ctrl_phospho_y(ii,:))
xlabel('Time (hour)');
ylabel(y_label_names{2});


hold off

ntitle('CTRL Phospho')

%plot ctrl RNA
subplot(4,3,3)
errorbar(x,ctrl_RNA_yexp(ii,:),ctrl_RNA_yneg(ii,:),ctrl_RNA_ypos(ii,:),'.', 'MarkerSize',20);

hold on
xlim([0 50])
plot(x,ctrl_RNA_y(ii,:))
xlabel('Time (hour)');
ylabel(y_label_names{3});


hold off

ntitle('CTRL RNA')

%plot gef protein
subplot(4,3,4)
errorbar(x,gef_prot_yexp(ii,:),gef_prot_yneg(ii,:),gef_prot_ypos(ii,:),'.', 'MarkerSize',20);

hold on
xlim([0 50])
plot(x,gef_prot_y(ii,:))
xlabel('Time (hour)');
ylabel(y_label_names{4});


hold off
ntitle('GEF Protein')

%plot GEF phospho
subplot(4,3,5)
errorbar(x,gef_phospho_yexp(ii,:),gef_phospho_yneg(ii,:),gef_phospho_ypos(ii,:),'.', 'MarkerSize',20);

hold on
xlim([0 50])
plot(x,gef_phospho_y(ii,:))
xlabel('Time (hour)');
ylabel(y_label_names{5});


hold off

ntitle('GEF Phospho')

%plot gef RNA
subplot(4,3,6)
errorbar(x,gef_RNA_yexp(ii,:),gef_RNA_yneg(ii,:),gef_RNA_ypos(ii,:),'.', 'MarkerSize',20);

hold on
xlim([0 50])
plot(x,gef_RNA_y(ii,:))
xlabel('Time (hour)');
ylabel(y_label_names{6});


hold off

ntitle('gef RNA')


%plot plx protein
subplot(4,3,7)
errorbar(x,plx_prot_yexp(ii,:),plx_prot_yneg(ii,:),plx_prot_ypos(ii,:),'.', 'MarkerSize',20);

hold on
xlim([0 50])
plot(x,plx_prot_y(ii,:))
xlabel('Time (hour)');
ylabel(y_label_names{7});


hold off
ntitle('PLX Protein')

%plot PLX phospho
subplot(4,3,8)
errorbar(x,plx_phospho_yexp(ii,:),plx_phospho_yneg(ii,:),plx_phospho_ypos(ii,:),'.', 'MarkerSize',20);

hold on
xlim([0 50])
plot(x,plx_phospho_y(ii,:))
xlabel('Time (hour)');
ylabel(y_label_names{8});



hold off

ntitle('PLX Phospho')

%plot PLX RNA
subplot(4,3,9)
errorbar(x,plx_RNA_yexp(ii,:),plx_RNA_yneg(ii,:),plx_RNA_ypos(ii,:),'.', 'MarkerSize',20);

hold on
xlim([0 50])
plot(x,plx_RNA_y(ii,:))
xlabel('Time (hour)');
ylabel(y_label_names{9});


hold off

ntitle('PLX RNA')


%plot gef_plx protein
subplot(4,3,10)
errorbar(x,gef_plx_prot_yexp(ii,:),gef_plx_prot_yneg(ii,:),gef_plx_prot_ypos(ii,:),'.', 'MarkerSize',20);

hold on
xlim([0 50])
plot(x,gef_plx_prot_y(ii,:))
xlabel('Time (hour)');
ylabel(y_label_names{10});


hold off
ntitle('GEF PLX Protein')

%plot GEF PLX phospho
subplot(4,3,11)
errorbar(x,gef_plx_phospho_yexp(ii,:),gef_plx_phospho_yneg(ii,:),gef_plx_phospho_ypos(ii,:),'.', 'MarkerSize',20);

hold on
xlim([0 50])
plot(x,gef_plx_phospho_y(ii,:))
xlabel('Time (hour)');
ylabel(y_label_names{11});


hold off

ntitle('GEF PLX Phospho')

%plot gef PLX RNA
subplot(4,3,12)
errorbar(x,gef_plx_RNA_yexp(ii,:),gef_plx_RNA_yneg(ii,:),gef_plx_RNA_ypos(ii,:),'.', 'MarkerSize',20);

hold on
xlim([0 50])
plot(x,gef_plx_RNA_y(ii,:))
xlabel('Time (hour)');
ylabel(y_label_names{12});


hold off

ntitle('GEF PLX RNA')

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
str = strcat('plots/Gene',num2str(ii));
print(h,str,'-dpdf','-r0')
end


%% to write the list of model reactions as an excel file
global rctns 
rctns = cell(size(reactions,3) ,1);


for i=1:1:length(eq_type)
    reac = reactions(find(reactions(:,1,i)),1,i)';
    prod = reactions(find(reactions(:,2,i)),2,i)';
    
    switch eq_type{i}
        case 'mass_r'
            write_reaction(reac,prod,i,1);
         
        case 'mass_f'
            write_reaction(reac,prod,i,0);
            
        case 'degradation'
            write_reaction(reac,prod,i,0);
            
        case 'production'
            write_reaction(reac,prod,i,0);
            
    end
end

xlswrite('data files\output data/reactions_list.xlsx',rctns);
