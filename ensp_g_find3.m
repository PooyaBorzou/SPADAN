function [] = ensp_g_find3()
%ensp_g_find finds ensp and ensg codes of prtns_data(1,:)
%   these codes are extracted from '../data files/biomart_uniprt_ensp_ensg_symbl.xlsx'
global prtns_data num_prtns mass_name
load('v6_biomart.mat')
prtns = prtns_data(1,1:num_prtns);

%prtns_data row 2 : gene symbol
for i=1:1:length(prtns)
    b_idx = find(strcmp(biomart(:,4) , prtns(i)))';
    [idx] = b_idx;
    a_idx = biomart(idx(find(idx)),1);
    if(find(idx))
        prtns_data(2,i)= a_idx(1);
    end
    clear idx
end

%prtns_data row 3 : ensg
for i=1:1:length(prtns)
    b_idx = find(strcmp(biomart(:,4) , prtns(i)))';
    [idx] = b_idx;
    a_idx = biomart(idx(find(idx)),2);
    if(find(idx))
        prtns_data(3,i)= a_idx(1);
    end
    clear idx
end

% prtns_data row 4-end  :ensps of proteins in the same group for mass spec
% data

%{
%list of family protein in signor families for each protein
[~,signor_families,~] = xlsread('../data files/signor_families.xlsx',1)
signor_families = signor_families(2:end,4:end);
b = size(signor_families,2);

for i=1:1:size(prtns,2)
    box = strfind(signor_families,prtns(i));
    box = not(cellfun('isempty',box));
    [r,~] = find(box);
    i
    if(r)
        same_family(i,1:b) = signor_families(r(1),1:b);
    end
end
%}

%list of homogroup protein in yasin mass data for each protein
a = size(mass_name,2);

for i=1:1:size(prtns,2)
    box = strfind(mass_name,prtns(i));
    box = not(cellfun('isempty',box));
    [r,~] = find(box);
    i
    if(r)
        same_group(i,1:a) = mass_name(r,1:a);
    end
end



   
for i=1:1:length(prtns)
        box = same_group(i,:);
    
%     [~,c,~]=unique(box);
%     tmp=true(size(box));
%     tmp(c)=false;
%     [box{tmp}] = deal('');
%     box  = box(c);
%     clear c

    id = find(not(cellfun('isempty',box)));
    list = box(id);
    box2 = ones(0,1);
    
    for k=1:1:length(list)
        
        box = find(strcmp(biomart(:,4) , list(k)))';
        [idx] = box;
        if(find(box))
            ensps =  cell2ensp(biomart(box,3));
        end
        box = biomart(idx(find(idx)),3);
        box2 = [box2;box];
    end
        
        if(length(box2))
            prtns_data(4:length(box2)+3,i)= box2;
        end
        
    end

end

