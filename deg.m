function [ ] = deg( reactnts,reaction_num)
%edited v6
%to add this reaction as degradation in the eqns(reaction_num)
%A >> null
global eqns kval_addr kval_cnt Vs degnum X_rel_param K

%forward equation
box = strcat(' K(',num2str(reaction_num),',',num2str(1),')');
kval_addr(reaction_num,1) = kval_cnt;
kval_cnt = kval_cnt+1;

box = strcat(box,' * X(',num2str(reactnts),')');


Vs(reaction_num) = strcat(Vs(reaction_num),'V(',num2str(reaction_num),') = ',box,';');

eqns(reactnts) = strcat(eqns(reactnts),' - ' , 'V(',num2str(reaction_num),')');



kk = zeros(size(K,1),size(K,2)); % kk is used to make X_rel_params for ILMA method 
kk(reaction_num,1) =1;
    
    X_rel_param(reactnts(1),1:length(kk(:))) =  ...
    bitor(X_rel_param(reactnts(1),1:length(kk(:))),transpose(kk(:)));





degnum = degnum+1;

end

