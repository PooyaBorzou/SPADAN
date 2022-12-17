function [ ] = mass2( reactnts,products,reaction_num ,reversible)
%edited v6
%to add this reaction as mass in the eqns(reaction_num)
%   A+B+ ... > C_prot+D + ..
%A_dot = A_dot - k*A*B
%B_bot =B_dot - k*A*B
%C_dot = C_dot + k*A*B
%D_dot = D_dot + k*A*B
% V2: 1)reversible option
%     2)write equations in V form (like my icee paper)

global eqns kval_addr kval_cnt Vs X_rel_param K

%forward equation
box_f = strcat(' K(',num2str(reaction_num),',',num2str(1),')');
kval_addr(reaction_num,1) = kval_cnt;
kval_cnt = kval_cnt+1;

for i=1:1:length(reactnts)
    
    box_f = strcat(box_f,' * X(',num2str(reactnts(i)),')');
end
box = box_f;

kk = zeros(size(K,1),size(K,2)); % kk is used to make X_rel_params for ILMA method 
%reverse equation
if(reversible)
box_r = strcat(' -K(',num2str(reaction_num),',',num2str(2),')');
kval_addr(reaction_num,2) = kval_cnt;
kval_cnt = kval_cnt+1;

  for i=1:1:length(products)
    
    box_r = strcat(box_r,' * X(',num2str(products(i)),')');
  end
box = strcat(box,box_r);
end

Vs(reaction_num) = strcat(Vs(reaction_num),'V(',num2str(reaction_num),') = ',box,';');

for i=1:1:length(reactnts)
    
    eqns(reactnts(i)) = strcat(eqns(reactnts(i)),' - ' , 'V(',num2str(reaction_num), ')');
    
    kk(reaction_num,1) =1;
    if(reversible)
        kk(reaction_num,2) =1;
    end
    
    X_rel_param(reactnts(i),1:length(kk(:))) =  ...
    bitor(X_rel_param(reactnts(i),1:length(kk(:))),transpose(kk(:)));

end

for i=1:1:length(products)
    
    eqns(products(i)) = strcat(eqns(products(i)),' + ' , 'V',num2str(reaction_num));

        kk(reaction_num,1) =1;
    if(reversible)
        kk(reaction_num,2) =1;
    end
    
    X_rel_param(products(i),1:length(kk(:))) =  ...
    bitor(X_rel_param(products(i),1:length(kk(:))),transpose(kk(:)));


end






end

