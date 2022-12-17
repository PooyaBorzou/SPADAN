function [ ] = write_reaction( reactnts,products,reaction_num ,reversible)

global  rctns nodes

%forward equation
if(length(reactnts)==0)
    box_f ='Null';
else
box_f = nodes{reactnts(1)};
end

if(length(reactnts)>1)
    
    for i=2:1:length(reactnts)
        
        box_f = strcat(box_f,' + ',nodes{reactnts(i)});
    end
end

switch reversible
    case 0
        box_f = strcat(box_f,' --> ');
    case 1
        box_f = strcat(box_f,' <--> ');
end
        


if(length(products)== 0)
    box_f = strcat(box_f,'Null');
else 
    box_f = strcat(box_f,nodes{products(1)});
end
    
if(length(products)>1)
    for i=2:1:length(products)
        box_f = strcat(box_f,' + ',nodes{products(i)});
    end
end



box_f = strcat(box_f,' : K(',num2str(reaction_num),',1) for forward direction ');
if (reversible)
    box_f = strcat(box_f,'and K(',num2str(reaction_num),',2) for reverse direction ');
end

rctns(reaction_num) = strcat(rctns(reaction_num),box_f);



end

