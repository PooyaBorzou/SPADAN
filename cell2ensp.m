function [ ensp_array ] = cell2ensp( cell_array )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    

    for i=1:1:length(cell_array)
        S  =  cell_array{i};
        S(isletter(S)) = [];
        ensp_array(i) = str2double(S);
    end
    


