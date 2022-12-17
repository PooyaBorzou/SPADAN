function sqr = eval_sqr_fmin_v20 (k_vect)
%for V7 and later
%eval_sqr used for iterative LMA method

%calculates the sum of sqares of differences between the experimental data
%and simulation based on input parameter vector
%this version is used for iterative method

global K_string check_params x_num elem_mass_ctrl A initial y paral  out kval_addr exprmnt_data exprmnt_time K  K_RNA K_phospho   exprmnt_prtns K_G nodes
global infi %it is added to fsolve to stop it
global K_id_total p_id_total n_id_total jpat exp_ctrl C_prot C_phos C_RNA
global exp_ctrl_max exp_ctrl_min exp_gef exp_gef_max exp_gef_min exp_plx exp_plx_max exp_plx_min K_id K_BRAF_id K_EGFR_id

last_K = K_string(check_params); %to solve the problem of sticking optimizer
infi =0;                          %to solve the problem of sticking optimizer

K_string(check_params) = abs(k_vect);

k_l = size(K,1) * size(K,2);
v1 = K_string(1:k_l);

%fh = @(t,X)  dynamics_auto(K,A,t,X);
fh = @(t,X)  dynamics_auto_v20(K,A,t,X,K_id_total,p_id_total,n_id_total);

K = reshape(v1,size(K,1),size(K,2));


%>>>>>>>>>>>>>>>>
%% Step1.1 calculate output of ctrl
opt = odeset('Jpattern',jpat);
fh = @(t,X) dynamics_auto_v20(K,A,t,X,K_id_total,p_id_total,n_id_total);
[~,y] = ode15s(fh,exprmnt_time ,initial,opt);

if(size(exp_ctrl,2) ~= size(y,1)) % if an eror occured in simulation
    infi = 1;
    sqr = 1e100;
    K_string(check_params) = last_K;
    ctrl_sqr = 1e100;
    
else
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
end

%% Step1.2 calculate output of gef (EGFRi) copied and uncomplete
% opt = odeset('Jpattern',jpat);
K1 = K; %K values which are related to EGFR reactions became zero
K1(K_EGFR_id,1:2) = K1(K_EGFR_id,1:2)*0;

if(infi)
    gef_sqr = 1e100;
else
    fh = @(t,X) dynamics_auto_v20(K1,A,t,X,K_id_total,p_id_total,n_id_total);
    [~,y] = ode15s(fh,exprmnt_time ,initial,opt);
    
    if(size(exp_gef,2) ~= size(y,1)) % if an eror occured in simulation
        infi = 1;
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
        
        clear K1 y model_prot model_phos model_RNA
    end
end

%% Step1.3 calculate output of plx (BRAFi) copied and uncomplete
% opt = odeset('Jpattern',jpat);
K2 = K; %K values which are related to BRAF reactions became zero
K2(K_BRAF_id,1:2) = K2(K_BRAF_id,1:2)*0;

if(infi)
    plx_sqr = 1e100;
else
    fh = @(t,X) dynamics_auto_v20(K2,A,t,X,K_id_total,p_id_total,n_id_total);
    [~,y] = ode15s(fh,exprmnt_time ,initial,opt);
    if(size(exp_plx,2) ~= size(y,1)) % if an eror occured in simulation
        infi = 1;
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
        
        clear K2 y model_prot model_phos model_RNA
    end
end


% to calculate the least square
%% Step1.2 calculate output of gef (EGFRi) copied and uncomplete
% opt = odeset('Jpattern',jpat);
K1 = K; %K values which are related to EGFR reactions became zero
K1(K_EGFR_id,1:2) = K1(K_EGFR_id,1:2)*0;

if(infi)
    gef_sqr = 1e100;
else
    fh = @(t,X) dynamics_auto_v20(K1,A,t,X,K_id_total,p_id_total,n_id_total);
    [~,y] = ode15s(fh,exprmnt_time ,initial,opt);
    
    if(size(exp_gef,2) ~= size(y,1)) % if an eror occured in simulation
        infi = 1;
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
        
        clear K1 y model_prot model_phos model_RNA
    end
end

%% Step1.3 calculate output of plx (BRAFi) copied and uncomplete
% opt = odeset('Jpattern',jpat);
K2 = K; %K values which are related to BRAF reactions became zero
K2(K_BRAF_id,1:2) = K2(K_BRAF_id,1:2)*0;

if(infi)
    plx_sqr = 1e100;
else
    fh = @(t,X) dynamics_auto_v20(K2,A,t,X,K_id_total,p_id_total,n_id_total);
    [~,y] = ode15s(fh,exprmnt_time ,initial,opt);
    if(size(exp_plx,2) ~= size(y,1)) % if an eror occured in simulation
        infi = 1;
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
        
        clear K2 y model_prot model_phos model_RNA
    end
    
    sqr = ctrl_sqr + gef_sqr + plx_sqr;
    %sprintf('ctrl: %1.1f\t gef: %1.1f\t plx: %1.1f\t total: %1.1f\n',ctrl_sqr,gef_sqr,plx_sqr,sqr)
    
    if(infi == 1)
        K_string(check_params) = last_K;
    end
end
