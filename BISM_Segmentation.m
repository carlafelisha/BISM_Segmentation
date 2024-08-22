%% Bayesian Inversion Stress Microscopy (Value Extraction)
% BISM Section is modified from Nier, V. et al., Inference of internal stress in a cell monolayer, Biophysical Journal, 110(7), 1625-1635, 2016
    
clear all;
close all;
    
%% Folder of JPG spheroid image sequence
dirPath = '/spheroid_img/';

%% Traction Force Microscopy Raw Data
DirectoryName = '';
Sessions = {''};

%% BISM Code
for session_num = 1:length(Sessions) 
    Session = Sessions{session_num} ;
 
    k0=1; % number of time steps
    ForceName = ['/Users/' DirectoryName '/' Session '/Traction_PIV' num2str(k0) '.txt'];
    e = dlmread(ForceName);
    
    R=sqrt(size(e,1)); % number of rows (y direction)
    C=sqrt(size(e,1)); % number of columns (x direction)
    N=R*C; % system size for traction force field
    Ninf=4*N+2*(C+R); % system size for stress tensor field
     
    %% Define spatial grid (rectangular)
    coeff=0.325; % 40x
    % coeff=0.216 ; % 60x 
    % PIV_step_size*Conv_pixtoum; PIV small step size * micromater per pixel 
    
    x =  e(1:71,1)*coeff;
    xmin = min(x);
    xmax = max(x);
    l = x(C)-x(C-1);
    y = x;
    ymin = xmin;
    ymax = xmax;
    [xg, yg]=meshgrid(x,y);
     
    %% Parameters
    BC=0; 
    meth_Lambda=3; 
    noise_value=0.; 
    fplot=1; 
    mult=10^(0); 
    meth_avg = 0;
    edgesize = 10;
    
    kM = 240;

    convertkPa = 1e-3;
    Pmin = -15;
    Pmax = 15;
    Tmin = -2;
    Tmax = 2;
    meanmin = -1200;
    meanmax = 1200;
    % graphics
    tmean = zeros(1,kM);
    dmean = zeros(1,kM);
    smean = zeros(1,kM);
    tstd = zeros(1,kM);
    dstd = zeros(1,kM);
    sstd = zeros(1,kM);
    % mean values + std dev

    %% Excel sheet columns
    timeframe = zeros(1, kM);
    a_average = zeros(1, kM);
    b_average = zeros(1, kM);
    c_average = zeros(1, kM);

    %% Building 2N*Ninf matrix A such that A*\sigma=T (discretized divergence operator)
    % Building matrix Ax (derivative wrt x)
    Ax1=diag(ones(C-1, 1), 1)-diag(ones(C, 1), 0); %bulk (order 2)
    Ax1=[Ax1,[zeros(C-1, 1); 1]];
    Ax=sparse(Ax1);
    for i=1:R-1
        Ax=blkdiag(Ax,Ax1); %iteration R times for all the rows
    end
    clear Ax1;
     
    % Building matrix Ay (derivative wrt y)
    Ay1=diag(ones((R-1)*C,1), C)-diag(ones(R*C,1 ), 0); %bulk (order 2)
    Ay=sparse([Ay1, [zeros((R-1)*C,C); eye(C)]]);
    clear Ay1;
     
    % Building A from Ax and Ay
    A=[Ax, sparse(N,N+C) Ay, sparse(N,N+R); sparse(N,N+R) Ay, sparse(N,N+C), Ax]/l;
    clear Ax Ay; 
     
    %% Building Ninf*Ninf B covariance matrix of the prior, S_0^{-1}=B/s_0^2
    %% Stress norm regularization
    B=speye(Ninf);
     
    %% Add to the prior a term enforcing the equality of shear components
    d1=diag(ones((R-1)*C, 1), C)+diag(ones(R*C,1), 0); %bulk (order 2)
    bd1=[d1,[zeros((R-1)*C, C); eye(C)]];
    clear d1;
    d2=-diag(ones(C-1, 1), 1)-diag(ones(C, 1),0); %bulk (order 2)
    d2=[d2,[zeros(C-1, 1); -1]];
    bd2=d2;
    % loop on rows
    for i=1:R-1
       bd2=blkdiag(bd2, d2); 
    end
    clear d2;
    Bdiff=sparse([zeros(N, 2*N+C+R), bd1, bd2]);
    clear bd1 bd2;
    alpha_xy=10^3; 
    % hyperparameter for the equality of shear components (has to be large)
    B=B+alpha_xy^2*sparse(Bdiff'*Bdiff);
    clear Bdiff;
     
    %% Add to the prior a term enforcing the free stress boundary conditions
    if BC==1
        Bbcx=zeros(2*C, C*(R+1));
        Bbcy=zeros(2*R, (C+1)*R);   
        Bbcx(1:C,1:C)=speye(C);
        Bbcx(C+1:2*C, 1+R*C:(R+1)*C)=speye(C);
        for i=1:R
             Bbcy(i, 1+(i-1)*(C+1))=1;
             Bbcy(i+R, C+1+(i-1)*(C+1))=1;
        end
        Bbc=sparse([Bbcy, sparse(2*R,(C+1)*R+2*C*(R+1)); sparse(2*C,(C+1)*R), Bbcx, sparse(2*C,(C+1)*R+C*(R+1)); sparse(2*C,(C+1)*R+C*(R+1)), Bbcx, sparse(2*C,(C+1)*R); sparse(2*R,(C+1)*R+2*C*(R+1)), Bbcy]); %xx, yy, yx and xy components
        clear Bbcx Bbcy; 
        % hyperparameter for the free stress boundary condition (has to be large)
        alpha_BC=10^3; 
        B=B+alpha_BC^2*sparse(Bbc'*Bbc);
        clear Bbc;
    end
     
    %% Traction force
    % loop on the time step
    for k0=1:kM
    
        ForceName = ['/Users/' DirectoryName '/' Session '/Traction_PIV' num2str(k0) '.txt'];
        e = dlmread(ForceName);
        vTx = e(:,3);
        vTy = e(:,4);
        clear e;

        T=[vTx; vTy]; % vector T (2N components)
        
        mTx=reshape(vTx, R, C);
        mTy=reshape(vTy, R, C);   
    
    %% Hyperparameters
    
    % Lambda
        if meth_Lambda==1
            stepM=10; 
            step=1;
            Lambda_MAP(step)=10^(-3);
     
            % iteration loop
            while step<stepM
                sigmaMAP=(Lambda_MAP(step)*B+l^2*A'*A)\(l^2*A'*T);
                sigmaMAP_norm(step)=sigmaMAP'*B*sigmaMAP;
                res_norm(step)=(T-A*sigmaMAP)'*(T-A*sigmaMAP);
                step=step+1;
                s0_2=sigmaMAP_norm(step-1)/(4*N+2*(C+R)+2);
                s_2=res_norm(step-1)/(2*N+2);
                Lambda_MAP(step)=l^2*s_2/s0_2;
            end
     
            % hyperparameter value   
            Lambda=Lambda_MAP(stepM);
            
            % MAP estimate of the noise amplitude
            noise_value_MAP=sqrt(s_2); 
    
        elseif meth_Lambda==2
        % (parametric) L-curve
            k=1;
            for lc=-10:0.1:0
                 Lambda_Lcurve(k)=10^lc;
                 Lsigma_post=(Lambda_Lcurve(k)*B+l^2*A'*A)\(l^2*A'*T);
                 sigma_norm(k)=Lsigma_post'*B*Lsigma_post;
                 res_norm(k)=(T-A*Lsigma_post)'*(T-A*Lsigma_post);
                 k=k+1;
            end
     
            % local slope of the L-curve in logarithmic scale
            dL=(log(sigma_norm(3:end))-log(sigma_norm(1:end-2)))./(log(res_norm(3:end))-log(res_norm(1:end-2)));
            % local curvature of the L-curve

            % hyperparameter
            % Lambda at the point of maximal curvature
            [max_curv,index_max]=max(ddL);
            Lambda=Lambda_Lcurve(index_max+2); 
    
            % L-curve estimate of the noise amplitude
            noise_value_Lcurve=sqrt(res_norm(index_max+2)/(2*N+2)); 
            
        elseif meth_Lambda==3
        % set the value of Lambda
            Lambda=10^(-6);
        end
    
    %% Stress and covariance matrix estimations
        sigma_post_obs{k0}=(Lambda*B+l^2*A'*A)\(l^2*A'*T);
    
        if noise_value ~= 0
            Spost{k0}=noise_value^2*l^2*((Lambda*B+l^2*A'*A)\speye(Ninf));
        end
        
    end
    clear B;
     
    %% Data interpolation and stress representation
    for k0=1:kM
    
        if noise_value ~= 0
            % Inferred variances of xx, yy, xy and yx stress components
            Spost_xx=reshape(diag(Spost{k0}(1:N+R, 1:N+R)), C+1, R)'; 
            Spost_yy=reshape(diag(Spost{k0}(N+R+1:2*N+C+R, N+R+1:2*N+C+R)), C, R+1)'; 
            Spost_xy=reshape(diag(Spost{k0}(2*N+C+R+1:3*N+2*C+R, 2*N+C+R+1:3*N+2*C+R)), C, R+1)'; 
            Spost_yx=reshape(diag(Spost{k0}(3*N+2*C+R+1:4*N+2*(C+R), 3*N+2*C+R+1:4*N+2*(C+R))), C+1, R)'; 
        end        
        vsigma_post_xx=sigma_post_obs{k0}(1:N+R); % inferred stress xx in N+R-vector form
        vsigma_post_yy=sigma_post_obs{k0}(N+R+1:2*N+C+R); % inferred stress yy in N+C-vector form
        vsigma_post_xy=sigma_post_obs{k0}(2*N+C+R+1:3*N+2*C+R); % inferred stress xy in N+C-vector form
        vsigma_post_yx=sigma_post_obs{k0}(3*N+2*C+R+1:4*N+2*(C+R)); % inferred stress yx in N+R-vector form
     
    % Interpolation on the grid of the force data
        for i=1:R
            for j=1:C       
                vsigma_xx((i-1)*C+j, 1)=(vsigma_post_xx((i-1)*(C+1)+j)+vsigma_post_xx((i-1)*(C+1)+j+1))/2;
                vsigma_yx((i-1)*C+j, 1)=(vsigma_post_yx((i-1)*(C+1)+j)+vsigma_post_yx((i-1)*(C+1)+j+1))/2;
                vsigma_yy((i-1)*C+j, 1)=(vsigma_post_yy((i-1)*C+j)+vsigma_post_yy((i-1)*C+j+C))/2;
                vsigma_xy((i-1)*C+j, 1)=(vsigma_post_xy((i-1)*C+j)+vsigma_post_xy((i-1)*C+j+C))/2; 
                if noise_value ~= 0
                    Ssigma_xx(i,j)=(Spost_xx(i, j)+Spost_xx(i, j+1))/4;
                    Ssigma_yx(i,j)=(Spost_yx(i, j)+Spost_yx(i, j+1))/4;
                    Ssigma_yy(i,j)=(Spost_yy(i, j)+Spost_yy(i+1, j))/4;
                    Ssigma_xy(i,j)=(Spost_xy(i, j)+Spost_xy(i+1, j))/4;
                end
            end
        end
        vsigma_post_xx=vsigma_xx; % inferred stress xx in N-vector form
        vsigma_post_yy=vsigma_yy; % inferred stress yy in N-vector form
        vsigma_post_xy=(vsigma_yx+vsigma_xy)/2; % inferred stress xy in N-vector form
        % Reshape the stress vector in a RxC-matrix form
        sigma_post_xx=reshape(vsigma_post_xx, C, R)'; 
        sigma_post_yy=reshape(vsigma_post_yy, C, R)'; 
        sigma_post_xy=reshape(vsigma_post_xy, C, R)'; 
    
        if noise_value ~= 0 % estimates of error bars
            Spost_xx=Ssigma_xx;
            Spost_yy=Ssigma_yy;
            Spost_xy=(Ssigma_xy+Ssigma_yx)/4;
            clear Ssigma_xx Ssigma_yy Ssigma_xy Ssigma_yx;
            spostxx=Spost_xx.^0.5; 
            spostyy=Spost_yy.^0.5;
            spostxy=Spost_xy.^0.5;
        end
     
        
        if meth_avg == 0
            sigma_mean_xx = zeros(R,C);
            sigma_mean_xy = zeros(R,C);
            sigma_mean_yy = zeros(R,C);
            sigma_mean_xx = sigma_post_xx;
            sigma_mean_xy = sigma_post_xy;
            sigma_mean_yy = sigma_post_yy;
            excluded_index = (xg < xmin+l*edgesize | xg > xmax-l*edgesize) | (yg < xmin+l*edgesize | yg > ymax-l*edgesize);
            sigma_mean_xx(excluded_index) = nan;
            sigma_mean_xy(excluded_index) = nan;
            sigma_mean_yy(excluded_index) = nan;
        elseif meth_avg == 1
            sigma_mean_xx = nan(R,C);
            sigma_mean_xy = nan(R,C);
            sigma_mean_yy = nan(R,C);
            included_index = (xg >= xmin+l*edgesize & xg < xmin+2*l*edgesize) |...
                (xg > xmax-2.*l*edgesize & xg <= xmax-l*edgesize) | ...
                (yg >= ymin+l*edgesize & yg < ymin+2*l*edgesize) | ...
                (yg > ymax-2*l*edgesize & yg <= ymax-l*edgesize) ;
            sigma_mean_xx(included_index) = sigma_post_xx(included_index);
            sigma_mean_xy(included_index) = sigma_post_xy(included_index);
            sigma_mean_yy(included_index) = sigma_post_yy(included_index);
            excluded_index = (xg < xmin+l*edgesize | xg > xmax-l*edgesize) | (yg < xmin+l*edgesize | yg > ymax-l*edgesize);
            sigma_mean_xx(excluded_index) = nan;
            sigma_mean_xy(excluded_index) = nan;
            sigma_mean_yy(excluded_index) = nan;
        end
        
        smean(k0) = mean(sigma_mean_xy, "all", "omitnan");
        dmean(k0) = mean((sigma_mean_xx - sigma_post_yy)./2., "all", "omitnan");
        tmean(k0) = mean((sigma_mean_xx + sigma_post_yy)./2., "all", "omitnan");
            
        sstd(k0) = std(sigma_post_xy, 0, "all", "omitnan");
        dstd(k0) = std((sigma_post_xx - sigma_post_yy)./2., 0, "all", "omitnan");
        tstd(k0) = std((sigma_post_xx + sigma_post_yy)./2., 0, "all", "omitnan");
     
    %% Accuracy of the inference
        Tinf=A*sigma_post_obs{k0}; 
        % divergence of the inferred stress Tinf =A*sigma_inf
        clear sigma_post_obs{k0};
        Tinfx=Tinf(1:N); 
        Tinfy=Tinf(N+1:2*N);
        clear Tinf;
        R2x=1-((vTx'-Tinfx')*(vTx-Tinfx))/((vTx'-mean(vTx))*(vTx-mean(vTx)));
        R2y=1-((vTy'-Tinfy')*(vTy-Tinfy))/((vTy'-mean(vTy))*(vTy-mean(vTy)));
        % R2: coefficient of determination between the data and the traction forces 
        % obtained by calculating the divergence of the inferred stress 
        % (should be close to 1)
        R2=(R2x+R2y)/2; 
    
    % Inferred stress mean values
        sxxm=sum(vsigma_post_xx)/N;
        syym=sum(vsigma_post_yy)/N;
        sxym=sum(vsigma_post_xy)/N;
     
    % Inferred hyperparameter
        stress.Lambda{k0}=Lambda;
        
    % Inferred stress values
        stress.sxx{k0}=sigma_post_xx;
        stress.syy{k0}=sigma_post_yy;
        stress.sxy{k0}=sigma_post_xy;
        stress.x{k0}=x;
        stress.y{k0}=y;
        stress.R2_T{k0}=R2;
        stress.mean_from_sigma{k0}=[sxxm,syym,sxym];   
    
    % Stress mean values from traction force
        if BC==1
            sxx_mean=-sum(sum(reshape(vTx, C, R)'.*(ones(R, 1)*(x-(xmax-xmin)/2))))/N;
            syy_mean=-sum(sum(reshape(vTy, C, R)'.*((y-(ymax-ymin)/2)'*ones(1, C))))/N;
            sxy_mean=-sum(sum(reshape(vTy, C, R)'.*(ones(R, 1)*(x-(xmax-xmin)/2))+reshape(vTx, C, R)'.*((y-(ymax-ymin)/2)'*ones(1,C))))/(2*N);
            stress.mean_from_t{k0}=[sxx_mean,syy_mean,sxy_mean];
        end
        
    % Values of error bars on the inferred stress
        if noise_value~=0
            stress.error_sxx{k0}=spostxx;
            stress.error_syy{k0}=spostyy;
            stress.error_sxy{k0}=spostxy;
        end
        
        %% Segmentation
        % Load spheroid image (jpg)
        fileName = sprintf('%s dic%04d.jpg', Session, k0);
        filePath = fullfile(dirPath, fileName);
        I = imread(filePath);
        % [I, m] = imread(filePath);

        % Resize spheroid image
        heatmap_height = round(length(y));
        heatmap_width = round(length(x));
        resized_image = imresize(I, [heatmap_width heatmap_height]);
       
        % Segmentation for OG size =========================
        % % Adjust data to span data range.
        I = imadjust(I);
        BW = imbinarize(im2gray(I), 'adaptive', 'Sensitivity', 0.0, 'ForegroundPolarity', 'bright'); % 0.6
         
        % Erode mask with default
        radius = 3; % change here
        decomposition = 0;
        se = strel('disk', radius, decomposition);
        BW = imerode(BW, se);
        
        % Dilate mask with disk
        radius = 34; %change here
        decomposition = 0;
        se = strel('disk', radius, decomposition);
        BW = imdilate(BW, se);
        
        % Fill holes
        BW = imfill(BW, 'holes');
        
        % Create masked image.
        maskedImage = I;
        maskedImage(~BW) = 0;
        
        %% Retaining largest component ============================
        % Find connected components
        CC = bwconncomp(BW);
        % Measure properties of image regions
        props = regionprops(CC, 'Area');
        % Find the largest connected component based on area
        [~, idx] = max([props.Area]);
        % Create a binary image with only the largest component
        largestComponent = false(size(BW));
        largestComponent(CC.PixelIdxList{idx}) = true;

        %% Dilation/Erosion Area Selection =========================

        pixel_per_micron = 0;
        if coeff == 0.325
            pixel_per_micron = 3.078156312;
        end
        if coeff == 0.216
            pixel_per_micron = 4.615;
        end

        % A. 10 microns around the spheroid 
        se_a = strel('disk', round(pixel_per_micron * 10));
        dilatedComponent_a = imdilate(largestComponent, se_a);
        dilation_a = dilatedComponent_a & ~largestComponent;
        dilation_a_resized = imresize(dilation_a, [heatmap_width heatmap_height]);

        % B. 10-20 microns border around spheroid
        se_b = strel('disk', round(pixel_per_micron * 20));
        dilatedComponent_b = imdilate(largestComponent, se_b);
        dilation_b = (dilatedComponent_b & ~largestComponent) & ~dilation_a;
        dilation_b_resized = imresize(dilation_b, [heatmap_width heatmap_height]);

        % C. 10 microns within spheroid
        erodedComponent_a = imerode(largestComponent, se_a);
        erosion_a = largestComponent & ~erodedComponent_a;
        erosion_a_resized = imresize(erosion_a, [heatmap_width heatmap_height]);

        %% Plotting 
        overlay_1 = imoverlay(I, dilation_a, 'red');
        overlay_2 = imoverlay(I, dilation_b, 'blue');
        overlay_3 = imoverlay(I, erosion_a, 'yellow');
        
        figure(k0);
        subplot(1,4,1); imshow(I); title('Original Grayscale Image');
        subplot(1,4,2); imshow(overlay_1); title('A. + 10 microns');
        subplot(1,4,3); imshow(overlay_2); title('B. + 10-20 microns');
        subplot(1,4,4); imshow(overlay_3); title('C. - 10 microns');
        saveas(k0, ['Segmentation' DirectoryName '_' Session '_' num2str((k0), '%d')  '_Lambda_' num2str((Lambda), '%.1e') '.jpg']  );
        close(k0);

        %% Multiplication of binary with heatmap
        heatmap = (convertkPa * (sigma_post_xx + sigma_post_yy) / 2.); %245.5920
      
        heatmap_a = (heatmap .* dilation_a_resized);
        heatmap_b = heatmap .* dilation_b_resized;
        heatmap_c = heatmap .* erosion_a_resized;

        selected_values_a = heatmap_a(dilation_a_resized == 1);
        selected_values_b = heatmap_b(dilation_b_resized == 1);
        selected_values_c = heatmap_c(erosion_a_resized == 1);
        

        timeframe(k0) = k0;
        a_average(k0) = mean(selected_values_a);
        b_average(k0) = mean(selected_values_b);
        c_average(k0) = mean(selected_values_c);
        
        % figure(5)
        % imagesc(x, y, heatmap_a);
        % colorbar;  % Add a colorbar to show the scale
        % axis tight;  % Adjust the axis limits to fit the data
        % set(gca, 'YDir', 'normal');  % Set y-axis to normal direction
        % xlabel('X-axis');
        % ylabel('Y-axis');
        % title('Heatmap A');
        % colormap(redblue);
        % caxis([-4 4]);


        %% Overlay heatmap and image 
        % figure(1000 + k0)
        % set(gcf,'Visible','off')
        % % Display spheroid image on axis1
        % ax1 = axes;
        % imshow(resized_image, 'Parent', ax1);
        % hold on;
        % 
        % % Display the heatmap overlay on axis2 
        % ax2 = axes;
        % imagesc(ax2, x, y, (convertkPa * (sigma_post_xx + sigma_post_yy) / 2.))
        % alpha(ax2, 0.6)
        % 
        % colormap(ax2, redblue);
        % colorbar(ax2);
        % set(ax2, 'YDir', 'normal');
        % set(ax2, 'FontSize', 17, 'FontName', 'Times');
        % xlabel(ax2, 'x (\mum)', 'FontSize', 17);
        % ylabel(ax2, 'y (\mum)', 'FontSize', 17);
        % title(ax2, 'Tension (kPa \mum)', 'FontSize', 17);
        % ax.Color = 'none';
        % set(ax2, 'Color', 'none');
        % caxis(ax2, [-4 4]);
        % axis(ax2, [xmin - 0.01 * (xmax - xmin), xmax + 0.01 * (xmax - xmin), ymin - 0.01 * (ymax - ymin), ymax + 0.01 * (ymax - ymin)]);
        % 
        % 
        % % Adjust position of ax2 to overlay on ax1
        % ax2.Position = ax1.Position;
        % 
        % % Link axes for zooming purposes, if needed
        % linkaxes([ax1, ax2], 'xy');
        % 
        % saveas(1000+k0, ['Tension_map_' DirectoryName '_' Session '_' num2str((k0), '%d')  '_Lambda_' num2str((Lambda), '%.1e') '.jpg']  );
        % close(1000+k0)

                        
    end
     
    T = table(timeframe', a_average', b_average', c_average', 'VariableNames', {'timeframe', 'a_average', 'b_average', 'c_average'});
    excelFileName = [DirectoryName, '_', Session, '.xlsx'];
    writetable(T, excelFileName, 'Sheet', 'Sheet1', 'Range', 'A1');

end

