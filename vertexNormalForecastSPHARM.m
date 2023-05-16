function [rot_inv , grown_rot_inv, lmcosi, grown_lmcosi] = vertexNormalForecastSPHARM(vert_struct, L)
% This function does the following 
% 1. rotates and normalizes each ooid band 
% 2. calculates their thickness
% 3. makes a forecast prediction of what each band should look like
% 4. calculate the spherical harmonic coefficients of both the band, and
% the forecasted band
% 5. calculates the rotationally-invariant coefficients (F(l) = Sum(A(l,m)^2)
%
%
% INPUT 
% vert_struct: a structure where vert_struct{ooid}{bands in ooid}
% L: the L_max of the spherical harmonic expansion
%
% OUTPUT
% rot_inv: rotational invariant coeff of the measured data
% grown_rot_inv: rotational invariant coeff of the forecasted data
% lmcosi: Full SPHARM coefficeints of measured data
% grown_lmcosi: Full SPHARM coefficeints of forecasted data
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % This are the indeces for the difference between bands
    diff_ind = 4:6;
    
    % Rotate the verteces, this is so that normalization works the same for
    % each sample
    for this_ooid = 1:length(vert_struct)
        for this_band = 1:length(vert_struct{this_ooid})

            vert_struct{this_ooid}{this_band} = pca_rotate(vert_struct{this_ooid}{this_band});

        end
    end

    
    % measure axes using PCA
    % measures all ooids in the structure within a loop
    axes = pcMeasureLoop(vert_struct);

    % Extract just the thickness of the growth interval
    for this_ooid = 1:length(vert_struct)
        for this_band = 1:length(vert_struct{this_ooid})-1
            thickness{this_ooid}{this_band} = axes{this_ooid}{this_band+1}(1:3) - axes{this_ooid}{this_band}(1:3);

        end
    end

    
    % Normalize Point Cloud, because SPHARM coefficients are size dependent
    % normalizes all ooids in the structure
    vert_nm = normalizeLoop(vert_struct);

    % Vertex Normal Growth 
    % Surface Normal Growth based on the normals of the vertex
    for this_ooid = 1:length(vert_struct)
        for this_band = 1:length(vert_struct{this_ooid})-1

            % size of growth increment is set by the largest increase in axis length
            increment = max(axes{this_ooid}{this_band+1}(diff_ind))/2; 

            % grow the ooid
            vert_grown{this_ooid}{this_band} = vertexNormalGrowth(vert_struct{this_ooid}{this_band},...
                increment);

        end
    end


    % Normalize grown ooids
    vert_grown_nm = normalizeLoop(vert_grown);
    
    % SPHARM Calculation
    % all ooids in the structure -- both the measured and grown data
    [L2err, lmcosi, ~, ~] = fitAllBands(vert_nm, L);
    [grown_L2err, grown_lmcosi, ~, ~] = fitAllBands(vert_grown_nm, L);
    
    
    % Calculate Rotational Invariance
    rot_inv = rotInvLoop(lmcosi);
    grown_rot_inv = rotInvLoop(grown_lmcosi);
    %keyboard
        
    
    
end


