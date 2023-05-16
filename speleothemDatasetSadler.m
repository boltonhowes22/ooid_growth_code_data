
function [age_diff, growth_rate] = speleothemDatasetSadler(path_to_data)
    spel = readtable(path_to_data);
    %%
    ID = spel.entity_id;                    % sample id
    depth = spel.depth_dating .* .1;      % convert from mm to microns
    age = spel.corr_age;                    % age in years

    % some weird negative ages, lets just avoid those
    ID = ID(age>0);
    depth = depth(age>0);
    age = age(age>0);

    %% Anytime the depth decreases, we will just make a new ID 
    % some things have multiple transects that is messing this up, so this
    % should solve that 
    % also, the ages are within uncertainty sometimes, so they arent always
    % monotonic. I will just start a new Id at those points and continue.
    

    id = 1;
    uniqueID(1) = 1;
    for x = 2:length(age)
        if depth(x)>depth(x-1) && ID(x) == ID(x-1) && age(x) > age(x-1)
           uniqueID(x) = id ;

        else

            % add one to id, and keep going
            id = id + 1;
            uniqueID(x) = id;
        end

    end





    thickness = [];
    age_diff = [];
    growth_rate = [];
    
    for x = 1:length(uniqueID)

        inds = uniqueID == uniqueID(x);
        depth_ID = depth(inds);
        age_ID = age(inds);

        % remove nans (hiatuses included in the dataset)
        depth_ID = depth_ID(~isnan(age_ID));
        age_ID = age_ID(~isnan(age_ID));


        % difference 
        thickness = [thickness; diff(depth_ID)];
        age_diff = [age_diff; diff(age_ID)];

        % growth rate
        growth_rate = [growth_rate; diff(depth_ID)./diff(age_ID)];
        
        if sum(growth_rate < 0 ) >0 
            keyboard
        end
        
    end


%{
figure(1); clf
hold on; box on; grid on;

sc = scatter(age_diff, growth_rate,15, 'filled', 'k');
sc.MarkerFaceAlpha = .01;

set(gca, 'Yscale', 'Log', 'Xscale', 'Log') 
%}