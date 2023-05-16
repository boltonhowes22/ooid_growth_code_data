function [var, lat, lon] = loadGLODAP(path, var_name)
    
    % Load in the GLODAP Data

    f = ncinfo(path);
    lat = ncread(path, 'lat');
    lon = ncread(path, 'lon');
    var = ncread(path,var_name);
end