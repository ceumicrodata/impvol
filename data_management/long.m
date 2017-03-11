function long_out = long(wide, parameters)
% This function creates a 2D matrix of size
% (n_countries × n_sectors) × n_years from a 3D matrix of size
% n_countries × n_sectors × n_years. The output matrix will have a
% structure of: within countries all the sectors are enumerated in rows,
% columns contain years.
% 
% e.g.
% a(:,:,1) = [111 121 131 141; 211 221 231 241; 311 321 331 341];
% a(:,:,2) = [112 122 132 142; 212 222 232 242; 312 322 332 342];
% parameters.n_countries = 3;
% parameters.n_sectors = 4;
% parameters.n_years = 2;
% wide(long(a,parameters), parameters)

    n_countries = parameters.n_countries;
    n_sectors = parameters.n_sectors;
    n_years = parameters.n_years;
    
    if size(wide) == [n_countries,n_sectors,n_years]
        long_out = reshape(permute(wide,[2,1,3]), n_countries * n_sectors, n_years);
    else
        disp('The size of the matrix does not meet the requirement');
    end
end

