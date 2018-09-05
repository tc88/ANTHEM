function dv = cv2dv(msh,cv)
% CV2DV brings vector fields from the control volumes
% Vhat to the dual volumes Vdual. There is no averaging involved here,
% only a relocation of quantities is done.
%
% Input:
% msh   struct as defined by src/msh.txt
%       required fields: np,Mx,My,Mz
% cv    vector field defined on the shifted volumes Vhat (3np-by-1)
%
% Output:
% dv    scalar field defined on dual volumes Vdual (np-by-1)
%
% See also pn2pe, ipe2ipn
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

% extract from msh struct
np = msh.np;

% half of each control volume contributes to one adjacent dual volume
cvHalf = 0.5*cv;

% initialize quantities on dual volume
dv = zeros(np,1);
% iterate over primary points
for idx = 1:np
    % if at lower x boundary
    if idx-msh.Mx < 1
        xback = 0;
    else
        xback = cvHalf(idx-msh.Mx);
    end
    % if at lower y boundary
    if idx-msh.My < 1
        yback = 0;
    else
        yback = cvHalf(np+idx-msh.My);
    end
    % if at lower z boundary
    if idx-msh.Mz < 1
        zback = 0;
    else
        zback = cvHalf(2*np+idx-msh.Mz);
    end
    % contributions from the front
    xfront = cvHalf(idx);
    yfront = cvHalf(np+idx);
    zfront = cvHalf(2*np+idx);
    % sum up the contributions
    dv(idx) = xback+xfront+yback+yfront+zback+zfront;
end

end