function ncm = noise_covariance_mtx(data)
% data should be samples x channels
ns=size(data,1);
nc=size(data,2);

% preallocate matrix
ncm=zeros(nc,nc);

% Loop over channels
for i=1:nc
    for j=1:nc
        ncm(i,j)=1/(ns)*sum(data(:,i).*conj(data(:,j)));
    end
end

% END
end