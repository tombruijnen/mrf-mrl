function [dict,idx] = extendedPhaseGraphs(T1list,T2list,B1list,SP,FA,TR,IS,MS,spoiled,EC)
% Function to do EPG simulations for all T1/T2/B1 and slice profile.

if nargin < 10
    EC = [];
end

% Easy dim handling
N=numel(FA);
Ne=numel(T1list)*numel(T2list)*numel(B1list);

% Check dimensionality of T1/T2 and B1list
if size(T1list,1)==1;T1list=permute(T1list,[2 1]);end
if size(T2list,1)==1;T2list=permute(T2list,[2 1]);end
if size(B1list,1)==1;B1list=permute(B1list,[2 1]);end

% Make the indexer idx [3 n_entries]
idx(1,:)=matrix_to_vec(repmat(T1list,[1 numel(T2list) numel(B1list)]));
idx(2,:)=matrix_to_vec(repmat(permute(T2list,[2 1 3]),[numel(T1list) 1 numel(B1list)]));
idx(3,:)=matrix_to_vec(repmat(permute(B1list,[3 2 1]),[numel(T1list) numel(T2list) 1]));

% Preallocate dictionairy
dict=single(zeros(Ne,N));

% If slice profile use conjugate symmetry
if numel(SP)>1
    hsp=floor(numel(SP)/2)+1;
else
    hsp=1;
end

% Loop over all properties and slice profile
tic
cnt=1; parfor_progress(Ne);% Counter for indexing
for b1=1:numel(B1list)
for t2=1:numel(T2list)
for t1=1:numel(T1list)
    
    % Create matrix for slice profile
    tmp=zeros(N,numel(SP));
    
    % Loop over slice profile
    for sp=1:hsp
        if T1list(t1) > 2 * T2list(t2)
            if isempty(EC)
                tmp(:,sp)=single(EPG_EC(1,B1list(b1)*FA,SP(sp),TR,T1list(t1),T2list(t2),IS,MS,spoiled));
            else
                tmp(:,sp)=single(EPG_EC(1,B1list(b1)*FA,SP(sp),TR,T1list(t1),T2list(t2),IS,MS,spoiled,EC));
            end
        end
    end
    
    % Conjugate symmetry for slice profile
    tmp(:,hsp+1:end)=-1*real(tmp(:,1:hsp-1))+1j*imag(tmp(:,1:hsp-1));
    parfor_progress;
    cnt=cnt+1;
    
    % Allocation
    dict(cnt,:)=sum(tmp,2);
    
end
end
end



% Reset progress
parfor_progress(0);

% Display time
disp(['Simulations took: ',num2str(toc),' s'])

% END
end