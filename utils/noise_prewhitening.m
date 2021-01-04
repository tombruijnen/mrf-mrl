function kspace_data = noise_prewhitening(kspace_data,noise_data,varargin)
%%NOISE PREWHITENING

if nargin < 3
	batches = 1;
else
    batches = varargin{1};
end

noise_covariance=noise_covariance_mtx(squeeze(noise_data));
noise_decorrelation=noise_decorrelation_mtx(noise_covariance);

batch_size = ceil(size(kspace_data,2) / batches);
for b = 1 : batches
	idx = 1 + (b-1) * batch_size : b * batch_size;
	idx(idx > size(kspace_data,2)) = [];
	kspace_data(:,idx,:,:,:,:,:,:,:,:,:,:) = permute(apply_noise_decorrelation_mtx(permute(kspace_data(:,idx,:,:,:,:,:,:,:,:,:,:),...
    		[1:3 5:12 4]),noise_decorrelation),[1:3 12 4:11]);
end

disp('+Receiver coils are whitened.')
% END
end