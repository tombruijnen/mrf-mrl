function kspace_data = radial_phase_correction_model(kspace_data,traj)
%%Do the model based radial phase correction based on the paper of Moussavi
%(2014). Currently only operatos in 2D space and therefore should be
% applied after the 1D FFT for stack-of-stars.
%
% Verbose mode: provide plot and estimate of residual phase at k0

kdim=size(kspace_data);
cp=kdim(1)/2+1;
rad_ang=mod(squeeze(angle(traj(1,1,:,1)+1j*traj(2,1,:,1)))+pi,2*pi);
model_pars=zeros([2 kdim(3:end)]);

% Loop over all partitions and estimate phase error model 
for p=1:prod(kdim(3:end))
    cph=angle(kspace_data(cp,:,p));
    model_pars(:,p)=radial_paramatrizephasemodel(cph',rad_ang);
end

% Loop over all partitions and apply phase correction
phi=@(theta,A)(A(1).*cos(theta)+A(2).*sin(theta)); 
phase_corr_mtx=zeros(size(kspace_data),'single');
for p=1:prod(kdim(3:end))
    cur_phase=phi(rad_ang',model_pars(:,p));
    phase_corr_mtx(:,:,p)=single(exp(-1j*repmat(cur_phase,[kdim(1) 1])));
end

% Apply phase correction by matrix multiplication
cph_pre=mean(matrix_to_vec(var(angle(kspace_data(cp,:,:)),[],2)));
kspace_data=kspace_data.*phase_corr_mtx;
cph_post=mean(matrix_to_vec(var(angle(kspace_data(cp,:,:)),[],2)));

disp(['>> Mean variance of k0 phase changed from: ',num2str(cph_pre),' -> ',num2str(cph_post)])
% END
end