function D = MRF_dictionary_umcu(idx, R, x)
% Calculates the dictionary structure required for Jakob his method, I
% modified it to work for FISP sequences.
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Jakob Asslaender, August 2016
% New York University School of Medicine, Center for Biomedical Imaging
% University Medical Center Freiburg, Medical Physics
% jakob.asslaender@nyumc.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get dictionary in adequate dimensions

x=x';
[u,s,~]=svd(x,'econ');
u = u(:,1:R);
x = u'*x;
D.u    = u;
D.s    = diag(s);

sos_dict = l2_norm(x, 1);
x = (x./repmat(sos_dict, [size(x,1) 1]));

D.magnetization = x;
D.normalization = sos_dict;
D.lookup_table  = idx;

D.parameter{1}  = 'T1';
D.parameter{2}  = 'T2';
D.parameter{3}  = 'B1';

disp(['Low-rank conserved ',num2str(sum(D.s(1:R)) / sum(D.s)),' of the energy.'])

end