%% this script calls drdq_corr() to get the WIMP differential rate
%% and saves the result to a text file.
A = 'Si';
M_GeV = 1;

%% good limits for 1 GeV
step_keV = 0.5E-4;
Erecoil_keV = [0:step_keV:500E-3]';

%% good limits for 0.5 GeV
%step_keV = 1E-5;
%Erecoil_keV = [0:step_keV:150E-3]';

%% good limits for 0.1 GeV
%step_keV = 1E-6;
%Erecoil_keV = [0:step_keV:50E-4]';

% sigma is required - although it doesn't matter since we're normalizing
sigma = 1E-41;
rate = drdq_corr(Erecoil_keV, A, M_GeV, sigma);

pdf = rate / trapz(Erecoil_keV, rate);
cdf = cumtrapz(Erecoil_keV, rate) / max(cumtrapz(Erecoil_keV, rate));

% create filename e.g.
% WIMP_xsection_1GeV_1E-3keV_steps.txt
filename = ['WIMP_Si_PDF_CDF_' num2str(M_GeV) 'GeV_' num2str(step_keV) 'keV_steps.txt'];

T = table(Erecoil_keV, pdf, cdf);
writetable(T, filename);

if true
   plot(Erecoil_keV, rate);
   set(gca, 'YScale', 'log');
   set(gca, 'XScale', 'log');
end

