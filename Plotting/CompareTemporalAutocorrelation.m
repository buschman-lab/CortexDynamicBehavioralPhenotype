data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\ALL';
fp = fig_params_vpa; 
savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\Figures\SupplementalFig3';
if ~exist(savedir)
    mkdir(savedir);
end
%% now load all the individual motifs 
motif_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution';
gp = general_params_vpa;
file_list = GrabFiles('\w*chunk\w*.mat', 0, {motif_dir});
[~, fn] = cellfun(@(x) fileparts(x), file_list, 'UniformOutput',0); 
mouse = cellfun(@(x) str2double(x(6:9)), fn,'UniformOutput',0);
mouse = [mouse{:}];
group = isVPA(mouse);

%load all the motifs
temp = cellfun(@(x) load(x,'H'),file_list,'UniformOutput',0);
H = cellfun(@(x) x.H, temp,'UniformOutput',0);
motif_id = arrayfun(@(n) repmat(mouse(n),size(H{n},1),1), 1:numel(H),'UniformOutput',0);%get the mouse id for each motif
motif_id = cat(1,motif_id{:});
H = cat(1,H{:});

%%
%Get the halflife of autocorrelation of Hs
dur = 13*2;
auto_rho = zeros(size(H,1),dur*2+1);
for i = 1:size(H,1)
   auto_rho(i,:)  = xcorr(H(i,:)-nanmean(H(i,:)),(dur),'normalized');     
end

%average per animal
auto_rho_avg = zeros(20,dur*2+1);
mouse_list = unique(mouse);
for i = 1:numel(mouse_list)
   auto_rho_avg(i,:) = nanmean(auto_rho(motif_id==mouse_list(i),:)); 
end

%compare between groups
temp_sal = auto_rho_avg(isVPA(mouse_list)==0,:);
temp_vpa = auto_rho_avg(isVPA(mouse_list)==1,:);
x = linspace((-1*dur)/13,(dur/13),(dur*2+1));

fp = fig_params_vpa;
figure; hold on;
shadedErrorBar(x,nanmean(temp_sal),std(temp_sal),'lineprops',{'color',fp.c_sal,'markerfacecolor','k'},'transparent',1)
shadedErrorBar(x,nanmean(temp_vpa),std(temp_vpa),'lineprops',{'color',fp.c_vpa,'markerfacecolor','k'},'transparent',1)
fp.FormatAxes(gca)
ylabel('Autocorrelation');
xlabel('time (s)');
title({'Group-Wise','Autocorrelation (H weights)'},'fontweight','normal')
set(gca,'units','centimeters','position',[3 3 4 4])

%% save off
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','TemporalAutocorrelation',savedir,1); close all




