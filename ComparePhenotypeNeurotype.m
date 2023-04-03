function ComparePhenotypeNeurotype(data_dir,varargin)
%Camden MacDowell - timeless
opts.file_string = 'chunk';
opts.behavior_col = 2:13; 
opts.neuro_col
% opts.nodynamics = 1; 

opts = ParseOptionalInputs(opts,varargin);

fp = fig_params_vpa;
%load mouse info
file_list = GrabFiles(['\w*',opts.file_string,'\w*'],0,{data_dir});
% file_list = substring_filename(file_list,'nodynamics',opts.nodynamics);

warning('off'); %excel loading creates silly warnings
[mouse, ~] = unique(MouseNumFromFileName(file_list));
grp = isVPA(mouse);

%get the neurotype and 
[pt,pt_null,~,~,~,pt_imaged,~] = Phenotype(data_dir,'behavior_col',opts.behavior_col,'file_string',opts.file_string);
[nt,nt_null] = Neurotype(data_dir,'file_string',opts.file_string,'neuro_col',opts.neuro_col); warning('on');

%Compare Neurotype and Phenotype axes
figure; hold on; 
lm = fitlm(pt,nt); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %poitns
p(2).Color = [0 0 0]; %fit
p(3).Color = [0 0 0]; %bounds lower
p(4).Color = [0 0 0]; %bound upper
plot(pt(grp==0),nt(grp==0),'.','markersize',15,'color',fp.c_sal)
plot(pt(grp==1),nt(grp==1),'.','markersize',15,'color',fp.c_vpa)
legend off
set(gca,'xtick',get(gca,'xlim'),'ytick',get(gca,'ylim'),...
    'xticklabel',{'SAL-like','VPA-like'},'yticklabel',{'SAL-like','VPA-like'});
ylabel('Neurotype')
xlabel('Phenotype')
[rho,p] = corr(pt(:),nt(:)); %correlation
fp.SetTitle(gca,{'Correlation Between Animal Phenotypes';sprintf('and Neurotypes rho=%0.2g p=%0.2g',rho,p)});
fp.FormatAxes(gca)

%Compare Neurotype and Phenotype axes Build Using Just the Imaged Animals
figure; hold on; 
lm = fitlm(pt_imaged,nt); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %poitns
p(2).Color = [0 0 0]; %fit
p(3).Color = [0 0 0]; %bounds lower
p(4).Color = [0 0 0]; %bound upper
plot(pt_imaged(grp==0),nt(grp==0),'.','markersize',15,'color',fp.c_sal)
plot(pt_imaged(grp==1),nt(grp==1),'.','markersize',15,'color',fp.c_vpa)
legend off
set(gca,'xtick',get(gca,'xlim'),'ytick',get(gca,'ylim'),...
    'xticklabel',{'SAL-like','VPA-like'},'yticklabel',{'SAL-like','VPA-like'});
ylabel('Neurotype')
xlabel('Phenotype')
[rho,p] = corr(pt_imaged(:),nt(:)); %correlation
fp.SetTitle(gca,{'Imaged Animal Phenotypes';sprintf('and Neurotypes rho=%0.2g p=%0.2g',rho,p)});
fp.FormatAxes(gca)


%Compare with the null neurotype
figure; hold on; 
lm = fitlm(pt,nt_null); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %poitns
p(2).Color = [0 0 0]; %fit
p(3).Color = [0 0 0]; %bounds lower
p(4).Color = [0 0 0]; %bound upper
plot(pt(grp==0),nt_null(grp==0),'.','markersize',15,'color',fp.c_sal)
plot(pt(grp==1),nt_null(grp==1),'.','markersize',15,'color',fp.c_vpa)
legend off
set(gca,'xtick',get(gca,'xlim'),'ytick',get(gca,'ylim'),...
    'xticklabel',{'SAL-like','VPA-like'},'yticklabel',{'SAL-like','VPA-like'});
ylabel('Null Neurotype')
xlabel('Phenotype')
[rho,p] = corr(pt(:),nt_null(:)); %correlation
fp.SetTitle(gca,{'Correlation Between Animal Phenotypes and';sprintf('NULL Neurotypes rho=%0.2g p=%0.2g',rho,p)});
fp.FormatAxes(gca)

%Compare with the null phenotype
figure; hold on; 
lm = fitlm(pt_null,nt); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %poitns
p(2).Color = [0 0 0]; %fit
p(3).Color = [0 0 0]; %bounds lower
p(4).Color = [0 0 0]; %bound upper
plot(pt_null(grp==0),nt(grp==0),'.','markersize',15,'color',fp.c_sal)
plot(pt_null(grp==1),nt(grp==1),'.','markersize',15,'color',fp.c_vpa)
legend off
set(gca,'xtick',get(gca,'xlim'),'ytick',get(gca,'ylim'),...
    'xticklabel',{'SAL-like','VPA-like'},'yticklabel',{'SAL-like','VPA-like'});
ylabel('Neurotype')
xlabel('Null Phenotype')
[rho,p] = corr(pt_null(:),nt(:)); %correlation
fp.SetTitle(gca,{'Correlation Between Animal NULL Phenotypes';sprintf('and Neurotypes rho=%0.2g p=%0.2g',rho,p)});
fp.FormatAxes(gca)


end %function end
