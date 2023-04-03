function ComparePhenotypeMotifs(data_dir,varargin)

opts.file_string = 'chunk';

opts = ParseOptionalInputs(opts,varargin);

%load mouse info
file_list = GrabFiles(['\w*',opts.file_string,'\w*'],0,{data_dir});
warning('off'); %excel loading creates silly warnings
[mouse, ~] = unique(MouseNumFromFileName(file_list));
grp = isVPA(mouse);

%get the phenotype axis
[pt,pt_null] = Phenotype(data_dir);
[motifs, temp_mouse] = CompareMotifLoadings(data_dir,'verbose',0,'type','box');
%confirm that temp_mouse and mouse are the same (e.g. same order for everything)
assert(sum(mouse == temp_mouse)==numel(mouse),'error there is an issue is the order of the mouse numbers');
motifs = cat(1,motifs{:})*100;

%loop through the motifs
for i = 1:size(motifs,2)
    figure; hold on; 
    lm = fitlm(pt,motifs(:,i)); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
    p(1).Color = [1 1 1]; p(2).Color = [0 0 0]; p(3).Color = [1 1 1]; p(4).Color = [1 1 1]; %data, line, lower bound, upper bound
    plot(pt(grp==0),motifs(grp==0,i),'.','markersize',15,'color',fp.c_sal)
    plot(pt(grp==1),motifs(grp==1,i),'.','markersize',15,'color',fp.c_vpa)
    legend off
    set(gca,'xtick',get(gca,'xlim'),'xticklabel',{'SAL-like','VPA-like'});
    ylabel('Relative PEV')
    xlabel('Phenotype')
    [rho,p] = corr(pt(:),motifs(:,i)); %correlation
    fp.SetTitle(gca,{'Correlation Between Animal Phenotypes and';sprintf('Motif %d rho=%0.2g p=%0.2g',i, rho,p)});
    fp.FormatAxes(gca)
    
%     %null axis
%     figure; hold on; 
%     lm = fitlm(pt_null,motifs(:,i)); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
%     p(1).Color = [1 1 1]; p(2).Color = [0 0 0]; p(3).Color = [1 1 1]; p(4).Color = [1 1 1]; %data, line, lower bound, upper bound
%     plot(pt_null(grp==0),motifs(grp==0,i),'.','markersize',15,'color',fp.c_sal)
%     plot(pt_null(grp==1),motifs(grp==1,i),'.','markersize',15,'color',fp.c_vpa)
%     legend off
%     set(gca,'xtick',get(gca,'xlim'),'xticklabel',{'SAL-like','VPA-like'});
%     ylabel('Relative PEV')
%     xlabel('NULL Phenotype')
%     [rho,p] = corr(pt_null(:),motifs(:,i)); %correlation
%     fp.SetTitle(gca,{'Correlation Between NULL Phenotypes and';sprintf('Motif %d rho=%0.2g p=%0.2g',i, rho,p)});
%     fp.FormatAxes(gca)   
end


end