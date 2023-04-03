function [data, mouse] = CompareMotifLoadings(data_dir,varargin)
%Camden MacDowell - timeless

opts.data_name = 'stats_refit';
opts.file_string = 'chunk';
opts.type = 'box';
opts.verbose = 1; 

opts = ParseOptionalInputs(opts,varargin);

file_list = GrabFiles(['\w*',opts.file_string,'\w*'],0,{data_dir});
[mouse, ~] = MouseNumFromFileName(file_list);
mouse_grp = isVPA(unique(mouse));

%load data
data_all = cellfun(@(x) load(x,opts.data_name),file_list,'UniformOutput',0);
data_all = cellfun(@(x) x.(opts.data_name),data_all,'UniformOutput',0);
        
fp = fig_params_vpa;

%Compare motif loadings
switch opts.type
    case 'cumsum' %when comparing motif discovery
        data = cellfun(@(x) cumsum(sort(x.loadings{1},'descend')), data_all,'UniformOutput',0);
        data = MakeCellsEqual(data,2,1);
        data = cat(1,data{:});
        data = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);

        vpa = cat(1,data{mouse_grp==1});
        vpa(:,isnan(nanmean(vpa)))=[];%remove any that are all NaN        
        vpa_avg = nanmean(vpa);
        vpa_ci = bootci(100,@nanmean,vpa);
        sal = cat(1,data{mouse_grp==0});
        sal(:,isnan(nanmean(sal)))=[];%remove any that are all NaN
        sal_avg = nanmean(sal);
        sal_ci = bootci(100,@nanmean,sal);
        

        if opts.verbose
            figure; hold on; 
            errorbar(1:numel(vpa_avg),vpa_avg,vpa_avg-vpa_ci(1,:),vpa_ci(2,:)-vpa_avg,'linewidth',2,'color',fp.c_vpa);
            errorbar(1:numel(sal_avg),sal_avg,sal_avg-sal_ci(1,:),sal_ci(2,:)-sal_avg,'linewidth',2,'color',fp.c_sal);
            min_length =min(size(vpa,2),size(sal,2));
            p = arrayfun(@(n) ranksum(vpa(:,n),sal(:,n)), 1:min_length,'UniformOutput',1);
            temp = min(cat(1,vpa_ci(:,min_length),sal_ci(:,min_length)));
            for i = 1:numel(temp); if p(i)<0.05; text(i,temp(i),sprintf('%0.2g',p(i))); end; end
            ylim([0 1])
            fp.SetTitle(gca,'Loadings');
            fp.FormatAxes(gca); 
            legend('VPA','SAL','Location','southeast')
        end
 
        
    case 'box' %when comparing same motifs fit to withheld data        
        data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
        data = cat(1,data{:});
        data = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);         
        vpa = cat(1,data{mouse_grp==1})*100;
        sal = cat(1,data{mouse_grp==0})*100;        
        %sort by expvar of saline
        [idx,sal] = SortMotifs(sal,'mean');        
        vpa = vpa(:,idx);
        
        %to return
        data = cellfun(@(x) x(:,idx), data, 'UniformOutput',0);        
        mouse = unique(mouse); 
        
        x = (0.5:1:size(vpa,2));
        
        if opts.verbose
            figure('position',[680   674   556   304]); hold on;     
            boxplot(gca,vpa,'Positions',x,'Widths',fp.b_width,'PlotStyle','compact','Colors',fp.c_vpa);                 
            boxplot(gca,sal,'Positions',x+0.25,'Widths',fp.b_width,'PlotStyle','compact','Colors',fp.c_sal);                
            %add significance
            p = arrayfun(@(n) ranksum(vpa(:,n),sal(:,n)), 1:size(vpa,2),'UniformOutput',1);
%             [~,p] = arrayfun(@(n) ttest2(vpa(:,n),sal(:,n)), 1:size(vpa,2),'UniformOutput',1);
            temp = max(cat(1,vpa,sal))+1;
            for i = 1:numel(temp); if p(i)<1; text(x(i),temp(i),sprintf('%0.2g',p(i)),'Rotation',60,'FontSize',fp.sig_fontsize); end; end

            %format plot
            ylim([0 max(temp)])
            set(gca,'xtick',x+0.125,'XTickLabel',1:size(x,2))
            xlabel('Motif'); 
            ylabel({'Percent Explained';'Variance'})
            fp.SetTitle(gca,'Motif Loadings');
            fp.FormatAxes(gca);         
        end
                
    otherwise
        error('unknown plot type');
end
        


end















