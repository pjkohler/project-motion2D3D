function ellipsePlot(rc_struct,sub_idx,rc_num,cond_colors,cond_labels,gca_opts,e_max,f_1,plot_ave)
    l_width = 1.5;
    if nargin < 9
        plot_ave = true;
    end
    if nargin < 8
        f_1 = 2; 
    else
    end
    if nargin < 7
        e_max = [];
    else
    end
    if nargin < 6
        gca_opts = {'tickdir','out','box','off','fontsize',12,'fontname','Helvetica','linewidth',l_width,'Color','none'};
    else
    end
    if nargin < 5
        cond_labels = {'h-in','h-anti','v-in','v-anti'};
    else
    end
    if nargin < 4
        c_brewer = load('colorBrewer.mat');
        cond_colors = c_brewer.rgb20(1:2:end,:);
        cond_colors([2,3],:) = cond_colors([3,2],:);
    else
    end
    
    n_conds = size(rc_struct.data,1);
    n_subs = size(rc_struct.data,2);
    
    [raw_real,raw_imag] = getRealImag(rc_struct.data);
    [comp_real,comp_imag] = getRealImag(rc_struct.comparisonData);
    raw_real = cellfun(@(x,y) cat(2,x,y), raw_real,comp_real,'uni',false);
    raw_imag = cellfun(@(x,y) cat(2,x,y), raw_imag,comp_imag,'uni',false);
    
    if nargin < 3
        rc_num = [1,6];
    else
    end
    
    if plot_ave
        n_bins = size(raw_real{1},1)+1;
    else
        n_bins = size(raw_real{1},1);
    end
    
    if nargin < 2
        sub_idx = true(1,n_subs);
    else
    end
    
    % loop over rc
    rc_count = size(raw_real{1},2); % rc and comparison
    for r = 1:rc_count
        rca_real = cellfun(@(x) squeeze(nanmean(x(:,r,:),3)),raw_real,'uni',false);
        rca_real = cell2mat(permute(rca_real,[3,2,1]));
        rca_imag = cellfun(@(x) squeeze(nanmean(x(:,r,:),3)),raw_imag,'uni',false);
        rca_imag = cell2mat(permute(rca_imag,[3,2,1]));
        rca_real(end+1,:,:) = squeeze(nanmean(rca_real));
        rca_imag(end+1,:,:) = squeeze(nanmean(rca_imag));
        rca_real = permute(rca_real,[2,3,1]);
        rca_imag = permute(rca_imag,[2,3,1]);
        
        for b = 1:n_bins
            for c = 1:n_conds
                xy_vals = cat(2,rca_real(sub_idx,c,b),rca_imag(sub_idx,c,b));
                nan_vals = sum(isnan(xy_vals),2)>0;
                [ampDiff,phaseDiff,~,error_ellipse] = fitErrorEllipse(xy_vals(~nan_vals,:),'95CI');
                ellipse_data(b,c,r).ellipse = error_ellipse;
                ellipse_data(b,c,r).means = nanmean(xy_vals);
                ellipse_data(b,c,r).phaseDiff = phaseDiff;
                ellipse_data(b,c,r).ampDiff = ampDiff;
                if b == 8 && c == 3
                    [ampDiff,phaseDiff,~,error_ellipse] = fitErrorEllipse(xy_vals(~nan_vals,:),'95CI',true);
                else
                end
            end
        end
        clear rca_real; clear rca_imag;
    end
    
    close all;
    for r = 1:length(rc_num) 
        ellipse_fig = figure;
        set(ellipse_fig,'units','centimeters');
        fig_pos = get(gcf,'pos');
        fig_pos(3) = 35;
        fig_pos(4) = 20;
        set(ellipse_fig,'pos',fig_pos);
        
        bin_vals = cell2mat(cellfun(@(x) str2num(x),rc_struct.settings.binLabels,'uni',false));
        bin_text = arrayfun(@(x) num2str(x,'%.1f'),bin_vals,'uni',false);
        log_step = diff(reallog(bin_vals(1:2))); % step size
        x_min = reallog(bin_vals(1))-log_step*.5;
        if plot_ave
            bin_text = [bin_text; 'ave'];
            x_max = reallog(bin_vals(end))+log_step*2.5; % add 2.5 steps
            extra_bins = arrayfun(@(x) reallog(bin_vals(end))+x, [log_step,log_step*2]);
        else
            x_max = reallog(bin_vals(end))+log_step*.5; % add 2.5 steps
            extra_bins = [];
        end
        
        
        for b = 1:n_bins
            if b == 1
                e_max = floor(max(max(abs(cat(1,ellipse_data(:,:,rc_num(r)).means))))*.75);
            else
            end
            if e_max > 1
                e_units = 1;
            else
                e_units = .5;
            end
            s_h(b) = subplot(2,n_bins,n_bins+b);
            hold on;
            grid_h(1) = plot(zeros(1,2),[-e_max,e_max*.75],'k-','linewidth',l_width);
            grid_h(2) = plot([-e_max,e_max],zeros(1,2),'k-','linewidth',l_width);
            for c = 1:n_conds
                cur_means = ellipse_data(b,c,rc_num(r)).means;
                cur_ellipse = ellipse_data(b,c,rc_num(r)).ellipse;
                
                temp_data = [cur_ellipse;cur_means];
                temp_phase = angle(complex(temp_data(:,1),temp_data(:,2)));
                temp_amp = abs(complex(temp_data(:,1),temp_data(:,2)));
                cycle_len = 1000/(f_1*str2num(rc_struct.settings.freqLabels{1}(1)));
                delay_ms = 180;
                delay_rad = 180/cycle_len*2*pi;
                plot_offset = pi/2;
                temp_phase = temp_phase+delay_rad+plot_offset;
                temp_complex = temp_amp.*exp(1i*temp_phase);
                temp_data = [real(temp_complex),imag(temp_complex)];
                cur_ellipse = temp_data(1:end-1,:);
                cur_means = temp_data(end,:);
                
                ellipse_data(b,c,rc_num(r)).means = cur_means;
                ellipse_data(b,c,rc_num(r)).ellipse = cur_ellipse;
                
                p_h(c) = plot([0 cur_means(1)],[0 cur_means(2)],'-','linewidth',l_width,'Color',cond_colors(c,:));
                plot(cur_ellipse(:,1), cur_ellipse(:,2),'-','linewidth',l_width,'Color',cond_colors(c,:));
            end
            text(0,e_max,bin_text(b),'fontsize',12,'fontname','Helvetica','horizontalalignment','center');
            xlim([-e_max,e_max]);
            ylim([-e_max,e_max]);
            axis square;
            if b == 1
                cycle_steps = linspace(0,2*pi-(pi/2),200)+plot_offset;
                [cycle_x,cycle_y] = pol2cart(cycle_steps,e_max*.75);
                pa_h = patch([e_max*.4,e_max*.8,e_max*.8,e_max*.4],[-e_max+.1 -e_max+.1 -.1 -.1],ones(1,3),'edgecolor','none');
                uistack(pa_h,'bottom');
                q_h = quiver(cycle_x(end),cycle_y(end)-e_max*.4,0,e_max*.8,50);
                set(q_h,'AutoScale','off','maxheadsize',20,'color',[200,200,200]./255,'linewidth',l_width)
                uistack(q_h,'bottom');
                arrayfun(@(x) uistack(x,'bottom'),grid_h);
                cycle_h = plot(cycle_x,cycle_y,'-','color',[200,200,200]./255,'linewidth',l_width);
                uistack(cycle_h,'top');
                set(gca,gca_opts{:},'ticklength',[0.03,0.03],'clipping','off','visible','on','xtick',-e_max:e_units:e_max,'ytick',-e_max:e_units:e_max);
                text(-e_max*.75,-e_max*.75,num2str(cycle_len,'%d ms'));
                xlabel('real','fontsize',12,'fontname','Helvetica');
                ylabel('imag','fontsize',12,'fontname','Helvetica');
            else
                set(gca,gca_opts{:},'clipping','off','visible','off');
                if b == 6
                    text(0,e_max*1.2,'displacement (arcmins)','fontsize',12,'fontname','Helvetica','horizontalalignment','center');
                else
                end
                
            end
           hold off
        end
        drawnow;
        for b = 1:n_bins
            sub_pos = get(s_h(b),'position');
            sub_pos([3,4]) = sub_pos([3,4])+sub_pos(3)*1.8;
            if b == 1
                sub_pos(1) = sub_pos(1)-sub_pos(3)*.5;
                old_pos = sub_pos;
                x_plotmin = sub_pos(1);
            else
                sub_pos(1) = old_pos(1)+old_pos(3)/2;
                old_pos = sub_pos;
                if b == n_bins
                    x_plotmax = sub_pos(1)+sub_pos(3);
                else
                end
            end
            set(s_h(b),'position',sub_pos);
        end
        for z = 1:2
            for c = 1:n_conds
                if z == 1
                    % plot amplitudes
                    big_plot(z) = subplot(2,n_bins,1:floor(n_bins/2));
                    big_min = 0; big_max = ceil(max(max(squeeze(rc_struct.stats.Amp(:,rc_num(r),:)))));
                    plot_vals = rc_struct.stats.Amp(1:n_bins,rc_num(r),c)';
                    plot_errs = cat(1,ellipse_data(1:n_bins,c,rc_num(r)).ampDiff)';
                    y_text = 'amplitude (\muV)';
                else
                     % plot phase
                    big_plot(z) = subplot(2,n_bins,(floor(n_bins/2)+1):n_bins);
                    big_min = 100; big_max = 200; 
                    temp_means = cat(1,ellipse_data(:,c,rc_num(r)).means);
                    temp_phases = angle(complex(temp_means(:,1),temp_means(:,2)))-(plot_offset);
                    temp_phases( temp_phases < 0 ) = 2*pi+temp_phases( temp_phases < 0 );
                    plot_vals = (temp_phases/(2*pi)*cycle_len)';
                    ph_error = cat(1,ellipse_data(1:n_bins,c,rc_num(r)).phaseDiff);
                    ph_error = ph_error .* pi / 180;
                    plot_errs = (ph_error ./(2*pi)*cycle_len)';
                    y_text = 'delay (ms)';
                end
                hold on
                if c > 2
                    cur_x = [reallog(bin_vals'),extra_bins(2)];
                else
                    cur_x = [reallog(bin_vals'),extra_bins(1)];
                end
                if plot_ave
                    if c == 1
                        plot(ones(1,2)*(reallog(bin_vals(end))+log_step*.5) ...
                            ,[big_min,big_max],'k-','linewidth',l_width);
                    else
                    end
                    p_h(c) = plot(cur_x(1:end-1),plot_vals(1:end-1),'-o','linewidth',l_width,'Color',cond_colors(c,:),'markersize',10,'markerfacecolor','w','markeredgecolor',cond_colors(c,:));
                    plot(cur_x(end),plot_vals(end),'o','linewidth',l_width,'Color',cond_colors(c,:),'markersize',10,'markerfacecolor','w','markeredgecolor',cond_colors(c,:));
                else
                    plot(cur_x,plot_vals,'-o','linewidth',l_width,'Color',cond_colors(c,:),'markersize',10,'markerfacecolor','w','markeredgecolor',cond_colors(c,:));
                end
                hold on
                h_e(c,:) = ErrorBars(cur_x,plot_vals,plot_errs,'type','bar','color',cond_colors(c,:),'cap',false);
                hold on
            end
            if z == 1
                legend(p_h,cond_labels,'location','northwest','box','off','fontsize',12,'fontname','Helvetica');
            end

            cellfun(@(x) uistack(x,'bottom'),h_e,'uni',false);
            xlim([x_min,x_max]);
            ylim([big_min,big_max]);
            
            xtick_vals = reallog([0.2,0.5,1,2,4,8,16]);
            xtick_labels = arrayfun(@(x) num2str(x),([0.2,0.5,1,2,4,8,16]),'uni',false);
            
            if plot_ave
                xtick_vals = [xtick_vals,xtick_vals(end) + log_step*1.5];
                xtick_labels = [xtick_labels,'ave'];
            else
            end
            set(gca,gca_opts{:},'ticklength',[0.015,0.015],'clipping','off','visible','on','xtick',xtick_vals,'xticklabels',xtick_labels);
            ylabel(y_text,'fontsize',12,'fontname','Helvetica');
            hold off
        end
        drawnow;
        amp_pos = get(big_plot(1),'position');
        plot_span = (x_plotmax - x_plotmin);
        amp_pos(1) = x_plotmin;
        amp_pos(3) = plot_span/2-plot_span/40;
        set(big_plot(1),'position',amp_pos);
        ph_pos = get(big_plot(2),'position');
        ph_pos(3) = plot_span/2-plot_span/40;
        ph_pos(1) = x_plotmax-ph_pos(3);
        set(big_plot(2),'position',ph_pos);
        if rc_num(r) == rc_count
            export_fig(ellipse_fig,'ellipseplot_comp.pdf','-transparent','-painters');
        else
            export_fig(ellipse_fig,sprintf('ellipseplot_rc%d.pdf',rc_num(r)),'-transparent','-painters');
        end
    end
end