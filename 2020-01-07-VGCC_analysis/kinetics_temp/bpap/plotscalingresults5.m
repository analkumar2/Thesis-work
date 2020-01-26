function[]=plotscalingresults5(casrimax,distances,synweight,trunksri)
    % replaces zeros (no calcium measured) in the original data with NaN
    casrimax=removezeros(casrimax);
    
    [nrruns,nrsegments]=size(synweight); % note nrruns is nrrun+1 because the synaptic weights are defined at the end of the run

    sortdistances(:,1)=distances;
    sortdistances(:,2)=1:size(distances,2);
    sorteddistances=sortrows(sortdistances,1);
    sorteddistances(:,3)=1:size(distances,2);
    sorteddistances2=sortrows(sorteddistances,2);
    colorlistjet=colormap(jet(size(distances,2)));
    
  
    %plot casrimax with subthreshold responses
    
    %thresholdca=0.0478391;
    thresholdca=0.052;
    nrskippedtr=6;
    
    figure()
    subplot('Position',[0.07 0.1 0.78 0.8])
    hold on
    for nrsegment=1:nrskippedtr:nrsegments
        plot(casrimax(1:nrruns-1,nrsegment)*1000,'Color',colorlistjet(sorteddistances2(nrsegment,3),:))
    end
    plot([1,nrruns-1],[thresholdca, thresholdca]*1000,'k--')
    xlabel('# run')
    ylabel('peak calcium')

    subplot('Position',[0.85 0.12 0.1 0.76])
    hold on
    for nrsegment=1:nrsegments
        plot(1,nrsegment, 'sq','MarkerFaceColor', colorlistjet(nrsegment,:),'MarkerEdgeColor',colorlistjet(nrsegment,:),'MarkerSize',14) 
    end
    ylim([0,size(distances,2)])
    axis('off')
    text(1.5,1,'0','FontSize',10)
    text(1.5,488,num2str(max(distances),4),'FontSize',10)      
    
    %plot casrimax as dots in stead of lines (incl subthreshold)
        figure()
    subplot('Position',[0.07 0.1 0.78 0.8])
    hold on
    for nrsegment=1:nrskippedtr:nrsegments
        plot(casrimax(1:nrruns-1,nrsegment)*1000,'.','MarkerSize',4,'Color',colorlistjet(sorteddistances2(nrsegment,3),:))
    end
    plot([1,nrruns-1],[thresholdca, thresholdca]*1000,'k--')
    xlabel('# run')
    ylabel('peak calcium')

    subplot('Position',[0.85 0.12 0.1 0.76])
    hold on
    for nrsegment=1:nrsegments
        plot(1,nrsegment, 'sq','MarkerFaceColor', colorlistjet(nrsegment,:),'MarkerEdgeColor',colorlistjet(nrsegment,:),'MarkerSize',14) 
    end
    ylim([0,size(distances,2)])
    axis('off')
    text(1.5,1,'0','FontSize',10)
    text(1.5,488,num2str(max(distances),4),'FontSize',10) 
    
    %plot casrimax without subthreshold responses. 
    figure()
    subplot('Position',[0.07 0.1 0.78 0.8])
    hold on
    
    casrimaxwithoutsub=casrimax;
    for run=1:nrruns-1
        if(min(casrimax(run,:))<0.030)
            casrimaxwithoutsub(run,:)=NaN;
            
        else
            casrimaxwithoutsub(run,:)=casrimax(run,:);
            plot([run,run],[105,108],'k-') % plot a line in case of suprathreshold responses
        end
    end
    
    for nrsegment=1:nrskippedtr:nrsegments
        plot(casrimaxwithoutsub(1:nrruns-1,nrsegment)*1000,'Color',colorlistjet(sorteddistances2(nrsegment,3),:))
    end
    plot([1,nrruns-1],[thresholdca, thresholdca]*1000,'k--')
    xlabel('# run')
    ylabel('peak calcium')
    
    subplot('Position',[0.85 0.12 0.1 0.76])
    hold on
    for nrsegment=1:nrsegments
        plot(1,nrsegment, 'sq','MarkerFaceColor', colorlistjet(nrsegment,:),'MarkerEdgeColor',colorlistjet(nrsegment,:),'MarkerSize',14) 
    end
    ylim([0,size(distances,2)])
    axis('off')
    text(1.5,1,'0','FontSize',10)
    text(1.5,488,num2str(max(distances),4),'FontSize',10)    
    
    % measure the amount of stimulations per segment
    stimnostim=casrimax>0;
    nrstimperseg=sum(stimnostim);
    figure()
    subplot(1,2,1)
    plot(nrstimperseg,'.')
    ylabel('nr stims')
    xlabel('segments nr')
    subplot(1,2,2)
    plot(distances,nrstimperseg,'.')
    xlabel('distances')
    
    % only plot synweight for synapses that are activated at least ones
    for nrsegment=1:nrsegments
        if(nrstimperseg(nrsegment)==0)
            synweight(:,nrsegment)=NaN;
        end
    end
   
    figure()
    subplot('Position',[0.07 0.1 0.78 0.8])
    hold on
    for nrsegment=1:nrskippedtr:nrsegments
        plot(synweight(1:nrruns-1,nrsegment)*10^6,'Color',colorlistjet(sorteddistances2(nrsegment,3),:))
    end
    xlabel('nr stim')
    ylabel('synaptic strength')

    subplot('Position',[0.85 0.12 0.1 0.76])
    hold on
    for nrsegment=1:nrsegments
        plot(1,nrsegment, 'sq','MarkerFaceColor', colorlistjet(nrsegment,:),'MarkerEdgeColor',colorlistjet(nrsegment,:),'MarkerSize',14) 
    end
    ylim([0,size(distances,2)])
    axis('off')
    text(1.5,1,'0','FontSize',10)
    text(1.5,488,num2str(max(distances),4),'FontSize',10)
    
    figure()
    hold on
    for nrsegment=1:nrsegments
        if(max(nrsegment==trunksri+1)==1)
            % trunk segment
            %plot(distances(nrsegment),synweight(1,nrsegment)*10^6,'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',8)
            plot(distances(nrsegment),synweight(200,nrsegment)*10^6,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8) 
        else
            % non trunk segment
            %plot(distances(nrsegment),synweight(1,nrsegment)*10^6,'o','MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor','w','MarkerSize',8)
            plot(distances(nrsegment),synweight(200,nrsegment)*10^6,'o','MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)            
        end
    end
    
    plot([0,max(distances)],[200,200],'-k')
    xlim([0 1000])
    %ylim([0 nrruns-1])
    xlabel('distance')
    ylabel('synaptic weight')
    title('synaptic weight at start (open) and end of scaling (closed)')
    
    figure()
    plot(synweight(1:200,:),casrimax(1:200,:),'.')
    ylabel('peak calcium')
    xlabel('synaptic strength')
    
    
end