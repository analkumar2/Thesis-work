function[]=plotscalingresults_pergroup1(casrimax,distances,synweight,trunksri)
    % replaces zeros (no calcium measured) in the original data with NaN
    casrimax=removezeros(casrimax);
    
    [nrruns,nrsegments]=size(synweight); % note nrruns is nrrun+1 because the synaptic weights are defined at the end of the run

    sortdistances(:,1)=distances;
    sortdistances(:,2)=1:size(distances,2);
    sorteddistances=sortrows(sortdistances,1);
    sorteddistances(:,3)=1:size(distances,2);
    sorteddistances2=sortrows(sorteddistances,2);
    colorlistjet=colormap(jet(size(distances,2)));
    
    thresholdca=0.0470;
    %thresholdca=0.052;
    nrskippedtr=40;
    
    total200min=0;
    total400min=0;
    total600min=0;
    total800min=0;
    total800plus=0;
    
    index200min=[];
    index400min=[];
    index600min=[];
    index800min=[];
    index800plus=[];
    
    % divide synapses in groups according to distance
    for nrsegment=1:nrsegments
        if(distances(nrsegment)<=200)
            total200min=total200min+1;
            index200min(total200min)=nrsegment;
        elseif(distances(nrsegment)<=400 && distances(nrsegment)>200)
            total400min=total400min+1;
            index400min(total400min)=nrsegment;
        elseif(distances(nrsegment)<=600 && distances(nrsegment)>400)
            total600min=total600min+1;
            index600min(total600min)=nrsegment;
        elseif(distances(nrsegment)<=800 && distances(nrsegment)>600)
            total800min=total800min+1;
            index800min(total800min)=nrsegment;
        elseif(distances(nrsegment)>800)
            total800plus=total800plus+1;
            index800plus(total800plus)=nrsegment;
        end
    end
    
    % determine mean in bins
    binsize=10;
    casrimaxbinned=binaverages(casrimax,binsize);
    size(casrimaxbinned)
    binvector=0.5*binsize:binsize:(nrruns-1);

    nrsyn200=sum(casrimaxbinned(:,index200min)>0,2);
    nrsyn400=sum(casrimaxbinned(:,index400min)>0,2);
    nrsyn600=sum(casrimaxbinned(:,index600min)>0,2);
    nrsyn800=sum(casrimaxbinned(:,index800min)>0,2);
    nrsyn800plus=sum(casrimaxbinned(:,index800plus)>0,2);

    %plot average casrimax for the different distances
    figure()
    hold on
    errorbar(binvector,nanmean(casrimaxbinned(:,index200min),2)*1000,nanstd(casrimaxbinned(:,index200min)')./sqrt(nrsyn200)'*1000,'r')
    errorbar(binvector,nanmean(casrimaxbinned(:,index400min),2)*1000,nanstd(casrimaxbinned(:,index400min)')./sqrt(nrsyn400)'*1000,'y')
    errorbar(binvector,nanmean(casrimaxbinned(:,index600min),2)*1000,nanstd(casrimaxbinned(:,index600min)')./sqrt(nrsyn600)'*1000,'g')
    errorbar(binvector,nanmean(casrimaxbinned(:,index800min),2)*1000,nanstd(casrimaxbinned(:,index800min)')./sqrt(nrsyn800)'*1000,'c')
    errorbar(binvector,nanmean(casrimaxbinned(:,index800plus),2)*1000,nanstd(casrimaxbinned(:,index800plus)')./sqrt(nrsyn800plus)'*1000,'b')
    %plot(binvector,binnedmean400*1000,'g')
    plot([1,nrruns-1],[thresholdca, thresholdca]*1000,'k--')
    set(gca,'FontName','Arial')
    xlabel('# run','Fontname','Arial')
    ylabel('peak calcium','Fontname','Arial')

    figure()
    hold on
    plot(binvector,nanmean(casrimaxbinned(:,index200min),2)*1000,'r')
    plot(binvector,nanmean(casrimaxbinned(:,index400min),2)*1000,'y')
    plot(binvector,nanmean(casrimaxbinned(:,index600min),2)*1000,'g')
    plot(binvector,nanmean(casrimaxbinned(:,index800min),2)*1000,'c')
    plot(binvector,nanmean(casrimaxbinned(:,index800plus),2)*1000,'b')
    %plot(binvector,binnedmean400*1000,'g')
    plot([1,nrruns-1],[thresholdca, thresholdca]*1000,'k--')
    set(gca,'FontName','Arial')
    xlabel('# run','Fontname','Arial')
    ylabel('peak calcium','Fontname','Arial')
    
    
    %plot casrimax with subthreshold responses
   
    figure()
    subplot('Position',[0.07 0.1 0.78 0.8])
    hold on
    for nrsegment=1:nrskippedtr:nrsegments
        plot(casrimax(1:nrruns-1,nrsegment)*1000,'Color',colorlistjet(sorteddistances2(nrsegment,3),:),'Linewidth',0.5)
    end
    plot([1,nrruns-1],[thresholdca, thresholdca]*1000,'k--')
    set(gca,'FontName','Arial')
    xlabel('# run','Fontname','Arial')
    ylabel('peak calcium','Fontname','Arial')

    subplot('Position',[0.85 0.12 0.1 0.76])
    hold on
    for nrsegment=1:nrsegments
        plot(1,nrsegment, 'sq','MarkerFaceColor', colorlistjet(nrsegment,:),'MarkerEdgeColor',colorlistjet(nrsegment,:),'MarkerSize',14) 
    end
    ylim([0,size(distances,2)])
    axis('off')
    text(1.5,1,'0','FontSize',10)
    text(1.5,488,num2str(max(distances),4),'FontSize',10)      
    
        %plot casrimax with subthreshold responses
   
    figure()
    subplot('Position',[0.07 0.1 0.78 0.8])
    hold on
    for nrsegment=1:nrskippedtr:nrsegments
        plot(casrimax(1:nrruns-1,nrsegment)*1000,'Color',colorlistjet(sorteddistances2(nrsegment,3),:),'Linewidth',0.5)
    end
    plot([1,nrruns-1],[thresholdca, thresholdca]*1000,'k--')
    set(gca,'FontName','Arial')
    xlabel('# run','Fontname','Arial')
    ylabel('peak calcium','Fontname','Arial')

    subplot('Position',[0.85 0.12 0.1 0.76])
    hold on
    for nrsegment=1:nrsegments
        plot(1,nrsegment, 'sq','MarkerFaceColor', colorlistjet(nrsegment,:),'MarkerEdgeColor',colorlistjet(nrsegment,:),'MarkerSize',14) 
    end
    ylim([0,size(distances,2)])
    axis('off')
    text(1.5,1,'0','FontSize',10)
    text(1.5,488,num2str(max(distances),4),'FontSize',10)   
 
    %plot casrimax only first 200 runs
    figure()
    subplot('Position',[0.07 0.1 0.78 0.8])
    hold on
    for nrsegment=1:nrskippedtr:nrsegments
        plot(casrimax(1:200,nrsegment)*1000,'Color',colorlistjet(sorteddistances2(nrsegment,3),:),'Linewidth',0.5)
    end
    plot([1,200],[thresholdca, thresholdca]*1000,'k--')
    set(gca,'FontName','Arial')
    xlabel('# run','Fontname','Arial')
    ylabel('peak calcium','Fontname','Arial')

    subplot('Position',[0.85 0.12 0.1 0.76])
    hold on
    for nrsegment=1:nrsegments
        plot(1,nrsegment, 'sq','MarkerFaceColor', colorlistjet(nrsegment,:),'MarkerEdgeColor',colorlistjet(nrsegment,:),'MarkerSize',14) 
    end
    ylim([0,size(distances,2)])
    axis('off')
    text(1.5,1,'0','FontSize',10)
    text(1.5,488,num2str(max(distances),4),'FontSize',10)  
    
    %plot casrimax as running average
        figure()
    subplot('Position',[0.07 0.1 0.78 0.8])
    hold on
    for nrsegment=1:nrskippedtr:nrsegments
        plot(binvector, casrimaxbinned(:,nrsegment)*1000,'-','Color',colorlistjet(sorteddistances2(nrsegment,3),:))
    end
    plot([1,nrruns-1],[thresholdca, thresholdca]*1000,'k--')
    set(gca,'FontName','Arial')
    xlabel('# run','Fontname','Arial')
    ylabel('peak calcium','Fontname','Arial')

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
    set(gca,'FontName','Arial')
    xlabel('# run','Fontname','Arial')
    ylabel('peak calcium','Fontname','Arial')
    
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
    set(gca,'FontName','Arial')
    ylabel('nr stims','Fontname','Arial')
    xlabel('segments nr','Fontname','Arial')
    subplot(1,2,2)
    plot(distances,nrstimperseg,'.')
    xlabel('distances','Fontname','Arial')
    
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
    set(gca,'FontName','Arial')
    xlabel('# stim','Fontname','Arial')
    ylabel('synaptic strength','Fontname','Arial')

    subplot('Position',[0.85 0.12 0.1 0.76])
    hold on
    for nrsegment=1:nrsegments
        plot(1,nrsegment, 'sq','MarkerFaceColor', colorlistjet(nrsegment,:),'MarkerEdgeColor',colorlistjet(nrsegment,:),'MarkerSize',14) 
    end
    ylim([0,size(distances,2)])
    axis('off')
    text(1.5,1,'0','FontSize',10)
    text(1.5,488,num2str(max(distances),4),'FontSize',10)
    
    % plot synweight only 200 runs
        figure()
    subplot('Position',[0.07 0.1 0.78 0.8])
    hold on
    for nrsegment=1:nrskippedtr:nrsegments
        plot(synweight(1:200,nrsegment)*10^6,'Color',colorlistjet(sorteddistances2(nrsegment,3),:))
    end
    set(gca,'FontName','Arial')
    xlabel('# stim','Fontname','Arial')
    ylabel('synaptic strength','Fontname','Arial')

    subplot('Position',[0.85 0.12 0.1 0.76])
    hold on
    for nrsegment=1:nrsegments
        plot(1,nrsegment, 'sq','MarkerFaceColor', colorlistjet(nrsegment,:),'MarkerEdgeColor',colorlistjet(nrsegment,:),'MarkerSize',14) 
    end
    ylim([0,size(distances,2)])
    axis('off')
    text(1.5,1,'0','FontSize',10)
    text(1.5,488,num2str(max(distances),4),'FontSize',10)
    
    % plot average synweights
    synweightsbinned=binaverages(synweight,binsize);
    nrsyn200_w=sum(synweightsbinned(:,index200min)>0,2);
    nrsyn400_w=sum(synweightsbinned(:,index400min)>0,2);
    nrsyn600_w=sum(synweightsbinned(:,index600min)>0,2);
    nrsyn800_w=sum(synweightsbinned(:,index800min)>0,2);
    nrsyn800plus_w=sum(synweightsbinned(:,index800plus)>0,2);
    
    figure()
    hold on
    errorbar(binvector,nanmean(synweightsbinned(:,index200min),2)*10^6,nanstd(synweightsbinned(:,index200min)')./sqrt(nrsyn200_w)'*10^6,'r')
    errorbar(binvector,nanmean(synweightsbinned(:,index400min),2)*10^6,nanstd(synweightsbinned(:,index400min)')./sqrt(nrsyn400_w)'*10^6,'y')
    errorbar(binvector,nanmean(synweightsbinned(:,index600min),2)*10^6,nanstd(synweightsbinned(:,index600min)')./sqrt(nrsyn600_w)'*10^6,'g')
    errorbar(binvector,nanmean(synweightsbinned(:,index800min),2)*10^6,nanstd(synweightsbinned(:,index800min)')./sqrt(nrsyn800_w)'*10^6,'c')
    errorbar(binvector,nanmean(synweightsbinned(:,index800plus),2)*10^6,nanstd(synweightsbinned(:,index800plus)')./sqrt(nrsyn800plus_w)'*10^6,'b')
    set(gca,'FontName','Arial')
    xlabel('# stim','Fontname','Arial')
    ylabel('synaptic strength','Fontname','Arial')
    ylim([0,500])
    
    % plot final synweights STILL TO CHANGE 200 VALUES
    figure()
    hold on
    for nrsegment=1:nrsegments
        if(max(nrsegment==trunksri+1)==1)
            % trunk segment
            plot(distances(nrsegment),synweight(1,nrsegment)*10^6,'xk')
            plot(distances(nrsegment),synweight(nrruns-1,nrsegment)*10^6,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4) 
        else
            % non trunk segment
            plot(distances(nrsegment),synweight(1,nrsegment)*10^6,'xk')
            plot(distances(nrsegment),synweight(nrruns-1,nrsegment)*10^6,'o','MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',4)            
        end
    end
    
    %plot([0,max(distances)],[200,200],'-k')
    xlim([0 1000])
    ylim([0,500])
    set(gca,'FontName','Arial')
    xlabel('distance','Fontname','Arial')
    ylabel('synaptic weight','Fontname','Arial')
    title('synaptic weight at start (open) and end of scaling (closed)')
    
    figure()
    set(gca,'FontName','Arial')
    plot(synweight(1:200,:),casrimax(1:200,:),'.')
    ylabel('peak calcium','Fontname','Arial')
    xlabel('synaptic strength','Fontname','Arial')
    
    
end