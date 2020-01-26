function[databinned]=binaverages(data,bininterval)
    nrbins=floor(size(data,1)/bininterval);
    for i=1:nrbins
        indexstart=bininterval*(i-1)+1;
        indexend=bininterval*i;
        databinned(i,:)=nanmean(data(indexstart:indexend,:));
    end
end