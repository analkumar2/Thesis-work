function[datamatrix]=removezeros(datamatrix)
    [nrruns,nrsegments]=size(datamatrix);
    for i=1:nrruns
        for m=1:nrsegments
            if(datamatrix(i,m)==0)
                datamatrix(i,m)=NaN;
            end
        end
    end
end