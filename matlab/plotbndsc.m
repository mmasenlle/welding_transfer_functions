function plotbndsc(bdb,ptype,phase)


pos=[1 1 1 1];
[bdb,ptype,phase,axs,pos,wbs,wbs2,coora,coorb]=qplotdef(bdb,ptype,phase,pos);


[rbdb,cbdb]=size(bdb);
lwbs2=length(wbs2);

if lwbs2 ~= cbdb,
    
    a=gca;
    apos=get(a,'pos');
    set(a,'box','on','xgrid','on','ygrid','on',...
        'gridlinestyle',':','drawmode','fast',...
        'nextplot','add','xlim',axs(1:2),'ylim',axs(3:4));
    
    bnd=qplotbd(phase,bdb,coora,coorb,axs);
    
else
    disp('Unplottable bounds.');
end
end