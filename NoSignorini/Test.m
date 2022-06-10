%%Test：将疏网格下的U值插值到密网格的UInterp值（正确）

X = [0,1];                         
Y = [0,1]; 
Time=1; 

intervals1=3;
intervals2=9;
Rate=round(intervals2/intervals1);

FNd1=intervals1*(intervals1+1);
FNd2=intervals2*(intervals2+1);
UInterp=zeros(FNd2,1);
U=1:1:FNd1;

for kx=1:(intervals1+1)
    for jx=1:intervals1
        if jx==1
            for ix=1:Rate
                UInterp( (kx-1)*intervals2*Rate+ix )=ix/Rate*U( (kx-1)*intervals1+jx);
            end
        else
            for ix=1:Rate
                UInterp( (kx-1)*intervals2*Rate+Rate*(jx-1)+ix ) =(1-ix/Rate) *U( (kx-1)*intervals1+jx-1 ) + ix/Rate *U( (kx-1)*intervals1+jx );
            end
        end
    end 
end

for kx=1:intervals1
    for jx=1:(Rate-1)
        for ix=1:intervals2
            UInterp( (kx-1)*intervals2*Rate+ jx*intervals2 +ix  )=(1-jx/Rate) *UInterp( (kx-1)*intervals2*Rate+ix ) +...
                jx/Rate *UInterp( kx*intervals2*Rate+ix );
        end
    end
end

