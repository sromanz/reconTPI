function dcf = dcftpi(crd,mrprot)

p = mrprot.MeasYaps.sWiPMemBlock.adFree{2}; % [%] 
ro       = mrprot.Meas.BaseResolution;

wi=zeros(size(crd,2)*size(crd,3)*size(crd,4),1);
temp =reshape(crd,3,[])';
w = sqrt(temp(:,1).^2+temp(:,2).^2+temp(:,3).^2);



for i=1:length(w)
    
    if (w(i) < 0.5*p)
        wi(i)=w(i)^2;
    else
        wi(i)=wi(i-1);
    end
end

dcf=reshape(wi,ro,1,[]); % data structure required by the gridding algorithm
end