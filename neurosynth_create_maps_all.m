%Author: Valter Lundegårdh

%Inspiration from Marguiles 2016 and Petri 2019

%Create maps as inputs for Neurosynth 

%% map back atlas based values to nifti file
function NSmaps=neurosynth_create_maps_all(vec_to_consider,refnii,allSDI)

% division in percentiles bins like Margulies 2016
for i=1:20
   thr1=percentile(allSDI,(i-1)*5);
   thr2=percentile(allSDI,i*5);

   map1=vec_to_consider>thr1;
   map2=vec_to_consider<=thr2;
   map=map1+map2;

   if(sum(map==2)==0)
       disp(sum(map==2))
       disp("No region found")
   end
   % Only where both thresholds are satisfied do we
   % want the region, i.e. map==2
   NSvecs(:,i)=map==2; 

   f=find(NSvecs(:,i)==1);
   NSmaps{i}=zeros(size(refnii));

   for j=1:size(f,1)
    NSmaps{i}(round(refnii)==f(j))=1; %Set those values to 1
   end
end

end