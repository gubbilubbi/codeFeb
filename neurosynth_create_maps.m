%Author: Valter Lundegårdh

%Inspiration from Marguiles 2016 and Petri 2019

%Create maps as inputs for Neurosynth 

%% map back atlas based values to nifti file
function NSmaps=neurosynth_create_maps(vec_to_consider,refnii)

% division in percentiles bins like Margulies 2016
for i=1:20
   thr1=percentile(vec_to_consider,(i-1)*5);
   thr2=percentile(vec_to_consider,i*5);

   map1=vec_to_consider>thr1;
   map2=vec_to_consider<=thr2;
   map=map1+map2;

   % Only where both thresholds are satisfied do we
   % want the region, i.e. map==2
   NSvecs(:,i)=map==2; 
   %disp(NSvecs)

   f=find(NSvecs(:,i)==1);
   %disp(f)
   NSmaps{i}=zeros(size(refnii));

   for j=1:size(f,1)
    NSmaps{i}(round(refnii)==f(j))=1; %Set those values to 1
   end
end

end