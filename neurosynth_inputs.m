%Author: Valter Lundegårdh

%Inspiration from Marguiles 2016 and Petri 2019

%% Generate maps for Neurosynth analysis (which will be in Python, code from Margulies 2016)
%The maps are saved in the folder data, inside the code folder

%select reference nifti map (MNI)
ref_nifti_path=which('HCP-MMP1_cortices_2mm.nii');
%ref_nifti_path=which('HCP-MMP1_onMNI152_2mm_Glasser360.nii');
%ref_nifti_path=which('MNI152_T1_2mm_brain.nii');
refhdr=spm_vol(ref_nifti_path);
refnii=spm_read_vols(refhdr);

% Adapted the image to have values 0-44 instead of 
% 0-22 and 101-122
refnii(refnii > 100)=refnii(refnii > 100)-100+22;

%select map to split into percentiles
map_to_consider=meanSDI(:,g);
map_to_consider(isnan(map_to_consider)==1)=0; %zero nans

NSmaps=neurosynth_create_maps(map_to_consider,refnii);

mkdir(strcat(mypath,'/results/my_masks_group',num2str(g)));
%save maps
for i=1:9
    refhdr.fname=strcat(mypath,'/results/my_masks_group',num2str(g),'/groups_SDI_0',num2str(i),'.nii');
    spm_write_vol(refhdr,NSmaps{i});
end
for i=10:20
    refhdr.fname=strcat(mypath,'/results/my_masks_group',num2str(g),'/groups_SDI_',num2str(i),'.nii');
    spm_write_vol(refhdr,NSmaps{i});
end
%%
% images using all SDI values for percentiles
NSmaps=neurosynth_create_maps_all(map_to_consider,refnii,reshape(meanSDI,1,[]));

mkdir(strcat(mypath,'/results/my_masks_group',num2str(g)));
%save maps
for i=1:9
    refhdr.fname=strcat(mypath,'/results/my_masks_group',num2str(g),'/all_SDI_0',num2str(i),'.nii');
    spm_write_vol(refhdr,NSmaps{i});
end
for i=10:20
    refhdr.fname=strcat(mypath,'/results/my_masks_group',num2str(g),'/all_SDI_',num2str(i),'.nii');
    spm_write_vol(refhdr,NSmaps{i});
end

%% for visualization of gradient chopped in 5th percentiles
sequence=[-1:0.1053:1];
sequence=[sequence,1];
clear meanSDIchopped

sequence_percentiles=prctile(meanSDI(:,g),5:5:100);
for i=1:nROI %for each node, which percentile it belongs to?
    [dist,index]=min(abs(meanSDI(i,g)-sequence_percentiles));
    if meanSDI(i,g)>sequence_percentiles(index)
        meanSDIchopped(i,g)=sequence(index+1);
    else
        meanSDIchopped(i,g)=sequence(index);
    end
end
