csv=readtable('/codeFeb/gravity_centers_HCP_atlas_sorted.csv');

%Change IDs to be 1-44
csv(23:end,2)=csv(23:end,2)-100+22;
id=table2array(csv(1:end,2));
num=44; %44 regions
name=cellstr(num2str(table2array(csv(:,2))));
sname=cellstr(num2str(table2array(csv(:,2))));

%Coordinates, hacky solution to align with brain figure
csv(:,3)=csv(:,3)-90; %x 90
csv(:,4)=csv(:,4)-125; %y 125
csv(:,5)=csv(:,5)-75; %z

% Swap x,y,z, left and right coordinates. 
% In the original data the left has a larger x than right.
for i=3:5
    temp=csv(1:22,i);
    csv(1:22,i)=csv(23:end,i);
    csv(23:end,i)=temp;
end

center=num2cell(table2array(csv(:,3:5)),2);

myCodeBook.id=id;
myCodeBook.num=num;
myCodeBook.name=name;
myCodeBook.sname=sname;
myCodeBook.center=center;

save('./gsp_data/dataFeb/misc/myCodeBook.mat','myCodeBook')


