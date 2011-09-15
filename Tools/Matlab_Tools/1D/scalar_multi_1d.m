% scalar_multi_1d('directory',filenumber,'variable')
function scalar_multi_1d(directory,variables,filenumber,filecount)

clf;
[m,n]=size(variables);
endfilenumber=filenumber+filecount-1;
nums=[filenumber:1:endfilenumber];
row=0;
for entry=variables
    column=0;
    for num=nums
        num=filecount*row+column+1;
        subplot(n,filecount,num);
        scalar_1d(directory,num,char(entry),0);
        column=column+1;
    end
    row=row+1;
end
