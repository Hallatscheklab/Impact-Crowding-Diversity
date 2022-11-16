%Prepare a thresholded two-color movie in tif format = filename

%Prepare the avi file of the tif file = videofilename.
%This is for checking trajectories of clusers.

%New files are produced in "Processed" directory.
%Open them with "virtual stack" to check your manual tracks.

currentpath=pwd;
filename = '0120-chamber-Otsu-Thresh.tif';
videofilename = '0120-chamber-Otsu-Thresh-jpeg.avi';
path = fullfile(currentpath,filename);
info = imfinfo(path);
MaxIndex = length(info);

dummy = struct('Area',{0},'Centroid',{[0.0,0.0]},'Extrema',{zeros(8,2,'double')},'time',{0});


% save new tiff images with circles
mkdir(fullfile(currentpath,'Processed'));
for i = 1 : MaxIndex
    im1_full = imread(path,'index',i);
    CC = bwconncomp(im1_full);
    Unassigned = regionprops(CC,'Centroid','Area','Extrema');
    [Unassigned(:).time] = deal(i);
    im1_inserted = im1_full;
    for j=1:length(Unassigned)
        x = Unassigned(j).Centroid(1);
        y = Unassigned(j).Centroid(2);
        rad = sqrt(Unassigned(j).Area);
        im1_inserted =insertShape(im1_inserted,'Circle',[x,y,rad],'Color','red');
    end
    newfilename = strcat('Track_',num2str(i,'%03d'),'.tif');
    imwrite(uint8(im1_inserted),fullfile(currentpath,'Processed',newfilename));
end



%open video for tracking
implay(videofilename);





%main semi-manual detection part
Alltracks = struct([]);
while 1
    %start new tracking
    prompt = 'add a new track? enter frame# to start. enter 0 to end  ';
    startframe = input(prompt);
    if startframe < 1
        break;
    end
    
    
    %Tracking part
    frame = startframe;
    track = dummy;
    while 1
        im_Ori=imread(path,'index',frame);
        CC = bwconncomp(im_Ori);
    
        Unassigned = regionprops(CC,'Centroid','Area','Extrema');
        [Unassigned(:).time] = deal(frame);
        
        newfilename = strcat('Track_',num2str(frame,'%03d'),'.tif');
        im_circ=imread(fullfile(currentpath,'Processed',newfilename));
        figure('Name',strcat('Ori_frame',int2str(frame))), imshow(im_circ);
        [Xclick,Yclick] = getpts;

        index = zeros(1,length(Xclick));
        for i = 1:length(Xclick)
            dist = zeros(1,length(Unassigned));
            for j = 1:length(Unassigned)
                dist(j) = sqrt((Xclick(i) - Unassigned(j).Centroid(1))^2 + (Yclick(i) - Unassigned(j).Centroid(2))^2);
            end
            index(i) = find(dist == min(dist));
        end
        
        im_inserted = im_circ;
        for i=index
            x = Unassigned(i).Centroid(1);
            y = Unassigned(i).Centroid(2);
            rad = sqrt(Unassigned(i).Area);
            im_inserted =insertShape(im_inserted,'Circle',[x,y,rad],'Color','yellow','LineWidth',2);
            display(i)
        end
        figure('Name','Confirmation'), imshow(im_inserted);
        prompt = 'OK? y[1]/no, try again[2]/end without this frame[3]/skip this frame[4]  ';
        confirmation = input(prompt);
        if confirmation == 1
            track2= Unassigned(index);
            for i = 1:length(track2)
                track(end+1) = track2(i);
            end
            newfilename = strcat('Track_',num2str(frame,'%03d'),'.tif');
            imwrite(uint8(im_inserted),fullfile(currentpath,'Processed',newfilename));
            frame = frame + 1;
            close all;
        elseif confirmation == 2
            close all;
        elseif confirmation == 3
            close all;
            break;
        elseif confirmation == 4
            frame = frame + 1;
            close all;
        end
    end
    Alltracks(end+1).track = track;
            
end

save('data.mat','Alltracks');