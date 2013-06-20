% This is a script to load the ideogram chromosome data.
clear

% Skip the header.
file = fopen('/tmp/pcarbo/ideogram','r');
fgetl(file);
data = textscan(file,'%s%s%s%d%d%d%d%s','Delimiter',',');
fclose(file);

% Convert the data to a struct.
ideogram.chr   = data{1};
ideogram.arm   = data{2};
ideogram.band  = data{3};
ideogram.start = data{6};
ideogram.stop  = data{7};
ideogram.stain = data{8};

% Save the ideogram data to a MAT file.
save('/tmp/pcarbo/ideogram.mat','ideogram','-v7.3');
