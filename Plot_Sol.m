myFolder = 'Output/Subdomains';
filepattern = fullfile(myFolder, '*.txt');
files = dir(filepattern);
labels = [];
for i=1:length(files)
  files(i).name;
  if contains(files(i).name, "Sol_")
    label = textscan(files(i).name,"Sol_%s");
    labels = [labels;label{1}{1}];
  end
end
[r,c] = size(labels);
omegax = 0.23; omegay = 0.49; omegapx = 0.331; omegapy = 0.667;
figure 
PDD_Sol = subplot(1,2,1);
Error = subplot(1,2,2);
hold(PDD_Sol,'on');
hold(Error,'on');
for i=1:r
    label = labels(i,:);
    filename = sprintf("Output/Subdomains/X_%s",labels(i,:));
    x = load(filename,'-ascii');
    len_x = length(x);
    filename = sprintf("Output/Subdomains/Y_%s",labels(i,:));
    y = load(filename,'-ascii');
    len_y = length(y);
    filename = sprintf("Output/Subdomains/Sol_%s",labels(i,:));
    u = load(filename,'-ascii');
    [xx,yy] = meshgrid(x,y);
    uu = reshape(u,len_x,len_y);
    subplot(PDD_Sol)
    colormap(PDD_Sol,'parula');
    s = surf(xx,yy,uu);
    s.EdgeColor = 'none';
    u_a = sin(omegax*pi*xx(:) + omegay*pi*yy(:)) + cos(omegapx*pi*xx(:) + omegapy*pi*yy(:));
    uu_a = reshape(u_a,len_x,len_y);
    subplot(Error)
    colormap(Error, 'gray');
    er = uu - uu_a;
    err = reshape(er,len_x,len_y);
    s = surf(xx,yy,err);
    s.EdgeColor = 'none';
end
subplot(PDD_Sol)
axis square
colorbar
title('PDD Sparse Solution')
subplot(Error)
title('PDD Sparse Error')
axis square 
colorbar