file = fopen('data_time0.25.txt','r');
maxtemp = str2double(fgetl(file));
timesteps = str2double(fgetl(file));
moves = str2double(fgetl(file));
x_v = str2num(fgetl(file));
n = size(x_v,2);
y_v = str2num(fgetl(file));
A = zeros(moves,5);
T_v = zeros(timesteps+1,n);
for i = 1:moves
    A(i,:) = str2num(fgetl(file));
end
for i = 1:(timesteps+1)
    T_v(i,:) = str2num(fgetl(file));
end
fclose(file);

W=sortrows(A,[1 3]);
style = cell(1, moves);
edgelabel = cell(1, moves);
nodelabel = cell(1,n);
T_e = zeros( moves,1);

for j=1:n
    nodelabel{j} = '';
end

h = figure;
axis tight;
filename = 'Welding.gif';
whitebg([0.5 0.5 0.5]);


for k = 1: (timesteps+1)
    for i = 1: moves
         if W(i,2) >= k
             style{i} = 'none';
             edgelabel{i} = '';
         elseif W(i,5)>0.5
             style{i} = '--';
             edgelabel{i} = int2str(W(i,2));
         else
             style{i} = '-';
             edgelabel{i} = int2str(W(i,2));
         end
     end

    for i = 1: moves
        if W(i,5)>0.5
            T_e(i)=0;
        else
            T_e(i)=(T_v(k,W(i,1))+T_v(k,W(i,3)))/2;
        end
     end

    G=digraph(W(:,1)',W(:,3)');

    p=plot(G,'XData',x_v,'YData',y_v);
    p.EdgeCData = T_e;
    p.EdgeColor = 'flat';
    p.LineStyle = style;
    p.EdgeLabel = edgelabel;
    p.MarkerSize = 6;
    p.NodeCData = T_v(k,:);
    p.NodeColor = 'flat';
    p.NodeLabel = nodelabel;
    colormap hot;
   caxis([0 maxtemp]);
    colorbar;
    view(2)

    frame = getframe(h);
    image = frame2im(frame);
    imind = image;
    [imind, cm] = rgb2ind(image,256);

    if k == 1
        imwrite(imind,cm,filename,'gif','DelayTime',1,'Loopcount', inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end
