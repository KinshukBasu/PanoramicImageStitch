data = xlsread('points.xlsx');

N = size(data,1);
data_left = data(:,1:2);
data_right = data(:,3:4);

% A = [0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0];
A = [];
D = [];

for i=1:1:N,
    tempmatrix = [0,0,0,-data(i,1),-data(i,2),-1,data(i,1)*data(i,4),data(i,2)*data(i,4);...
                data(i,1),data(i,2),1,0,0,0,-data(i,1)*data(i,3),-data(i,2)*data(i,3)];
    
    A = [A;tempmatrix];
    
    tempvector = [-data(i,4);data(i,3)];
    D = [D;tempvector];
end

h = inv(A'*A)*A'*D;
H = [h(1),h(2),h(3);h(4),h(5),h(6);h(7),h(8),1];

clearvars tempmatrix tempvector A D;

pic1 = imread('leftpic.jpg');
pic2 = imread('rightpic.jpg');
[maxrow,maxcol,dimen] = size(pic1);

rmax_result = maxrow+100;
cmax_result = 2*maxcol;
result = uint8(zeros(rmax_result,cmax_result,dimen));

origin_x = maxcol;
origin_y = 100;

% %---------Forward mapping method-----------------------
% for i=1:1:maxrow,
%     for j=1:1:maxcol,
%         result(origin_y+i,origin_x+j,:) = pic2(i,j,:);
%         
%         old_coods = [j;i;1];        %Bcz j is columns -> X dirn, i is rows -> Y dirn
%         wvector = (H*old_coods)';
%         wx = wvector(1)/wvector(3);
%         wy = wvector(2)/wvector(3);
% 
%         
%         xnew = int64(origin_x+ wx);
%         ynew = int64(origin_y+ wy);        
%         result(ynew,xnew,:) = pic1(i,j,:);
%     end
% end
% imshow(result);

%-----------Reverse Mapping Method----------------
% %Say I fix the origin at (50,maxcol)

for i=1:1:maxrow,
    for j=1:1:maxcol,
        result(origin_y+i,origin_x+j,:) = pic2(i,j,:);
    end
end

for i=1:1:rmax_result,
    for j=1:1:cmax_result,
        xdash = j-origin_x;
        ydash = i-origin_y;
        
        wvector_old = inv(H)*[xdash;ydash;1];
        
        xold = int64(wvector_old(1)/wvector_old(3));
        yold = int64(wvector_old(2)/wvector_old(3));
        
        if( (xold>0 && xold<=maxcol) && (yold>0 && yold<=maxrow) ),
            result(i,j,:) = pic1(yold,xold,:);
        end
    end
end
imshow(result);
        
        
