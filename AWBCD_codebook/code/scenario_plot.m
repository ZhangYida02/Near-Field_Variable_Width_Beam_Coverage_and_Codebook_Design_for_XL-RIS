figure();  
hold on; 

%
edges = [  
    1, 2;  
    2, 3; 
    3, 4;  
    4, 1;  
      
    5, 6;  
    6, 7;
    7, 8;   
    8, 5;
      
    1, 5;  
    2, 6; 
    3, 7;   
    4, 8; 
];  
for i = 1:size(edges, 1)  
    line([x_room(edges(i, 1)), x_room(edges(i, 2))], [y_room(edges(i, 1)), y_room(edges(i, 2))], [z_room(edges(i, 1)), z_room(edges(i, 2))], 'Color', 'k');  
end  



%
scatter3(AN_x_rotated,AN_y_rotated,AN_z_rotated, 'ro'); 
% plot3([x_BS,x_BS+AN_v(1)],[y_BS,y_BS+AN_v(2)],[z_BS,z_BS+AN_v(3)],'b')
scatter3(x_BS,y_BS, z_BS, 'g*'); 
plot3([x_BS,x_BS],[y_BS,y_BS],[0,z_BS],'k')

%  
scatter3(x_cell,y_cell,z_cell, 'ro'); 
scatter3(x_RIS,y_RIS, z_RIS, 'g*'); 
axis([0 40 0 40 -1 10]) 

%
scatter3(x_aim,y_aim,z_aim,'y'); 
% scatter3(aim_x,aim_y,aim_z, 'g*'); 

% scatter3(mother_x, mother_y,mother_z, 'ro'); % 母点
% scatter3(x_user, y_user,z_user ,'b.');                % 子点

hold off; 

xlabel('X');  
ylabel('Y');  
zlabel('Z');  
grid on
view(3)
daspect([1 1 1])