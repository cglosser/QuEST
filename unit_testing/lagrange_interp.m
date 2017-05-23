pnts_count = 6;
scale = 10;
syms x
order = 3;
%random points
x_pnts = rand(1,pnts_count).*scale;
y_pnts = rand(1,pnts_count).*scale;
%non-random points
x_pnts = [0 1 2 3 4 5];
y_pnts = [0 1 2 3 4 5];

poly_mat = sym(zeros(pnts_count,1));

%polynomial calculations
for j=1:pnts_count;
    poly = 1;
    for i=1:pnts_count;
        if i~=j;
            poly = poly*(x-x_pnts(i))/(x_pnts(j)-x_pnts(i));
        end
    end
    poly_mat(j) = poly;
end

% derivative calculations 
poly_mat_d1 = diff(poly_mat,x,1);
poly_mat_d2 = diff(poly_mat,x,2);
poly_mat_d22 = diff(poly_mat_d1,x,1);

x = x_pnts;

d0_values = sum(subs(poly_mat));
d1_values = sum(subs(poly_mat_d1));
d2_values = sum(subs(poly_mat_d2));

storage_matrix = [d0_values;d1_values;d2_values]

%plotting
x = min(x_pnts):.01:max(x_pnts);
%hold on
scatter(x_pnts,d2_values)
scatter(x_pnts,d1_values)
scatter(x_pnts,d0_values)
plot(x,sum(subs(poly_mat)))
plot(x,sum(subs(poly_mat_d1)))
plot(x,sum(subs(poly_mat_d2)))
%hold off

x = 0:0.1:5;
%plot(x,subs(poly_mat))
hold on 
plot(x,subs(poly_mat_d1))

hold off