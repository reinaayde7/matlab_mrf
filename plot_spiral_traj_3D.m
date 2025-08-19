
x = kxall(adcpad+1:uplimit,1);
y = kyall(adcpad+1:uplimit,1);
z_levels = linspace(0,15,15);

figure, hold on
view(3)

for i=1:length(z_levels)
    if mod(i,2)
        z = z_levels(i) *ones(size(x));
        plot3(x,y,z,'b', 'LineWidth',2)
    else
        z = z_levels(i) *ones(size(x));
        plot3(x,y,z,'r--')
    end
end
axis off

% axis([-5,5,-5,5]);axis square;
% title('1 spiral arm'); xlabel('kx (1/cm)');ylabel('ky (1/cm)');