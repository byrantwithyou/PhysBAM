function x0 = data_zeros(x,y)

    % Indices of Approximate Zero-Crossings
    % (you can also use your own 'find' method here, although it has 
    %  this pesky difference of 1-missing-element because of diff...)
    dy = find(y(:).*circshift(y(:), [-1 0]) <= 0);   

    % Do linear interpolation of near-zero-crossings
    x0 = NaN(size(dy,1)-1,1);    
    for k1 = 1:size(dy,1)-1

        b = [[1;1] [x(dy(k1)); x(dy(k1)+1)]] \ ...   
            [y(dy(k1)); y(dy(k1)+1)]; 

        x0(k1) = -b(1)/b(2);   

    end

end

