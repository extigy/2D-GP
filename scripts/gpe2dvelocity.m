function [mgx,mgy,velx,vely] = gpe2dvelocity(gridx,gridy,phase)
        dspace=(gridx(2)-gridx(1));
        dims = size(phase);
        for i = 1:dims(1)
        for j = 1:dims(2)
            temp1 = anglediff(phase(bc(i+1,dims(1)),j),phase(bc(i-1,dims(1)),j));
            vely(i,j) = real(temp1)/(2*dspace);
            mgx(i,j) = gridx(j);
            mgy(i,j) = gridy(i);
        end
        end
        
        for i = 1:dims(1)
        for j = 1:dims(2)
            temp1 = anglediff(phase(i,bc(j+1,dims(2))),phase(i,bc(j-1,dims(2))));
            velx(i,j) = real(temp1)/(2*dspace);
        end
        end
end

function r = bc(i,dim)
    if(i == 0)
        r = dim;
    elseif(i == dim+1)
        r = 1;
    else
        r = i;
    end
end
    
function d = anglediff(th1, th2)
    d = atan2(sin(th1-th2), cos(th1-th2));
end   