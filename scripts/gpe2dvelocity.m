function [mgx,mgy,velx,vely] = gpe2dvelocity(gridx,gridy,ophase)
    
        dspace=(gridx(2)-gridx(1));
        phase = unwrap(ophase);
        dims = size(ophase);
        for i = 2:dims(1)-1
        for j = 2:dims(2)-1
            if (phase(i+1,j)-phase(i-1,j)<-(pi/2.0d0))
                temp1 = phase(i+1,j)-(phase(i-1,j) - pi);
            elseif (phase(i+1,j)-phase(i-1,j)>(pi/2.0d0))
                temp1 = phase(i+1,j)-(phase(i-1,j) + pi);
            else
                temp1 = phase(i+1,j)-phase(i-1,j);
            end
            mgx(i,j) = gridx(j);
            mgy(i,j) = gridy(i);
            vely(i,j) = real(temp1)/(2*dspace);
        end
            vely(i,dims(2)) = vely(i,dims(2)-1);
            vely(i,1) = vely(i,2);
        end
        for j = 1:dims(2)
            vely(dims(1),j) = vely(dims(1)-1,j);
            vely(1,j) = vely(2,j);
        end


        phase = unwrap(ophase,[],2);
        for i = 2:dims(1)-1
        for j = 2:dims(2)-1
            if (phase(i,j+1)-phase(i,j-1)<-(pi/2.0d0))
                temp1 = phase(i,j+1)-(phase(i,j-1) - pi);
            elseif (phase(i,j+1)-phase(i,j-1)>(pi/2.0d0))
                temp1 = phase(i,j+1)-(phase(i,j-1) + pi);
            else
                temp1 = phase(i,j+1)-phase(i,j-1);
            end
            velx(i,j) = real(temp1)/(2*dspace);
        end
            velx(i,dims(2)) = velx(i,dims(2)-1);
            velx(i,1) = velx(i,2);
        end
        for j = 1:dims(2)
            velx(dims(1),j) = velx(dims(1)-1,j);
            velx(1,j) = velx(2,j);
        end

end
