function [mgx,mgy,velx,vely] = gpe2dvelocity(gridx,gridy,ophase)
    
        dspace=(gridx(2)-gridx(1));
        phase = unwrap(ophase);
        dims = size(ophase);
        for i = 1:dims(1)-1
        for j = 1:dims(2)-1
            if (phase(i+1,j)-phase(i,j)<-(pi/2.0d0))
                temp1 = phase(i+1,j)-(phase(i,j) - pi);
            elseif (phase(i+1,j)-phase(i,j)>(pi/2.0d0))
                temp1 = phase(i+1,j)-(phase(i,j) + pi);
            else
                temp1 = phase(i+1,j)-phase(i,j);
            end
            mgx(i,j) = gridx(j);
            mgy(i,j) = gridy(i);
            vely(i,j) = real(temp1)/dspace;
        end
            vely(i,dims(2)) = vely(i,dims(2)-1);
        end
        for j = 1:dims(2)
            vely(dims(1),j) = vely(dims(1)-1,j);
        end
        phase = unwrap(ophase,[],2);
        for i = 1:dims(1)-1
        for j = 1:dims(2)-1
            if (phase(i,j+1)-phase(i,j)<-(pi/2.0d0))
                temp1 = phase(i,j+1)-(phase(i,j) - pi);
            elseif (phase(i,j+1)-phase(i,j)>(pi/2.0d0))
                temp1 = phase(i,j+1)-(phase(i,j) + pi);
            else
                temp1 = phase(i,j+1)-phase(i,j);
            end
            velx(i,j) = real(temp1)/dspace;
        end
            velx(i,dims(2)) = velx(i,dims(2)-1);
        end
        for j = 1:dims(2)
            velx(dims(1),j) = velx(dims(1)-1,j);
        end

end
