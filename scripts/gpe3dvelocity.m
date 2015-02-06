function [mgx,mgy,mgz,velx,vely,velz] = gpe3dvelocity(gridx,gridy,gridz,ophase)
    
        dspace=(gridx(2)-gridx(1));
        phase = unwrap(ophase);
        dims = size(ophase);
        for i = 2:dims(1)-1
        for j = 2:dims(2)-1
        for k = 2:dims(3)-1
            if (phase(i+1,j,k)-phase(i-1,j,k)<-(pi/2.0d0))
                temp1 = phase(i+1,j,k)-(phase(i-1,j,k) - pi);
            elseif (phase(i+1,j,k)-phase(i-1,j,k)>(pi/2.0d0))
                temp1 = phase(i+1,j,k)-(phase(i-1,j,k) + pi);
            else
                temp1 = phase(i+1,j,k)-phase(i-1,j,k);
            end
            mgx(i-1,j-1,k-1) = gridx(i);
            mgy(i-1,j-1,k-1) = gridy(j);
            mgz(i-1,j-1,k-1) = gridz(k);
            velx(i-1,j-1,k-1) = real(temp1)/dspace;
        end
        end
        end
        phase = unwrap(ophase,[],2);
        for i = 2:dims(1)-1
        for j = 2:dims(2)-1
        for k = 2:dims(3)-1
            if (phase(i,j+1,k)-phase(i,j-1,k)<-(pi/2.0d0))
                temp1 = phase(i,j+1,k)-(phase(i,j-1,k) - pi);
            elseif (phase(i,j+1,k)-phase(i,j-1,k)>(pi/2.0d0))
                temp1 = phase(i,j+1,k)-(phase(i,j-1,k) + pi);
            else
                temp1 = phase(i,j+1,k)-phase(i,j-1,k);
            end
            vely(i-1,j-1,k-1) = real(temp1)/dspace;
        end
        end
        end
        phase = unwrap(ophase,[],3);
        for i = 2:dims(1)-1
        for j = 2:dims(2)-1
        for k = 2:dims(3)-1
            if (phase(i,j,k+1)-phase(i,j,k-1)<-(pi/2.0d0))
                temp1 = phase(i,j,k+1)-(phase(i,j,k-1) - pi);
            elseif (phase(i,j,k+1)-phase(i,j,k-1)>(pi/2.0d0))
                temp1 = phase(i,j,k+1)-(phase(i,j,k-1) + pi);
            else
                temp1 = phase(i,j,k+1)-phase(i,j,k-1);
            end
            velz(i-1,j-1,k-1) = real(temp1)/dspace;
        end
        end
        end
        

end
