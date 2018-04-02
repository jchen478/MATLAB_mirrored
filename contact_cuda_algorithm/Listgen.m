function  fiberstatus=Listgen(fiber,minthread,fiberarrays,fiberstatus)
%%min thread is is the identity of the smallest thread for which there is overlapping work    
callfiber=fiber;
N = size(fiberarrays(fiber).list,2);
    for n=1:N
        fiber=fiberarrays(callfiber).list(n);
        minthread=min(fiber,minthread);
        
        if fiber~=minthread && fiberstatus(fiber)==0
            fiberstatus(fiber)=2*((1+sign(fiber-minthread))/2-1/2)
            fiberstatus=Listgen(fiber,minthread,fiberarrays,fiberstatus);
        else
            %fiber
        end
        if callfiber~=minthread
            fiberstatus(callfiber)=1;
        end
    end
end