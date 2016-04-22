module Procrustes

#define the fundimental type
type Coords
    rawCoords::Array
    alignCoords::Array
    centroidSize::Vector
end

#main function
function procAlign(data::Array,scale=true,tol=10e-5)
    
    centSize=centriodSize(data)
    rawCoords=copy(data)
    procAlign!(data,scale,tol)
    return(Coords(rawCoords,data,centSize));

end

function procAlign!(x::Array, scale, tol)

    centerShape!(x)
    
    if scale
        scaleShape!(x)
    end

    #pick random sample for frist iterative alignment
    ind=rand(1:(size(x)[3]));
    meanShapeIn=x[:,:,ind];
    meanShapeOut=zeros(size(meanShapeIn));
    procError=100;

    #iterate until meanshape converges
    while procError > tol
    
        for sample in 1:size(x)[3]
            x[:,:,sample]=align_two_shapes(meanShapeIn[:,:,1],x[:,:,sample])
        end

        meanShapeOut=meanShape(x)
        procError=procDist(meanShapeIn,meanShapeOut)

        meanShapeIn=meanShapeOut
    end
  
end

function align_two_shapes(refShape::Array, movingShape::Array)
    H=transpose(refShape)*movingShape
    
    U,S,V=svd(H)
    R=V*transpose(U)
    
    return(movingShape*R)
end

function centerShape!(x::Array)
    #for each sample
    for k in 1:size(x)[3];
        for j in 1:size(x)[2];
            
            colMean=mean(x[:,j,k])
            for i in 1:size(x)[1]
                x[i,j,k]=x[i,j,k]-colMean
            end     
        end
    end
end


function scaleShape!(x::Array)
    for k in 1:size(x)[3]
        x[:,:,k]=x[:,:,k]/sqrt(sum(var(x[:,:,k], 1)))
    end  
end

function meanShape(x::Array)
    return(mean(x,3))
end



function procDist(x::Array, y::Array)
    return(sqrt(sum((x-y).^2)))
end


function centriodSize(x::Array)
    
    centSize=zeros(size(x)[3])
    
    #for each sample
    for k in 1:size(x)[3];
        
        colMean=mean(x[:,:,k],1)
        for i in 1:size(x)[1]
            centSize[k]+=(sum((x[i,:,k]-colMean).^2))
        end
    end
    
    map!(sqrt,centSize)
    return(centSize)
end

end # module
