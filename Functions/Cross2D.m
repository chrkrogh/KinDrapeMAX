function result = Cross2D(a,b)
    if length(a) > 2 || length(b) > 2
        error('Vectors must have length 2')
    end
   result = a(1)*b(2) - a(2)*b(1); 
end