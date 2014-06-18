#miscellenious helper functions and constants

import Base.@math_const


Q95 = sqrt(2.)*erfinv(0.95)





range(x) = (minimum(x), maximum(x))



#95% normal quantile





eps2 = sqrt(eps())









intervalgaussian(U, a, b) = broadcast(intervalgaussian, U, a, b)

function intervalgaussian (U::Real,a::Real, b::Real)

    a, b = min(a,b), max(a,b)

    c =  erf(a/sqrt(2.))   

    d =  erf(b/sqrt(2.))  

    sqrt(2.)*erfinv(c + U.*(d-c))  

end







function scalar(x)

assert(length(x) == 1)

x[1]

end



function norma!(x)

 minx = min(x)

 x = (x - minx) / (max(x)-minx)

 x

end





function cut(x, a, b)

    warn("cut(): better use clamp()")

    clamp(x, a,b)

end

#compute normal confidence interval for Monte Carlo estimate, precision aware rounding



function mc(M)

  m = mean(M)

  ste = std(M)/sqrt(length(M))

  round([m, Q95*ste], int(2-log(10,Q95*ste)))

end



#for sequential estimates, it suffices to provide the running sum Z, running sum of squares Z2 and the number of observations k





va(k, Z, Z2) =  Z2/k-(Z/k)^2
ste(k, Z, Z2) = sqrt(Z2/k^2-(Z/k)^2/k)



function mc2(k, Z::Real, Z2::Real, mode = :ci)
    m = Z/k
    if mode ==: ci 
        F = Q95 else F = 1. 
    end

    try
        ste = sqrt(Z2/k-m.^2)/sqrt(k)
        res = (m, F*ste)
        return res
    catch
        return (m, NaN.*m)
    end
end


function mc2(k, Z::Vector, Z2::Vector, mode = :ci)
   map((z,z2) -> mc2(k, z, z2, mode), Z, Z2)
end
function mc2str(k, Z, Z2)
    m, ste = mc2(k, Z, Z2, :se)
    m, ci = round([m, Q95*ste], int(2-log(10,Q95*ste)))
    "$m+-$(ci)ci "
end
function mc2str(k, Z::Vector, Z2::Vector)
   reduce(*,map((z,z2) -> mc2str(k, z, z2), Z, Z2))
end




function isignif_og(x, digits, base)

    if base == 10

        ifloor(log10(abs(x)) - digits + 1.)

    elseif base == 2

        ifloor(log2(abs(x)) - digits + 1.)

    else

        ifloor(log2(abs(x))/log2(base) - digits + 1.)

    end

end





function roundste(out::IOBuffer, m, ste::FloatingPoint)

    if 

        !isfinite(m) print(out, m)

        return out

    end



    if isnan(ste) 

        print(out, m, " ± NaN (se)")

        return out

    elseif isinf(ste)

        print(out, m, " ± Inf (se)")

        return out

    elseif ste < eps(m)

        print(out, m)

        return out

    end



    assert(ste >= 0.) 

    og = max(isignif_og(m, 3, 10) - isignif_og(ste, 3, 10)+3,0)

    grisu(out, m, Base.Grisu.PRECISION, og)

    print(out, " ± ")

    grisu(out, ste, Base.Grisu.PRECISION, 3)

    print(out, " (se)")

out

end



function roundste(out::IOBuffer, r)

    m, ste = r

    roundste(out, m[1], ste[1]) 

    for i in 2:length(m)

        print(", ")

        roundste(out, m[i], ste[i]) 

    end

end





function strmc2(k, Z, Z2) 

    s = IOBuffer(true, true)

    truncate(s,0)

    roundste(s, mc2(k, Z, Z2, false))

    takebuf_string(s)

end





#stdv(k, Z, Z2) = va(k, Z, Z2)/k

# X = [ZL, L]

# X2 = X * X'   



function selfadjmc(k, X, X2)

    m = X[1]/X[2]

    v = (va(k, X[1], X2[1,1]) - 2m*(X2[1,2]/k - X[1]*X[2]/k^2) + m^2*va(k, X[2], X2[2,2]))/k

#   println([va(k, X[1], X2[1,1]), - 2m*(X2[1,2]/k - X[1]*X[2]/k^2),m^2*va(k, X[2], X2[2,2])])

    try

        ste = sqrt(v)

        res = [m, ste]

        return res

    catch

        return [m, NaN.*m]

    end

    

end





function grisu(io::IO, x::FloatingPoint, mode::Int32, n::Int)

    if isnan(x) return write(io, "NaN"); end

    if isinf(x)

        if x < 0 write(io,'-') end

        write(io, "Inf")

        return

    end

    Base.Grisu.@grisu_ccall(x, mode, n)

    pdigits = pointer(Base.Grisu.DIGITS)

    neg = Base.Grisu.NEG[1]

    len = int(Base.Grisu.LEN[1])

    pt  = int(Base.Grisu.POINT[1])

#    if mode == PRECISION

#        while len > 1 && DIGITS[len] == '0'

#            len -= 1

#        end

#    end

    if neg write(io,'-') end

    if pt <= -4 || pt > 6 # .00001 to 100000.

        # => #.#######e###

        write(io, pdigits, 1)

        write(io, '.')

        if len > 1

            write(io, pdigits+1, len-1)

        else

            write(io, '0')

        end

        write(io, 'e')

        write(io, dec(pt-1))

        return

    elseif pt <= 0

        # => 0.00########

        write(io, "0.")

        while pt < 0

            write(io, '0')

            pt += 1

        end

        write(io, pdigits, len)

    elseif pt >= len

        # => ########00.0

        write(io, pdigits, len)

        while pt > len

            write(io, '0')

            len += 1

        end

        write(io, ".0")

    else # => ####.####

        write(io, pdigits, pt)

        write(io, '.')

        write(io, pdigits+pt, len-pt)

    end

    nothing

end



#plot (2d) sample path



function pl(x::Array{ FloatingPoint,2})

    p = FramedPlot()

    #compute range

    R1 = range(x[1,:]) 

    R2 = range(x[2,:])



    #widen a bit    

    R1 = (R1[1] -0.1*(R1[2]-R1[1]), R1[2] + 0.1*(R1[2]-R1[1]))

    R2 = (R2[1] -0.1*(R2[2]-R2[1]), R2[2] + 0.1*(R2[2]-R2[1]))



    setattr(p, "xrange", R1)

    setattr(p, "yrange", R2)



    add(p, Curve(x[1,:],x[2,:]))

    Winston.display(p)



    p   

end
