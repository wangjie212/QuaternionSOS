function standardterm(a::Vector{Vector{UInt16}},n)
    if length(a[3])<=1
        a[1]=sort(a[1])
        a[2]=sort(a[2])
        return a
    else
        Ind=[]
        i=1
        while i >= 1
            if a[3][i]==a[3][i+1]-UInt16(n)||a[3][i]==a[3][i+1]+UInt16(n)
                push!(a[1],a[3][i]<a[3][i+1] ? a[3][i] : a[3][i+1])
                push!(a[2],a[3][i]<a[3][i+1] ? a[3][i+1] : a[3][i])
                # push!(Ind,i,i+1)
                deleteat!(a[3],[i,i+1])
                a[1]=sort(a[1])
                a[2]=sort(a[2])
                if length(a[3])<=1
                    return a
                    break
                else
                    i=1
                end
                # if i+3<=length(a[3])
                #     i=i+2
                # else
                #     deleteat!(a[3],Ind)
                #     a[1]=sort(a[1])
                #     a[2]=sort(a[2])
                #     return a
                #     break
                # end
            else
                if i+2<=length(a[3])
                    i=i+1
                else
                    a[1]=sort(a[1])
                    a[2]=sort(a[2])
                    return a
                    break
                end
            end
        end
    end  
end

function qtermadd(a::Vector{Vector{UInt16}},b::Vector{Vector{UInt16}},n)
    return standardterm([append!(star(a,n)[1],b[1]),append!(star(a,n)[2],b[2]),append!(star(a,n)[3],b[3])],n)
end

function qtermadd3(a::Vector{Vector{UInt16}},b::Vector{Vector{UInt16}},c::Vector{Vector{UInt16}},n)
    return standardterm([append!(star(a,n)[1],b[1],c[1]),append!(star(a,n)[2],b[2],c[2]),append!(star(a,n)[3],b[3],c[3])],n)
end

function _qcyclic_canon(a::Vector{Vector{UInt16}},n)
    if length(a[3])<=1
        return a
    else
        a=[a[1],a[2],minimum([[a[3][i+1:length(a[3])];a[3][1:i]] for i=0:length(a[3])-1])]
        a=standardterm(a,n)
        return a
        # standardterm([a[1],a[2],minimum([[a[3][i+1:length(a[3])];a[3][1:i]] for i=0:length(a[3])-1])],n)
    end
end

function _qcyclic_canon(w::Monomial{false})
    ind = w.z .> 0
    wz = w.z[ind]
    wv = w.vars[ind]
    lw = length(wz)
    if lw == 0
        return w
    else
        return minimum([prod([wv[i+1:lw]; wv[1:i]] .^ [wz[i+1:lw]; wz[1:i]]) for i=0:lw-1])
    end
end

function qcyclic_canon(supp, coe, n; type=QuaternionF64)
    nsupp = [min(_qcyclic_canon(word,n), _qcyclic_canon([word[1],word[2],reverse(word[3])],n)) for word in supp]
    sort!(nsupp)
    unique!(nsupp)
    l = length(nsupp)
    ncoe = zeros(type, l)
    for (i,item) in enumerate(supp)
        bi = min(_qcyclic_canon(item,n), _qcyclic_canon([item[1],item[2],reverse(item[3])],n))
        Locb = bfind(nsupp, l, bi)
        ncoe[Locb] += coe[i]
    end
    return nsupp,ncoe
end

function randomsymfunc(q,n,d,rng;conjugates=false,coelimit=false)
    mon=Monomial{false}[1]
    for j=1:d
        if conjugates!=false
            append!(mon,monomials(q, j))
        else
            append!(mon,monomials(q[1:n], j))
        end
    end
    monc=Monomial{false}[]
    for i=1:length(mon)
        temp=prod(reverse(mon[i].vars).^reverse(mon[i].z))
        push!(monc,temp(q[1:n]=>q[n+1:2n],q[n+1:2n]=>q[1:n]))
        # push!(monc,temp(q[1]=>q[3],q[2]=>q[4],q[3]=>q[1],q[4]=>q[2]))
    end
    n=length(mon)
    if coelimit!=false
        A = 2 .* rand(rng,n, n) .- 1  # 生成范围在-1到1的随机矩阵
        A_symmetric = (A + A') / 2
    else
        A_symmetric=Symmetric(rand(rng,n,n))
    end
    return transpose(monc)*A_symmetric*mon
end

function qrandomsymfunc(q,n,d,rng;conjugates=false)
    mon=Monomial{false}[1]
    for j=1:d
        if conjugates!=false
            append!(mon,monomials(q, j))
        else
            append!(mon,monomials(q[1:n], j))
        end
    end
    monc=Monomial{false}[]
    for i=1:length(mon)
        temp=prod(reverse(mon[i].vars).^reverse(mon[i].z))
        push!(monc,temp(q[1:n]=>q[n+1:2n],q[n+1:2n]=>q[1:n]))
    end
    t=length(mon)
    println(t)
        A = rand(rng,QuaternionF64,t,t) 
        A_symmetric = (A + adjoint(A)) / 2
    # to real
    B = conj.(A_symmetric)
    return transpose(monc)*A_symmetric*mon, transpose(mon)*B*monc
end

function bfind(A, l, a)
    if l == 0
        return 0
    end
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if ndims(A) == 2
            temp = A[:, mid]
        else
            temp = A[mid]
        end
        if temp == a
            return mid
        elseif temp < a
            low = mid + 1
        else
            high = mid - 1
        end
    end
    return nothing
end


function ncbfind(A, l, a)
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if A[mid] == a
           return mid
        elseif A[mid] < a
            high = mid - 1
        else
            low = mid + 1
        end
    end
    return nothing
end

function star(a::Vector{Vector{UInt16}},n)
    if length(a[3])==0
        return a
    else
        a[3]=a[3][end:-1:1]
        for i=1:length(a[3])
            if a[3][i]<=UInt16(n)
                a[3][i]=a[3][i]+UInt16(n)
            else
                a[3][i]=a[3][i]-UInt16(n)
            end
        end
        return standardterm(a,n)
    end
end

function star(w::Monomial{false})
    return prod(reverse(w.vars).^reverse(w.z))
end

function star(p::Polynomial{false})
    return coefficients(p)'*star.(monomials(p))
end

function get_qncbasis(n, d; ind=Vector{UInt16}(1:2n), binary=false,conjubasis=false)
    basis=[[UInt16[],UInt16[],UInt16[]]]
    for i = 1:d
        if conjubasis!=false
        append!(basis, _get_qncbasis_deg2(n, i, ind=ind, binary=binary))
        else
        append!(basis, _get_qncbasis_deg(n, i, ind=ind, binary=binary))
        end
    end
    return basis
end

function bget_qncbasis(n, d; ind=Vector{UInt16}(1:2n), binary=false,conjubasis=false)
    basis=[[UInt16[],UInt16[],UInt16[]]]
    for i = 1:d
        if conjubasis!=false
        append!(basis, _get_qncbasis_deg2(n, i, ind=ind, binary=binary))
        else
        append!(basis, _get_qncbasis_deg(n, i, ind=ind, binary=binary))
        end
    end
    basistemp=deepcopy(basis)
    if !conjubasis
        for i=1:n
            for j=1:length(basistemp)
                a=deepcopy(basistemp[j])
                ltemp=qtermadd([UInt16[],UInt16[],[UInt16(i)]],a,n)
                push!(basis,ltemp)
            end
        end
    end
    unique!(basis)
    return basis
end

function _get_qncbasis_deg(n, d; ind=Vector{UInt16}(1:n), binary=false)
    if d > 0
        basis=[[UInt16[],UInt16[],UInt16[]]]
        for i = 1:n
            temp = _get_qncbasis_deg(n, d-1, ind=ind, binary=binary)
            if binary == false || d == 1
                lm=length(temp)
                for j =1:lm
                    push!(temp[j][3],ind[i])
                    append!(basis,temp)
                end
            else
                for item in temp
                    if item[3][end] != ind[i]
                        push!(basis, [item[1],item[2],[item[3];ind[i]]])
                    end
                end
            end
        end
        de=[]
        for i=1:length(basis)
            if length(basis[i][1])+length(basis[i][2])+length(basis[i][3])<d
                push!(de,i)
            end
        end
        deleteat!(basis,de)
        unique!(basis)
        return basis
    else
        return [[UInt16[],UInt16[],UInt16[]]]
    end
end

function _get_qncbasis_deg2(n, d; ind=Vector{UInt16}(1:2n), binary=false)
    if d > 0
        basis=[[UInt16[],UInt16[],UInt16[]]]
        for i = 1:2n
            temp = _get_qncbasis_deg2(n, d-1, ind=ind, binary=binary)
            if binary == false || d == 1
                lm=length(temp)
                for j =1:lm
                    if  length(temp[j][3])!=0
                        if temp[j][3][end]!= ind[i]-ind[n] && temp[j][3][end]!= ind[i]+ind[n]
                        push!(temp[j][3],ind[i])
                        else
                        ab=(temp[j][3][end]<ind[i]) ? temp[j][3][end] : ind[i]
                        push!(temp[j][1],ab)
                        push!(temp[j][2],ab+ind[n])
                        deleteat!(temp[j][3],length(temp[j][3]))
                        end
                    else
                    push!(temp[j][3],ind[i])
                    end
                temp[j][1]=sort(temp[j][1])
                temp[j][2]=sort(temp[j][2])
                end
            end
            append!(basis,temp)
        end       
        de=[]
        for i=1:length(basis)
            if length(basis[i][1])+length(basis[i][2])+length(basis[i][3])<d
                push!(de,i)
            end
        end
        deleteat!(basis,de)
        unique!(basis)
        return basis
    else
        return [[UInt16[],UInt16[],UInt16[]]]
    end
end

function qreduce_unitnorm(a; nb=0)
    a = [UInt16[],UInt16[],copy(a[3])]
    return a
end

function qresort(supp, coe; nb=0)
    if nb > 0
        supp = qreduce_unitnorm.(supp, nb=nb)
    end
    nsupp = deepcopy(supp)
    sort!(nsupp)
    unique!(nsupp)
    l = length(nsupp)
    ncoe = zeros(typeof(coe[1]), l)
    for i in eachindex(supp)
        locb = bfind(nsupp, l, supp[i])
        ncoe[locb] += coe[i]
    end
    return nsupp,ncoe
end
