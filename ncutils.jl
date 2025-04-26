function get_basis(n, d)
    lb = binomial(n+d, d)
    basis = zeros(UInt8, n, lb)
    i = 0
    t = 1
    while i < d+1
        t += 1
        if basis[n,t-1] == i
           if i < d
              basis[1,t] = i+1
           end
           i += 1
        else
            j = findfirst(x->basis[x,t-1]!=0, 1:n)
            basis[:,t] = basis[:,t-1]
            if j == 1
               basis[1,t] -= 1
               basis[2,t] += 1
            else
               basis[1,t] = basis[j,t] - 1
               basis[j,t] = 0
               basis[j+1,t] += 1
            end
        end
    end
    return basis
end

function _cyclic_canon(a::Vector{UInt16})
    if isempty(a)
        return a
    else
        return minimum([[a[i+1:length(a)]; a[1:i]] for i=0:length(a)-1])
    end
end

function _cyclic_canon(w::Monomial{false})
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

function cyclic_canon(supp, coe; type=Float64)
    nsupp = [min(_cyclic_canon(word), _cyclic_canon(reverse(word))) for word in supp]
    sort!(nsupp)
    unique!(nsupp)
    l = length(nsupp)
    ncoe = zeros(type, l)
    for (i,item) in enumerate(supp)
        bi = min(_cyclic_canon(item), _cyclic_canon(reverse(item)))
        Locb = bfind(nsupp, l, bi)
        ncoe[Locb] += coe[i]
    end
    return nsupp,ncoe
end

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
# function _qcyclic_canon(a::Vector{Vector{UInt16}})
#     if isempty(a)
#         return a
#     else
#         return [a[1],a[2],minimum([[a[3][i+1:length(a[3])];a[3][1:i]] for i=0:length(a[3])-1])]
#     end
# end

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

function _sym_canon(a::Vector{UInt16})
    i = 1
    while i <= Int(ceil((length(a)-1)/2))
        if a[i] < a[end+1-i]
            return a
        elseif a[i] > a[end+1-i]
            return reverse(a)
        else
            i += 1
        end
    end
    return a
end

function _sym_canon(w::Monomial{false})
    return min(w, star(w))
end

function is_sym(a::Vector{UInt16})
    l = Int(ceil((length(a)-1)/2))
    return isequal(a[1:l], a[end:-1:end-l+1])
end

function sym_canon(supp, coe; type=Float64)
    nsupp = copy(supp)
    nsupp = _sym_canon.(nsupp)
    sort!(nsupp)
    unique!(nsupp)
    l = length(nsupp)
    ncoe = zeros(type, l)
    for (i,item) in enumerate(supp)
        Locb = bfind(nsupp, l, _sym_canon(item))
        ncoe[Locb] += coe[i]
    end
    return nsupp,ncoe
end

function get_ncbasis(n, d; ind=Vector{UInt16}(1:n), binary=false)
    basis = [UInt16[]]
    for i = 1:d
        append!(basis, _get_ncbasis_deg(n, i, ind=ind, binary=binary))
    end
    return basis
end
function _get_ncbasis_deg(n, d; ind=Vector{UInt16}(1:n), binary=false)
    if d > 0
        basis = Vector{UInt16}[]
        for i = 1:n
            temp = _get_ncbasis_deg(n, d-1, ind=ind, binary=binary)
            if binary == false || d == 1
                push!.(temp, ind[i])
                # println(temp)
                append!(basis, temp)
                # println(basis)
            else
                for item in temp
                    # printlln(item[end])
                    if item[end] != ind[i]
                        push!(basis, [item;ind[i]])
                        # println(basis)
                    end
                end
            end
        end
        return basis
    else
        return [UInt16[]]
    end
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
# function _get_qncbasis_deg(n, d; ind=Vector{UInt16}(1:2n), binary=false)
#     if d > 0
#         # basis = Vector{Vector{Vector{UInt16}}}[[[],[],[]]]
#         basis=[[UInt16[],UInt16[],UInt16[]]]
#         for i = 1:2n
#             temp = _get_qncbasis_deg(n, d-1, ind=ind, binary=binary)
#             # println(temp)
#             # println(temp,length(temp))
#             if binary == false || d == 1
#                 lm=length(temp)
#                 for j =1:lm
#                     if length(temp[j][3])!=0
#                         if temp[j][3][end]!= ind[i]-ind[n] && temp[j][3][end]!= ind[i]+ind[n]
#                         push!(temp[j][3],ind[i])
#                         else
#                         ab=(temp[j][3][end]<ind[i]) ? temp[j][3][end] : ind[i]
#                         push!(temp[j][1],ab)
#                         push!(temp[j][2],ab+ind[n])
#                         deleteat!(temp[j][3],length(temp[j][3]))
#                         end
#                     else
#                         # println(ind[i])
#                         # println(typeof(temp[j][3]))
#                         push!(temp[j][3],ind[i])
#                         # println(temp[j][3])
#                     end
#                 end
#             end
#                 # println(temp)
#                 append!(basis, temp)
#             # else
#             #     for item in temp
#             #         # printlln(item[end])
#             #         if item[end] != ind[i]
#             #             push!(basis, [item;ind[i]])
#             #             # println(basis)
#             #         end
#             #     end
#             # end
#         end
#         de=[]
#         # for i=1:length(basis)
#         #     temp3=basis[i][3][basis[i][3].<=ind[n]]
#         #     if length(basis[i][1])+length(temp3)<d && length(basis[i][2])+length(basis[i][3])-length(temp3)<d
#         #         push!(de,i)
#         #     end
#         # end
#         for i=1:length(basis)
#             if length(basis[i][1])+length(basis[i][2])+length(basis[i][3])<d
#                 push!(de,i)
#             end
#         end
#         deleteat!(basis,de)
#         unique!(basis)
#         return basis
#     else
#         return [[UInt16[],UInt16[],UInt16[]]]
#     end
# end
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
# function _get_qncbasis_deg2(n, d; ind=Vector{UInt16}(1:2n), binary=false)
#     if d > 0
#         # basis = Vector{Vector{Vector{UInt16}}}[[[],[],[]]]
#         basis=[[UInt16[],UInt16[],UInt16[]]]
#         Temp = _get_qncbasis_deg2(n, d-1, ind=ind, binary=binary)
#         # println(Temp)
#         if binary == false || d == 1
#             lm=length(Temp)
#             for j =1:lm
#                 temp=deepcopy(Temp)
#                 # println(temp)
#                 temp3=temp[j][3][temp[j][3].<=ind[n]]
#                 # println(temp3)
#                 lengthq = length(temp[j][1])+length(temp3)
#                 lengthqbar = length(temp[j][2])+length(temp[j][3])-length(temp3)
#                 # println(lengthq)
#                 # Ind=Vector{UInt16}()
#                 for k = 0:d-lengthq                 
#                     for s = 0:d-lengthqbar
#                         Ind=Vector{UInt16}()
#                         if k!=0
#                             for kk =1:k
#                             append!(Ind,Vector{UInt16}(1:n))
#                             end
#                         end   
#                         if s!=0
#                             for ss=1:s
#                                 append!(Ind,Vector{UInt16}(n+1:2n))
#                             end
#                         end
#                         # println(Ind)
#                         if length(Ind)!=0
#                             TTemp=[[[UInt16[],UInt16[],UInt16[]]]]
#                             TTemp[1][1]=deepcopy(temp[j])
#                             for i = 1:s+k
#                                 Indtemp=Ind[1+(i-1)*n:i*n]
#                                 ttemp=[[UInt16[],UInt16[],UInt16[]]]
#                                 for t = 1:length(Indtemp)
#                                     for ii= 1:length(TTemp[i])
#                                         # println(TTemp)
#                                         sTemp=deepcopy(TTemp[i][ii])
#                                         if  length(sTemp[3])!=0
#                                             if sTemp[3][end]!= Indtemp[t]-ind[n] && sTemp[3][end]!= Indtemp[t]+ind[n]
#                                             push!(sTemp[3],Indtemp[t])
#                                             else
#                                             ab=(sTemp[3][end]<Indtemp[t]) ? sTemp[3][end] : Indtemp[t]
#                                             push!(sTemp[1],ab)
#                                             push!(sTemp[2],ab+ind[n])
#                                             deleteat!(sTemp[3],length(sTemp[3]))
#                                             end
#                                         else
#                                         # println(sTemp)
#                                         # println(TTemp)
#                                         push!(sTemp[3],Indtemp[t])
#                                         # println(a,b,c)
#                                         # println(sTemp)
#                                         # println(TTemp)
#                                         end
#                                         # println(ttemp,TTemp[i][ii])
#                                         sTemp[1]=sort(sTemp[1])
#                                         sTemp[2]=sort(sTemp[2])
#                                         push!(ttemp,sTemp)
#                                         # deleteat!(ttemp,1)
#                                         # println(ttemp)
#                                         # println(TTemp)
#                                         # append!(TTemp[i+1],TTemp[i][ii])
#                                     end
#                                 end
#                                 # println(TTemp,ttemp)
#                                 push!(TTemp,ttemp)
#                                 # println(TTemp)
#                             end
#                             append!(basis, TTemp[1+s+k])
#                             # println(basis)
#                         end
#                     end                    
#                 end
#                 # for k = 1:d-lengthq
#                 #     append!(Ind,Vector{UInt16}(n+1:2n))
#                 # end
#                 # for i = 1:length(Ind)
#                 #     if  length(temp[j][3])!=0
#                 #         if temp[j][3][end]!= Ind[i]-ind[n] && temp[j][3][end]!= Ind[i]+ind[n]
#                 #         push!(temp[j][3],Ind[i])
#                 #         else
#                 #         ab=(temp[j][3][end]<Ind[i]) ? temp[j][3][end] : Ind[i]
#                 #         push!(temp[j][1],ab)
#                 #         push!(temp[j][2],ab+Ind[n])
#                 #         deleteat!(temp[j][3],length(temp[j][3]))
#                 #         end
#                 #     else
#                 #     # println(ind[i])
#                 #     # println(typeof(temp[j][3]))
#                 #     push!(temp[j][3],Ind[i])
#                 #     # println(temp[j][3])
#                 #     end
#                 #     println(temp)
#                 #     append!(basis, temp)
#                 #     println(basis)
#                 #   end
#             end
#         end
#             # else
#             #     for item in temp
#             #         # printlln(item[end])
#             #         if item[end] != ind[i]
#             #             push!(basis, [item;ind[i]])
#             #             # println(basis)
#             #         end
#             #     end
#             # end
#         de=[]
#         for i=1:length(basis)
#             temp3=basis[i][3][basis[i][3].<=ind[n]]
#             if length(basis[i][1])+length(temp3)<d && length(basis[i][2])+length(basis[i][3])-length(temp3)<d
#                 push!(de,i)
#             end
#         end
#         # println(de)
#         # println(basis)
#         deleteat!(basis,de)
#         unique!(basis)
#         return basis
#     else
#         return [[UInt16[],UInt16[],UInt16[]]]
#     end
# end
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
function reduce_cons!(word::Vector{UInt16}; constraint="unipotent")
    i = 1
    while i < length(word)
        if word[i] == word[i+1]
            deleteat!(word, i)
            if constraint == "unipotent"
                deleteat!(word, i)
            end
        else
            i += 1
        end
    end
    return word
end

function reduce_cons(w::Monomial{false}; constraint="unipotent")
    ind = w.z .> 1
    if constraint == "unipotent"
        w.z[ind] .= 0
    else
        w.z[ind] .= 1
    end
    return prod(w.vars .^ w.z)
end

function reduce!(word::Vector{UInt16}; obj="eigen", partition=0, constraint=nothing)
    if obj == "trace"
        word = min(_cyclic_canon(word), _cyclic_canon(reverse(word)))
    else
        if partition > 0 && constraint === nothing
            word = min(_comm(word, partition), _comm(reverse(word), partition))
        elseif partition == 0 && constraint !== nothing
            cword = copy(word)
            word = min(reduce_cons!(word, constraint = constraint), reduce_cons!(reverse(cword), constraint = constraint))
        elseif partition > 0 && constraint !== nothing
            word = min(reduce_cons!(_comm(word, partition), constraint = constraint), reduce_cons!(_comm(reverse(word), partition), constraint = constraint))
        else
            word = _sym_canon(word)
        end
    end
    return word
end

function reduce(word::Monomial{false}, x; obj="eigen", partition=0, constraint=nothing)
    if obj == "trace"
        word = min(_cyclic_canon(word), _cyclic_canon(star(word)))
    else
        if partition > 0 && constraint === nothing
            word = min(_comm(word, x, partition), _comm(star(word), x, partition))
        elseif partition == 0 && constraint !== nothing
            word = min(reduce_cons(word, constraint = constraint), reduce_cons(star(word), constraint = constraint))
        elseif partition > 0 && constraint !== nothing
            word = min(reduce_cons(_comm(word, x, partition), constraint = constraint), reduce_cons(_comm(star(word), x, partition), constraint = constraint))
        else
            word = _sym_canon(word)
        end
    end
    return word
end

function _comm(word::Vector{UInt16}, partition)
    ind1 = word .<= partition
    ind2 = word .> partition
    return [word[ind1]; word[ind2]]
end

function _comm(w::Monomial{false}, x, partition)
    ind1 = w.vars .>= x[partition]
    ind2 = w.vars .< x[partition]
    return prod([w.vars[ind1]; w.vars[ind2]] .^ [w.z[ind1]; w.z[ind2]])
end

function bfind(A, l, a; lt=isless, rev=false)
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if isequal(A[mid], a)
           return mid
        elseif lt(A[mid], a)
            if rev == false
                low = mid+1
            else
                high = mid-1
            end
        else
            if rev == false
                high = mid-1
            else
                low = mid+1
            end
        end
    end
    return nothing
end

function permutation(a)
    b = sparse(a)
    ua = convert(Vector{UInt16}, b.nzind)
    na = convert(Vector{UInt16}, b.nzval)
    return _permutation(ua, na)
end

function _permutation(ua, na)
    if !isempty(ua)
        perm = Vector{UInt16}[]
        for i = 1:length(ua)
            nua = copy(ua)
            nna = copy(na)
            if na[i] == 1
                deleteat!(nua, i)
                deleteat!(nna, i)
            else
                nna[i] -= 1
            end
            temp = _permutation(nua, nna)
            push!.(temp, ua[i])
            append!(perm, temp)
        end
        return perm
    else
        return [UInt16[]]
    end
end

function polys_info(pop, x)
    n = length(x)
    m = length(pop)-1
    coe = Vector{Vector{Float64}}(undef, m+1)
    supp = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    for k = 1:m+1
        mon = monomials(pop[k])
        coe[k] = coefficients(pop[k])
        supp[k] = [UInt16[] for i=1:length(mon)]
        for i = 1:length(mon)
            ind = mon[i].z .> 0
            vars = mon[i].vars[ind]
            exp = mon[i].z[ind]
            for j = 1:length(vars)
                l = bfind(x, n, vars[j], rev=true)
                append!(supp[k][i], l*ones(UInt16, exp[j]))
            end
        end
    end
    return n,supp,coe
end

function poly_info(f, x)
    n = length(x)
    mon = monomials(f)
    coe = coefficients(f)
    lm = length(mon)
    supp = [UInt16[] for i=1:lm]
    for i = 1:lm
        ind = mon[i].z .> 0
        vars = mon[i].vars[ind]
        exp = mon[i].z[ind]
        for j = 1:length(vars)
            k = bfind(x, n, vars[j], rev=true)
            append!(supp[i], k*ones(UInt16, exp[j]))
        end
    end
    return n,supp,coe
end

function isless_td(a, b)
    if length(a) < length(b)
        return true
    elseif length(a) > length(b)
        return false
    else
        return a < b
    end
end

function sym_cyclic(word)
    return min(_cyclic_canon(word), _cyclic_canon(reverse(word)))
end

function _cyclic_basis(var)
    basis = _permutation(var, ones(length(var)))
    return unique(_cyclic_canon.(basis))
end

function sym(word, vargroup)
    cword = copy(word)
    ind = [gind(word[i], vargroup) for i = 1:length(word)]
    uind = unique(ind)
    nind = [count(ind .== uind[i]) for i = 1:length(uind)]
    k = 0
    for i = 1:length(uind)
        cword[k+1:k+nind[i]] = reverse(cword[k+1:k+nind[i]])
        k += nind[i]
    end
    return min(word, cword)
end

function iscomm(a, vargroup)
    for i = 1:length(a)-1
        if a[i] > a[i+1] && gind(a[i], vargroup) != gind(a[i+1], vargroup)
            return false
        end
    end
    return true
end

function gind(k, vargroup)
    return findfirst(i -> k <= sum(vargroup[1:i]), 1:length(vargroup))
end

function res_comm!(a, vargroup)
    i = 1
    while i < length(a)
        if a[i] > a[i+1] && gind(a[i], vargroup) != gind(a[i+1], vargroup)
            temp = a[i]
            a[i] = a[i+1]
            a[i+1] = temp
            if i > 1
                i -= 1
            end
        else
            i += 1
        end
    end
    return a
end

function issym(word, vargroup)
    ind = [gind(word[i], vargroup) for i = 1:length(word)]
    uind = unique(ind)
    nind = [count(ind .== uind[i]) for i = 1:length(uind)]
    k = 0
    for i = 1:length(uind)
        temp = word[k+1:k+nind[i]]
        if reverse(temp) != temp
            return false
        end
        k += nind[i]
    end
    return true
end

function star(w::Monomial{false})
    return prod(reverse(w.vars).^reverse(w.z))
end

function star(p::Polynomial{false})
    return coefficients(p)'*star.(monomials(p))
end

# generate an SOHS polynomial with variables vars and degree 2d
function add_SOHS!(model, vars, d; obj="eigen", partition=0, constraint=nothing)
    basis = vcat([MultivariatePolynomials.monomials(vars, i) for i = 0:d]...)
    if constraint !== nothing
        basis = basis[[all(item.z .< 2) for item in basis]]
    end
    if partition > 0
        ind = Int[]
        for (i,item) in enumerate(basis)
            vs = item.vars[item.z .> 0]
            if findfirst(j -> vs[j] < vars[partition] && vs[j+1] >= vars[partition], 1:length(vs)-1) === nothing
                push!(ind, i)
            end
        end
        basis = basis[ind]
    end
    sohs = 0
    pos = @variable(model, [1:length(basis), 1:length(basis)], PSD)
    for j = 1:length(basis), k = j:length(basis)
        word = reduce(star(basis[j])*basis[k], vars, obj=obj, partition=partition, constraint=constraint)
        if j == k
            @inbounds sohs += pos[j,k]*word
        else
            @inbounds sohs += 2*pos[j,k]*word
        end
    end
    return sohs
end

# generate a polynomial with variables vars and degree d
function add_poly!(model, vars, d; obj="eigen", partition=0, constraint=nothing)
    basis = vcat([MultivariatePolynomials.monomials(vars, i) for i = 0:d]...)
    if constraint !== nothing
        basis = basis[[all(item.z .< 2) for item in basis]]
    end
    if partition > 0
        ind = Int[]
        for (i,item) in enumerate(basis)
            vs = item.vars[item.z .> 0]
            if findfirst(j -> vs[j] < vars[partition] && vs[j+1] >= vars[partition], 1:length(vs)-1) === nothing
                push!(ind, i)
            end
        end
        basis = basis[ind]
    end
    free = @variable(model, [1:length(basis)])
    poly = free'*basis
    return poly
end

function arrange(p, vars; obj="eigen", partition=0, constraint=nothing)
    mons = monomials(p)
    coe = coefficients(p)
    mons = [reduce(mon, vars, obj=obj, partition=partition, constraint=constraint) for mon in mons]
    nmons = unique(sort(mons))
    ncoe = zeros(typeof(coe[1]), length(nmons))
    for (i,item) in enumerate(coe)
        Locb = bfind(nmons, length(nmons), mons[i])
        ncoe[Locb] += coe[i]
    end
    return nmons,ncoe
end
