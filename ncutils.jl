function standardterm(b::Vector{Vector{UInt16}}, n)
    a = deepcopy(b)
    if length(a[3]) <= 1
        a[1] = sort(a[1])
        a[2] = sort(a[2])
        return a
    else
        Ind = []
        i = 1
        while i >= 1
            if a[3][i] == a[3][i+1] - UInt16(n) || a[3][i] == a[3][i+1] + UInt16(n)
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

function qtermaddleft(a::Vector{Vector{UInt16}},b::Vector{Vector{UInt16}},n)
    return standardterm([append!(star(a,n)[1],b[1]),append!(star(a,n)[2],b[2]),append!(star(a,n)[3],b[3])],n)
end
function qtermadd(a::Vector{Vector{UInt16}},b::Vector{Vector{UInt16}},n)
    return standardterm([append!(a[1],star(b,n)[1]),append!(a[2],star(b,n)[2]),append!(a[3],star(b,n)[3])],n)
end

#= function qtermadd3(a::Vector{Vector{UInt16}},b::Vector{Vector{UInt16}},c::Vector{Vector{UInt16}},n)
    return standardterm([append!(star(a,n)[1],b[1],c[1]),append!(star(a,n)[2],b[2],c[2]),append!(star(a,n)[3],b[3],c[3])],n)
end =#
function qtermadd3(a::Vector{Vector{UInt16}},b::Vector{Vector{UInt16}},c::Vector{Vector{UInt16}},n)
    return standardterm([append!(a[1],b[1],star(c,n)[1]),append!(a[2],b[2],star(c,n)[2]),append!(a[3],b[3],star(c,n)[3])],n)
end
function qtermadd4(a::Vector{Vector{UInt16}},b::Vector{Vector{UInt16}},c::Vector{Vector{UInt16}},d::Vector{Vector{UInt16}},n)
    return standardterm([append!(a[1],b[1],c[1],d[1]),append!(a[2],b[2],c[2],d[2]),append!(a[3],b[3],c[3],d[3])],n)
end

function qtermadd3left(a1::Vector{Vector{UInt16}},b1::Vector{Vector{UInt16}},c1::Vector{Vector{UInt16}},n)
    a = deepcopy(a1)
    b = deepcopy(b1)
    c = deepcopy(c1)
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

function _qcyclic_canon(w)
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

function mono_to_term(m, q, n)
    # if m == 1
    #     return standardterm([UInt16[], UInt16[], UInt16[]], n)
    # end
    varmap = Dict(q[i] => UInt16(i) for i in 1:length(q))
    t = [UInt16[], UInt16[], UInt16[]]
    for (v, z) in zip(m.vars, m.z)
        idx = varmap[v]
        append!(t[3], fill(idx, z))
    end
    return standardterm(t, n)
end

function randomsymfunc(q,n,d;conjugates=false,coelimit=false)
    mon = NCMono[1]
    for j=1:d
        if conjugates!=false
            append!(mon,monomials(q, j))
        else
            append!(mon,monomials(q[1:n], j))
        end
    end
    monc = NCMono[]
    for i=1:length(mon)
        temp=prod(reverse(mon[i].vars).^reverse(mon[i].z))
        push!(monc,temp(q[1:n]=>q[n+1:2n],q[n+1:2n]=>q[1:n]))
        # push!(monc,temp(q[1]=>q[3],q[2]=>q[4],q[3]=>q[1],q[4]=>q[2]))
    end
    r=length(mon)
    if coelimit!=false
        A = 2 .* rand(r, r) .- 1  # 生成范围在-1到1的随机矩阵
        A_symmetric = (A + A') / 2
    else
        A_symmetric=Symmetric(rand(r,r))
    end
    #return transpose(monc)*A_symmetric*mon
    #new
    f = transpose(monc)*A_symmetric*mon
    # 把 mon 转成内部 term
    monterm = [mono_to_term(mon[i], q, n) for i in 1:r]

    fsupp = Vector{Vector{UInt16}}[]
    fcoe = QuaternionF64[]

    # 目标 support 按 w_j * w_i^* 的顺序生成
    for i = 1:r, j = 1:r
        if A_symmetric[i, j]!=0
            a = deepcopy(monterm[i])
            b = deepcopy(monterm[j])
            bi = qtermadd(b, a, n)   # w_j * w_i^*
            push!(fsupp, bi)
            push!(fcoe, A_symmetric[i, j])
        end
    end

    return f, fsupp, fcoe
end
function randomsymfunc2(q,n,d;conjugates=false,coelimit=false)
    mon = NCMono[1]
    for j=1:d
        if conjugates!=false
            append!(mon,monomials(q, j))
        else
            append!(mon,monomials(q[1:n], j))
        end
    end
    monc = NCMono[]
    for i=1:length(mon)
        temp=prod(reverse(mon[i].vars).^reverse(mon[i].z))
        push!(monc,temp(q[1:n]=>q[n+1:2n],q[n+1:2n]=>q[1:n]))
        # push!(monc,temp(q[1]=>q[3],q[2]=>q[4],q[3]=>q[1],q[4]=>q[2]))
    end
    r=length(mon)
    # if coelimit!=false
    #     A = 2 .* rand(r, r) .- 1  # 生成范围在-1到1的随机矩阵
    #     A_symmetric = (A + A') / 2
    # else
    #     A_symmetric=Symmetric(rand(r,r))
    # end
    A_symmetric = sparse_symmetric_pm1(r, density=0.05, zero_diag=true)
    #return transpose(monc)*A_symmetric*mon
    #new
    f = transpose(monc)*A_symmetric*mon
    # 把 mon 转成内部 term
    monterm = [mono_to_term(mon[i], q, n) for i in 1:r]

    fsupp = Vector{Vector{UInt16}}[]
    fcoe = QuaternionF64[]

    # 目标 support 按 w_j * w_i^* 的顺序生成
    for i = 1:r, j = 1:r
        if A_symmetric[i, j]!=0
            a = deepcopy(monterm[i])
            b = deepcopy(monterm[j])
            bi = qtermadd(b, a, n)   # w_j * w_i^*
            push!(fsupp, bi)
            push!(fcoe, A_symmetric[i, j])
        end
    end

    return f, fsupp, fcoe
end

function cliques_randomsymfunc(q,n,cn,size,d;conjugates=false,coelimit=false)
    f = 0*q[1]^0
    g = DynamicPolynomials.Polynomial{DynamicPolynomials.NonCommutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, Int64}[]
    # s = max(1, (n - size) ÷ (cn - 1)) 
    fsupp = Vector{Vector{UInt16}}[]
    fcoe = QuaternionF64[]
    for k = 1:cn
        mon = NCMono[1]
        for j=1:d
            if conjugates!=false
                append!(mon,monomials(vcat(q[1 + (k - 1) * (size - 1) : min(1 + k * (size - 1), n)],q[1 + (k - 1) * (size - 1)+n : min(1 + k * (size - 1), n)+n]), j))
            else
                # append!(mon,monomials(q[1 + (k-1)*s : min(1 + (k-1)*s + size - 1, n)], j))
                append!(mon,monomials(q[1 + (k - 1) * (size - 1) : min(1 + k * (size - 1), n)], j))
            end
        end
        monc = NCMono[]
        for i = 1:length(mon)
            temp = prod(reverse(mon[i].vars).^reverse(mon[i].z))
            push!(monc,temp(q[1:n]=>q[n+1:2n],q[n+1:2n]=>q[1:n]))
        end
        n_mon=length(mon)
        if coelimit!=false
            A = 2 .* rand(n_mon, n_mon) .- 1  # 生成范围在-1到1的随机矩阵
            A_symmetric = (A + A') / 2
        else
            A_symmetric=Symmetric(rand(n_mon, n_mon))
        end
        f = f + transpose(monc)*A_symmetric*mon
        # push!(g,1-sum(q[i]*q[i+n] for i = 1 + (k-1)*s : min(1 + (k-1)*s + size - 1, n)))
        push!(g,1-sum(q[i]*q[i+n] for i = 1 + (k - 1) * (size - 1) : min(1 + k * (size - 1), n)))
        monterm = [mono_to_term(mon[i], q, n) for i in 1:n_mon]
        # 目标 support 按 w_j * w_i^* 的顺序生成
        for i = 1:n_mon, j = 1:n_mon
            if A_symmetric[i, j]!=0
                a = deepcopy(monterm[i])
                b = deepcopy(monterm[j])
                bi = qtermadd(b, a, n)   # w_j * w_i^*
                push!(fsupp, bi)
                push!(fcoe, A_symmetric[i, j])
            end
        end
    end
    return f,g,fsupp,fcoe
end

using SparseArrays

function sparserandomsymfunc(q,n,d,sparsity;conjugates=false,coelimit=false)
    mon=NCMono[1]
    for j=1:d
        if conjugates!=false
            append!(mon, MP.monomials(q, j))
        else
            append!(mon, MP.monomials(q[1:n], j))
        end
    end
    monc=NCMono[]
    for i=1:length(mon)
        temp=prod(reverse(mon[i].vars).^reverse(mon[i].z))
        push!(monc,temp(q[1:n]=>q[n+1:2n],q[n+1:2n]=>q[1:n]))
        # push!(monc,temp(q[1]=>q[3],q[2]=>q[4],q[3]=>q[1],q[4]=>q[2]))
    end
    n=length(mon)
    if coelimit!=false
        A = 2 .* rand(n, n) .- 1  # 生成范围在-1到1的随机矩阵
        A_symmetric = (A + A') / 2
    else
        A_symmetric=Symmetric(rand(n,n))
    end
    B = sparse_symmetric_binary_matrix(n,sparsity)
    return transpose(monc)*(A_symmetric.*B)*mon
end

function cliques_sparserandomsymfunc(q,n,cn,size,d,sparsity;conjugates=false,coelimit=false)
    f = 0*q[1]^0
    g = DynamicPolynomials.Polynomial{DynamicPolynomials.NonCommutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, Int64}[]
    s = max(1, (n - size) ÷ (cn - 1)) 
    for k = 1:cn
        mon = NCMono[1]
        for j=1:d
            if conjugates!=false
                append!(mon,monomials(vcat(q[1 + (k-1)*s : min(1 + (k-1)*s + size - 1, n)],q[1 + (k-1)*s+n : min(1 + (k-1)*s + size - 1, n)+n]), j))
            else
                append!(mon,monomials(q[1 + (k-1)*s : min(1 + (k-1)*s + size - 1, n)], j))
            end
        end
        monc = NCMono[]
        for i = 1:length(mon)
            temp = prod(reverse(mon[i].vars).^reverse(mon[i].z))
            push!(monc,temp(q[1:n]=>q[n+1:2n],q[n+1:2n]=>q[1:n]))
        end
        n_mon=length(mon)
        if coelimit!=false
            A = 2 .* rand(n_mon, n_mon) .- 1  # 生成范围在-1到1的随机矩阵
            A_symmetric = (A + A') / 2
        else
            A_symmetric=Symmetric(rand(n_mon, n_mon))
        end
        B = sparse_symmetric_binary_matrix(n_mon,sparsity)
        f = f + transpose(monc)*(A_symmetric.*B)*mon
        push!(g,1-sum(q[i]*q[i+n] for i = 1 + (k-1)*s : min(1 + (k-1)*s + size - 1, n)))
    end
    monterm = [mono_to_term(mon[i], q, n) for i in 1:n_mon]

    fsupp = Vector{Vector{UInt16}}[]
    fcoe = QuaternionF64[]
    A_sparsity = A_symmetric.*B
    # 目标 support 按 w_j * w_i^* 的顺序生成
    for i = 1:n_mon, j = 1:n_mon
        if A_sparsity[i, j]!=0
            a = deepcopy(monterm[i])
            b = deepcopy(monterm[j])
            bi = qtermadd(b, a, n)   # w_j * w_i^*
            push!(fsupp, bi)
            push!(fcoe, A_sparsity[i, j])
        end
    end
    return f,g,fsupp,fcoe
end

function sparse_symmetric_binary_matrix(n::Int, sparsity::Float64)
    I = Int[]
    J = Int[]
    
    # 遍历上三角（含对角线）生成索引
    for i in 1:n
        # 处理对角线元素
        if rand(1) < sparsity
            push!(I, i)
            push!(J, i)
        end
        
        # 处理上三角非对角元素
        for j in (i+1):n
            if rand(1) < sparsity
                # 添加对称位置索引
                push!(I, i); push!(J, j)
                push!(I, j); push!(J, i)
            end
        end
    end
    
    # 创建稀疏矩阵（自动合并重复项）
    sparse(I, J, 1, n, n) .!= 0  # 转换为布尔型稀疏矩阵
end

function qrandomsymfunc(q,n,d;conjugates=false, sparsity = 0.0)
    mon=NCMono[1]
    for j=1:d
        if conjugates!=false
            append!(mon,monomials(q, j))
        else
            append!(mon,monomials(q[1:n], j))
        end
    end
    monc=NCMono[]
    for i=1:length(mon)
        temp=prod(reverse(mon[i].vars).^reverse(mon[i].z))
        push!(monc,temp(q[1:n]=>q[n+1:2n],q[n+1:2n]=>q[1:n]))
    end
    t=length(mon)
    A = rand(QuaternionF64,t,t) 
    A_symmetric = (A + adjoint(A)) / 2
        # 把 mon 转成内部 term
    monterm = [mono_to_term(mon[i], q, n) for i in 1:t]

    fsupp = Vector{Vector{UInt16}}[]
    fcoe = QuaternionF64[]

    # 目标 support 按 w_j * w_i^* 的顺序生成
    for i = 1:t, j = 1:t
        if A_symmetric[i, j]!=0
            a = deepcopy(monterm[i])
            b = deepcopy(monterm[j])
            bi = qtermadd(b, a, n) 
            # bi = qtermaddleft(a, b, n) 
            idx = findfirst(x -> x == bi, fsupp)

            if idx === nothing
                push!(fsupp, bi)
                push!(fcoe, A_symmetric[i,j])
            else
                fcoe[idx] += A_symmetric[i,j]
            end  # w_j * w_i^*
            # push!(fsupp, bi)
            # push!(fcoe, A_symmetric[i, j])
        end
    end
    if sparsity > 0.0
        C = sparse_symmetric_binary_matrix(t,sparsity)
        B = conj.(A_symmetric.*C)
        return transpose(monc)*(A_symmetric.*C)*mon, transpose(mon)*B*monc,fsupp,fcoe
    else
        # to real
        # B = conj.(A_symmetric)
        # return transpose(monc)*(A_symmetric)*mon, transpose(mon)*B*monc,fsupp,fcoe
        # for i = 1:t, j = 1:t
        #     if A_symmetric[i, j]!=0
        #         a = deepcopy(monterm[i])
        #         b = deepcopy(monterm[j])
        #         bi = qtermadd(b, a, n)   # w_j * w_i^*
        #         push!(fsupp, bi)
        #         push!(fcoe, A_symmetric[i, j])
        #     end
        # end
        # fr = deepcopy(monterm[1])
        # for i = 1:t, j = 1:t
        #     if A_symmetric[i, j]!=0
        #         a = deepcopy(monterm[i])
        #         b = deepcopy(monterm[j])
        #         bi = qtermadd(b, a, n) 
        #         fr += A_symmetric[i, j]*bi
        #     end
        # end
        B = transpose(A_symmetric)
        return transpose(monc)*(A_symmetric)*mon, transpose(mon)*B*monc,fsupp,fcoe

    end
end

function qrandomsymfunc2(q,n,d,Q;conjugates=false, sparsity = 0.0,given=false)
    mon=NCMono[1]
    for j=1:d
        if conjugates!=false
            append!(mon,monomials(q, j))
        else
            append!(mon,monomials(q[1:n], j))
        end
    end
    monc=NCMono[]
    for i=1:length(mon)
        temp=prod(reverse(mon[i].vars).^reverse(mon[i].z))
        push!(monc,temp(q[1:n]=>q[n+1:2n],q[n+1:2n]=>q[1:n]))
    end
    t=length(mon)
    if given == false
        # A = rand(QuaternionF64,t,t) 
        # A_symmetric = (A + adjoint(A)) / 2
        # A_symmetric = rand_hermitian_int_quat(t, 10)
        A_symmetric = sparse_hermitian_int_quat(t, 0.2, 10)
    else
        A_symmetric = Q
    end
        # 把 mon 转成内部 term
    monterm = [mono_to_term(mon[i], q, n) for i in 1:t]

    fsupp = Vector{Vector{UInt16}}[]
    fcoe = QuaternionF64[]

    # 目标 support 按 w_j * w_i^* 的顺序生成
    for i = 1:t, j = 1:t
        if A_symmetric[i, j]!=0
            a = deepcopy(monterm[i])
            b = deepcopy(monterm[j])
            println(i,j,a,b)
            bi = qtermadd(b, a, n) 
            # bi = qtermaddleft(a, b, n) 
            idx = findfirst(x -> x == bi, fsupp)

            if idx === nothing
                push!(fsupp, bi)
                push!(fcoe, A_symmetric[i,j])
            else
                fcoe[idx] += A_symmetric[i,j]
            end  # w_j * w_i^*
            # push!(fsupp, bi)
            # push!(fcoe, A_symmetric[i, j])
        end
    end
    if sparsity > 0.0
        C = sparse_symmetric_binary_matrix(t,sparsity)
        B = conj.(A_symmetric.*C)
        return transpose(monc)*(A_symmetric.*C)*mon, transpose(mon)*B*monc,fsupp,fcoe
    else
        # to real
        # B = conj.(A_symmetric)
        # return transpose(monc)*(A_symmetric)*mon, transpose(mon)*B*monc,fsupp,fcoe
        # for i = 1:t, j = 1:t
        #     if A_symmetric[i, j]!=0
        #         a = deepcopy(monterm[i])
        #         b = deepcopy(monterm[j])
        #         bi = qtermadd(b, a, n)   # w_j * w_i^*
        #         push!(fsupp, bi)
        #         push!(fcoe, A_symmetric[i, j])
        #     end
        # end
        # fr = deepcopy(monterm[1])
        # for i = 1:t, j = 1:t
        #     if A_symmetric[i, j]!=0
        #         a = deepcopy(monterm[i])
        #         b = deepcopy(monterm[j])
        #         bi = qtermadd(b, a, n) 
        #         fr += A_symmetric[i, j]*bi
        #     end
        # end
        B = transpose(A_symmetric)
        return transpose(monc)*(A_symmetric)*mon, transpose(mon)*B*monc,fsupp,fcoe

    end
end
# using Quaternions, LinearAlgebra, SparseArrays

# """
# 生成稀疏 Hermitian 四元数矩阵，元素为整数四元数。

# 参数:
# - t: 矩阵维度
# - density: 非零元密度（0~1），实际控制上三角（含对角）的填充率
# - bound: 整数分量的取值范围（-bound 到 bound）
# - zero_diag: 是否允许对角线为零（默认 false，对角线至少非零以保持可逆性？可根据需求调整）
# """
function sparse_hermitian_int_quat(t, density=0.2, bound=10, zero_diag=false)
    # 预分配稀疏矩阵存储（COO 格式）
    I = Int[]
    J = Int[]
    V = Quaternion{Int}[]
    
    for i in 1:t
        # 对角线：至少一个非零（除非 zero_diag=true）
        if zero_diag
            if rand(1) < density
                a = rand(-bound:bound)
                push!(I, i); push!(J, i); push!(V, quat(a, 0, 0, 0))
            end
        else
            # 对角线总是非零（实数整数）
            a = rand(-bound:bound)
            # 可选：避免全零矩阵，确保至少一个非零
            push!(I, i); push!(J, i); push!(V, quat(a, 0, 0, 0))
        end
        
        # 上三角部分 (i < j)
        for j in i+1:t
            if rand(1) < density
                a = rand(-bound:bound)
                b = rand(-bound:bound)
                c = rand(-bound:bound)
                d = rand(-bound:bound)
                q = quat(a, b, c, d)
                println(i,j)
                push!(I, i); push!(J, j); push!(V, q)
                # 下三角共轭对称将自动由 sparse 构造时处理？不，我们需要显式添加下三角以得到完整矩阵
                # 但最终若需完整矩阵，可先构建上三角，再镜像共轭。
            end
        end
    end
    
    # 构建上三角稀疏矩阵
    Q_tri = sparse(I, J, V, t, t)
    # 补全为 Hermitian：下三角 = adjoint(上三角)
    Q = Q_tri + adjoint(Q_tri) - Diagonal(diag(Q_tri))
    # for i = 1:7
    #     Q[i,i]= quat(0,0,0,0)
    # end
    # Q[1,6] = 0
    # Q[6,1] = 0
    # Q[2,4] = quat(1,1,0,0)
    # Q[4,2] = quat(1,-1,0,0)
    # Q[2,6] = quat(1,0,0,0)
    # Q[6,2] = quat(1,0,0,0)
    # Q[5,7] = quat(1,0,0,0)
    # Q[7,5] = quat(1,0,0,0)
    # # Q[2,4] = 0
    # # Q[4,2] = 0
    # # Q[2,6] = 0
    # # Q[6,2] = 0
    # # Q[5,7] = 0
    # # Q[7,5] = 0
    # # 减去重复的对角线
    return Q
end

function sparse_symmetric_pm1(t; density=0.2, zero_diag=false)
    I = Int[]
    J = Int[]
    V = Int[]

    for i in 1:t

        # 对角线
        if !zero_diag
            v = rand((-1, 1))
            push!(I, i)
            push!(J, i)
            push!(V, v)
        end

        # 上三角
        for j in i+1:t
            if rand(1) < density
                v = rand((-1, 1))

                # (i,j)
                push!(I, i)
                push!(J, j)
                push!(V, v)

                # (j,i)
                push!(I, j)
                push!(J, i)
                push!(V, v)
            end
        end
    end

    return sparse(I, J, V, t, t)
end

function rand_hermitian_int_quat(t, bound=10)
    Q = zeros(Quaternion{Int}, t, t)
    for i in 1:t
        for j in 1:i
            if i == j
                # 对角线：实数整数
                Q[i,i] = quat(rand(-bound:bound), 0, 0, 0)
            else
                # 下三角：随机整数四元数
                a = rand(-bound:bound)
                b = rand(-bound:bound)
                c = rand(-bound:bound)
                d = rand(-bound:bound)
                q = quat(a, b, c, d)
                Q[i,j] = q
                Q[j,i] = conj(q)   # 共轭填充上三角
            end
        end
    end
    return Q
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

function star(b::Vector{Vector{UInt16}},n)
    a = deepcopy(b)
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

function star(w)
    return prod(reverse(w.vars).^reverse(w.z))
end

function star(p)
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

function bget_qncbasis(n, d, k; ind=Vector{UInt16}(1:2n), binary=false)
    basis=[[UInt16[],UInt16[],UInt16[]]]
    for i = 1:d
        # if conjubasis!=false
        # append!(basis, _get_qncbasis_deg2(n, i, ind=ind, binary=binary))
        # else
        append!(basis, _get_qncbasis_deg(n, i, ind=ind, binary=binary))
        # end
    end
    basistemp=deepcopy(basis)
    # if !conjubasis
        # for i=1:n
            for i=1:length(basistemp)
                a=deepcopy(basistemp[i])
                ltemp=qtermadd([UInt16[],UInt16[],[UInt16(k)]],a,n)
                push!(basis,ltemp)
            end
        # end
    # end
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

#generate the standard monomial basis in the sparse form
function get_qncbasis(var::Vector{T}, n, d; ind=Vector{UInt16}(1:2n), binary=false,conjubasis=false) where T <: Union{UInt16, Int}
    basis=[[UInt16[],UInt16[],UInt16[]]]
    temp = copy(var)
    var_bar = unique([var;temp.+ind[n]])
    for i = 1:d
        if conjubasis!=false
            append!(basis, _get_qncbasis_deg2(var_bar, n, i, ind=ind, binary=binary))
        else
            append!(basis, _get_qncbasis_deg(var, n, i, ind=ind, binary=binary))
        end
    end
    return basis
end

function _get_qncbasis_deg(var::Vector{T}, n, d; ind=Vector{UInt16}(1:n), binary=false) where T <: Union{UInt16, Int}
    var = var[var.<= ind[n]]
    if d > 0
        basis=[[UInt16[],UInt16[],UInt16[]]]
        for i in var
            temp = _get_qncbasis_deg(var,n, d-1, ind=ind, binary=binary)
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
        for i = 1:length(basis)
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

function _get_qncbasis_deg2(var::Vector{T}, n, d; ind=Vector{UInt16}(1:2n), binary=false) where T <: Union{UInt16, Int}
    if d > 0
        basis=[[UInt16[],UInt16[],UInt16[]]]
        for i in var
            temp = _get_qncbasis_deg2(var,n, d-1, ind=ind, binary=binary)
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
        for i = 1:length(basis)
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

function qreduce_unitnorm(b,n; nb=0)
    b = standardterm(b, n)
    a = [UInt16[],UInt16[],copy(b[3])]
    return a
end

function qresort(supp, coe, n; nb=0)
    if nb > 0
        supp = qreduce_unitnorm.(supp,n,nb=nb)
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

"""
Canonical form for quaternion monomials under real-part equivalence:
Re(m1) = Re(m2)

Rule:
- ignore ordering under real part symmetry
- treat (q_i q_j) ~ (q_j q_i)
- treat conjugate-reversed forms as identical
"""
function canonical_qmono(m::Vector{Vector{UInt16}}, n::Int)

    # scalar parts
    r  = sort(deepcopy(m[1]))
    im = sort(deepcopy(m[2]))

    v = deepcopy(m[3])

    if isempty(v)
        return [r, im, v]
    end
    if v[1] == v[end] + UInt16(n)
        push!(r,v[end])
        push!(im,v[1])
        v = v[2:end-1]
    elseif v[1] == v[end] - UInt16(n)
        push!(im,v[end])
        push!(r,v[1])
        v = v[2:end-1]
    end

    ################################################
    # Step 1: cyclic canonicalization of full word
    ################################################

    rotations = Vector{Vector{UInt16}}()

    d = length(v)

    for k in 0:d-1
        rot = vcat(v[k+1:end], v[1:k])
        push!(rotations, rot)
    end

    v_cyclic = minimum(rotations)
    # if v_cyclic[1] == v_cyclic[end] + UInt16(n)
    #     push!(r,v_cyclic[end])
    #     push!(im,v_cyclic[1])
    #     v_cyclic = v_cyclic[2:end-1]
    # elseif v_cyclic[1] == v_cyclic[end] - UInt16(n)
    #     push!(im,v_cyclic[end])
    #     push!(r,v_cyclic[1])
    #     v_cyclic = v_cyclic[2:end-1]
    # end
    ################################################
    # Step 2: standardterm
    ################################################

    m1 = standardterm([r, im, v_cyclic], n)

    ################################################
    # Step 3: canonicalize remaining nc-part again
    ################################################

    v2 = m1[3]

    if !isempty(v2)

        rotations2 = Vector{Vector{UInt16}}()

        d2 = length(v2)

        for k in 0:d2-1
            rot2 = vcat(v2[k+1:end], v2[1:k])
            push!(rotations2, rot2)
        end

        v2 = minimum(rotations2)
    end

    ################################################
    # Final standard form
    ################################################

    return standardterm([m1[1], m1[2], v2],n)
    # return m1
end

function skew_entry(X,t,r,bs)

    if t < r
        idx = Int((2bs-t)*(t-1)/2) + (r-t)
        return X[idx]

    elseif t > r
        idx = Int((2bs-r)*(r-1)/2) + (t-r)
        return -X[idx]

    else
        return 0.0
    end
end


function generate_comm_constraints(n)

    Qi = QuaternionF64(0,1,0,0)
    Qj = QuaternionF64(0,0,1,0)
    Qk = QuaternionF64(0,0,0,1)

    pop_num = n*(n-1)

    supp = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef,pop_num)
    coe  = Vector{Vector{QuaternionF64}}(undef,3)

    cnt = 1

    for i in 1:n
        for j in 1:n

            i == j && continue

            supp[cnt] = [

                # -q_i * i * q_j
                # [mono_to_term(q[i],q,n),mono_to_term(q[j],q,n)],
                [[UInt16[], UInt16[], [UInt16(i)]],[UInt16[], UInt16[], [UInt16(j)]]],

                # +i * conj(q_i) * q_j
                # [mono_to_term(1,q,n),mono_to_term(q[i+n]*q[j],q,n)],
                [[UInt16[], UInt16[], UInt16[]],[UInt16[], UInt16[], [UInt16(i+n), UInt16(j)]]],

                # +q_j*q_i*i
                # [mono_to_term(q[j]*q[i],q,n), mono_to_term(1,q,n)],
                [[UInt16[], UInt16[], [UInt16(j), UInt16(i)]],[UInt16[], UInt16[], UInt16[]]],

                # -q_j*i*conj(q_i)
                # [mono_to_term(q[j],q,n), mono_to_term(q[i+n],q,n)]
                [[UInt16[], UInt16[], [UInt16(j)]],[UInt16[], UInt16[], [UInt16(i+n)]]]
            ]
            # supp[cnt][1] = [[UInt16[], UInt16[], [UInt16(i)]],[UInt16[], UInt16[], [UInt16(j)]]]
            # supp[cnt][2] = [[UInt16[], UInt16[], UInt16[]],[UInt16[], UInt16[], [UInt16(i+n), UInt16(j)]]]
            # supp[cnt][3] = [[UInt16[], UInt16[], [UInt16(j), UInt16(i)]],[UInt16[], UInt16[], UInt16[]]]
            # supp[cnt][4] = [[UInt16[], UInt16[], [UInt16(j)]],[UInt16[], UInt16[], [UInt16(i+n)]]]
            cnt += 1

        end
    end
    coe[1] = QuaternionF64[-Qi,Qi,Qi,-Qi]
    coe[2] = QuaternionF64[-Qj,Qj,Qj,-Qj]
    coe[3] = QuaternionF64[-Qk,Qk,Qk,-Qk]
    return supp, coe

end

