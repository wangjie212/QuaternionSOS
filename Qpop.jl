const Poly = DP.Polynomial{DP.Commutative{DP.CreationOrder}, Graded{LexOrder}, Float64}
const NCPoly = DP.Polynomial{DP.NonCommutative{DP.CreationOrder}, Graded{LexOrder}, Float64}
const NCMono = DP.Monomial{DP.NonCommutative{DP.CreationOrder}, Graded{LexOrder}}

function qtssos(pop, z, n::Int, d; fsupp=nothing, fcoe=nothing, numeq=0, rncnumeq=0, qncnumeq=0, RemSig=false, nb=0, CS="MF",cliques=[],TS="block",
    merge=false, md=3, solver="Mosek", QUIET=false, solve=true, solution=false, ipart=true,
    mosek_setting=mosek_para(), normality=0, conjubasis=true,addrcons=true,addicons=false)
    if addrcons
        rcomm_constraints = generate_rcomm_constraints(q, n)
        supp,coe = QPolys_info(vcat(pop,rcomm_constraints), z, n)
        numeq = numeq + length(rcomm_constraints)
        rncnumeq = rncnumeq + length(rcomm_constraints)
    else
        supp,coe = QPolys_info(pop, z, n)
    end
    for i = 1:length(supp)
        supp[i] = standardterm.(supp[i],n)
    end
    if ipart && fsupp !== nothing && fcoe !== nothing
        supp[1] = fsupp
        coe[1] = fcoe
    elseif (fsupp === nothing) ⊻ (fcoe === nothing)
        error("fsupp and fcoe must be provided together.")
    end
    println("*********************************** QSOS ***********************************")
    println("QSOS is launching...")
    if nb > 0
        supp[1],coe[1] = qresort(supp[1],coe[1],n;nb=nb)
    end
    m = length(supp) - 1
    dc = zeros(Int, m)
    for i = 1:m
        dc[i] = maximum([length(supp[i+1][j][1])+length(supp[i+1][j][2])+length(supp[i+1][j][3]) for j=1:length(supp[i+1])])
    end
    if cliques != []
        cql = length(cliques)
        cliquesize = length.(cliques)
    else
        time = @elapsed begin
        CS = CS == true ? "MF" : CS
        cliques,cql,cliquesize = clique_decomp(n, m, dc, rncnumeq, supp, order=d, alg=CS)
        end
        if CS != false && QUIET == false
            mc = maximum(cliquesize)
            println("Obtained the variable cliques in $time seconds.\nThe maximal size of cliques is $mc.")
        end
    end
    I,J,ncc = assign_constraint(n, m, numeq, supp, cliques, cql)
    rlorder = d*ones(Int, cql)
    ebasis = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, cql)
    basis = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, cql)
    for i = 1:cql
        ebasis[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(J[i]))
        if normality > 0
                basis[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(I[i])+1+cliquesize[i])
                basis[i][1] = get_qncbasis(cliques[i], n, rlorder[i]; conjubasis=conjubasis)
                for s = 1:cliquesize[i]
                    basis[i][s+1] = get_qncbasis(cliques[i], n, normality; conjubasis=false)
                    temp = deepcopy(basis[i][s+1])
                    for l=1:length(temp)
                        a=deepcopy(temp[l])
                        # ltemp=qtermadd([UInt16[],UInt16[],[cliques[i][s]]],a,n)
                        if ipart
                            ltemp=qtermadd(a,[UInt16[],UInt16[],[cliques[i][s]]],n)
                        else
                            ltemp=qtermaddleft([UInt16[],UInt16[],[cliques[i][s]]],a,n)
                        end
                        push!(basis[i][s+1],ltemp)
                    end
                    if nb > 0
                        basis[i][s+1] = qreduce_unitnorm.(basis[i][s+1], n,nb=nb)
                        unique!(basis[i][s+1])
                    end
                    unique!(basis[i][s+1])
                end
                for s = 1:length(I[i])
                    basis[i][s+1+cliquesize[i]] = get_qncbasis(cliques[i], n, rlorder[i]-Int(ceil(dc[I[i][s]]/2)); conjubasis=conjubasis)
                    if nb > 0
                        basis[i][s+1] = qreduce_unitnorm.(basis[i][s+1], n, nb=nb)
                        unique!(basis[i][s+1])
                    end
                end
        else
            basis[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(I[i])+1)
            basis[i][1] = get_qncbasis(cliques[i], n, d; conjubasis=conjubasis)
            if nb > 0
                basis[i][1] = qreduce_unitnorm.(basis[i][1], n, nb=nb)
                unique!(basis[i][1])
            end
            for s = 1:length(I[i])
                basis[i][s+1] = get_qncbasis(cliques[i], n, rlorder[i]-Int(ceil(dc[I[i][s]]/2)))
                if nb > 0
                    basis[i][s+1] = qreduce_unitnorm.(basis[i][s+1], n, nb=nb)
                    unique!(basis[i][s+1])
                end
            end
        end        
        for s = 1:length(J[i])
            ebasis[i][s] = get_qncbasis(cliques[i], n, rlorder[i]-Int(ceil(dc[J[i][s]]/2)))
            if nb > 0
                ebasis[i][s] = qreduce_unitnorm.(ebasis[i][s], n, nb=nb)
                unique!(ebasis[i][s])
            end
        end
    end
    tsupp = deepcopy(supp[1])
    for i = 2:m+1, j = 1:length(supp[i])
        push!(tsupp, supp[i][j])
    end
    sort!(tsupp)
    unique!(tsupp)
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    blocks,eblocks,cl,blocksize = get_blocks(n,rlorder, I, J, supp, cliques, cliquesize, cql, tsupp, basis, ebasis, TS=TS, ConjugateBasis=conjubasis, normality=normality, merge=merge, md=md, nb=nb)
    if QUIET == false
        mb = maximum(maximum.([maximum.(blocksize[i]) for i = 1:cql]))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    opt, ksupp, SDP_status =  qsolvesdp(n, m, rncnumeq, qncnumeq, rlorder, supp, coe, basis, ebasis, cliques, cql, cliquesize, I, J, ncc, blocks, eblocks, cl, blocksize, numeq=numeq, QUIET=QUIET, TS=TS, solver=solver, solve=solve, solution=solution, ipart=ipart,
    nb=nb, mosek_setting=mosek_setting, normality=normality, conjubasis=conjubasis,addicons=addicons)
    return opt
end

function QPolys_info(pop, z, n)
    coe = Vector{Vector{QuaternionF64}}(undef, length(pop))
    supp = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(pop))
    for k in eachindex(pop)
        mon = MultivariatePolynomials.monomials(pop[k])
        coe[k] = MultivariatePolynomials.coefficients(pop[k])
        lm = length(mon)
        supp[k] = [[[], [], []] for i=1:lm]
        temp=Vector{Vector{UInt16}}[]
        Ind=[]
        for i = 1:lm
            ind = mon[i].z.> 0
            vars = mon[i].vars[ind]
            exp = mon[i].z[ind]
            for j in eachindex(vars)
                l = ncbfind(z, 2n, vars[j])
                for t= 1: exp[j]
                    if length(supp[k][i][3])!=0
                        if supp[k][i][3][end]!= l-n && supp[k][i][3][end]!= l+n
                        push!(supp[k][i][3],l)
                        else
                        ab=(supp[k][i][3][end]<l) ? supp[k][i][3][end] : l
                        push!(supp[k][i][1],ab)
                        push!(supp[k][i][2],ab+UInt16(n))
                        deleteat!(supp[k][i][3],length(supp[k][i][3]))
                        end
                    else
                        push!(supp[k][i][3],l)
                    end
                end
            end
            supp[k][i][1]=sort(supp[k][i][1])
            supp[k][i][2]=sort(supp[k][i][2])
        end
    end
    return supp,coe
end

function qsolvesdp(n, m, rncnumeq, qncnumeq, rlorder, supp::Vector{Vector{Vector{Vector{UInt16}}}}, coe, basis, ebasis, cliques, cql, cliquesize, I, J, ncc, blocks, eblocks, cl, blocksize; 
    numeq=0, nb=0, QUIET=false, TS=false, solver="Mosek", solve=true,  solution=false, ipart=true, 
    mosek_setting=mosek_para(), normality=0, conjubasis=false,addicons=false)
    tsupp = Vector{Vector{UInt16}}[]
    ttsupp = Vector{Vector{UInt16}}[]
    if addicons
        addsupp,addcoe = generate_icomm_constraints(n)
    end
    for i = 1:cql
        a = normality > 0 ? cliquesize[i] + 1 : 1
        for s = 1 : a, j = 1:cl[i][s], k = 1:blocksize[i][s][j], r = 1:blocksize[i][s][j]
            at=deepcopy(basis[i][s][blocks[i][s][j][k]])
            b=deepcopy(basis[i][s][blocks[i][s][j][r]])
            if ipart
                @inbounds bi = qtermadd(b,at,n)
            else
                @inbounds bi = qtermaddleft(at,b,n)
            end
            if nb > 0
                bi = qreduce_unitnorm(bi,n, nb=nb)
            end

            if ipart != true
                bi= canonical_qmono(bi,n)
            else
                bi_cano = canonical_qmono(bi,n)
            end
            if nb > 0
                bi = qreduce_unitnorm(bi,n,nb=nb)
                if ipart 
                    bi_cano= qreduce_unitnorm(bi_cano,n, nb=nb)
                end
            end
            push!(tsupp, bi)
            if ipart
                push!(ttsupp,bi_cano)
            end
        end
        for (j, k) in enumerate(J[i])
            if m - rncnumeq < J[i][j] <= m - qncnumeq
                bs = length(eblocks[i][j])
                for t = 1:bs, r = 1:bs
                    for s = 1:length(supp[k+1])
                        a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                        b=deepcopy(supp[k+1][s])
                        c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                        if ipart
                            @inbounds bi = qtermadd3(b,c,a,n)
                        else
                            @inbounds bi = qtermadd3left(a,b,c,n)
                        end
                        if nb > 0
                            bi = qreduce_unitnorm(bi,n;nb=nb)
                        end
                        if ipart != true
                            bi= canonical_qmono(bi,n)
                        else
                            bi_cano = canonical_qmono(bi,n)
                        end
                        if nb > 0
                            bi= qreduce_unitnorm(bi,n, nb=nb)
                            if ipart 
                                bi_cano= qreduce_unitnorm(bi_cano,n, nb=nb)
                            end
                        end
                        push!(tsupp, bi)
                        if ipart
                            push!(ttsupp,bi_cano)
                        end
                    end
                    for s = 1:length(supp[k+1])
                        a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                        b=deepcopy(supp[k+1][s])
                        c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                        a_star = star(deepcopy(a),n)
                        b_star = star(deepcopy(b),n)
                        if ipart
                            @inbounds bi = qtermadd3(c,a_star,b,n)
                        else
                            @inbounds bi = qtermadd3left(a,b_star,c,n)
                        end
                        if nb > 0
                            bi = qreduce_unitnorm(bi,n;nb=nb)
                        end
                        if ipart != true
                            bi= canonical_qmono(bi,n)
                        else
                            bi_cano = canonical_qmono(bi,n)
                        end
                        if nb > 0
                            bi= qreduce_unitnorm(bi,n, nb=nb)
                            if ipart 
                                bi_cano= qreduce_unitnorm(bi_cano,n, nb=nb)
                            end
                        end
                        push!(tsupp, bi)
                        if ipart
                            push!(ttsupp,bi_cano)
                        end
                    end
                    if addicons
                        for ss = 1:length(addsupp)
                            for tt = 1:4
                                a=deepcopy(ebasis[i][j][t])
                                b=deepcopy(ebasis[i][j][r])
                                c=deepcopy(addsupp[ss][tt][1])
                                d=deepcopy(addsupp[ss][tt][2])
                                a_star = star(deepcopy(a),n)
                                @inbounds bi = qtermadd4(d,b,a_star,c,n)
                                if nb > 0
                                    bi = qreduce_unitnorm(bi,n, nb=nb)
                                end
                                if ipart != true
                                    bi= canonical_qmono(bi,n)
                                else
                                    bi_cano = canonical_qmono(bi,n)
                                end
                                if nb > 0
                                    bi = qreduce_unitnorm(bi,n, nb=nb)
                                    if ipart 
                                        bi_cano= qreduce_unitnorm(bi_cano,n, nb=nb)
                                    end
                                end
                                push!(tsupp, bi)
                                if ipart
                                    push!(ttsupp,bi_cano)
                                end
                            end
                        end
                        for ss = 1:length(addsupp)
                            for tt = 1:4
                                a=deepcopy(ebasis[i][j][t])
                                b=deepcopy(ebasis[i][j][r])
                                c=deepcopy(addsupp[ss][tt][1])
                                d=deepcopy(addsupp[ss][tt][2])
                                a_star = star(deepcopy(a),n)
                                c_star = star(deepcopy(c),n)
                                d_star = star(deepcopy(d),n)
                                @inbounds bi = qtermadd4(c_star,b,a_star,d_star,n)
                                if nb > 0
                                    bi = qreduce_unitnorm(bi,n, nb=nb)
                                end
                                if ipart != true
                                    bi= canonical_qmono(bi,n)
                                else
                                    bi_cano = canonical_qmono(bi,n)
                                end
                                if nb > 0
                                    bi = qreduce_unitnorm(bi,n, nb=nb)
                                    if ipart 
                                        bi_cano= qreduce_unitnorm(bi_cano,n, nb=nb)
                                    end
                                end
                                push!(tsupp, bi)
                                if ipart
                                    push!(ttsupp,bi_cano)
                                end
                            end
                        end
                    end
                end
            elseif J[i][j] > m - qncnumeq
                bs = length(eblocks[i][j])
                for t = 1:bs, r = 1:bs
                    for s = 1:length(supp[k+1])
                        a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                        b=deepcopy(supp[k+1][s])
                        c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                        @inbounds bi = qtermadd3(b,c,a,n)
                        if nb > 0
                            bi = qreduce_unitnorm(bi,n;nb=nb)
                        end
                        if ipart != true
                            bi= canonical_qmono(bi,n)
                        end
                        if nb > 0
                            bi = qreduce_unitnorm(bi,n, nb=nb)
                        end
                        push!(tsupp, bi)
                    end
                    for s = 1:length(supp[k+1])
                        a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                        b=deepcopy(supp[k+1][s])
                        c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                        a_star = star(a,n)
                        @inbounds bi = qtermadd3(c,a_star,b,n)
                        if nb > 0
                            bi = qreduce_unitnorm(bi,n;nb=nb)
                        end
                        if ipart != true
                            bi= canonical_qmono(bi,n)
                        end
                        if nb > 0
                            bi= qreduce_unitnorm(bi,n, nb=nb)
                        end
                        push!(tsupp, bi)
                    end
                end
            end
        end
    end
    if TS != false
        gsupp = get_gsupp(rlorder, basis, ebasis, supp, cql, I, J, ncc, blocks, eblocks, cl, blocksize, cliquesize, ConjugateBasis=conjubasis, nb=nb, normality=normality)
        append!(tsupp, gsupp)
    end
    sort!(tsupp)
    unique!(tsupp)
    sort!(ttsupp)
    unique!(ttsupp)
    ksupp = tsupp
    #initial set
    objv = SDP_status= nothing
    if solve == true
        ltsupp = length(tsupp)
        lttsupp = length(ttsupp)
        if QUIET == false
            println("Assembling the SDP...")
        end
        model = Model(Mosek.Optimizer)
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
        rcons = [AffExpr(0) for i=1:ltsupp]
        if ipart==true || qncnumeq > 0
            icons = [AffExpr(0) for i=1:ltsupp]
            jcons = [AffExpr(0) for i=1:ltsupp]
            kcons = [AffExpr(0) for i=1:ltsupp]
        end
        hnom = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, n)
        pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, cql)
        for i = 1:cql
            a = normality > 0 ? cliquesize[i] + 1 : 1
            pos[i] = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, length(I[i])+a)
            for j = 1:a
                pos[i][j] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[i][j])
                for l = 1:cl[i][j]
                    @inbounds bs = blocksize[i][j][l]
                    if bs == 1
                        pos[i][j][l] = @variable(model, lower_bound=0)
                        at=deepcopy(basis[i][j][blocks[i][j][l][1]])
                        b=deepcopy(basis[i][j][blocks[i][j][l][1]])
                        if ipart
                            @inbounds bi = qtermadd(b,at,n)
                        else
                            @inbounds bi = qtermaddleft(at,b,n)
                        end
                        if nb > 0
                            bi = qreduce_unitnorm(bi,n;nb=nb)
                        end
                        bi_cano= canonical_qmono(bi,n)
                        if bi_cano != canonical_qmono(deepcopy(bi_cano),n)
                            println(bi_cano)
                        end 
                        if nb > 0
                            bi = qreduce_unitnorm(bi, n, nb=nb)
                            if ipart
                                bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                            end
                        end
                        Locb = bfind(tsupp, ltsupp, bi)
                        if ipart == true
                            Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                            @inbounds add_to_expression!(rcons[Locb_cano], pos[i][j][l])
                            @inbounds add_to_expression!(icons[Locb], pos[i][j][l])
                            @inbounds add_to_expression!(jcons[Locb], pos[i][j][l])
                            @inbounds add_to_expression!(kcons[Locb], pos[i][j][l])
                        else
                            @inbounds add_to_expression!(rcons[Locb], pos[i][j][l])
                        end
                    else
                        if ipart == true
                            pos[i][j][l] = @variable(model, [1:4bs, 1:4bs], PSD)
                        else
                            pos[i][j][l] = @variable(model, [1:bs, 1:bs], PSD)
                        end
                        for t = 1:bs, r = 1:bs
                            at=deepcopy(basis[i][j][blocks[i][j][l][t]])
                            b=deepcopy(basis[i][j][blocks[i][j][l][r]])
                            if ipart
                                @inbounds bi = qtermadd(b,at,n)
                            else
                                @inbounds bi = qtermaddleft(at,b,n)
                            end
                            if nb > 0
                                bi = qreduce_unitnorm(bi,n;nb=nb)
                            end
                            if ipart != true
                                bi= canonical_qmono(bi,n)
                            else
                                bi_cano= canonical_qmono(bi,n)
                            end
                            if nb > 0
                                bi = qreduce_unitnorm(bi, n, nb=nb)
                                if ipart
                                    bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                end
                            end
                            Locb = bfind(tsupp, ltsupp, bi)
                            if ipart == true
                                Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                                @inbounds add_to_expression!(rcons[Locb_cano], pos[i][j][l][t,r]+pos[i][j][l][t+bs,r+bs]+pos[i][j][l][t+2*bs,r+2*bs]+pos[i][j][l][t+3*bs,r+3*bs])
                                @inbounds add_to_expression!(icons[Locb], pos[i][j][l][t+bs,r]-pos[i][j][l][t,r+bs]+pos[i][j][l][t+3*bs,r+2*bs]-pos[i][j][l][t+2*bs,r+3*bs]) 
                                @inbounds add_to_expression!(jcons[Locb], pos[i][j][l][t+2*bs,r]-pos[i][j][l][t,r+2*bs]-pos[i][j][l][t+3*bs,r+bs]+pos[i][j][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], pos[i][j][l][t+3*bs,r]-pos[i][j][l][t,r+3*bs]+pos[i][j][l][t+2*bs,r+bs]-pos[i][j][l][t+bs,r+2*bs])
                            else
                                @inbounds add_to_expression!(rcons[Locb], pos[i][j][l][t,r])
                            end            
                        end
                    end
                end
            end
        end
        for i = 1:cql, (j, k) in enumerate(I[i])
            a = normality > 0 ? cliquesize[i] + 1 : 1
            pos[i][j+a] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[i][j+a])
            for l = 1: cl[i][j+a]
                bs = blocksize[i][j+a][l]
                if bs == 1
                    pos[i][j+a][l] = @variable(model, lower_bound=0)
                    for s = 1:length(supp[k+1])
                        at=deepcopy(basis[i][j+a][blocks[i][j+a][l][1]])
                        b=deepcopy(supp[k+1][s])
                        c=deepcopy(basis[i][j+a][blocks[i][j+a][l][1]])
                        if ipart
                            @inbounds bi = qtermadd3(b,c,at,n)
                        else
                            @inbounds bi = qtermadd3left(at,b,c,n)
                        end
                        if nb > 0
                            bi = qreduce_unitnorm(bi,n;nb=nb)
                        end
                        if ipart != true
                                bi= canonical_qmono(bi,n)
                        else
                            bi_cano= canonical_qmono(bi,n)
                        end
                        if nb > 0
                            bi = qreduce_unitnorm(bi, n, nb=nb)
                            if ipart
                                bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                            end
                        end
                        Locb = bfind(tsupp, ltsupp, bi)
                        if ipart == true
                            Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                            @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), pos[i][j+a][l])
                            @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[1], pos[i][j+a][l])
                            @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[2], pos[i][j+a][l])
                            @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], pos[i][j+a][l])
                        else
                            @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[i][j+a][l])
                        end
                    end 
                else
                    if ipart == true
                        pos[i][j+a][l]= @variable(model, [1:4bs, 1:4bs], PSD)
                    else
                        pos[i][j+a][l]= @variable(model, [1:bs, 1:bs], PSD)
                    end
                    for t = 1:bs, r = 1:bs
                        for s = 1:length(supp[k+1])
                            at=deepcopy(basis[i][j+a][blocks[i][j+a][l][t]])
                            b=deepcopy(supp[k+1][s])
                            c=deepcopy(basis[i][j+a][blocks[i][j+a][l][r]])
                            if ipart
                                @inbounds bi = qtermadd3(b,c,at,n)
                            else
                                @inbounds bi = qtermadd3left(at,b,c,n)
                            end
                            if nb > 0
                                bi = qreduce_unitnorm(bi,n;nb=nb)
                            end
                            if ipart != true
                                bi= canonical_qmono(bi,n)
                            else
                                bi_cano= canonical_qmono(bi,n)
                            end
                            if nb > 0
                                bi = qreduce_unitnorm(bi, n, nb=nb)
                                if ipart
                                    bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                end
                            end
                            Locb = bfind(tsupp, ltsupp, bi)
                            if ipart == true
                                Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                                @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), pos[i][j+a][l][t,r]+pos[i][j+a][l][t+bs,r+bs]+pos[i][j+a][l][t+2*bs,r+2*bs]+pos[i][j+a][l][t+3*bs,r+3*bs])
                                @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[1], pos[i][j+a][l][t+bs,r]-pos[i][j+a][l][t,r+bs]+pos[i][j+a][l][t+3*bs,r+2*bs]-pos[i][j+a][l][t+2*bs,r+3*bs])
                                @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[2], pos[i][j+a][l][t+2*bs,r]-pos[i][j+a][l][t,r+2*bs]-pos[i][j+a][l][t+3*bs,r+bs]+pos[i][j+a][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[3], pos[i][j+a][l][t+3*bs,r]-pos[i][j+a][l][t,r+3*bs]+pos[i][j+a][l][t+2*bs,r+bs]-pos[i][j+a][l][t+bs,r+2*bs])
    
                                @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], pos[i][j+a][l][t,r]+pos[i][j+a][l][t+bs,r+bs]+pos[i][j+a][l][t+2*bs,r+2*bs]+pos[i][j+a][l][t+3*bs,r+3*bs])
                                @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), pos[i][j+a][l][t+bs,r]-pos[i][j+a][l][t,r+bs]+pos[i][j+a][l][t+3*bs,r+2*bs]-pos[i][j+a][l][t+2*bs,r+3*bs])
                                #+A_J-A_K->-A_J+A_K/4.12
                                @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[3], pos[i][j+a][l][t+2*bs,r]-pos[i][j+a][l][t,r+2*bs]-pos[i][j+a][l][t+3*bs,r+bs]+pos[i][j+a][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[2], pos[i][j+a][l][t+3*bs,r]-pos[i][j+a][l][t,r+3*bs]+pos[i][j+a][l][t+2*bs,r+bs]-pos[i][j+a][l][t+bs,r+2*bs])
                                            
                                @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[2], pos[i][j+a][l][t,r]+pos[i][j+a][l][t+bs,r+bs]+pos[i][j+a][l][t+2*bs,r+2*bs]+pos[i][j+a][l][t+3*bs,r+3*bs])
                                #\A_R(J_X)-\A_I(K_X)-\A_J(R_X)+\A_K(I_X)=b_J->-\A_J(R_X)-\A_K(I_X)+\A_R(J_X)+\A_I(K_X)=b_J
                                @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[3], pos[i][j+a][l][t+bs,r]-pos[i][j+a][l][t,r+bs]+pos[i][j+a][l][t+3*bs,r+2*bs]-pos[i][j+a][l][t+2*bs,r+3*bs])
                                @inbounds add_to_expression!(jcons[Locb], real(coe[k+1][s]), pos[i][j+a][l][t+2*bs,r]-pos[i][j+a][l][t,r+2*bs]-pos[i][j+a][l][t+3*bs,r+bs]+pos[i][j+a][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[1], pos[i][j+a][l][t+3*bs,r]-pos[i][j+a][l][t,r+3*bs]+pos[i][j+a][l][t+2*bs,r+bs]-pos[i][j+a][l][t+bs,r+2*bs])
                                #\A_R(K_X)+\A_I(J_X)-\A_J(I_X)-\A_K(R_X)=b_K->-\A_K(R_X)+\A_J(I_X)-\A_I(J_X)+\A_R(K_X)=b_K                
                                @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], pos[i][j+a][l][t,r]+pos[i][j+a][l][t+bs,r+bs]+pos[i][j+a][l][t+2*bs,r+2*bs]+pos[i][j+a][l][t+3*bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[2], pos[i][j+a][l][t+bs,r]-pos[i][j+a][l][t,r+bs]+pos[i][j+a][l][t+3*bs,r+2*bs]-pos[i][j+a][l][t+2*bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], pos[i][j+a][l][t+2*bs,r]-pos[i][j+a][l][t,r+2*bs]-pos[i][j+a][l][t+3*bs,r+bs]+pos[i][j+a][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), pos[i][j+a][l][t+3*bs,r]-pos[i][j+a][l][t,r+3*bs]+pos[i][j+a][l][t+2*bs,r+bs]-pos[i][j+a][l][t+bs,r+2*bs])
                            else
                                @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[i][j+a][l][t,r])
                            end                   
                        end
                    end
                end
            end
        end

        ## equality constraints

        epos = Vector{Vector{Any}}(undef, cql)
        if addicons
            addepos = Vector{Vector{Any}}(undef, cql)
        end
        for i = 1:cql
            if !isempty(J[i])
                epos[i] = Vector{Any}(undef, length(J[i]))
                # println(cliques[i])
                if addicons
                    if cql > 1
                        addepos[i] = Vector{Any}(undef, 3*n*(n-1))
                    else
                        addepos[i] = Vector{Any}(undef, 3*length(J[i]))
                    end
                end
                for (j, k) in enumerate(J[i])
                    if J[i][j] <= m - rncnumeq - qncnumeq
                        bs = length(eblocks[i][j])
                        if bs == 1
                            epos[i][j] = @variable(model)
                            for s = 1:length(supp[k+1])
                                a=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                                b=deepcopy(supp[k+1][s])
                                c=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                                if ipart
                                    @inbounds bi = qtermadd3(b,c,a,n)
                                else
                                    @inbounds bi = qtermadd3left(a,b,c,n)
                                end
                                if nb > 0
                                    bi = qreduce_unitnorm(bi,n;nb=nb)
                                end
                                if ipart != true
                                    bi= canonical_qmono(bi,n)
                                else
                                    bi_cano= canonical_qmono(bi,n)
                                end
                                if nb > 0
                                    bi = qreduce_unitnorm(bi, n, nb=nb)
                                    if ipart
                                        bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                    end
                                end
                                Locb = bfind(tsupp, ltsupp, bi)
                                if ipart == true
                                    Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                                    @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j])
                                    @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[1], epos[i][j])
                                    @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j])
                                    @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], epos[i][j])
                                else
                                    @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j])
                                end
                            end
                        else
                            if ipart == true
                                epos[i][j]= @variable(model, [1:4bs, 1:4bs], Symmetric)
                            else
                                epos[i][j]= @variable(model, [1:bs, 1:bs], Symmetric)
                            end
                            for t = 1:bs, r = 1:bs
                                for s = 1:length(supp[k+1])
                                    a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                                    b=deepcopy(supp[k+1][s])
                                    c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                                    if ipart
                                        @inbounds bi = qtermadd3(b,c,a,n)
                                    else
                                        @inbounds bi = qtermadd3left(a,b,c,n)
                                    end
                                    if nb > 0
                                        bi = qreduce_unitnorm(bi,n;nb=nb)
                                    end
                                    if ipart != true
                                        bi= canonical_qmono(bi,n)
                                    else
                                        bi_cano= canonical_qmono(bi,n)
                                    end
                                    if nb > 0
                                        bi = qreduce_unitnorm(bi, n, nb=nb)
                                        if ipart
                                            bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                        end
                                    end
                                    Locb = bfind(tsupp, ltsupp, bi)
                                    if ipart == true
                                        Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                                        @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[1], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[2], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[3], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
            
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        #+A_J-A_K->-A_J+A_K/4.12
                                        @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[3], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        #\A_R(J_X)-\A_I(K_X)-\A_J(R_X)+\A_K(I_X)=b_J->-\A_J(R_X)-\A_K(I_X)+\A_R(J_X)+\A_I(K_X)=b_J            
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], real(coe[k+1][s]), epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[1], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        #\A_R(K_X)+\A_I(J_X)-\A_J(I_X)-\A_K(R_X)=b_K->-\A_K(R_X)+\A_J(I_X)-\A_I(J_X)+\A_R(K_X)=b_K               
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                    else
                                        @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][t,r])
                                    end                   
                                end
                            end
                        end
                    ### asymmetric equality constraints with real coefficients
                    elseif J[i][j] <= m - qncnumeq
                        bs = length(eblocks[i][j])
                        if bs == 1
                            epos[i][j] = @variable(model)
                            for s = 1:length(supp[k+1])
                                a=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                                b=deepcopy(supp[k+1][s])
                                c=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                                if ipart
                                    @inbounds bi = qtermadd3(b,c,a,n)
                                else
                                    @inbounds bi = qtermadd3left(a,b,c,n)
                                end
                                if nb > 0
                                    bi = qreduce_unitnorm(bi,n;nb=nb)
                                end
                                if ipart != true
                                    bi= canonical_qmono(bi,n)
                                else
                                    bi_cano= canonical_qmono(bi,n)
                                end
                                if nb > 0
                                    bi = qreduce_unitnorm(bi, n, nb=nb)
                                    if ipart
                                        bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                    end
                                end
                                Locb = bfind(tsupp, ltsupp, bi)
                                if ipart == true
                                    Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                                    @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j])
                                    @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[1], epos[i][j])
                                    @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j])
                                    @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], epos[i][j])
                                else
                                    @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j])
                                end
                            end
                            for s = 1:length(supp[k+1])
                                a=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                                b=deepcopy(supp[k+1][s])
                                c=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                                a_star = star(a,n)
                                b_star = star(b,n)
                                if ipart
                                    @inbounds bi = qtermadd3(c,a_star,b,n)
                                else
                                    @inbounds bi = qtermadd3left(a,b_star,c,n)
                                end
                                if nb > 0
                                    bi = qreduce_unitnorm(bi,n;nb=nb)
                                end
                                if ipart != true
                                    bi= canonical_qmono(bi,n)
                                else
                                    bi_cano= canonical_qmono(bi,n)
                                end
                                if nb > 0
                                    bi = qreduce_unitnorm(bi, n, nb=nb)
                                    if ipart
                                        bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                    end
                                end
                                Locb = bfind(tsupp, ltsupp, bi)
                                if ipart == true
                                    Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                                    @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j])
                                    @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[1], epos[i][j])
                                    @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j])
                                    @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], epos[i][j])
                                else
                                    @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j])
                                end
                            end
                        else
                            if ipart == true
                                epos[i][j]= @variable(model, [1:4bs, 1:4bs])
                            else
                                epos[i][j]= @variable(model, [1:bs, 1:bs])
                            end
                            for t = 1:bs, r = 1:bs
                                for s = 1:length(supp[k+1])
                                    a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                                    b=deepcopy(supp[k+1][s])
                                    c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                                    if ipart
                                        @inbounds bi = qtermadd3(b,c,a,n)
                                    else
                                        @inbounds bi = qtermadd3left(a,b,c,n)
                                    end
                                    if nb > 0
                                        bi = qreduce_unitnorm(bi,n;nb=nb)
                                    end
                                    if ipart != true
                                        bi= canonical_qmono(bi,n)
                                    else
                                        bi_cano= canonical_qmono(bi,n)
                                    end
                                    if nb > 0
                                        bi = qreduce_unitnorm(bi, n, nb=nb)
                                        if ipart
                                            bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                        end
                                    end
                                    Locb = bfind(tsupp, ltsupp, bi)
                                    if ipart == true
                                        Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                                        @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[1], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[2], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[3], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
            
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        #+A_J-A_K->-A_J+A_K/4.12
                                        @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[3], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        #\A_R(J_X)-\A_I(K_X)-\A_J(R_X)+\A_K(I_X)=b_J->-\A_J(R_X)-\A_K(I_X)+\A_R(J_X)+\A_I(K_X)=b_J            
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], real(coe[k+1][s]), epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[1], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        #\A_R(K_X)+\A_I(J_X)-\A_J(I_X)-\A_K(R_X)=b_K->-\A_K(R_X)+\A_J(I_X)-\A_I(J_X)+\A_R(K_X)=b_K               
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                    else
                                        @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][t,r])
                                    end               
                                end
                                for s = 1:length(supp[k+1])
                                    a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                                    b=deepcopy(supp[k+1][s])
                                    c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                                    a_star = star(a,n)
                                    b_star = star(b,n)
                                    if ipart
                                        @inbounds bi = qtermadd3(c,a_star,b,n)
                                    else
                                        @inbounds bi = qtermadd3left(a,b_star,c,n)
                                    end
                                    if nb > 0
                                        bi = qreduce_unitnorm(bi, n;nb=nb)
                                    end
                                    if ipart != true
                                        bi= canonical_qmono(bi,n)
                                    else
                                        bi_cano= canonical_qmono(bi,n)
                                    end
                                    if nb > 0
                                        bi = qreduce_unitnorm(bi, n, nb=nb)
                                        if ipart
                                            bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                        end
                                    end
                                    Locb = bfind(tsupp, ltsupp, bi)
                                    if ipart == true
                                        Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                                        @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j][r,t]+epos[i][j][r+bs,t+bs]+epos[i][j][r+2*bs,t+2*bs]+epos[i][j][r+3*bs,t+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[1], -epos[i][j][r+bs,t]+epos[i][j][r,t+bs]-epos[i][j][r+3*bs,t+2*bs]+epos[i][j][r+2*bs,t+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[2], -epos[i][j][r+2*bs,t]+epos[i][j][r,t+2*bs]+epos[i][j][r+3*bs,t+bs]-epos[i][j][r+bs,t+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[3], -epos[i][j][r+3*bs,t]+epos[i][j][r,t+3*bs]-epos[i][j][r+2*bs,t+bs]+epos[i][j][r+bs,t+2*bs])
            
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][r,t]+epos[i][j][r+bs,t+bs]+epos[i][j][r+2*bs,t+2*bs]+epos[i][j][r+3*bs,t+3*bs])
                                        #+A_J-A_K->-A_J+A_K/4.12
                                        @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), -epos[i][j][r+bs,t]+epos[i][j][r,t+bs]-epos[i][j][r+3*bs,t+2*bs]+epos[i][j][r+2*bs,t+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[3], -epos[i][j][r+2*bs,t]+epos[i][j][r,t+2*bs]+epos[i][j][r+3*bs,t+bs]-epos[i][j][r+bs,t+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[2], -epos[i][j][r+3*bs,t]+epos[i][j][r,t+3*bs]-epos[i][j][r+2*bs,t+bs]+epos[i][j][r+bs,t+2*bs])
                                        #\A_R(J_X)-\A_I(K_X)-\A_J(R_X)+\A_K(I_X)=b_J->-\A_J(R_X)-\A_K(I_X)+\A_R(J_X)+\A_I(K_X)=b_J            
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][r,t]+epos[i][j][r+bs,t+bs]+epos[i][j][r+2*bs,t+2*bs]+epos[i][j][r+3*bs,t+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[3], -epos[i][j][r+bs,t]+epos[i][j][r,t+bs]-epos[i][j][r+3*bs,t+2*bs]+epos[i][j][r+2*bs,t+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], real(coe[k+1][s]), -epos[i][j][r+2*bs,t]+epos[i][j][r,t+2*bs]+epos[i][j][r+3*bs,t+bs]-epos[i][j][r+bs,t+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[1], -epos[i][j][r+3*bs,t]+epos[i][j][r,t+3*bs]-epos[i][j][r+2*bs,t+bs]+epos[i][j][r+bs,t+2*bs])
                                        #\A_R(K_X)+\A_I(J_X)-\A_J(I_X)-\A_K(R_X)=b_K->-\A_K(R_X)+\A_J(I_X)-\A_I(J_X)+\A_R(K_X)=b_K               
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][r,t]+epos[i][j][r+bs,t+bs]+epos[i][j][r+2*bs,t+2*bs]+epos[i][j][r+3*bs,t+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[2], -epos[i][j][r+bs,t]+epos[i][j][r,t+bs]-epos[i][j][r+3*bs,t+2*bs]+epos[i][j][r+2*bs,t+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], -epos[i][j][r+2*bs,t]+epos[i][j][r,t+2*bs]+epos[i][j][r+3*bs,t+bs]-epos[i][j][r+bs,t+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), -epos[i][j][r+3*bs,t]+epos[i][j][r,t+3*bs]-epos[i][j][r+2*bs,t+bs]+epos[i][j][r+bs,t+2*bs])
                                        
                                    else
                                        @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][r,t])
                                    end         
                                end
                            end
                        end
                    ### asymmetric equality constraints with left quaternion coefficients
                    else
                        bs = length(eblocks[i][j])
                        if bs == 1
                            epos[i][j] = @variable(model)
                            for s = 1:length(supp[k+1])
                                a=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                                b=deepcopy(supp[k+1][s])
                                c=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                                # if ipart
                                @inbounds bi = qtermadd3(b,c,a,n)
                                if nb > 0
                                    bi = qreduce_unitnorm(bi, n;nb=nb)
                                end
                                if ipart != true
                                    bi= canonical_qmono(bi,n)
                                else
                                    bi_cano= canonical_qmono(bi,n)
                                end
                                if nb > 0
                                    bi = qreduce_unitnorm(bi, n, nb=nb)
                                    if ipart
                                        bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                    end
                                end
                                Locb = bfind(tsupp, ltsupp, bi)
                                Locb_cano = bfind(tsupp, ltsupp, bi_cano)
                                @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j])
                                if qncnumeq > 0
                                    @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[1], epos[i][j])
                                    @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j])
                                    @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], epos[i][j])
                                end
                            end
                            for s = 1:length(supp[k+1])
                                a=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                                b=deepcopy(supp[k+1][s])
                                c=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                                a_star = star(a,n)
                                @inbounds bi = qtermadd3(c,a_star,b,n)
                                if nb > 0
                                    bi = qreduce_unitnorm(bi, n;nb=nb)
                                end
                                if ipart != true
                                    bi= canonical_qmono(bi,n)
                                else
                                    bi_cano= canonical_qmono(bi,n)
                                end
                                if nb > 0
                                    bi = qreduce_unitnorm(bi, n, nb=nb)
                                    if ipart
                                        bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                    end
                                end
                                Locb = bfind(tsupp, ltsupp, bi)
                                Locb_cano = bfind(tsupp, ltsupp, bi_cano)
                                @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j])
                                if qncnumeq > 0
                                    @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[1], epos[i][j])
                                    @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j])
                                    @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], epos[i][j])
                                end
                            end
                        else
                            epos[i][j]= @variable(model, [1:4bs, 1:4bs])
                            for t = 1:bs, r = 1:bs
                                for s = 1:length(supp[k+1])
                                    a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                                    b=deepcopy(supp[k+1][s])
                                    c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                                    @inbounds bi = qtermadd3(b,c,a,n)
                                    if nb > 0
                                        bi = qreduce_unitnorm(bi, n;nb=nb)
                                    end
                                    if ipart != true
                                        bi= canonical_qmono(bi,n)
                                    else
                                        bi_cano= canonical_qmono(bi,n)
                                    end
                                    if nb > 0
                                        bi = qreduce_unitnorm(bi, n, nb=nb)
                                        if ipart
                                            bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                        end
                                    end
                                    Locb = bfind(tsupp, ltsupp, bi)
                                    Locb_cano = bfind(tsupp, ltsupp, bi_cano)
                                    if qncnumeq > 0
                                    # if ipart == true
                                        @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[1], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[2], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[3], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
            
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        #+A_J-A_K->-A_J+A_K/4.12
                                        @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[3], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        #\A_R(J_X)-\A_I(K_X)-\A_J(R_X)+\A_K(I_X)=b_J->-\A_J(R_X)-\A_K(I_X)+\A_R(J_X)+\A_I(K_X)=b_J            
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], real(coe[k+1][s]), epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[1], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        #\A_R(K_X)+\A_I(J_X)-\A_J(I_X)-\A_K(R_X)=b_K->-\A_K(R_X)+\A_J(I_X)-\A_I(J_X)+\A_R(K_X)=b_K               
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                    end                 
                                end
                                for s = 1:length(supp[k+1])
                                    a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                                    b=deepcopy(supp[k+1][s])
                                    c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                                    a_star = star(a,n)
                                    @inbounds bi = qtermadd3(c,a_star,b,n)
                                    if nb > 0
                                        bi = qreduce_unitnorm(bi, n;nb=nb)
                                    end
                                    if ipart != true
                                        bi= canonical_qmono(bi,n)
                                    else
                                        bi_cano= canonical_qmono(bi,n)
                                    end
                                    if nb > 0
                                        bi = qreduce_unitnorm(bi, n, nb=nb)
                                        if ipart
                                            bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                        end
                                    end
                                    Locb = bfind(tsupp, ltsupp, bi)
                                    Locb_cano = bfind(tsupp, ltsupp, bi_cano)
                                    if qncnumeq > 0
                                        @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[1], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[2], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[3], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
            
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        #+A_J-A_K->-A_J+A_K/4.12
                                        @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[3], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        #\A_R(J_X)-\A_I(K_X)-\A_J(R_X)+\A_K(I_X)=b_J->-\A_J(R_X)-\A_K(I_X)+\A_R(J_X)+\A_I(K_X)=b_J            
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], real(coe[k+1][s]), epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[1], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        #\A_R(K_X)+\A_I(J_X)-\A_J(I_X)-\A_K(R_X)=b_K->-\A_K(R_X)+\A_J(I_X)-\A_I(J_X)+\A_R(K_X)=b_K               
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                    end                 
                                end
                            end
                        end
                    end
                end
                ### im(q_i)*q_j = q_j*im(q_i), i!=j
                if addicons
                    addebasis = ebasis[i][end]
                    bs = length(addebasis)
                    if bs == 1
                        for rr = 1:3
                            for ss = 1:length(addsupp)
                                addepos[i][(rr-1)*length(addsupp)+ss] = @variable(model)
                                for tt= 1:4
                                    a = deepcopy(addebasis[1])
                                    b = deepcopy(addebasis[1])
                                    c = deepcopy(addsupp[ss][tt][1])
                                    d = deepcopy(addsupp[ss][tt][2])
                                    a_star = star(deepcopy(a),n)
                                    @inbounds bi = qtermadd4(d,b,a_star,c,n)
                                    if nb > 0
                                        bi = qreduce_unitnorm(bi, n;nb=nb)
                                    end
                                    if ipart != true
                                        bi= canonical_qmono(bi,n)
                                    else
                                        bi_cano= canonical_qmono(bi,n)
                                    end
                                    if nb > 0
                                        bi = qreduce_unitnorm(bi, n, nb=nb)
                                        if ipart
                                            bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                        end
                                    end
                                    Locb = bfind(tsupp, ltsupp, bi)
                                    if ipart
                                        Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                                        @inbounds add_to_expression!(rcons[Locb_cano], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss])
                                        @inbounds add_to_expression!(icons[Locb], imag_part(addcoe[rr][tt])[1], addepos[i][(rr-1)*length(addsupp)+ss])
                                        @inbounds add_to_expression!(jcons[Locb], imag_part(addcoe[rr][tt])[2], addepos[i][(rr-1)*length(addsupp)+ss])
                                        @inbounds add_to_expression!(kcons[Locb], imag_part(addcoe[rr][tt])[3], addepos[i][(rr-1)*length(addsupp)+ss])
                                    else
                                        @inbounds add_to_expression!(rcons[Locb], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss])
                                    end
                                end
                                for tt= 1:4
                                    a = deepcopy(addebasis[1])
                                    b = deepcopy(addebasis[1])
                                    c = deepcopy(addsupp[ss][tt][1])
                                    d = deepcopy(addsupp[ss][tt][2])
                                    a_star = star(deepcopy(a),n)
                                    c_star = star(deepcopy(c),n)
                                    d_star = star(deepcopy(d),n)
                                    @inbounds bi = qtermadd4(c_star,b,a_star,d_star,n)
                                    # println(addcoe[rr][tt],bi)
                                    if nb > 0
                                        bi = qreduce_unitnorm(bi, n;nb=nb)
                                    end
                                    if ipart != true
                                        bi= canonical_qmono(bi,n)
                                    else
                                        bi_cano= canonical_qmono(bi,n)
                                    end
                                    if nb > 0
                                        bi = qreduce_unitnorm(bi, n, nb=nb)
                                        if ipart
                                            bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                        end
                                    end
                                    Locb = bfind(tsupp, ltsupp, bi)
                                    if ipart
                                        Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                                        @inbounds add_to_expression!(rcons[Locb_cano], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss])
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(addcoe[rr][tt])[1], addepos[i][(rr-1)*length(addsupp)+ss])
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(addcoe[rr][tt])[2], addepos[i][(rr-1)*length(addsupp)+ss])
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(addcoe[rr][tt])[3], addepos[i][(rr-1)*length(addsupp)+ss])
                                    else
                                        @inbounds add_to_expression!(rcons[Locb], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss])
                                    end
                                end
                            end
                        end
                    else
                        for rr = 1:3
                            for ss = 1:length(addsupp)
                                addepos[i][(rr-1)*length(addsupp)+ss] = @variable(model, [1:bs, 1:bs])
                                for t = 1:bs, r = 1:bs
                                    for tt= 1:4
                                        a = deepcopy(addebasis[t])
                                        b = deepcopy(addebasis[r])
                                        c = deepcopy(addsupp[ss][tt][1])
                                        d = deepcopy(addsupp[ss][tt][2])
                                        a_star = star(deepcopy(a),n)
                                        @inbounds bi = qtermadd4(d,b,a_star,c,n)
                                        if nb > 0
                                            bi = qreduce_unitnorm(bi,n;nb=nb)
                                        end
                                        if ipart != true
                                            bi= canonical_qmono(bi,n)
                                        else
                                            bi_cano= canonical_qmono(bi,n)
                                        end
                                        if nb > 0
                                            bi = qreduce_unitnorm(bi, n, nb=nb)
                                            if ipart
                                                bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                            end
                                        end
                                        Locb = bfind(tsupp, ltsupp, bi)
                                        if ipart == true
                                            Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                                            @inbounds add_to_expression!(rcons[Locb_cano], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss][t,r])
                
                                            @inbounds add_to_expression!(icons[Locb], imag_part(addcoe[rr][tt])[1], addepos[i][(rr-1)*length(addsupp)+ss][t,r])
                                
                                            @inbounds add_to_expression!(jcons[Locb], imag_part(addcoe[rr][tt])[2], addepos[i][(rr-1)*length(addsupp)+ss][t,r])
                                        
                                            @inbounds add_to_expression!(kcons[Locb], imag_part(addcoe[rr][tt])[3], addepos[i][(rr-1)*length(addsupp)+ss][t,r])
                        
                                        else
                                            @inbounds add_to_expression!(rcons[Locb], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss][t,r])
                                        end
                                    end
                                end               
                                for t = 1:bs, r = 1:bs
                                    for tt = 1:4
                                        a = deepcopy(addebasis[t])
                                        b = deepcopy(addebasis[r])
                                        c = deepcopy(addsupp[ss][tt][1])
                                        d = deepcopy(addsupp[ss][tt][2])
                                        a_star = star(deepcopy(a),n)
                                        c_star = star(deepcopy(c),n)
                                        d_star = star(deepcopy(d),n)
                                        @inbounds bi = qtermadd4(c_star,b,a_star,d_star,n)
                                        if nb > 0
                                            bi = qreduce_unitnorm(bi,n;nb=nb)
                                        end
                                        if ipart != true
                                            bi= canonical_qmono(bi,n)
                                        else
                                            bi_cano= canonical_qmono(bi,n)
                                        end
                                        if nb > 0
                                            bi = qreduce_unitnorm(bi, n, nb=nb)
                                            if ipart
                                                bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                            end
                                        end
                                        Locb = bfind(tsupp, ltsupp, bi)
                                        if ipart == true
                                            Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                                            @inbounds add_to_expression!(rcons[Locb_cano], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss][r,t])
                
                                            @inbounds add_to_expression!(icons[Locb], -imag_part(addcoe[rr][tt])[1], addepos[i][(rr-1)*length(addsupp)+ss][r,t])
                                
                                            @inbounds add_to_expression!(jcons[Locb], -imag_part(addcoe[rr][tt])[2], addepos[i][(rr-1)*length(addsupp)+ss][r,t])
                                        
                                            @inbounds add_to_expression!(kcons[Locb], -imag_part(addcoe[rr][tt])[3], addepos[i][(rr-1)*length(addsupp)+ss][r,t])
                        
                                        else
                                            @inbounds add_to_expression!(rcons[Locb], real(addcoe[rr][tt]), addepos[i][(rr-1)*length(addsupp)+ss][r,t])
                                        end
                                    end  
                                end             
                            end
                        end
                    end
                end
            end
        end
        rbc = zeros(ltsupp)
        if ipart == true
            ibc = zeros(ltsupp)
            jbc = zeros(ltsupp)
            kbc = zeros(ltsupp)
        end
        ncons = ltsupp
        if QUIET == false
            println("There are $ncons affine constraints.")
        end
        for i = 1:length(supp[1])
            mono = supp[1][i]
            mono_cano = canonical_qmono(mono,n)
            if ipart != true
                mono = canonical_qmono(deepcopy(mono),n)
            end

            if nb > 0
                mono = qreduce_unitnorm(mono,n, nb=nb)
                mono_cano = qreduce_unitnorm(mono_cano,n, nb=nb)
            end

            Locb = bfind(tsupp, ltsupp, mono)
            if ipart == true
                Locb_cano = bfind(ttsupp, lttsupp, mono_cano)
            end

            if Locb === nothing || (ipart == true && Locb_cano === nothing)
                println(mono)
                @error "The monomial basis is not enough!"
                return nothing,ksupp,nothing,nothing,nothing
            else
                if ipart == true
                    rbc[Locb_cano] += real(coe[1][i])
                    ibc[Locb] += imag_part(coe[1][i])[1]
                    jbc[Locb] += imag_part(coe[1][i])[2]
                    kbc[Locb] += imag_part(coe[1][i])[3]
                else
                    rbc[Locb] += real(coe[1][i])
                end
            end
        end
        @variable(model, lower)
        rcons[1] += lower
        ###remove conjugate pairs
        if ipart == true
            keep = trues(ltsupp)
            visited = Set{Int}()

            for i = 1:ltsupp

                if i in visited
                    continue
                end

                mono = standardterm(tsupp[i], n)

                mono_star = star(deepcopy(mono), n)

                j = bfind(tsupp, ltsupp, mono_star)

                # 没找到共轭项
                if j === nothing
                    push!(visited, i)
                    continue
                end

                # self-adjoint monomial
                if j == i
                    push!(visited, i)
                    continue
                end

                push!(visited, i)
                push!(visited, j)

                keep[max(i,j)] = false
            end
            keep1 = trues(lttsupp)
            visited1 = Set{Int}()

            for i = 1:lttsupp

                if i in visited1
                    continue
                end

                mono = standardterm(ttsupp[i], n)

                mono_star = canonical_qmono(star(deepcopy(mono), n), n)

                j = bfind(ttsupp, lttsupp, mono_star)
                
                # 没找到共轭项
                if j === nothing
                    push!(visited1, i)
                    continue
                end

                # self-adjoint monomial
                if j == i
                    push!(visited1, i)
                    continue
                end

                push!(visited1, i)
                push!(visited1, j)

                # 两个都出现
                keep1[max(i,j)] = false

            end
            inds = findall(keep)
            inds1 = findall(keep1)
        else
            keep = trues(ltsupp)
            visited = Set{Int}()

            for i = 1:ltsupp

                if i in visited
                    continue
                end

                mono = standardterm(tsupp[i], n)

                mono_star = canonical_qmono(star(deepcopy(mono), n), n)

                j = bfind(tsupp, ltsupp, mono_star)
                # 没找到共轭项
                if j === nothing
                    push!(visited, i)
                    continue
                end

                # self-adjoint monomial
                if j == i
                    push!(visited, i)
                    continue
                end

                push!(visited, i)
                push!(visited, j)

                keep[max(i,j)] = false

            end
            inds = findall(keep)
        end
        if ipart == true || qncnumeq > 0
            @constraint(model,rcon[k=1:length(inds1)],rcons[inds1[k]] == rbc[inds1[k]])  
            @constraint(model, icon[k=1:length(inds)], icons[inds[k]] == ibc[inds[k]])
            @constraint(model, jcon[k=1:length(inds)], jcons[inds[k]] == jbc[inds[k]])
            @constraint(model, kcon[k=1:length(inds)], kcons[inds[k]] == kbc[inds[k]])
        else
            @constraint(model,rcon[k=1:length(inds)],rcons[inds[k]] == rbc[inds[k]])
        end
        
        @objective(model, Max, lower)
        end
        if QUIET == false
            println("SDP assembling time: $time seconds.")
            println("Solving the SDP...")
        end
        time = @elapsed begin
        optimize!(model)
        end
        if QUIET == false
            println("SDP solving time: $time seconds.")
        end
        SDP_status = termination_status(model)
        cons=all_constraints(model; include_variable_in_set_constraints = false)
        objv = objective_value(model)
        if SDP_status != MOI.OPTIMAL
            println("termination status: $SDP_status")
            status = primal_status(model)
            println("solution status: $status")
        end
        println("optimum = $objv")
    end
    return objv,ksupp,SDP_status
end

function add_clique!(G, nodes)
    for i in 1:length(nodes)-1, j in i+1:length(nodes)
        add_edge!(G, nodes[i], nodes[j])
    end
end

function max_cliques(G)
    cliques = convert(Vector{Vector{UInt16}}, maximal_cliques(G))
    sort!.(cliques)
    cliquesize = length.(cliques)
    cql = length(cliquesize)
    return cliques,cql,cliquesize
end

function chordal_cliques!(G; method="MF", minimize=false)
    # choose algorithm
    alg = method == "MF" && minimize ? CliqueTrees.MinimalChordal(CliqueTrees.MF())  :
          method == "MD" && minimize ? CliqueTrees.MinimalChordal(CliqueTrees.MMD()) :
          method == "MF" && !minimize ? CliqueTrees.MF()                              :
          method == "MD" && !minimize ? CliqueTrees.MMD()                             :
          error()

    # compute maximal cliques
    label, tree = CliqueTrees.cliquetree(G; alg)
    
    # triangulate graph
    F = CliqueTrees.FilledGraph(tree)
    
    for edge in edges(F)
        add_edge!(G, label[src(edge)], label[dst(edge)])
    end
    
    # return maximal cliques
    maximal_cliques = Vector{Vector{UInt16}}(undef, length(tree))
    
    for (i, clique) in enumerate(tree)
        maximal_cliques[i] = sort!(label[clique])
    end
    
    cliquesize = length.(maximal_cliques)
    cql = length(cliquesize)
    return maximal_cliques, cql, cliquesize
end

function quaternion_to_real(qpop,q)
    n = Int(length(q)/2)
    @polyvar x[1:4n]
    pop = Vector{Poly}(undef, length(qpop))
    I=quat(0,1,0,0)
    J=quat(0,0,1,0)
    K=quat(0,0,0,1)
    for (i,qp) in enumerate(qpop)
        temp = qp(q[1:n]=>x[1:n]+I*x[n+1:2n]+J*x[2n+1:3n]+K*x[3n+1:4n], q[n+1:2n]=>x[1:n]-I*x[n+1:2n]-J*x[2n+1:3n]-K*x[3n+1:4n])
        # println(temp)
        pop[i] = real.(MultivariatePolynomials.coefficients(temp))'*MultivariatePolynomials.monomials(temp)
    end
    return pop,x
end

function clique_decomp(n, m, dc, rncnumeq, supp::Vector{Vector{Vector{Vector{UInt16}}}}; order="min", alg="MF")
    if alg == false
        cliques = [UInt16[i for i=1:n]]
        cql = 1
        cliquesize = [n]
    else
        G = SimpleGraph(n)
        for i = 1:m+1-rncnumeq
            if i == 1 || order == Int(ceil(dc[i-1]/2))
            # if i == 1
                for j = 1:length(supp[i])
                    # println(unique([supp[i][j][1];supp[i][j][2];supp[i][j][3]]))
                    temp1 = copy(supp[i][j][3])
                    temp2 = copy(supp[i][j][3])
                    # println(temp)
                    # temp1 .=  [ x > UInt16(n/2) ? x - UInt16(n/2) : x for x in temp1]
                    temp1 .=  [ x > UInt16(n) ? x - UInt16(n) : x for x in temp1]
                    # temp2 .=  [ x <= UInt16(n/2) ? x + UInt16(n/2) : x for x in temp2]
                    # println(temp)
                    # add_clique!(G, unique([supp[i][j][1];supp[i][j][2];temp1;temp2;supp[i][j][3]]))
                    add_clique!(G, unique([supp[i][j][1];temp1]))
                    # add_clique!(G, unique([supp[i][j][1];supp[i][j][2];supp[i][j][3]]))
                end
            else
                temp1 = copy(supp[i][1][3])
                temp1 .=  [ x > UInt16(n) ? x - UInt16(n) : x for x in temp1]
                # temp = copy([supp[i][1][1];supp[i][1][2];supp[i][1][3]])
                temp = copy([supp[i][1][1];temp1])
                for j = 2:length(supp[i])
                    temp2 = copy(supp[i][j][3])
                    temp2 .=  [ x > UInt16(n) ? x - UInt16(n) : x for x in temp2]
                    # append!(temp, [supp[i][j][1];supp[i][j][2];supp[i][j][3]])
                    append!(temp, [supp[i][j][1];temp2])
                end
                add_clique!(G, unique(temp))
            end
        end
        if alg == "NC"
            cliques,cql,cliquesize = max_cliques(G)
        else
            cliques,cql,cliquesize = chordal_cliques!(G, method=alg, minimize=true)
        end
    end
    uc = unique(cliquesize)
    sizes = [sum(cliquesize.== i) for i in uc]
    println("-----------------------------------------------------------------------------")
    println("The clique sizes of varibles:\n$uc\n$sizes")
    println("-----------------------------------------------------------------------------")
    return cliques,cql,cliquesize
end

function get_blocks(m, l, d, tsupp, supp::Vector{Vector{Vector{Vector{UInt16}}}}, basis, ebasis; nb=0, normality=0, nvar=0, TS="block", ConjugateBasis=false, merge=false, md=3)
    # if (ConjugateBasis == false && normality > 0) || (ConjugateBasis == true && normality >= d)
    if normality > 0
        uk = m + 1 + nvar
    else
        uk = m + 1 
    end
    blocks = Vector{Vector{Vector{Int}}}(undef, uk)
    blocksize = Vector{Vector{Int}}(undef, uk)
    cl = Vector{Int}(undef, uk)
    eblocks = Vector{Vector{Int}}(undef, l)
    if TS == false
        for k = 1:uk
            lb = length(basis[k])
            blocks[k],blocksize[k],cl[k] = [Vector(1:lb)],[lb],1
        end
        for k = 1:l
            eblocks[k] = Vector(1:length(ebasis[k]))
        end
    else
        for k = 1:uk
            if k == 1
                G = get_graph(tsupp, basis[1], nb=nb, ConjugateBasis=ConjugateBasis)
            # elseif ((ConjugateBasis == false && normality > 0) || (ConjugateBasis == true && normality >= d)) && k <= 1 + nvar
            elseif  normality > 0 && k <= 1 + nvar
                G = get_graph(tsupp, basis[k], nb=nb)
            elseif normality > 0 && k > 1 + nvar
                G = get_graph(tsupp, supp[k-1-nvar], basis[k], nb=nb)
            else
                G = get_graph(tsupp, supp[k-1], basis[k], nb=nb)
            end
            if TS == "block"
                blocks[k] = connected_components(G)
                blocksize[k] = length.(blocks[k])
                cl[k] = length(blocksize[k])
            else
                blocks[k],cl[k],blocksize[k] = chordal_cliques!(G, method=TS)
                if merge == true
                    blocks[k],cl[k],blocksize[k] = clique_merge!(blocks[k], d=md, QUIET=true)
                end
            end
        end
        for k = 1:l
            eblocks[k] = get_eblock(tsupp, supp[k+m], ebasis[k], nb=nb)
        end
    end
    return blocks,eblocks,cl,blocksize
end

function get_blocks(n, rlorder, I, J, supp::Vector{Vector{Vector{Vector{UInt16}}}}, cliques, cliquesize, cql, tsupp, basis, ebasis; TS="block", ConjugateBasis=false, nb=0, normality=1, merge=false, md=3)
    blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, cql)
    eblocks = Vector{Vector{Vector{Int}}}(undef, cql)
    cl = Vector{Vector{Int}}(undef, cql)
    blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
    for i = 1:cql
        # temp_tsupp = deepcopy(tsupp)
        # for j in eachindex(temp_tsupp)
        #     temp_tsupp[j][3] = [x > UInt16(n) ? x - UInt16(n) : x for x in temp_tsupp[j][3]]
        # end
        # ksupp = TS ? tsupp[[issubset(union(tsupp[j][1], tsupp[j][2], temp_tsupp[j][3]), cliques[i]) for j = 1: length(tsupp)]] : nothing
        ksupp = TS == false ? nothing : tsupp[[
            let 
                modified_3rd_l = [x > UInt16(n) ? x - UInt16(n) : x for x in tsupp[j][3]]
                # modified_3rd_g = [x <= UInt16(n) ? x + UInt16(n) : x for x in tsupp[j][3]]
                # issubset(union(tsupp[j][1], tsupp[j][2], modified_3rd_l, modified_3rd_g), cliques[i])
                issubset(union(tsupp[j][1], modified_3rd_l), cliques[i])
            end 
            for j in eachindex(tsupp)
        ]]
        # ksupp = TS == false ? nothing : tsupp[[issubset(union(tsupp[j][1], tsupp[j][2], tsupp[j][3].=[x <= UInt16(n) ? x - UInt16(n) : x for x in tsupp[j][3]]), cliques[i]) for j in eachindex(tsupp)]]
        # println(ksupp)
        blocks[i],eblocks[i],cl[i],blocksize[i] = get_blocks(length(I[i]), length(J[i]), rlorder[i], ksupp, supp[[I[i]; J[i]].+1], basis[i], 
        ebasis[i], TS=TS, ConjugateBasis=ConjugateBasis, merge=merge, md=md, nb=nb, normality=normality, nvar=cliquesize[i])
    end
    return blocks,eblocks,cl,blocksize
end

function assign_constraint(n,m, numeq, supp::Vector{Vector{Vector{Vector{UInt16}}}}, cliques, cql)
    I = [Int[] for i=1:cql]
    J = [Int[] for i=1:cql]
    ncc = Int[]
    for i = 1:m
        temp1 = copy(supp[i+1][1][3])
        temp1 .=  [ x > UInt16(n) ? x - UInt16(n) : x for x in temp1]
        temp = copy([supp[i][1][1];temp1])
        for j = 2:length(supp[i+1])
            temp2 = copy(supp[i+1][j][3])
            temp2 .=  [ x > UInt16(n) ? x - UInt16(n) : x for x in temp2]
            append!(temp, [supp[i+1][j][1];temp2])
        end
        ind = findall(k->issubset(unique(temp), cliques[k]), 1:cql)
        # ind = findall(k->issubset(unique(reduce(vcat, [[item[1];item[2];item[3]] for item in supp[i+1]])), cliques[k]), 1:cql)
        if isempty(ind)
            push!(ncc, i)
        elseif i <= m - numeq
            push!.(I[ind], i)
        else
            push!.(J[ind], i)
        end
    end
    return I,J,ncc
end

function get_graph(tsupp::Vector{Vector{Vector{UInt16}}}, basis; nb=0, ConjugateBasis=false)
    lb = length(basis)
    G = SimpleGraph(lb)
    ltsupp = length(tsupp)
    for i = 1:lb, j = 1:lb
        a=deepcopy(basis[i])
        b=deepcopy(basis[j])
        bi = qtermadd(a,b,n)
        if nb > 0
            bi = qreduce_unitnorm(bi, nb=nb)
        end
        if bfind(tsupp, ltsupp, bi) !== nothing
            add_edge!(G, i, j)
        end
    end
    return G
end

function get_graph(tsupp::Vector{Vector{Vector{UInt16}}}, supp, basis; nb=0, ConjugateBasis=false)
    lb = length(basis)
    ltsupp = length(tsupp)
    G = SimpleGraph(lb)
    for i = 1:lb, j = 1:lb
        r = 1
        while r <= length(supp)
                a=deepcopy(basis[i])
                b=deepcopy(supp[r])
                c=deepcopy(basis[j])
                bi = qtermadd3(a,b,c,n)
            if nb > 0
                bi = qreduce_unitnorm(bi, nb=nb)
            end
            if bfind(tsupp, ltsupp, bi) !== nothing
                break
            else
                r += 1
            end
        end
        if r <= length(supp)
            add_edge!(G, i, j)
        end
    end
    return G
end

function get_eblock(tsupp::Vector{Vector{Vector{UInt16}}}, hsupp::Vector{Vector{Vector{UInt16}}}, basis::Vector{Vector{Vector{UInt16}}}; nb=nb)
    ltsupp = length(tsupp)
    eblock = Int[]
    lb = length(basis)
    for i = 1:lb, j = 1:lb
        flag = 0
        for temp in hsupp
            a=deepcopy(basis[i])
            b=deepcopy(temp)
            c=deepcopy(basis[j])
            bi = qtermadd3(a,b,c,n)
            if nb > 0
                bi = qreduce_unitnorm(bi, nb=nb)
            end
            if bfind(tsupp, ltsupp, bi) !== nothing
                flag = 1
                break
            end
        end
        if flag == 1
            if i == j
                push!(eblock,i)
            else
                append!(eblock, [i,j])
            end
        end
    end
    return eblock
end

function get_gsupp(rlorder, basis, ebasis, supp, cql, I, J, ncc, blocks, eblocks, cl, blocksize, cliquesize; norm=false, nb=0, ConjugateBasis=false, normality=1)
    if norm == true
        gsupp = Vector{UInt16}[]
    else
        gsupp = Vector{Vector{UInt16}}[]
    end
    for i = 1:cql
        for (j, w) in enumerate(I[i]), l = 1:cl[i][j+1], t = 1:blocksize[i][j+1][l], r = 1:blocksize[i][j+1][l], item in supp[w+1]
            ind1 = blocks[i][j+1][l][t]
            ind2 = blocks[i][j+1][l][r]
            a=deepcopy(basis[i][j+1][ind1])
            b=deepcopy(item)
            c=deepcopy(basis[i][j+1][ind2])
            @inbounds bi = qtermadd3(a,b,c,n)
            if nb > 0
                bi = qreduce_unitnorm(bi;nb=nb)
            end
            push!(gsupp, bi)
        end
        for (j, w) in enumerate(J[i]), k = 1: length(eblocks[i][j]),r = 1: length(eblocks[i][j]), item in supp[w+1]
            a=deepcopy(ebasis[i][j][eblocks[i][j][k]])
            b=deepcopy(item)
            c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
            @inbounds bi = qtermadd3(a,b,c,n)
            if nb > 0
                bi = qreduce_unitnorm(bi, nb=nb)
            end
            push!(gsupp, bi)
        end
    end
    for i in ncc, j = 1:length(supp[i+1])
        @inbounds bi = supp[i+1][j]
        if nb > 0
            bi = qreduce_unitnorm(bi, nb=nb)
        end
        push!(gsupp, bi)
    end
    return gsupp
end
