const Poly = DP.Polynomial{DP.Commutative{DP.CreationOrder}, Graded{LexOrder}, Float64}
const NCPoly = DP.Polynomial{DP.NonCommutative{DP.CreationOrder}, Graded{LexOrder}, Float64}
const NCMono = DP.Monomial{DP.NonCommutative{DP.CreationOrder}, Graded{LexOrder}}
using MathOptInterface
# const MOI = MathOptInterface

function qs_tssos_first(pop, z, n::Int, d; fsupp=nothing, fcoe=nothing, numeq=0, rncnumeq=0, qncnumeq=0, RemSig=false, nb=0,CS="MF",cliques=[],TS="block",
    merge=false, md=3, solver="Mosek", reducebasis=false, QUIET=false, solve=true, solution=false, ipart=true,
    dualize=false, balanced=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), 
    writetofile=false, normality=0, NormalSparse=false, conjubasis=true)
    supp,coe = QPolys_info(pop, z, n)
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
        cliques,cql,cliquesize = clique_decomp(n, m, dc, supp, order=d, alg=CS)
        end
        if CS != false && QUIET == false
            mc = maximum(cliquesize)
            println("Obtained the variable cliques in $time seconds.\nThe maximal size of cliques is $mc.")
        end
    end
    I,J,ncc = assign_constraint(n, m, numeq, supp, cliques, cql)
    # I,J,NJ,ncc = assign_constraint2(n, m, m1, numeq, supp, cliques, cql)
    rlorder = d*ones(Int, cql)
    # basis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, m+1)
    # basis[1] = get_qncbasis(n, d; conjubasis=conjubasis)
    ebasis = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, cql)
    # nebasis = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, cql)
    basis = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, cql)
    # if QUIET == false
    #     lb = length(basis[1])
    #     println("Maximal PSD block size: $lb.")
    # end
    # for s = 1:m
    #     basis[s+1] = get_qncbasis(n, d-Int(ceil(dc[s]/2));conjubasis=conjubasis)
    # end
        # if nb > 0
    #     for s = 1:m+1
    #         basis[s] = qreduce_unitnorm.(basis[s], nb=nb)
    #         unique!(basis[s])
    #     end
    # end
    # println(supp[1])
    for i = 1:cql
        ebasis[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(J[i]))
        # nebasis[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(NJ[i]))
        if normality > 0
            # if normality < d
            #     basis[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(I[i])+1)
            #     basis[i][1] = get_qncbasis(cliques[i], n, d; conjubasis=conjubasis)
            #     if nb > 0
            #         basis[i][1] = qreduce_unitnorm.(basis[i][s+1], nb=nb)
            #         unique!(basis[i][1])
            #     end
            #     for s = 1:length(I[i])
            #         basis[i][s+1] = get_qncbasis(cliques[i], n, rlorder[i]-Int(ceil(dc[I[i][s]]/2)))
            #         if nb > 0
            #             basis[i][s+1] = qreduce_unitnorm.(basis[i][s+1], nb=nb)
            #             unique!(basis[i][s+1])
            #         end
            #     end
            # else
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
            # end
        else
            basis[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(I[i])+1)
            basis[i][1] = get_qncbasis(cliques[i], n, d; conjubasis=conjubasis)
            # println(length(basis[i][1]))
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
        # for s = 1:length(NJ[i])
        #     nebasis[i][s] = get_qncbasis(cliques[i], n, rlorder[i]-Int(ceil(dc[NJ[i][s]]/2)))
        #     if nb > 0
        #         nebasis[i][s] = qreduce_unitnorm.(nebasis[i][s], n, nb=nb)
        #         unique!(nebasis[i][s])
        #     end
        # end
    end
    #     basis[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(I[i])+1)
    #     # println(cliques[i])
    #     basis[i][1] = get_qncbasis(cliques[i], n, d; conjubasis=conjubasis)
    #     # println(basis[i][1])
    #     for s = 1:length(I[i])
    #         basis[i][s+1] = get_qncbasis(cliques[i], n, rlorder[i]-Int(ceil(dc[I[i][s]]/2)))
    #         if nb > 0
    #             basis[i][s+1] = qreduce_unitnorm.(basis[i][s+1], nb=nb)
    #             unique!(basis[i][s+1])
    #         end
    #     end
    #     for s = 1:length(J[i])
    #         ebasis[i][s] = get_qncbasis(cliques[i], n, rlorder[i]-Int(ceil(dc[J[i][s]]/2)))
    #         if nb > 0
    #             ebasis[i][s] = qreduce_unitnorm.(ebasis[i][s+1], nb=nb)
    #             unique!(ebasis[i][s])
    #         end
    #     end
    # end
    # tsupp = canonical_qmono.(deepcopy(supp[1]),n)
    tsupp = deepcopy(supp[1])
    # println(tsupp)
    # for k = 1:length(tsupp)
    #     if tsupp[k] != canonical_qmono(tsupp[k],n)
    #         println("add 1st:",tsupp[k],canonical_qmono(tsupp[k],n))
    #     end
    # end
    for i = 2:m+1, j = 1:length(supp[i])
        # push!(tsupp, canonical_qmono(supp[i][j],n))
        push!(tsupp, supp[i][j])
    end
    sort!(tsupp)
    unique!(tsupp)
    # for k = 1:length(tsupp)
    #     if tsupp[k] != canonical_qmono(tsupp[k],n)
    #         println("add 1st:",tsupp[k],canonical_qmono(tsupp[k],n))
    #     end
    # end
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,eblocks,cl,blocksize = get_blocks(n,rlorder, I, J, supp, cliques, cliquesize, cql, tsupp, basis, ebasis, TS=TS, ConjugateBasis=conjubasis, normality=normality, merge=merge, md=md, nb=nb)
    if QUIET == false
        mb = maximum(maximum.([maximum.(blocksize[i]) for i = 1:cql]))
        println("Obtained the block structure in $time seconds.\nThe maximal size of blocks is $mb.")
    end
    if solution
        opt,ksupp,SDP_status,sol,moment = qsolvesdp(n, m, rncnumeq, qncnumeq, rlorder, supp, coe, basis, ebasis, cliques, cql, cliquesize, I, J, ncc, blocks, eblocks, cl, blocksize, numeq=numeq, QUIET=QUIET, TS=TS, solver=solver, solve=solve, solution=solution, ipart=ipart, balanced=balanced,
        nb=nb, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, dualize=dualize, writetofile=writetofile, normality=normality, NormalSparse=NormalSparse,conjubasis=conjubasis)
        return opt,sol,moment
    else
        opt, ksupp, SDP_status =  qsolvesdp(n, m, rncnumeq, qncnumeq, rlorder, supp, coe, basis, ebasis, cliques, cql, cliquesize, I, J, ncc, blocks, eblocks, cl, blocksize, numeq=numeq, QUIET=QUIET, TS=TS, solver=solver, solve=solve, solution=solution, ipart=ipart, balanced=balanced,
        nb=nb, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, dualize=dualize, writetofile=writetofile, normality=normality, NormalSparse=NormalSparse,conjubasis=conjubasis)
        return opt
    end
    # opt,ksupp,SDP_status,sol,moment = qsolvesdp(n, m, rlorder, supp, coe, basis, ebasis, cliques, cql, cliquesize, I, J, ncc, blocks, eblocks, cl, blocksize, numeq=numeq, QUIET=QUIET, TS=TS, solver=solver, solve=solve, solution=solution, ipart=ipart, balanced=balanced,
    # nb=nb, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, dualize=dualize, writetofile=writetofile, normality=normality, NormalSparse=NormalSparse,conjubasis=conjubasis)
    # return opt,sol,moment
end
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
    numeq=0, nb=0, QUIET=false, TS=false, solver="Mosek", solve=true, dualize=false, solution=false, ipart=true, 
    cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false, balanced=false, normality=0, NormalSparse=false,conjubasis=false)
    tsupp = Vector{Vector{UInt16}}[]
    ttsupp = Vector{Vector{UInt16}}[]
    checktsupp = Vector{Vector{UInt16}}[]
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
            # @inbounds bi = qtermadd(b,at,n)
            if nb > 0
                bi = qreduce_unitnorm(bi,n, nb=nb)
            end
            push!(checktsupp,bi)
            # ====== NEW: canonical mapping ======
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
            # println("ebasis",ebasis[i][j])
            if m - rncnumeq < J[i][j] <= m - qncnumeq
                bs = length(eblocks[i][j])
                # println(bs)
                for t = 1:bs, r = 1:bs
                # for t = 1:bs-1, r = t+1:bs
                    for s = 1:length(supp[k+1])
                        a=deepcopy(ebasis[i][j][t])
                        b=deepcopy(supp[k+1][s])
                        # c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                        c=deepcopy(ebasis[i][j][r])
                        if ipart
                            @inbounds bi = qtermadd3(b,c,a,n)
                        else
                            @inbounds bi = qtermadd3left(a,b,c,n)
                            # @inbounds bi = qtermadd3left(a,c,b,n)
                        end
                                    #@inbounds bi = qtermadd3(a,b,c,n)
                                    # @inbounds bi = qtermadd3(b,c,a,n)
                        if nb > 0
                            bi = qreduce_unitnorm(bi,n;nb=nb)
                        end
                        # if bi == Vector{UInt16}[[], [], [0x0001, 0x0001, 0x0004]]
                        #     println("a",ebasis[i][j][eblocks[i][j][t]],"b",b,"c",c)
                        # end
                        push!(checktsupp,bi)
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
                            # push!(tsupp,bi_cano)
                        end
                    end
                    for s = 1:length(supp[k+1])
                        # a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                        a=deepcopy(ebasis[i][j][t])
                        b=deepcopy(supp[k+1][s])
                        # c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                        c=deepcopy(ebasis[i][j][r])
                        a_star = star(deepcopy(a),n)
                        b_star = star(deepcopy(b),n)
                        if ipart
                            # @inbounds bi = qtermadd3(b_star,c,a,n)
                            @inbounds bi = qtermadd3(c,a_star,b,n)
                        else
                            @inbounds bi = qtermadd3left(a,b_star,c,n)
                            # @inbounds bi = qtermadd3left(a,c,b,n)
                        end
                                    #@inbounds bi = qtermadd3(a,b,c,n)
                                    # @inbounds bi = qtermadd3(b,c,a,n)
                        if nb > 0
                            bi = qreduce_unitnorm(bi,n;nb=nb)
                        end
                        # if bi == Vector{UInt16}[[], [], [0x0001, 0x0001, 0x0004]]
                        #     println("a",a,"b",b_star,"c",c)
                        # end
                        push!(checktsupp,bi)
                        if ipart != true
                            bi= canonical_qmono(bi,n)
                        else
                            bi_cano = canonical_qmono(bi,n)
                            if bi_cano != canonical_qmono(deepcopy(bi_cano),n)
                                println(bi_cano)
                            end 
                        #     bi_star = star(deepcopy(bi),n)
                        #     if bi_star == bi
                        #         bi = canonical_qmono(bi,n)
                        #     end
                        end
                        # bi_canon = canonical_qmono(bi,n)
                        if nb > 0
                            bi= qreduce_unitnorm(bi,n, nb=nb)
                            if ipart 
                                bi_cano= qreduce_unitnorm(bi_cano,n, nb=nb)
                            end
                        end
                        push!(tsupp, bi)
                        if ipart
                            push!(ttsupp,bi_cano)
                            # push!(tsupp,bi_cano)
                        end
                    end
                end
            elseif J[i][j] > m - qncnumeq
                bs = length(eblocks[i][j])
                for t = 1:bs, r = 1:bs
                # for t = 1:bs-1, r = t+1:bs
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
                        # else
                        #     bi_star = star(deepcopy(bi),n)
                        #     if bi_star == bi
                        #         bi = canonical_qmono(bi,n)
                        #     end
                        end
                        # bi_canon = canonical_qmono(bi,n)
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
                        # else
                        #     bi_star = star(deepcopy(bi),n)
                        #     if bi_star == bi
                        #         bi = canonical_qmono(bi,n)
                        #     end
                        end
                        # bi_canon = canonical_qmono(bi,n)
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
    # if solution == true && TS != false
    #     ksupp = deepcopy(tsupp)
    #     for i = 1:cql, j = 1:cliquesize[i]
    #         push!(tsupp, [UInt16[], UInt16[], UInt16[cliques[i][j]]])
    #         for k = 1:cliquesize[i]
    #             if ipart
    #                 bi = qtermadd([UInt16[], UInt16[], UInt16[cliques[i][k]]],[UInt16[], UInt16[], UInt16[cliques[i][j]]],n)
    #             else
    #                 bi = qtermaddleft([UInt16[], UInt16[], UInt16[cliques[i][j]]],[UInt16[], UInt16[], UInt16[cliques[i][k]]],n)
    #             end
    #             #bi = qtermadd([UInt16[], UInt16[], UInt16[cliques[i][j]]],[UInt16[], UInt16[], UInt16[cliques[i][k]]],n)
    #             # bi = qtermadd([UInt16[], UInt16[], UInt16[cliques[i][k]]],[UInt16[], UInt16[], UInt16[cliques[i][j]]],n)
    #             if nb > 0
    #                 bi = qreduce_unitnorm(bi, n, nb=nb)
    #             end
    #             push!(tsupp, bi)
    #         end
    #     end
    # end
    # tsupp = canonical_qmono.(tsupp,n)
    sort!(tsupp)
    unique!(tsupp)
    sort!(ttsupp)
    unique!(ttsupp)
    if tsupp == canonical_qmono.(tsupp,n)
        println(" tsupp all cano")
    end
    # sort!(checktsupp)
    # unique!(checktsupp)
    # found = true
    # k=1
    # for i = 1:length(checktsupp)
    #     term = checktsupp[i]
    #     term_star = star(deepcopy(term), n)
    #     if term_star in checktsupp
    #         k+=1
    #     else
    #         println(checktsupp[i])
    #     end
    #     if k == length(checktsupp)
    #         println("checktsupp true")
    #     end
    # end
    # checkttsupp = ttsupp[setdiff(1:end, [4,8,9,10,30,33,41,42])]
    # for i = 1:length(checkttsupp)
    #     term = checkttsupp[i]
    #     term_star = canonical_qmono(star(deepcopy(term), n),n)
    #     if term_star in checkttsupp
    #         k+=1
    #     else
    #         println(checkttsupp[i])
    #     end
    #     println("term",term,"canoconj",term_star)
    #     if k == length(checkttsupp)
    #         println("checkttsupp true")
    #     end
    # end
    # for k = 1:length(ttsupp)
    #     if ttsupp[k] == star(deepcopy(ttsupp[k]), n)
    #         println("自共轭",ttsupp[k], " at index ", k)
    #     end
    # end

        # found = false

        # for i = 1:length(useless)

        #     mono = standardterm(tsupp[useless[i]], n)

        #     # 共轭
        #     mono_star = canonical_qmono(star(deepcopy(mono), n),n)
        #     # mono_star = star(deepcopy(mono), n)

        #     # 在 tsupp[inds] 里找
        #     for j = 1:length(useless)

        #         if i == j
        #             continue
        #         end

        #         if tsupp[useless[j]]== mono_star

        #             println("conjugate pair found:")
        #             println("i = ", useless[i], " : ", mono)
        #             println("j = ", useless[j], " : ", tsupp[useless[j]])
        #             println()

        #             found = true
        #         end
        #     end
        # end

        # if !found
        #     println("in remove No conjugate pairs remain.")
        # end
        
    # if !ipart
    #     tsupp = canonical_qmono.(tsupp,n)
    # end
    # for k = 1:length(tsupp)
    #     if tsupp[k] != canonical_qmono(tsupp[k],n)
    #         println(tsupp[k])
    #         println(canonical_qmono(tsupp[k],n))
    #     end
    # end
    # end
    if solution == true && TS != false
        sort!(ksupp)
        unique!(ksupp)
    else
        ksupp = tsupp
    end
    # ksupp = deepcopy(tsupp)
    # sort!(tsupp)
    # unique!(tsupp)
    #initial set
    objv = SDP_status= nothing
    if solve == true
        ltsupp = length(tsupp)
        lttsupp = length(ttsupp)
        if QUIET == false
            println("Assembling the SDP...")
        end
        if solver == "Mosek"
            if dualize == false
                model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => mosek_setting.tol_pfeas, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => mosek_setting.tol_dfeas, 
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => mosek_setting.tol_relgap, "MSK_DPAR_OPTIMIZER_MAX_TIME" => mosek_setting.time_limit, "MSK_IPAR_NUM_THREADS" => mosek_setting.num_threads))
            else
                model = Model(dual_optimizer(Mosek.Optimizer))
            end
        elseif solver == "COSMO"
            model = Model(optimizer_with_attributes(COSMO.Optimizer, "eps_abs" => cosmo_setting.eps_abs, "eps_rel" => cosmo_setting.eps_rel, "max_iter" => cosmo_setting.max_iter, "time_limit" => cosmo_setting.time_limit))
        elseif solver == "SDPT3"
            model = Model(optimizer_with_attributes(SDPT3.Optimizer))
        elseif solver == "SDPNAL"
            model = Model(optimizer_with_attributes(SDPNAL.Optimizer))
        else
            @error "The solver is currently not supported!"
            return nothing,nothing,nothing,nothing,nothing
        end
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        time = @elapsed begin
        # rcons = [AffExpr(0) for i=1:lttsupp]
        rcons = [AffExpr(0) for i=1:ltsupp]
        if ipart==true || qncnumeq > 0
            icons = [AffExpr(0) for i=1:ltsupp]
            jcons = [AffExpr(0) for i=1:ltsupp]
            kcons = [AffExpr(0) for i=1:ltsupp]
        end
        hnom = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, n)
        pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, cql)
        for i = 1:cql
            # if solution == true && TS != false
            #     bs = cliquesize[i] + 1
            #     if ipart == true
            #         pos0 = @variable(model, [1:4bs, 1:4bs], PSD)
            #     else
            #         pos0 = @variable(model, [1:bs, 1:bs], PSD)
            #     end
            #     for t = 1:bs, r = 1:bs
            #         if t == 1 && r == 1
            #             bi = [UInt16[], UInt16[],UInt16[]]
            #         elseif t == 1 && r > 1
            #             bi = [UInt16[],UInt16[], UInt16[cliques[i][r-1]]]
            #         elseif t > 1 && r == 1
            #             bi = [UInt16[],UInt16[], UInt16[cliques[i][t-1]]]
            #         else
            #             if ipart
            #                 bi = qtermadd([UInt16[],UInt16[], UInt16[cliques[i][r-1]]],[UInt16[],UInt16[], UInt16[cliques[i][t-1]]],n)
            #             else
            #                 bi = qtermaddleft([UInt16[],UInt16[], UInt16[cliques[i][t-1]]], [UInt16[],UInt16[], UInt16[cliques[i][r-1]]],n)
            #             end
            #             #bi = qtermadd([UInt16[],UInt16[], UInt16[cliques[i][t-1]]], [UInt16[],UInt16[], UInt16[cliques[i][r-1]]],n)
            #             # bi = qtermadd([UInt16[],UInt16[], UInt16[cliques[i][r-1]]],[UInt16[],UInt16[], UInt16[cliques[i][t-1]]],n)
            #             if nb > 0
            #                 bi = qreduce_unitnorm(bi, n, nb=nb)
            #             end
            #         end
            #         # Locb = bfind(tsupp, ltsupp, bi)
            #         # if ipart != true
            #         #     bi= canonical_qmono(bi,n)
            #         # else
            #         #     bi_star = star(deepcopy(bi),n)
            #         #     if bi_star == bi
            #         #         bi = canonical_qmono(bi,n)
            #         #     end
            #         # end
            #         bi_cano = canonical_qmono(bi,n)
            #         if nb > 0
            #             bi = qreduce_unitnorm(bi, n, nb=nb)
            #             bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
            #         end
            #         Locb = bfind(tsupp, ltsupp, bi)
            #         Locb_cano = bfind(tsupp, ltsupp, bi_cano)
            #         if ipart == true
            #             # @inbounds add_to_expression!(rcons[Locb], pos0[t,r]+pos0[t+bs,r+bs]+pos0[t+2*bs,r+2*bs]+pos0[t+3*bs,r+3*bs])
            #             @inbounds add_to_expression!(rcons[Locb_cano], pos0[t,r]+pos0[t+bs,r+bs]+pos0[t+2*bs,r+2*bs]+pos0[t+3*bs,r+3*bs])
            #             @inbounds add_to_expression!(icons[Locb], pos0[t+bs,r]-pos0[t,r+bs]+pos0[t+3*bs,r+2*bs]-pos0[t+2*bs,r+3*bs]) 
            #             @inbounds add_to_expression!(jcons[Locb], pos0[t+2*bs,r]-pos0[t,r+2*bs]-pos0[t+3*bs,r+bs]+pos0[t+bs,r+3*bs])
            #             @inbounds add_to_expression!(kcons[Locb], pos0[t+3*bs,r]-pos0[t,r+3*bs]+pos0[t+2*bs,r+bs]-pos0[t+bs,r+2*bs])
            #         else
            #             # @inbounds add_to_expression!(rcons[Locb], pos0[t,r])
            #             @inbounds add_to_expression!(rcons[Locb_cano], pos0[t,r])

            #         end            
            #     end
            # end
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
                        #@inbounds bi = qtermadd(at, b,n)
                        # @inbounds bi = qtermadd(b,at,n)
                        if nb > 0
                            bi = qreduce_unitnorm(bi,n;nb=nb)
                        end
                        # if ipart != true
                        #     bi= canonical_qmono(bi,n)
                        # else
                        #     bi_star = star(deepcopy(bi),n)
                        #     if bi_star == bi
                        #         bi = canonical_qmono(bi,n)
                        #     end
                        # end
                        bi_cano= canonical_qmono(bi,n)
                        if bi_cano != canonical_qmono(deepcopy(bi_cano),n)
                            println(bi_cano)
                        end 
                        if nb > 0
                            bi = qreduce_unitnorm(bi, n, nb=nb)
                            bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                        end
                        # bi = canonical_qmono(bi,n)
                        Locb = bfind(tsupp, ltsupp, bi)
                        # @inbounds add_to_expression!(rcons[Locb], pos[i][j][l])
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
                            #@inbounds bi = qtermadd(at,b,n)
                            # @inbounds bi = qtermadd(b,at,n)
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
                                bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                            end
                            # bi = canonical_qmono(bi,n)
                            Locb = bfind(tsupp, ltsupp, bi)
                            if ipart == true
                                Locb_cano = bfind(ttsupp, lttsupp, bi_cano)
                                @inbounds add_to_expression!(rcons[Locb_cano], pos[i][j][l][t,r]+pos[i][j][l][t+bs,r+bs]+pos[i][j][l][t+2*bs,r+2*bs]+pos[i][j][l][t+3*bs,r+3*bs])
                                @inbounds add_to_expression!(icons[Locb], pos[i][j][l][t+bs,r]-pos[i][j][l][t,r+bs]+pos[i][j][l][t+3*bs,r+2*bs]-pos[i][j][l][t+2*bs,r+3*bs]) 
                                @inbounds add_to_expression!(jcons[Locb], pos[i][j][l][t+2*bs,r]-pos[i][j][l][t,r+2*bs]-pos[i][j][l][t+3*bs,r+bs]+pos[i][j][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], pos[i][j][l][t+3*bs,r]-pos[i][j][l][t,r+3*bs]+pos[i][j][l][t+2*bs,r+bs]-pos[i][j][l][t+bs,r+2*bs])
                            else
                                # @inbounds add_to_expression!(rcons[Locb], pos[i][j][l][t,r])
                                @inbounds add_to_expression!(rcons[Locb], pos[i][j][l][t,r])
                            end            
                        end
                    end
                end
            end
        end
        for i = 1:cql, (j, k) in enumerate(I[i])
            # a = normality >= rlorder[i] ? 1 + cliquesize[i] : 1
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
                            # @inbounds bi = qtermadd3left(at,c,b,n)
                        end
                        #@inbounds bi = qtermadd3(at,b,c,n)
                        # @inbounds bi = qtermadd3(b,c,at,n)
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
                            bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
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
                                # @inbounds bi = qtermadd3left(at,c,b,n)
                            end
                            #@inbounds bi = qtermadd3(at,b,c,n)
                            # @inbounds bi = qtermadd3(b,c,at,n)
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
                                bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
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
        # epos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, cql)
        epos = Vector{Vector{Any}}(undef, cql)
        for i = 1:cql
            # println(cliques[i],J[i])
            if !isempty(J[i])
                # epos[i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, length(J[i]))
                epos[i] = Vector{Any}(undef, length(J[i]))
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
                                    # @inbounds bi = qtermadd3left(a,c,b,n)
                                end
                                #@inbounds bi = qtermadd3(a,b,c,n)
                                # @inbounds bi = qtermadd3(b,c,a,n)
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
                                    bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
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
                                # epos[i][j] = @variable(model,[1:Int(4bs*(4bs-1)/2)])
                            else
                                epos[i][j]= @variable(model, [1:bs, 1:bs], Symmetric)
                                # epos[i][j] = @variable(model,[1:Int(bs*(bs-1)/2)])
                                # epos[i][j] = @variable(model, [1:bs, 1:bs])
                                # for t in 1:bs
                                #     @constraint(model, epos[i][j][t,t] == 0)
                                #     for r in t+1:bs
                                #         @constraint(model, epos[i][j][r,t] == -epos[i][j][t,r])
                                #     end
                                # end
                            end
                            for t = 1:bs, r = 1:bs
                            # for t = 1:bs-1, r = t+1:bs
                                for s = 1:length(supp[k+1])
                                    a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                                    b=deepcopy(supp[k+1][s])
                                    c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                                    if ipart
                                        @inbounds bi = qtermadd3(b,c,a,n)
                                    else
                                        @inbounds bi = qtermadd3left(a,b,c,n)
                                        # @inbounds bi = qtermadd3left(a,c,b,n)
                                    end
                                    #@inbounds bi = qtermadd3(a,b,c,n)
                                    # @inbounds bi = qtermadd3(b,c,a,n)
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
                                        bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
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
                                        # @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][Int((2bs-t)*(t-1)/2) + (r-t)])
                                    end                   
                                end
                            end
                        end
                    elseif J[i][j] <= m - qncnumeq
                        bs = length(eblocks[i][j])
                        if bs == 1
                            epos[i][j] = @variable(model)
                            # epos[i][j] = 0.0
                            for s = 1:length(supp[k+1])
                                a=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                                b=deepcopy(supp[k+1][s])
                                c=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                                if ipart
                                    @inbounds bi = qtermadd3(b,c,a,n)
                                else
                                    @inbounds bi = qtermadd3left(a,b,c,n)
                                    # @inbounds bi = qtermadd3left(a,c,b,n)
                                end
                                #@inbounds bi = qtermadd3(a,b,c,n)
                                # @inbounds bi = qtermadd3(b,c,a,n)
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
                                    bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
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
                                    # @inbounds bi = qtermadd3(b_star,c,a,n)
                                else
                                    @inbounds bi = qtermadd3left(a,b_star,c,n)
                                    # @inbounds bi = qtermadd3left(a,c,b,n)
                                end
                                #@inbounds bi = qtermadd3(a,b,c,n)
                                # @inbounds bi = qtermadd3(b,c,a,n)
                                if nb > 0
                                    bi = qreduce_unitnorm(bi,n;nb=nb)
                                end
                                # if ipart != true
                                #     bi= canonical_qmono(bi,n)
                                # else
                                #     bi_star = star(deepcopy(bi),n)
                                #     if bi_star == bi
                                #         bi = canonical_qmono(bi,n)
                                #     end
                                # end
                                bi_cano= canonical_qmono(bi,n)
                                if bi_cano != canonical_qmono(deepcopy(bi_cano),n)
                                    println(bi_cano)
                                end 
                                if nb > 0
                                    bi = qreduce_unitnorm(bi, n, nb=nb)
                                    bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
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
                                # epos[i][j] = @variable(model,[1:Int(4bs*(4bs-1)/2)])
                            else
                                epos[i][j]= @variable(model, [1:bs, 1:bs])
                                # epos[i][j] = @variable(model,[1:Int(bs*(bs-1)/2)])
                                # epos[i][j] = @variable(model, [1:bs, 1:bs])
                            end
                            for t = 1:bs, r = 1:bs
                            # for t = 1:bs-1, r = t+1:bs
                                for s = 1:length(supp[k+1])
                                    a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                                    b=deepcopy(supp[k+1][s])
                                    c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                                    if ipart
                                        @inbounds bi = qtermadd3(b,c,a,n)
                                    else
                                        @inbounds bi = qtermadd3left(a,b,c,n)
                                        # @inbounds bi = qtermadd3left(a,c,b,n)
                                    end
                                    #@inbounds bi = qtermadd3(a,b,c,n)
                                    # @inbounds bi = qtermadd3(b,c,a,n)
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
                                        bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
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
                                        # @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][Int((2bs-t)*(t-1)/2) + (r-t)])
                                    end               
                                end
                                for s = 1:length(supp[k+1])
                                    a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                                    b=deepcopy(supp[k+1][s])
                                    c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                                    a_star = star(a,n)
                                    b_star = star(b,n)
                                    if ipart
                                        # @inbounds bi = qtermadd3(b_star,c,a,n)
                                        @inbounds bi = qtermadd3(c,a_star,b,n)
                                        # @inbounds bi = qtermadd3(c,a_star,b,n)
                                    else
                                        @inbounds bi = qtermadd3left(a,b_star,c,n)
                                    end
                                    #@inbounds bi = qtermadd3(a,b,c,n)
                                    # @inbounds bi = qtermadd3(b,c,a,n)
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
                                        bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                    end
                                    Locb = bfind(tsupp, ltsupp, bi)
                                    if ipart == true
                                        # @inbounds add_to_expression!(rcons[Locb_cano], real(coe[k+1][s]), epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        # @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[1], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        # @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[2], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        # @inbounds add_to_expression!(rcons[Locb_cano], imag_part(coe[k+1][s])[3], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
            
                                        # @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        # #+A_J-A_K->-A_J+A_K/4.12
                                        # @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        # @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[3], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        # @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        # #\A_R(J_X)-\A_I(K_X)-\A_J(R_X)+\A_K(I_X)=b_J->-\A_J(R_X)-\A_K(I_X)+\A_R(J_X)+\A_I(K_X)=b_J            
                                        # @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        # @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        # @inbounds add_to_expression!(jcons[Locb], real(coe[k+1][s]), epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        # @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[1], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        # #\A_R(K_X)+\A_I(J_X)-\A_J(I_X)-\A_K(R_X)=b_K->-\A_K(R_X)+\A_J(I_X)-\A_I(J_X)+\A_R(K_X)=b_K               
                                        # @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                        # @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                        # @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                        # @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                        # 
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
                                        # @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][Int((2bs-t)*(t-1)/2) + (r-t)])
                                    end
                                    # if bi == Vector{UInt16}[[0x0002], [0x0004], [0x0004, 0x0001]]
                                    #     println("t", t, "r", r, ", Value: ", -epos[i][j][r,t+bs], "Icon: ", icons[Locb])
                                    # end               
                                end
                            end
                        end
                    else
                        bs = length(eblocks[i][j])
                        if bs == 1
                            epos[i][j] = @variable(model)
                            # epos[i][j] = 0.0
                            for s = 1:length(supp[k+1])
                                a=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                                b=deepcopy(supp[k+1][s])
                                c=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                                # if ipart
                                @inbounds bi = qtermadd3(b,c,a,n)
                                # else
                                #     @inbounds bi = qtermadd3left(a,b,c,n)
                                    # @inbounds bi = qtermadd3left(a,c,b,n)
                                # end
                                #@inbounds bi = qtermadd3(a,b,c,n)
                                # @inbounds bi = qtermadd3(b,c,a,n)
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
                                    bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
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
                                # if ipart
                                @inbounds bi = qtermadd3(c,a_star,b,n)
                                # else
                                #     @inbounds bi = qtermadd3left(a,b_star,c,n)
                                #     # @inbounds bi = qtermadd3left(a,c,b,n)
                                # end
                                #@inbounds bi = qtermadd3(a,b,c,n)
                                # @inbounds bi = qtermadd3(b,c,a,n)
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
                                    bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
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
                            # if ipart == true
                            epos[i][j]= @variable(model, [1:4bs, 1:4bs])
                                # epos[i][j] = @variable(model,[1:Int(4bs*(4bs-1)/2)])
                            # else
                            #     epos[i][j]= @variable(model, [1:bs, 1:bs])
                                # epos[i][j] = @variable(model,[1:Int(bs*(bs-1)/2)])
                                # epos[i][j] = @variable(model, [1:bs, 1:bs])
                            # end
                            for t = 1:bs, r = 1:bs
                            # for t = 1:bs-1, r = t+1:bs
                                for s = 1:length(supp[k+1])
                                    a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                                    b=deepcopy(supp[k+1][s])
                                    c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                                    # if ipart
                                    @inbounds bi = qtermadd3(b,c,a,n)
                                    # else
                                    # @inbounds bi = qtermadd3left(a,b,c,n)
                                        # @inbounds bi = qtermadd3left(a,c,b,n)
                                    # end
                                    #@inbounds bi = qtermadd3(a,b,c,n)
                                    # @inbounds bi = qtermadd3(b,c,a,n)
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
                                        bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
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
                                    # else
                                    #     @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][t,r])
                                    #     # @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][Int((2bs-t)*(t-1)/2) + (r-t)])
                                    # end                   
                                end
                                for s = 1:length(supp[k+1])
                                    a=deepcopy(ebasis[i][j][eblocks[i][j][t]])
                                    b=deepcopy(supp[k+1][s])
                                    c=deepcopy(ebasis[i][j][eblocks[i][j][r]])
                                    # b_star = star(b,n)
                                    a_star = star(a,n)
                                    # if ipart
                                    @inbounds bi = qtermadd3(c,a_star,b,n)
                                    # else
                                    #     @inbounds bi = qtermadd3left(a,b_star,c,n)
                                    #     # @inbounds bi = qtermadd3left(a,c,b,n)
                                    # end
                                    #@inbounds bi = qtermadd3(a,b,c,n)
                                    # @inbounds bi = qtermadd3(b,c,a,n)
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
                                        bi_cano = qreduce_unitnorm(bi_cano, n, nb=nb)
                                    end
                                    Locb = bfind(tsupp, ltsupp, bi)
                                    Locb_cano = bfind(tsupp, ltsupp, bi_cano)
                                    # if ipart == true
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
                                    # else
                                        # @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][r,t])
                                        # @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][Int((2bs-t)*(t-1)/2) + (r-t)])
                                    # end                   
                                end
                            end
                        end
                    end
                end
            end
        end
        # println(length(ncc))
        for i in ncc
            if i <= m - numeq
                pos0 = @variable(model, lower_bound=0)
            else
                pos0 = @variable(model)
            end
            for j = 1:length(supp[i+1])
                Locb = bfind(tsupp, ltsupp, supp[i+1][j])
                if ipart == true
                    @inbounds add_to_expression!(icons[Locb], imag_part(coe[i+1][j])[1], pos0)
                    @inbounds add_to_expression!(jcons[Locb], imag_part(coe[i+1][j])[2], pos0)
                    @inbounds add_to_expression!(kcons[Locb], imag_part(coe[i+1][j])[3], pos0)
                end
                @inbounds add_to_expression!(rcons[Locb], real(coe[i+1][j]), pos0)
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
        # println(supp[1])
        # visited = Set{Any}()

        # for i = 1:length(supp[1])

        #     mono = canonical_qmono(supp[1][i], n)

        #     if nb > 0
        #         mono = qreduce_unitnorm(mono, nb=nb)
        #     end

        #     mono_star = canonical_qmono(star(deepcopy(mono), n),n)

        #     if nb > 0
        #         mono_star = qreduce_unitnorm(mono_star, nb=nb)
        #     end

        #     # 如果共轭已经出现过，则跳过
        #     if mono_star in visited
        #         continue
        #     end

        #     # 记录当前项
        #     push!(visited, mono)

        #     Locb = bfind(tsupp, ltsupp, mono)

        #     if Locb === nothing
        #         println(mono)
        #         @error "The monomial basis is not enough!"
        #         return nothing,ksupp,nothing,nothing,nothing
        #     else
        #         rbc[Locb] += real(coe[1][i])

        #         if ipart == true ||qncnumeq > 0
        #             ibc[Locb] += imag_part(coe[1][i])[1]
        #             jbc[Locb] += imag_part(coe[1][i])[2]
        #             kbc[Locb] += imag_part(coe[1][i])[3]
        #         end
        #     end
        # end
        for i = 1:length(supp[1])
            # mono = standardterm(supp[1][i],n)
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
                # rbc[Locb_cano] += real(coe[1][i])
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

                # mono_star = canonical_qmono(star(deepcopy(mono), n), n)
                mono_star = star(deepcopy(mono), n)
                # if mono == canonical_qmono(deepcopy(mono),n)
                #     mono_star = canonical_qmono(star(deepcopy(mono), n),n)
                # else
                #     mono_star = star(deepcopy(mono), n)
                # end
                # mono_star = star(standardterm(deepcopy(mono), n), n)
                # if mono != canonical_qmono(deepcopy(mono),n)
                #     # mono_star = canonical_qmono(star(deepcopy(mono), n),n)
                #     mono_star = star(deepcopy(mono), n)
                #     j = bfind(tsupp, ltsupp, mono_star)
                # else
                #     j = nothing
                # end

                j = bfind(tsupp, ltsupp, mono_star)
                # j1 = bfind(tsupp, ltsupp, mono_star1)

                # 没找到共轭项
                if j === nothing
                    push!(visited, i)
                    println("No conjugate found for ", mono, " at index ", i)
                    continue
                end

                # self-adjoint monomial
                if j == i
                    push!(visited, i)
                    continue
                end

                push!(visited, i)
                push!(visited, j)


                # if abs(rbc[i]) > 1e-12 || (ipart && (abs(ibc[i]) > 1e-12 || abs(jbc[i]) > 1e-12 || abs(kbc[i]) > 1e-12))
                #     keep[j] = false

                # elseif abs(rbc[j]) > 1e-12 || (ipart && (abs(ibc[j]) > 1e-12 || abs(jbc[j]) > 1e-12 || abs(kbc[j]) > 1e-12))
                #     keep[i] = false

                # else
                    # 都是0，默认保留较小index
                keep[max(i,j)] = false
                # keep[j] = false
                # end
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
                # j1 = bfind(tsupp, ltsupp, mono_star1)

                # 没找到共轭项
                if j === nothing
                    push!(visited1, i)
                    nonc+=1
                    # println("No conjugate found for ", mono, " at index ", i)
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

                if abs(rbc[i]) > 1e-12 || (ipart && (abs(ibc[i]) > 1e-12 || abs(jbc[i]) > 1e-12 || abs(kbc[i]) > 1e-12))
                    keep1[j] = false

                elseif abs(rbc[j]) > 1e-12 || (ipart && (abs(ibc[j]) > 1e-12 || abs(jbc[j]) > 1e-12 || abs(kbc[j]) > 1e-12))
                    keep1[i] = false

                else
                    # 都是0，默认保留较小index
                keep1[max(i,j)] = false
                # keep[j] = false
                end
            end
            # @constraint(model, rcon[i=1:ltsupp], rcons[i]==rbc[i])
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
                # j1 = bfind(tsupp, ltsupp, mono_star1)

                # 没找到共轭项
                if j === nothing
                    push!(visited, i)
                    println("No conjugate found for ", mono, " at index ", i)
                    continue
                end

                # self-adjoint monomial
                if j == i
                    push!(visited, i)
                    continue
                end

                push!(visited, i)
                push!(visited, j)


                # if abs(rbc[i]) > 1e-12 || (ipart && (abs(ibc[i]) > 1e-12 || abs(jbc[i]) > 1e-12 || abs(kbc[i]) > 1e-12))
                #     keep[j] = false

                # elseif abs(rbc[j]) > 1e-12 || (ipart && (abs(ibc[j]) > 1e-12 || abs(jbc[j]) > 1e-12 || abs(kbc[j]) > 1e-12))
                #     keep[i] = false

                # else
                    # 都是0，默认保留较小index
                keep[max(i,j)] = false
                # keep[j] = false
                # end
            end
            inds = findall(keep)
        end
        # println(keep[80],keep[87])
        # println("有",nonc,"项找不到共轭的轮转")
        # println("no cano cons number:",length(inds),"cano cons number:",length(inds1))
        # useless = setdiff(1:length(ttsupp), inds1)
        # if Vector{UInt16}[[], [], [0x0001, 0x0004,0x0004]] in ttsupp[inds1]
        #     println("term 1 in useless")
        # end

        # println("===== checking keep conjugate pairs =====")

        # found = false

        # for i = 1:length(tsupp[inds])

        #     mono = standardterm(tsupp[inds[i]], n)

        #     # 共轭
        #     mono_star = canonical_qmono(star(deepcopy(mono), n),n)

        #     # 在 tsupp[inds] 里找
        #     for j = i:length(tsupp[inds])

        #         if i == j
        #             continue
        #         end

        #         if tsupp[inds[j]] == mono_star

        #             println("conjugate pair found:")
        #             println("i = ", i, " : ", mono)
        #             println("j = ", inds[j], " : ", tsupp[inds[j]])
        #             println()

        #             found = true
        #         end
        #     end
        # end
        # if !found
        #     println("in keep No conjugate pairs remain.")
        # end
        # println("===== checking remove conjugate pairs =====")

        # found = false

        # for i = 1:length(useless)

        #     mono = standardterm(ttsupp[useless[i]], n)

        #     # 共轭
        #     if mono == canonical_qmono(mono,n)
        #         mono_star = canonical_qmono(star(deepcopy(mono), n),n)
        #     else
        #         mono_star = star(deepcopy(mono), n)
        #     end
        #     # mono_star = star(deepcopy(mono), n)

        #     # 在 tsupp[inds] 里找
        #     for j = 1:length(useless)

        #         if i == j
        #             continue
        #         end

        #         if ttsupp[useless[j]]== mono_star

        #             println("conjugate pair found:")
        #             println("i = ", useless[i], " : ", mono)
        #             println("j = ", useless[j], " : ", ttsupp[useless[j]])
        #             println()

        #             found = true
        #         end
        #     end
        # end

        # if !found
        #     println("in remove No conjugate pairs remain.")
        # end
        # for i = 1:length(inds)
        #     println([Int.(tsupp[inds[i]][1]);Int.(tsupp[inds[i]][2]);Int.(tsupp[inds[i]][3])])
        # end
        # for k = 1:length(inds1)
        #     if ttsupp[inds1[k]] == star(deepcopy(ttsupp[inds1[k]]), n)
        #         println(ttsupp[inds1[k]], " at index ", inds1[k])
        #     end
        # end
        # for i = 1:length(inds1)
        #     term = ttsupp[inds1[i]]
        #     term_star = star(deepcopy(term), n)
        #     if term != term_star
        #         println("keep term",i,":",term)
        #     end
        # end
        # @constraint(model,rcon[k=1:length(inds1)],rcons[inds1[k]] == rbc[inds1[k]])       
        if ipart == true || qncnumeq > 0
            # @constraint(model, icon[i=1:ltsupp], icons[i]==ibc[i])
            # @constraint(model, jcon[i=1:ltsupp], jcons[i]==jbc[i])
            # @constraint(model, kcon[i=1:ltsupp], kcons[i]==kbc[i])
            @constraint(model,rcon[k=1:length(inds1)],rcons[inds1[k]] == rbc[inds1[k]])  
            @constraint(model, icon[k=1:length(inds)], icons[inds[k]] == ibc[inds[k]])
            @constraint(model, jcon[k=1:length(inds)], jcons[inds[k]] == jbc[inds[k]])
            @constraint(model, kcon[k=1:length(inds)], kcons[inds[k]] == kbc[inds[k]])
        else
            @constraint(model,rcon[k=1:length(inds)],rcons[inds[k]] == rbc[inds[k]])
        end
        # visited = Set{Int}()
        # for i = 1:ltsupp
        #     if i in visited
        #         continue
        #     end
        #     mono = standardterm(tsupp[i], n)
        #     mono_star = star(deepcopy(mono), n)
        #     j = bfind(tsupp, ltsupp, mono_star)
        #     println(mono,"at index", i, " icons ", icons[inds[i]],"value ", ibc[inds[i]])
        #     println(mono_star,"at index", j," conj icons", icons[inds[j]] ,"value ", ibc[inds[j]])
        #     push!(visited, i)
        #     push!(visited, j)
        # end
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
        if writetofile != false
            write_to_file(dualize(model), writetofile)
        end
        SDP_status = termination_status(model)
        cons=all_constraints(model; include_variable_in_set_constraints = false)
        # println("cons number:",length(rcons[inds]))
        # for (i,a) in in enumerate(tsupp[inds])
        #     print("term:", [Int.(c) for (j,c) in enumerate(a)])
        # end
        # for i = 1:length(inds)
        #     println( [Int.(tsupp[i][1]);Int.(tsupp[i][2]);Int.(tsupp[i][3])])
        # end
        # println("term:",[Int.(a) for a in tsupp[inds]])
        # for (i, c) in enumerate(cons)
        #     println("Constraint $i: ")
        #     println("func = ", JuMP.constraint_object(c).func)
        #     println("set  = ", JuMP.constraint_object(c).set)
        # end
        objv = objective_value(model)
        if SDP_status != MOI.OPTIMAL
            println("termination status: $SDP_status")
            status = primal_status(model)
            println("solution status: $status")
        end
        println("optimum = $objv")
        # println(blocksize[1])
        # rmeasure = [-dual(c) for c in rcon]
        # imeasure = nothing
        # jmeasure = nothing
        # kmeasure = nothing
        # if ipart == true
        #     imeasure = [-dual(c) for c in icon]
        #     jmeasure = [-dual(c) for c in jcon]
        #     kmeasure = [-dual(c) for c in kcon]
        # end
        if solution == true
            rmeasure = [-dual(c) for c in rcon]
            imeasure = nothing
            jmeasure = nothing
            kmeasure = nothing
            if ipart == true
                imeasure = [-dual(c) for c in icon]
                jmeasure = [-dual(c) for c in jcon]
                kmeasure = [-dual(c) for c in kcon]
            end
            moment = get_qmoment(n,rmeasure, imeasure, jmeasure, kmeasure, tsupp, cql, blocks, cl, blocksize, basis;ipart=ipart,nb = nb)
            # A = qmat_to_realmat(moment[1][1])
            # evals, evecs = eigen(Symmetric(A))
            # idx_max = argmax(evals)
            # v = evecs[:, idx_max]

            # sol = extract_first_order_quaternion_solution(v, tsupp, n)

            # # 打印结果
            # println("提取的一阶 quaternion 变量解:")
            # for (i,q) in enumerate(sol)
            #     println("q[$i] = $q, 模长 = ", norm(q))
            # end
            # sol = extract_quaternion_solution(moment)
            # linear_monos = find_all_linear_monomials(tsupp, ltsupp)
            # sol = Quaternion{Float64}[]
            # for i in 1:n
            #     idx = bfind(tsupp, ltsupp, [UInt16[], UInt16[], [UInt16(i)]])
            #     println("i=$i, idx=$idx")
            #     println("  r = ", rmeasure[idx])
            #     println("  i = ", imeasure[idx])
            #     println("  j = ", jmeasure[idx])
            #     println("  k = ", kmeasure[idx])
            # end
            # for i in 1:n
            #     if haskey(linear_monos, i)
            #         idx = linear_monos[i]
            #         a = rmeasure[idx]
            #         b = imeasure !== nothing ? -imeasure[idx] : 0.0
            #         c = jmeasure !== nothing ? -jmeasure[idx] : 0.0
            #         d = kmeasure !== nothing ? -kmeasure[idx] : 0.0
            #         push!(sol, normalize(Quaternion(a, b, c, d)))
            #     else
            #         println("Warning: no linear moment for q_$i")
            #     end
            # end
            # println(tsupp)
            # one_order_monomials = [[UInt16[], UInt16[], UInt16[i]] for i in 1:n]
            # for (i, mono) in enumerate(one_order_monomials)
            #     idx = findfirst(m -> m == mono, tsupp)
            #     if idx !== nothing
            #         println("mono $mono found at index $idx, value = ", rmeasure[idx])
            #     else
            #         println("mono $mono NOT found in tsupp")
            #     end
            # end
            # for i in 1:n
            #     println("q[$i] in tsupp? ", any(m -> m == q[i], tsupp))
            #     println("conj(q[$i]) in tsupp? ", any(m -> m == q[i+n], tsupp))
            # end
            # if TS != false || CS != false
            #     sol = extract_solution_from_sparse_moment(moment, cliques, n)
            # else
            sol = Quaternion{Float64}[]
            for i = 1:n
                a = rmeasure[bfind(tsupp, ltsupp, [UInt16[], UInt16[], UInt16[i]])]
                b = ipart ? -imeasure[bfind(tsupp, ltsupp, [UInt16[], UInt16[], UInt16[i]])] : 0.0
                c = ipart ? -jmeasure[bfind(tsupp, ltsupp, [UInt16[], UInt16[], UInt16[i]])] : 0.0
                d = ipart ? -kmeasure[bfind(tsupp, ltsupp, [UInt16[], UInt16[], UInt16[i]])] : 0.0
                push!(sol, Quaternion(a, b, c, d))
            end
        end
    end
    if solution
        return objv,ksupp,SDP_status,sol,moment
    else
        return objv,ksupp,SDP_status
    end
    # return objv,ksupp,SDP_status,sol,moment
end
function find_all_linear_monomials(tsupp, ltsupp)
    linear_monomials = Dict{Int, Int}()  # 变量索引 => tsupp 索引
    # 遍历所有单项式
    for idx in 1:length(tsupp)
        bi = tsupp[idx]
        total_degree = sum(length(b) for b in bi)
        if total_degree == 1
            # 找出是哪一个变量，假设只在其中一个部分有1个元素
            for part in bi
                if length(part) == 1
                    var_idx = Int(part[1])
                    linear_monomials[var_idx] = idx
                end
            end
        end
    end
    return linear_monomials
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
function clique_decomp(n, m, dc, supp::Vector{Vector{Vector{Vector{UInt16}}}}; order="min", alg="MF")
    if alg == false
        cliques = [UInt16[i for i=1:n]]
        cql = 1
        cliquesize = [n]
    else
        G = SimpleGraph(n)
        for i = 1:m+1
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
# function get_blocks2(n, rlorder, I, J, NJ, supp::Vector{Vector{Vector{Vector{UInt16}}}}, cliques, cliquesize, cql, tsupp, basis, ebasis,nebasis; TS="block", ConjugateBasis=false, nb=0, normality=1, merge=false, md=3)
#     blocks = Vector{Vector{Vector{Vector{Int}}}}(undef, cql)
#     eblocks = Vector{Vector{Vector{Int}}}(undef, cql)
#     neblocks = Vector{Vector{Vector{Int}}}(undef, cql)
#     cl = Vector{Vector{Int}}(undef, cql)
#     blocksize = Vector{Vector{Vector{Int}}}(undef, cql)
#     for i = 1:cql
#         # temp_tsupp = deepcopy(tsupp)
#         # for j in eachindex(temp_tsupp)
#         #     temp_tsupp[j][3] = [x > UInt16(n) ? x - UInt16(n) : x for x in temp_tsupp[j][3]]
#         # end
#         # ksupp = TS ? tsupp[[issubset(union(tsupp[j][1], tsupp[j][2], temp_tsupp[j][3]), cliques[i]) for j = 1: length(tsupp)]] : nothing
#         ksupp = TS == false ? nothing : tsupp[[
#             let 
#                 modified_3rd_l = [x > UInt16(n) ? x - UInt16(n) : x for x in tsupp[j][3]]
#                 # modified_3rd_g = [x <= UInt16(n) ? x + UInt16(n) : x for x in tsupp[j][3]]
#                 # issubset(union(tsupp[j][1], tsupp[j][2], modified_3rd_l, modified_3rd_g), cliques[i])
#                 issubset(union(tsupp[j][1], modified_3rd_l), cliques[i])
#             end 
#             for j in eachindex(tsupp)
#         ]]
#         # ksupp = TS == false ? nothing : tsupp[[issubset(union(tsupp[j][1], tsupp[j][2], tsupp[j][3].=[x <= UInt16(n) ? x - UInt16(n) : x for x in tsupp[j][3]]), cliques[i]) for j in eachindex(tsupp)]]
#         # println(ksupp)
#         blocks[i],eblocks[i],cl[i],blocksize[i] = get_blocks(length(I[i]), length(J[i]), rlorder[i], ksupp, supp[[I[i]; J[i]].+1], basis[i], 
#         ebasis[i], TS=TS, ConjugateBasis=ConjugateBasis, merge=merge, md=md, nb=nb, normality=normality, nvar=cliquesize[i])
#     end
#     return blocks,eblocks,cl,blocksize
# end
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
function assign_constraint2(n,m,m1,numeq, supp::Vector{Vector{Vector{Vector{UInt16}}}}, cliques, cql)
    I = [Int[] for i=1:cql]
    J = [Int[] for i=1:cql]
    NJ = [Int[] for i=1:cql]
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
        elseif i <= m - m1
            push!.(J[ind], i)
        else
            push!.(NJ[ind], i)
        end
    end
    return I,J,NJ,ncc
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
function get_qmoment(n,rmeasure, imeasure, jmeasure, kmeasure, tsupp, cql, blocks, cl, blocksize, basis; ipart = true, nb=0)
    ctype = ipart == true ? QuaternionF64 : Float64
    moment = Vector{Vector{Matrix{ctype}}}(undef, cql)
    for i = 1:cql
        moment[i] = Vector{Matrix{Quaternion{Float64}}}(undef, cl[i][1])
        for l = 1:cl[i][1]
            bs = blocksize[i][1][l]
            rtemp = zeros(Float64, bs, bs)
            if ipart == true
                itemp = zeros(Float64, bs, bs)
                jtemp = zeros(Float64, bs, bs)
                ktemp = zeros(Float64, bs, bs)
            end
            for t = 1:bs, r = 1:bs
                a = deepcopy(basis[i][1][blocks[i][1][l][t]])
                b = deepcopy(basis[i][1][blocks[i][1][l][r]])
                bi = qtermadd(a, b, n)

                if nb > 0
                    bi = qreduce_unitnorm(bi, n,nb=nb)
                end

                Locb = bfind(tsupp, length(tsupp), bi)
                rtemp[t,r] = rmeasure[Locb]
                if ipart == true
                    itemp[t,r] = imeasure[Locb]
                    jtemp[t,r] = jmeasure[Locb]
                    ktemp[t,r] = kmeasure[Locb]
                end
            end

            rtemp = (rtemp + rtemp') / 2
            if ipart == true
                itemp = (itemp - itemp') / 2
                jtemp = (jtemp - jtemp') / 2
                ktemp = (ktemp - ktemp') / 2
            end

            qmat = Matrix{ctype}(undef, bs, bs)
            for t = 1:bs, r = 1:bs
                if ipart == true
                    qmat[t, r] = Quaternion(rtemp[t,r], itemp[t,r], jtemp[t,r], ktemp[t,r])
                else
                    qmat[t, r] = rtemp[t,r]
                end
            end
            moment[i][l] = qmat
        end
    end
    return moment
end
function qmat_to_realmat(Q::Matrix{QuaternionF64})
    n = size(Q, 1)
    A = zeros(Float64, 4n, 4n)

    for r = 1:n, c = 1:n
        q = Q[r,c]
        a, b, c_, d = real(q), imag_part(q)[1], imag_part(q)[2], imag_part(q)[3]
        base = 4*(r-1)+1
        colbase = 4*(c-1)+1

        A[base,     colbase]     = a
        A[base,     colbase+1]   = -b
        A[base,     colbase+2]   = -c_
        A[base,     colbase+3]   = -d

        A[base+1,   colbase]     = b
        A[base+1,   colbase+1]   = a
        A[base+1,   colbase+2]   = -d
        A[base+1,   colbase+3]   = c_

        A[base+2,   colbase]     = c_
        A[base+2,   colbase+1]   = d
        A[base+2,   colbase+2]   = a
        A[base+2,   colbase+3]   = -b

        A[base+3,   colbase]     = d
        A[base+3,   colbase+1]   = -c_
        A[base+3,   colbase+2]   = b
        A[base+3,   colbase+3]   = a
    end

    return A
end

function extract_first_order_quaternion_solution(v::AbstractVector{Float64}, tsupp::Vector, qcount::Int)
    qsol = QuaternionF64[]

    # 找出 tsupp 中一阶单变量的索引：形如 [[], [], [i]]，且 i 只有一个元素
    one_order_indices = []
    for (idx, mono) in enumerate(tsupp)
        # 这里根据你tsupp格式，判断是否为一阶单变量
        if length(mono) == 3 && length(mono[3]) == 1 && isempty(mono[1]) && isempty(mono[2])
            push!(one_order_indices, idx)
        end
    end

    # 检查数量是否匹配
    if length(one_order_indices) < qcount
        error("找到的一阶单变量数量 $(length(one_order_indices)) 小于 qcount = $qcount,请确认 tsupp 格式和 qcount")
    elseif length(one_order_indices) > qcount
        @warn "找到的一阶单变量数量 $(length(one_order_indices)) 多于 qcount = $qcount,只提取前 $qcount 个变量"
        one_order_indices = one_order_indices[1:qcount]
    end

    # 按 tsupp 中顺序依次提取对应的 quaternion
    for idx in one_order_indices
        base = 4*(idx - 1) + 1
        q = QuaternionF64(v[base], v[base+1], v[base+2], v[base+3])
        push!(qsol, q)
    end

    return qsol
end
function extract_clique_quaternions(v::Vector{Float64}, clique_vars::Vector{UInt16})
    qsol = QuaternionF64[]
    for (idx, var_idx) in enumerate(clique_vars)
        base = 4*(idx - 1) + 1
        q = QuaternionF64(v[base], v[base+1], v[base+2], v[base+3])
        push!(qsol, q)
    end
    return qsol
end

function extract_solution_from_sparse_moment(moment::Vector{Vector{Matrix{QuaternionF64}}}, cliques::Vector{Vector{UInt16}}, n::Int)
    # cliques 是一个列表，每个元素是该 clique 涉及的变量索引向量
    var_candidates = Dict{Int, Vector{QuaternionF64}}()
    for var in 1:n
        var_candidates[var] = QuaternionF64[]
    end

    # 遍历所有 clique
    for (block_idx, clique_vars) in enumerate(cliques)
        Q = moment[1][block_idx]  # 取第一个 cql 下、第 block_idx 个子矩阵
        A = qmat_to_realmat(Q)
        evals, evecs = eigen(Symmetric(A))
        v = evecs[:, argmax(evals)]
        # 提取该 clique 中的解
        qsol_local = extract_clique_quaternions(v, clique_vars)
        for (var_idx, q) in zip(clique_vars, qsol_local)
            push!(var_candidates[var_idx], q)
        end
    end

    # 融合每个变量的候选解
    qsol_global = Vector{QuaternionF64}(undef, n)
    for var in 1:n
        qs = var_candidates[var]
        if isempty(qs)
            # 若某变量在所有 clique 中均未出现（不常见），可设为默认单位四元数
            qsol_global[var] = QuaternionF64(1.0, 0.0, 0.0, 0.0)
        else
            s_avg = mean(q.s for q in qs)
            v1_avg = mean(q.v1 for q in qs)
            v2_avg = mean(q.v2 for q in qs)
            v3_avg = mean(q.v3 for q in qs)
            q_avg = QuaternionF64(s_avg, v1_avg, v2_avg, v3_avg)
            qsol_global[var] = q_avg / norm(q_avg)
        end
    end

    return qsol_global
end