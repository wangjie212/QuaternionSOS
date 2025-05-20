function qs_tssos_first(pop::Vector{Polynomial{false,T}}, z, n::Int, d; numeq=0, RemSig=false, nb=0,CS="MF",cliques=[],TS="block",
    merge=false, md=3, solver="Mosek", reducebasis=false, QUIET=false, solve=true, solution=false, ipart=true, 
    dualize=false, balanced=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), 
    writetofile=false, normality=0, NormalSparse=false, conjubasis=true) where {T<:Number}
    supp,coe = Qpolys_info(pop, z, n)
    println("*********************************** QSSOS ***********************************")
    println("QSSOS is launching...")
    if nb > 0
        supp[1],coe[1] = qresort(supp[1],coe[1];nb=nb)
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
        cliques,cql,cliquesize = clique_decomp(2n, m, dc, supp, order=d, alg=CS)
        end
        if CS != false && QUIET == false
            mc = maximum(cliquesize)
            println("Obtained the variable cliques in $time seconds.\nThe maximal size of cliques is $mc.")
        end
    end
    I,J,ncc = assign_constraint(m, numeq, supp, cliques, cql)
    rlorder = d*ones(Int, cql)
    # basis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, m+1)
    # basis[1] = get_qncbasis(n, d; conjubasis=conjubasis)
    ebasis = Vector{Vector{Vector{Vector{Vector{UInt16}}}}}(undef, cql)
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
                        ltemp=qtermadd([UInt16[],UInt16[],[cliques[i][s]]],a,n)
                        push!(basis[i][s+1],ltemp)
                    end
                    if nb > 0
                        basis[i][s+1] = qreduce_unitnorm.(basis[i][s+1], nb=nb)
                        unique!(basis[i][s+1])
                    end
                end
                for s = 1:length(I[i])
                    basis[i][s+1+cliquesize[i]] = get_qncbasis(cliques[i], n, rlorder[i]-Int(ceil(dc[I[i][s]]/2)))
                    if nb > 0
                        basis[i][s+1] = qreduce_unitnorm.(basis[i][s+1], nb=nb)
                        unique!(basis[i][s+1])
                    end
                end
            # end
        else
            basis[i] = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(I[i])+1)
            basis[i][1] = get_qncbasis(cliques[i], n, d; conjubasis=conjubasis)
            if nb > 0
                basis[i][1] = qreduce_unitnorm.(basis[i][1], nb=nb)
                unique!(basis[i][1])
            end
            for s = 1:length(I[i])
                basis[i][s+1] = get_qncbasis(cliques[i], n, rlorder[i]-Int(ceil(dc[I[i][s]]/2)))
                if nb > 0
                    basis[i][s+1] = qreduce_unitnorm.(basis[i][s+1], nb=nb)
                    unique!(basis[i][s+1])
                end
            end
        end        
        for s = 1:length(J[i])
            ebasis[i][s] = get_qncbasis(cliques[i], n, rlorder[i]-Int(ceil(dc[J[i][s]]/2)))
            if nb > 0
                ebasis[i][s] = qreduce_unitnorm.(ebasis[i][s+1], nb=nb)
                unique!(ebasis[i][s])
            end
        end
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
    tsupp = deepcopy(supp[1])
    for i = 2:m+1, j = 1:length(supp[i])
        push!(tsupp, supp[i][j])
    end
    sort!(tsupp)
    unique!(tsupp)
    if TS != false && QUIET == false
        println("Starting to compute the block structure...")
    end
    time = @elapsed begin
    blocks,eblocks,cl,blocksize = get_blocks(n,rlorder, I, J, supp, cliques, cliquesize, cql, tsupp, basis, ebasis, TS=TS, ConjugateBasis=conjubasis, normality=normality, merge=merge, md=md, nb=nb)
    opt,ksupp,SDP_status = qsolvesdp(n, m, rlorder, supp, coe, basis, ebasis, cliques, cql, cliquesize, I, J, ncc, blocks, eblocks, cl, blocksize, numeq=numeq, QUIET=QUIET, TS=TS, solver=solver, solve=solve, solution=solution, ipart=ipart, balanced=balanced,
    nb=nb, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, dualize=dualize, writetofile=writetofile, normality=normality, NormalSparse=NormalSparse,conjubasis=conjubasis)
    return opt
end
end

function Qpolys_info(pop, z, n)
    coe = Vector{Vector{QuaternionF64}}(undef, length(pop))
    supp = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(pop))
    for k in eachindex(pop)
        # mon = monomials(pop[k])
        # coe[k] = coefficients(pop[k])
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
                # l = bfind(z, 2n, vars[j])
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
            # if !(supp[k][i] in temp)
            #     push!(temp,supp[k][i])
            # else
            #     push!(Ind,i)
            #     a=findall(isequal(supp[k][i]),temp)
            #     coe[k][a[1]]=coe[k][a[1]]+coe[k][i]
            # end
        end
        # supp[k]=deleteat!(supp[k],Ind)
        # coe[k]=deleteat!(coe[k],Ind)
    end
    return supp,coe
end
function qsolvesdp(n, m, rlorder, supp::Vector{Vector{Vector{Vector{UInt16}}}}, coe, basis, ebasis, cliques, cql, cliquesize, I, J, ncc, blocks, eblocks, cl, blocksize; 
    numeq=0, nb=0, QUIET=false, TS=false, solver="Mosek", solve=true, dualize=false, solution=false, ipart=true, 
    cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false, balanced=false, normality=0, NormalSparse=false,conjubasis=false)
    tsupp = Vector{Vector{UInt16}}[]
    for i = 1:cql
        # a = normality >= rlorder[i] ? cliquesize[i] + 1 : 1
        a = normality > 0 ? cliquesize[i] + 1 : 1
        for s = 1 : a, j = 1:cl[i][s], k = 1:blocksize[i][s][j], r = 1:blocksize[i][s][j]
            at=deepcopy(basis[i][s][blocks[i][s][j][k]])
            b=deepcopy(basis[i][s][blocks[i][s][j][r]])
            @inbounds bi = qtermadd(at,b,n)
            if nb > 0
                bi = qreduce_unitnorm(bi, nb=nb)
            end
            push!(tsupp,bi)
        end
    end
    if TS != false
        gsupp = get_gsupp(rlorder, basis, ebasis, supp, cql, I, J, ncc, blocks, eblocks, cl, blocksize, cliquesize, ConjugateBasis=conjubasis, nb=nb, normality=normality)
        append!(tsupp, gsupp)
    end
    # for k = 1:length(basis[1]), r = 1:length(basis[1])
    #     a=deepcopy(basis[1][k])
    #     b=deepcopy(basis[1][r])
    #     @inbounds bi = qtermadd(a,b,n)
    #     if nb > 0
    #         bi = qreduce_unitnorm(bi, nb=nb)
    #     end
    #     push!(tsupp,bi)
    # end
    ksupp = deepcopy(tsupp)
    # if normality > 0
    #     wbasis = [[UInt16[],UInt16[],UInt16[]]]
    #     for t=1:n
    #         wbasis= bget_qncbasis(n, normality,t)
    #         if nb > 0
    #             wbasis= qreduce_unitnorm.(wbasis, nb=nb)
    #             unique!(wbasis)
    #         end
    #         ws = length(wbasis)
    #         for k = 1:ws, r = 1:ws
    #             a=deepcopy(wbasis[k])
    #             b=deepcopy(wbasis[r])
    #             @inbounds bi = qtermadd(a,b,n)
    #             if nb > 0
    #                 bi = qreduce_unitnorm(bi, nb=nb)
    #             end
    #             push!(tsupp,bi)
    #         end
    #     end
    # end
    sort!(tsupp)
    unique!(tsupp)
    #initial set
    objv = SDP_status= nothing
    if solve == true
        ltsupp = length(tsupp)
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
        rcons = [AffExpr(0) for i=1:ltsupp]
        if ipart==true
            icons = [AffExpr(0) for i=1:ltsupp]
            jcons = [AffExpr(0) for i=1:ltsupp]
            kcons = [AffExpr(0) for i=1:ltsupp]
        end
        hnom = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, n)
        # if normality > 0
        #     for t=1:n
        #         wbasis= bget_qncbasis(n, normality,t)
        #         if nb > 0
        #             wbasis= qreduce_unitnorm.(wbasis, nb=nb)
        #             unique!(wbasis)
        #         end
        #         ws = length(wbasis)
        #         if NormalSparse == false
        #             hnom[t] = @variable(model, [1:ws, 1:ws], PSD)
        #             for j = 1:ws, k = 1:ws
        #                 a=deepcopy(wbasis[j])
        #                 b=deepcopy(wbasis[k])
        #                 bi =qtermadd(a,b,n)
        #                 if nb > 0
        #                     bi = qreduce_unitnorm(bi, nb=nb)
        #                 end
        #                 Locb = bfind(tsupp, ltsupp, bi)
        #                 @inbounds add_to_expression!(rcons[Locb], hnom[t][j,k])
        #             end
        #         end
        #     end
        # end
        # pos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m-numeq+1)
        # pos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m+1)
        pos = Vector{Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}}(undef, cql)
        for i = 1:cql
            # a = normality >= rlorder[i] ? 1 + cliquesize[i] : 1
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
                        @inbounds bi = qtermadd(at, b,n)
                        if nb > 0
                            bi = qreduce_unitnorm(bi;nb=nb)
                        end
                        Locb = bfind(tsupp, ltsupp, bi)
                        @inbounds add_to_expression!(rcons[Locb], pos[i][j][l])
                    else
                        if ipart == true
                            pos[i][j][l] = @variable(model, [1:4bs, 1:4bs], PSD)
                        else
                            pos[i][j][l] = @variable(model, [1:bs, 1:bs], PSD)
                        end
                        for t = 1:bs, r = 1:bs
                            at=deepcopy(basis[i][j][blocks[i][j][l][t]])
                            b=deepcopy(basis[i][j][blocks[i][j][l][r]])
                            @inbounds bi = qtermadd(at,b,n)
                            if nb > 0
                                bi = qreduce_unitnorm(bi;nb=nb)
                            end
                            Locb = bfind(tsupp, ltsupp, bi)
                            if ipart == true
                                @inbounds add_to_expression!(rcons[Locb], pos[i][j][l][t,r]+pos[i][j][l][t+bs,r+bs]+pos[i][j][l][t+2*bs,r+2*bs]+pos[i][j][l][t+3*bs,r+3*bs])
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
            # pos[i][1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[i][1])
            # for l = 1:cl[i][1]
            #     @inbounds bs = blocksize[i][1][l]
            #     if bs == 1
            #         pos[i][1][l] = @variable(model, lower_bound=0)
            #         a=deepcopy(basis[i][1][blocks[i][1][l][1]])
            #         b=deepcopy(basis[i][1][blocks[i][1][l][1]])
            #         @inbounds bi = qtermadd(a, b,n)
            #         if nb > 0
            #             bi = qreduce_unitnorm(bi;nb=nb)
            #         end
            #         Locb = bfind(tsupp, ltsupp, bi)
            #         @inbounds add_to_expression!(rcons[Locb], pos[i][1][l])
            #     else
            #         if ipart == true
            #             pos[i][1][l] = @variable(model, [1:4bs, 1:4bs], PSD)
            #         else
            #             pos[i][1][l] = @variable(model, [1:bs, 1:bs], PSD)
            #         end
            #         for t = 1:bs, r = 1:bs
            #             a=deepcopy(basis[i][1][blocks[i][1][l][t]])
            #             b=deepcopy(basis[i][1][blocks[i][1][l][r]])
            #             @inbounds bi = qtermadd(a,b,n)
            #             if nb > 0
            #                 bi = qreduce_unitnorm(bi;nb=nb)
            #             end
            #             Locb = bfind(tsupp, ltsupp, bi)
            #             if ipart == true
            #                 @inbounds add_to_expression!(rcons[Locb], pos[i][1][l][t,r]+pos[i][1][l][t+bs,r+bs]+pos[i][1][l][t+2*bs,r+2*bs]+pos[i][1][l][t+3*bs,r+3*bs])
            #                 @inbounds add_to_expression!(icons[Locb], pos[i][1][l][t+bs,r]-pos[i][1][l][t,r+bs]+pos[i][1][l][t+3*bs,r+2*bs]-pos[i][1][l][t+2*bs,r+3*bs]) 
            #                 @inbounds add_to_expression!(jcons[Locb], pos[i][1][l][t+2*bs,r]-pos[i][1][l][t,r+2*bs]-pos[i][1][l][t+3*bs,r+bs]+pos[i][1][l][t+bs,r+3*bs])
            #                 @inbounds add_to_expression!(kcons[Locb], pos[i][1][l][t+3*bs,r]-pos[i][1][l][t,r+3*bs]+pos[i][1][l][t+2*bs,r+bs]-pos[i][1][l][t+bs,r+2*bs])
            #             else
            #                 @inbounds add_to_expression!(rcons[Locb], pos[i][1][l][t,r])
            #             end            
            #         end
            #     end
            # end
        end
        # pos[1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, 1)
        # @inbounds bs = length(basis[1])
        # if bs == 1
        #     pos[1][1] = @variable(model, lower_bound=0,base_name="Q")
        #     a=deepcopy(basis[1][1])
        #     b=deepcopy(basis[1][1])
        #     @inbounds bi = qtermadd(a, b,n)
        #     if nb > 0
        #         bi = qreduce_unitnorm(bi;nb=nb)
        #     end
        #     Locb = bfind(tsupp, ltsupp, bi)
        #     @inbounds add_to_expression!(rcons[Locb], pos[1][1])
        # else
        #     if ipart == true
        #         pos[1][1] = @variable(model, [1:4bs, 1:4bs], PSD,base_name="Q")
        #     else
        #         pos[1][1] = @variable(model, [1:bs, 1:bs], PSD,base_name="Q")
        #     end
        #     for t = 1:bs, r = 1:bs
        #         a=deepcopy(basis[1][t])
        #         b=deepcopy(basis[1][r])
        #         @inbounds bi = qtermadd(a,b,n)
        #         if nb > 0
        #             bi = qreduce_unitnorm(bi;nb=nb)
        #         end
        #         Locb = bfind(tsupp, ltsupp, bi)
        #         if ipart == true
        #             @inbounds add_to_expression!(rcons[Locb], pos[1][1][t,r]+pos[1][1][t+bs,r+bs]+pos[1][1][t+2*bs,r+2*bs]+pos[1][1][t+3*bs,r+3*bs])
        #             @inbounds add_to_expression!(icons[Locb], pos[1][1][t+bs,r]-pos[1][1][t,r+bs]+pos[1][1][t+3*bs,r+2*bs]-pos[1][1][t+2*bs,r+3*bs]) 
        #             @inbounds add_to_expression!(jcons[Locb], pos[1][1][t+2*bs,r]-pos[1][1][t,r+2*bs]-pos[1][1][t+3*bs,r+bs]+pos[1][1][t+bs,r+3*bs])
        #             @inbounds add_to_expression!(kcons[Locb], pos[1][1][t+3*bs,r]-pos[1][1][t,r+3*bs]+pos[1][1][t+2*bs,r+bs]-pos[1][1][t+bs,r+2*bs])
        #         else
        #             @inbounds add_to_expression!(rcons[Locb], pos[1][1][t,r])
        #         end            
        #     end
        # end
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
                        @inbounds bi = qtermadd3(at,b,c,n)
                        if nb > 0
                            bi = qreduce_unitnorm(bi;nb=nb)
                        end
                        Locb = bfind(tsupp, ltsupp, bi)
                        @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[i][j+a][l])
                        if ipart == true
                            @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[1], pos[i][j+a][l])
                            @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[2], pos[i][j+a][l])
                            @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], pos[i][j+a][l])
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
                            @inbounds bi = qtermadd3(at,b,c,n)
                            if nb > 0
                                bi = qreduce_unitnorm(bi;nb=nb)
                            end
                            Locb = bfind(tsupp, ltsupp, bi)
                            if ipart == true
                                @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[i][j+a][l][t,r]+pos[i][j+a][l][t+bs,r+bs]+pos[i][j+a][l][t+2*bs,r+2*bs]+pos[i][j+a][l][t+3*bs,r+3*bs])
                                @inbounds add_to_expression!(rcons[Locb], imag_part(coe[k+1][s])[1], pos[i][j+a][l][t+bs,r]-pos[i][j+a][l][t,r+bs]+pos[i][j+a][l][t+3*bs,r+2*bs]-pos[i][j+a][l][t+2*bs,r+3*bs])
                                @inbounds add_to_expression!(rcons[Locb], imag_part(coe[k+1][s])[2], pos[i][j+a][l][t+2*bs,r]-pos[i][j+a][l][t,r+2*bs]-pos[i][j+a][l][t+3*bs,r+bs]+pos[i][j+a][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(rcons[Locb], imag_part(coe[k+1][s])[3], pos[i][j+a][l][t+3*bs,r]-pos[i][j+a][l][t,r+3*bs]+pos[i][j+a][l][t+2*bs,r+bs]-pos[i][j+a][l][t+bs,r+2*bs])
    
                                @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], pos[i][j+a][l][t,r]+pos[i][j+a][l][t+bs,r+bs]+pos[i][j+a][l][t+2*bs,r+2*bs]+pos[i][j+a][l][t+3*bs,r+3*bs])
                                @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), pos[i][j+a][l][t+bs,r]-pos[i][j+a][l][t,r+bs]+pos[i][j+a][l][t+3*bs,r+2*bs]-pos[i][j+a][l][t+2*bs,r+3*bs])
                                @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[3], pos[i][j+a][l][t+2*bs,r]-pos[i][j+a][l][t,r+2*bs]-pos[i][j+a][l][t+3*bs,r+bs]+pos[i][j+a][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[2], pos[i][j+a][l][t+3*bs,r]-pos[i][j+a][l][t,r+3*bs]+pos[i][j+a][l][t+2*bs,r+bs]-pos[i][j+a][l][t+bs,r+2*bs])
                                            
                                @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[2], pos[i][j+a][l][t,r]+pos[i][j+a][l][t+bs,r+bs]+pos[i][j+a][l][t+2*bs,r+2*bs]+pos[i][j+a][l][t+3*bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], pos[i][j+a][l][t+bs,r]-pos[i][j+a][l][t,r+bs]+pos[i][j+a][l][t+3*bs,r+2*bs]-pos[i][j+a][l][t+2*bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), pos[i][j+a][l][t+2*bs,r]-pos[i][j+a][l][t,r+2*bs]-pos[i][j+a][l][t+3*bs,r+bs]+pos[i][j+a][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], pos[i][j+a][l][t+3*bs,r]-pos[i][j+a][l][t,r+3*bs]+pos[i][j+a][l][t+2*bs,r+bs]-pos[i][j+a][l][t+bs,r+2*bs])
                                                
                                @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], pos[i][j+a][l][t,r]+pos[i][j+a][l][t+bs,r+bs]+pos[i][j+a][l][t+2*bs,r+2*bs]+pos[i][j+a][l][t+3*bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[2], pos[i][j+a][l][t+bs,r]-pos[i][j+a][l][t,r+bs]+pos[i][j+a][l][t+3*bs,r+2*bs]-pos[i][j+a][l][t+2*bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[1], pos[i][j+a][l][t+2*bs,r]-pos[i][j+a][l][t,r+2*bs]-pos[i][j+a][l][t+3*bs,r+bs]+pos[i][j+a][l][t+bs,r+3*bs])
                                @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), pos[i][j+a][l][t+3*bs,r]-pos[i][j+a][l][t,r+3*bs]+pos[i][j+a][l][t+2*bs,r+bs]-pos[i][j+a][l][t+bs,r+2*bs])
                            else
                                @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[i][j+a][l][t,r])
                            end                   
                        end
                    end
                end
            end
        end
        epos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, cql)
        for i = 1:cql
            if !isempty(J[i])
                epos[i] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, length(J[i]))
                for (j, k) in enumerate(J[i])
                    bs = length(eblocks[i][j])
                    if bs == 1
                        epos[i][j] = @variable(model)
                        for s = 1:length(supp[k+1])
                            a=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                            b=deepcopy(supp[k+1][s])
                            c=deepcopy(ebasis[i][j][eblocks[i][j][1]])
                            @inbounds bi = qtermadd3(a,b,c,n)
                            if nb > 0
                                bi = qreduce_unitnorm(bi;nb=nb)
                            end
                            Locb = bfind(tsupp, ltsupp, bi)
                            @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j])
                            if ipart == true
                                @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[1], epos[i][j])
                                @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j])
                                @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], epos[i][j])
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
                                @inbounds bi = qtermadd3(a,b,c,n)
                                if nb > 0
                                    bi = qreduce_unitnorm(bi;nb=nb)
                                end
                                Locb = bfind(tsupp, ltsupp, bi)
                                if ipart == true
                                    @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                    @inbounds add_to_expression!(rcons[Locb], imag_part(coe[k+1][s])[1], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                    @inbounds add_to_expression!(rcons[Locb], imag_part(coe[k+1][s])[2], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                    @inbounds add_to_expression!(rcons[Locb], imag_part(coe[k+1][s])[3], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
        
                                    @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                    @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                    @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                    @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[2], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                                
                                    @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                    @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                    @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                    @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                                    
                                    @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], epos[i][j][t,r]+epos[i][j][t+bs,r+bs]+epos[i][j][t+2*bs,r+2*bs]+epos[i][j][t+3*bs,r+3*bs])
                                    @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[2], epos[i][j][t+bs,r]-epos[i][j][t,r+bs]+epos[i][j][t+3*bs,r+2*bs]-epos[i][j][t+2*bs,r+3*bs])
                                    @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[1], epos[i][j][t+2*bs,r]-epos[i][j][t,r+2*bs]-epos[i][j][t+3*bs,r+bs]+epos[i][j][t+bs,r+3*bs])
                                    @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), epos[i][j][t+3*bs,r]-epos[i][j][t,r+3*bs]+epos[i][j][t+2*bs,r+bs]-epos[i][j][t+bs,r+2*bs])
                                else
                                    @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), epos[i][j][t,r])
                                end                   
                            end
                        end
                    end
                end
            end
        end
        # for k=1:m
        #     pos[k+1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, 1)
        #     bs = length(basis[k+1])
        #     if bs == 1
        #         if k <= m-numeq
        #             pos[k+1][1] = @variable(model, lower_bound=0)
        #         else
        #             pos[k+1][1] = @variable(model)
        #         end
        #         for s = 1:length(supp[k+1])
        #             a=deepcopy(basis[k+1][1])
        #             b=deepcopy(supp[k+1][s])
        #             c=deepcopy(basis[k+1][1])
        #             @inbounds bi = qtermadd3(a,b,c,n)
        #             if nb > 0
        #                 bi = qreduce_unitnorm(bi;nb=nb)
        #             end
        #             Locb = bfind(tsupp, ltsupp, bi)
        #             @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[k+1][1])
        #             if ipart == true
        #                 @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[1], pos[k+1][1])
        #                 @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[2], pos[k+1][1])
        #                 @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], pos[k+1][1])
        #             end
        #         end 
        #     else
        #         if ipart == true
        #             if k <= m-numeq
        #                 pos[k+1][1] = @variable(model, [1:4bs, 1:4bs], PSD)
        #             else
        #                 pos[k+1][1] = @variable(model, [1:4bs, 1:4bs], Symmetric)
        #             end
        #         else
        #             if k <= m-numeq
        #                 pos[k+1][1] = @variable(model, [1:bs, 1:bs], PSD)
        #             else
        #                 pos[k+1][1] = @variable(model, [1:bs, 1:bs], Symmetric)
        #             end
        #         end
        #         for t = 1:bs, r = 1:bs
        #             for s = 1:length(supp[k+1])
        #                 a=deepcopy(basis[k+1][t])
        #                 b=deepcopy(supp[k+1][s])
        #                 c=deepcopy(basis[k+1][r])
        #                 @inbounds bi = qtermadd3(a,b,c,n)
        #                 if nb > 0
        #                     bi = qreduce_unitnorm(bi;nb=nb)
        #                 end
        #                 Locb = bfind(tsupp, ltsupp, bi)
        #                 if ipart == true
        #                     @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[k+1][1][t,r]+pos[k+1][1][t+bs,r+bs]+pos[k+1][1][t+2*bs,r+2*bs]+pos[k+1][1][t+3*bs,r+3*bs])
        #                     @inbounds add_to_expression!(rcons[Locb], imag_part(coe[k+1][s])[1], pos[k+1][1][t+bs,r]-pos[k+1][1][t,r+bs]+pos[k+1][1][t+3*bs,r+2*bs]-pos[k+1][1][t+2*bs,r+3*bs])
        #                     @inbounds add_to_expression!(rcons[Locb], imag_part(coe[k+1][s])[2], pos[k+1][1][t+2*bs,r]-pos[k+1][1][t,r+2*bs]-pos[k+1][1][t+3*bs,r+bs]+pos[k+1][1][t+bs,r+3*bs])
        #                     @inbounds add_to_expression!(rcons[Locb], imag_part(coe[k+1][s])[3], pos[k+1][1][t+3*bs,r]-pos[k+1][1][t,r+3*bs]+pos[k+1][1][t+2*bs,r+bs]-pos[k+1][1][t+bs,r+2*bs])

        #                     @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], pos[k+1][1][t,r]+pos[k+1][1][t+bs,r+bs]+pos[k+1][1][t+2*bs,r+2*bs]+pos[k+1][1][t+3*bs,r+3*bs])
        #                     @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), pos[k+1][1][t+bs,r]-pos[k+1][1][t,r+bs]+pos[k+1][1][t+3*bs,r+2*bs]-pos[k+1][1][t+2*bs,r+3*bs])
        #                     @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[3], pos[k+1][1][t+2*bs,r]-pos[k+1][1][t,r+2*bs]-pos[k+1][1][t+3*bs,r+bs]+pos[k+1][1][t+bs,r+3*bs])
        #                     @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[2], pos[k+1][1][t+3*bs,r]-pos[k+1][1][t,r+3*bs]+pos[k+1][1][t+2*bs,r+bs]-pos[k+1][1][t+bs,r+2*bs])
                                        
        #                     @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[2], pos[k+1][1][t,r]+pos[k+1][1][t+bs,r+bs]+pos[k+1][1][t+2*bs,r+2*bs]+pos[k+1][1][t+3*bs,r+3*bs])
        #                     @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], pos[k+1][1][t+bs,r]-pos[k+1][1][t,r+bs]+pos[k+1][1][t+3*bs,r+2*bs]-pos[k+1][1][t+2*bs,r+3*bs])
        #                     @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), pos[k+1][1][t+2*bs,r]-pos[k+1][1][t,r+2*bs]-pos[k+1][1][t+3*bs,r+bs]+pos[k+1][1][t+bs,r+3*bs])
        #                     @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], pos[k+1][1][t+3*bs,r]-pos[k+1][1][t,r+3*bs]+pos[k+1][1][t+2*bs,r+bs]-pos[k+1][1][t+bs,r+2*bs])
                                            
        #                     @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], pos[k+1][1][t,r]+pos[k+1][1][t+bs,r+bs]+pos[k+1][1][t+2*bs,r+2*bs]+pos[k+1][1][t+3*bs,r+3*bs])
        #                     @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[2], pos[k+1][1][t+bs,r]-pos[k+1][1][t,r+bs]+pos[k+1][1][t+3*bs,r+2*bs]-pos[k+1][1][t+2*bs,r+3*bs])
        #                     @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[1], pos[k+1][1][t+2*bs,r]-pos[k+1][1][t,r+2*bs]-pos[k+1][1][t+3*bs,r+bs]+pos[k+1][1][t+bs,r+3*bs])
        #                     @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), pos[k+1][1][t+3*bs,r]-pos[k+1][1][t,r+3*bs]+pos[k+1][1][t+2*bs,r+bs]-pos[k+1][1][t+bs,r+2*bs])
        #                 else
        #                     @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[k+1][1][t,r])
        #                 end                   
        #             end
        #         end
        #     end
        # end
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
        # itsupp = nothing
        # ind = [item[1] != item[2] for item in tsupp]
        # itsupp = tsupp[ind]
        # icons = icons[ind]
        # jcons = jcons[ind]
        # kcons = kcons[ind]
        # ibc = zeros(length(itsupp))
        # jbc = zeros(length(itsupp))
        # kbc = zeros(length(itsupp))
        # ncons += length(itsupp)
        if QUIET == false
            println("There are $ncons affine constraints.")
        end
        for i = 1:length(supp[1])
            Locb = bfind(tsupp, ltsupp, supp[1][i])
            if Locb === nothing
                # println(supp[1][i])
                @error "The monomial basis is not enough!"
                return nothing,ksupp,nothing,nothing,nothing
            else
                rbc[Locb] = rbc[Locb]+real(coe[1][i])
                if ipart == true
                    ibc[Locb] = imag_part(coe[1][i])[1]
                    jbc[Locb] = imag_part(coe[1][i])[2]
                    kbc[Locb] = imag_part(coe[1][i])[3]
                end
            end
        end
        @variable(model, lower)
        rcons[1] += lower
        @constraint(model, rcon[i=1:ltsupp], rcons[i]==rbc[i])
        if ipart == true
            @constraint(model, icon[i=1:ltsupp], icons[i]==ibc[i])
            @constraint(model, jcon[i=1:ltsupp], jcons[i]==jbc[i])
            @constraint(model, kcon[i=1:ltsupp], kcons[i]==kbc[i])
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
        if writetofile != false
            write_to_file(dualize(model), writetofile)
        end
        SDP_status = termination_status(model)
        cons=all_constraints(model; include_variable_in_set_constraints = false)
        # println(cons)
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
    pop = Vector{Polynomial{true,Float64}}(undef, length(qpop))
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
                    temp1 .=  [ x > UInt16(n/2) ? x - UInt16(n/2) : x for x in temp1]
                    temp2 .=  [ x <= UInt16(n/2) ? x + UInt16(n/2) : x for x in temp2]
                    # println(temp)
                    add_clique!(G, unique([supp[i][j][1];supp[i][j][2];temp1;temp2;supp[i][j][3]]))
                    # add_clique!(G, unique([supp[i][j][1];supp[i][j][2];supp[i][j][3]]))
                end
            else
                temp = copy([supp[i][1][1];supp[i][1][2];supp[i][1][3]])
                for j = 2:length(supp[i])
                    append!(temp, [supp[i][j][1];supp[i][j][2];supp[i][j][3]])
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
function get_blocks(m, l, d, tsupp, supp::Vector{Vector{Vector{Vector{UInt16}}}}, basis, ebasis; nb=0, normality=1, nvar=0, TS="block", ConjugateBasis=false, merge=false, md=3)
    if (ConjugateBasis == false && normality > 0) || (ConjugateBasis == true && normality >= d)
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
            elseif ((ConjugateBasis == false && normality > 0) || (ConjugateBasis == true && normality >= d)) && k <= 1 + nvar
                G = get_graph(tsupp, basis[k], nb=nb, ConjugateBasis=true)
            elseif ((ConjugateBasis == false && normality > 0) || (ConjugateBasis == true && normality >= d)) && k > 1 + nvar
                G = get_graph(tsupp, supp[k-1-nvar], basis[k], nb=nb, ConjugateBasis=ConjugateBasis)
            else
                G = get_graph(tsupp, supp[k-1], basis[k], nb=nb, ConjugateBasis=ConjugateBasis)
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
                # 临时计算修改后的第三元素（不实际修改原数据）
                modified_3rd_l = [x > UInt16(n) ? x - UInt16(n) : x for x in tsupp[j][3]]
                modified_3rd_g = [x <= UInt16(n) ? x + UInt16(n) : x for x in tsupp[j][3]]
                issubset(union(tsupp[j][1], tsupp[j][2], modified_3rd_l, modified_3rd_g), cliques[i])
            end 
            for j in eachindex(tsupp)
        ]]
        # ksupp = TS == false ? nothing : tsupp[[issubset(union(tsupp[j][1], tsupp[j][2], tsupp[j][3].=[x > UInt16(n) ? x - UInt16(n) : x for x in tsupp[j][3]]), cliques[i]) for j in eachindex(tsupp)]]
        # println(ksupp)
        blocks[i],eblocks[i],cl[i],blocksize[i] = get_blocks(length(I[i]), length(J[i]), rlorder[i], ksupp, supp[[I[i]; J[i]].+1], basis[i], 
        ebasis[i], TS=TS, ConjugateBasis=ConjugateBasis, merge=merge, md=md, nb=nb, normality=normality, nvar=cliquesize[i])
    end
    return blocks,eblocks,cl,blocksize
end
function assign_constraint(m, numeq, supp::Vector{Vector{Vector{Vector{UInt16}}}}, cliques, cql)
    I = [Int[] for i=1:cql]
    J = [Int[] for i=1:cql]
    ncc = Int[]
    for i = 1:m
        ind = findall(k->issubset(unique(reduce(vcat, [[item[1];item[2];item[3]] for item in supp[i+1]])), cliques[k]), 1:cql)
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
