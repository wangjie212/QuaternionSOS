using TSSOS
using DynamicPolynomials
using Quaternions
using MosekTools
using JuMP
using LinearAlgebra
include("function.jl")
include("ncutils.jl")

function qs_tssos_first(pop::Vector{Polynomial{false,T}}, z, n, d; numeq=0, RemSig=false, nb=0, CS="MF", cliques=[], minimize=false, 
    TS=false, merge=false, md=3, solver="Mosek", reducebasis=false, QUIET=false, solve=true, tune=false, solution=false, ipart=true, 
    dualize=false, balanced=false, MomentOne=false, Gram=false, Mommat=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), 
    writetofile=false, normality=0, NormalSparse=false,conjubasis=true) where {T<:Number}
    supp,coe = Qpolys_info(pop, z, n)
    opt= qs_tssos_first(supp, coe, n, d, numeq=numeq, RemSig=RemSig, nb=nb, CS=CS, cliques=cliques, minimize=minimize,
    TS=TS, merge=merge, md=md, solver=solver, reducebasis=reducebasis, QUIET=QUIET, solve=solve, tune=tune, 
    solution=solution, ipart=ipart, dualize=dualize, balanced=balanced, MomentOne=MomentOne, Gram=Gram, Mommat=Mommat, 
    cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, writetofile=writetofile, normality=normality, NormalSparse=NormalSparse,conjubasis=conjubasis,cpop=pop, z=z)
    return opt
end

function qs_tssos_first(supp::Vector{Vector{Vector{Vector{UInt16}}}}, coe, n, d; numeq=0, RemSig=false, nb=0, CS="MF", cliques=[], 
    minimize=false, TS=false, merge=false, md=3, solver="Mosek", reducebasis=false, QUIET=false, solve=true, tune=false, solution=false, 
    ipart=true, dualize=false, balanced=false, MomentOne=false, Gram=false, Mommat=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), 
    writetofile=false, normality=0, NormalSparse=false,conjubasis=true,cpop=nothing, z=nothing)
    println("*********************************** TSSOS ***********************************")
    println("NCTSSOS is launching...")
    if nb > 0
        supp[1],coe[1] = qresort(supp[1],coe[1];nb=nb)
    end
    m = length(supp) - 1
    dc = zeros(Int, m)
    for i = 1:m       
        # dc[i] = maximum([max(length(supp[i+1][j][1])+length(supp[i+1][j][3][supp[i+1][j][3].<=UInt16(n)]), length(supp[i+1][j][2])+length(supp[i+1][j][3])-length(supp[i+1][j][3][supp[i+1][j][3].<=UInt16(n)])) for j=1:length(supp[i+1])])
        dc[i] = maximum([length(supp[i+1][j][1])+length(supp[i+1][j][2])+length(supp[i+1][j][3]) for j=1:length(supp[i+1])])
    end
    # cliques = [UInt16[i for i=1:n]]
    # cql=1
    # Ic,Jc,ncc = assign_constraint2(m, numeq, supp, cliques, cql)
    # println(I)
    # if d == "min"
    #     rlorder = isempty(Ic) && isempty(Jc) ? 1 : maximum(dc[[Ic; Jc]])
    # else
    #     rlorder = d
    # end
    rlorder=d
    basis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, m+1)
    # hbasis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, numeq)
    basis[1]=get_qncbasis(n, rlorder;conjubasis=conjubasis)
    for s = 1:m
            basis[s+1] = get_qncbasis(n, rlorder-Int(ceil(dc[s]/2));conjubasis=conjubasis)
    end
    if nb > 0
        for s=1:m+1
        basis[s]= qreduce_unitnorm.(basis[s], nb=nb)
        unique!(basis[s])
        end
    end
    # for s = 1:length(Jc)
    #         # temp = get_basis(cliques[1], rlorder-dc[J[1][s]])
    #         # temp = get_qncbasis(n, rlorder-Int(ceil(dc[Jc[s]]/2)))
    #         # hbasis[s] = vec([[item1, item2] for item1 in temp, item2 in temp])
    #         hbasis[s] = get_qncbasis(n, rlorder-Int(ceil(dc[Jc[s]]/2)))
    #         # sort!(hbasis[s])
    # end  
    # tsupp = copy(supp[1])
    # for i = 2:m+1, j = 1:length(supp[i])
    #     if supp[i][j][1] <= supp[i][j][2]
    #         push!(tsupp, supp[i][j])
    #     end
    # end
    # sort!(tsupp)
    # unique!(tsupp)
    opt,ksupp,SDP_status= qsolvesdp(n, m, rlorder, supp, coe, basis, numeq=numeq, QUIET=QUIET, TS=TS, solver=solver, solve=solve, tune=tune, solution=solution, ipart=ipart, MomentOne=MomentOne, balanced=balanced,
    Gram=Gram, Mommat=Mommat, nb=nb, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, dualize=dualize, writetofile=writetofile, normality=normality, NormalSparse=NormalSparse,conjubasis=conjubasis)
    return opt
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
function qsolvesdp(n, m, rlorder, supp::Vector{Vector{Vector{Vector{UInt16}}}}, coe, basis; 
    numeq=0, nb=0, QUIET=false, TS=false, solver="Mosek", tune=false, solve=true, dualize=false, solution=false, Gram=false, MomentOne=false, ipart=true, Mommat=false, 
    cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false,balanced=false, normality=0, NormalSparse=false,conjubasis=false)
    tsupp = Vector{Vector{UInt16}}[]
    # println(basis[1])
    # for k = 1:length(basis[1]), r = k:length(basis[1])
    for k = 1:length(basis[1]), r = 1:length(basis[1])
        # @inbounds bi = [basis[1][k], basis[1][r]]
        a=deepcopy(basis[1][k])
        b=deepcopy(basis[1][r])
        @inbounds bi = qtermadd(a,b,n)
        if nb > 0
            bi = qreduce_unitnorm(bi, nb=nb)
        end
        push!(tsupp,bi)
    end
    # println(tsupp)
    ksupp = deepcopy(tsupp)
    if normality > 0
        wbasis = [[UInt16[],UInt16[],UInt16[]]]
        wbasis= get_qncbasis(n, normality)
        temp =get_qncbasis(n, 1)
        bs = length(wbasis)
        for i = 2:n+1
            if normality >= rlorder
                for j = 1:bs, k = 1:bs
                    a=deepcopy(temp[i])
                    b=deepcopy(temp[i])
                    c=deepcopy(wbasis[j])
                    d=deepcopy(wbasis[k])
                    bi =qtermadd3(c,qtermadd(a,b,n),d,n)
                    push!(tsupp, bi)
                end
                for j = 1:bs, k = 1:bs
                    a=deepcopy(temp[i])
                    b=deepcopy(wbasis[j])
                    c=deepcopy(wbasis[k])
                    bi = qtermadd3(b,a,c,n)
                    push!(tsupp, bi)
                end
                for j = 1:bs, k = 1:bs
                    a=deepcopy(temp[i+n])
                    b=deepcopy(wbasis[j])
                    c=deepcopy(wbasis[k])
                    bi = qtermadd3(b,a,c,n)
                    push!(tsupp, bi)
                end
            end
        end
    end
    sort!(tsupp)
    unique!(tsupp)
    #initial set
    objv = moment = GramMat = SDP_status= nothing
    if solve == true
        ltsupp = length(tsupp)
        # println(ltsupp)
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
        icons = [AffExpr(0) for i=1:ltsupp]
        jcons = [AffExpr(0) for i=1:ltsupp]
        kcons = [AffExpr(0) for i=1:ltsupp]
        if normality > 0
            bs = length(wbasis)
            for i = 2:n+1
                if NormalSparse == false
                        # if ipart == true
                        #     hnom = @variable(model, [1:4bs, 1:4bs], PSD)
                        # else
                    hnom = @variable(model, [1:2bs, 1:2bs], PSD)
                        # end
                    for j = 1:bs, k = j:bs
                        a=deepcopy(wbasis[j])
                        b=deepcopy(wbasis[k])
                        bi =qtermadd(a,b,n)
                            # if bi[1] <= bi[2]
                        Locb = bfind(tsupp, ltsupp, bi)
                                # if bi[1] != bi[2] && ipart == true
                                #     @inbounds add_to_expression!(icons[Locb], hnom[j,k+2bs]-hnom[k,j+2bs])
                                # end
                            # else
                                # Locb = bfind(tsupp, ltsupp, bi[2:-1:1])
                                # if ipart == true
                                #     @inbounds add_to_expression!(icons[Locb], -1, hnom[j,k+2bs]-hnom[k,j+2bs])
                                # end
                            # end
                            # if ipart == true
                            #     @inbounds add_to_expression!(rcons[Locb], hnom[j,k]+hnom[j+2bs,k+2bs])
                            # else
                        if j==k
                            @inbounds add_to_expression!(rcons[Locb], hnom[j,k])
                        else
                            @inbounds add_to_expression!(rcons[Locb], 2, hnom[j,k])
                        end
                            # end
                        a=deepcopy(temp[i])
                        b=deepcopy(temp[i])
                        c=deepcopy(wbasis[j])
                        d=deepcopy(wbasis[k])
                        bi =qtermadd3(c,qtermadd(a,b,n),d,n)
                            # a=deepcopy(temp[i])
                            # b=deepcopy(wbasis[j])
                            # c=deepcopy(wbasis[k])
                            # bi = qtermadd3(a,b,c,n)
                            # bi = [sadd(wbasis[s][j], [cliques[s][i]]), sadd(wbasis[s][k], [cliques[s][i]])]
                            # if bi[1] <= bi[2]
                        Locb = bfind(tsupp, ltsupp, bi)
                        # @inbounds add_to_expression!(rcons[Locb], hnom[j+bs,k+bs])
                        if j==k
                            @inbounds add_to_expression!(rcons[Locb], hnom[j+bs,k+bs])
                        else
                            @inbounds add_to_expression!(rcons[Locb], 2, hnom[j+bs,k+bs])
                        end
                    end
                    # for j = 1:bs, k = 1:bs
                    #     a=deepcopy(temp[i+n])
                    #     b=deepcopy(wbasis[j])
                    #     c=deepcopy(wbasis[k])
                    #     bi = qtermadd3(a,b,c,n)
                    #         # bi = [sadd(wbasis[s][j], [cliques[s][i]]), wbasis[s][k]]
                    #         # if bi[1] <= bi[2]
                    #     Locb = bfind(tsupp, ltsupp, bi)
                    #             # if bi[1] != bi[2]
                    #                 # if ipart == true
                    #                 #     @inbounds add_to_expression!(icons[Locb], hnom[j,k+3bs]-hnom[k+bs,j+2bs])
                    #                 #     @inbounds add_to_expression!(rcons[Locb], hnom[j,k+bs]+hnom[j+2bs,k+3bs])
                    #                 # else
                    #     @inbounds add_to_expression!(rcons[Locb], hnom[j,k+bs])
                    #                 # end
                    #             # else
                    #                 # if ipart == true
                    #                     # @inbounds add_to_expression!(rcons[Locb], 2, hnom[j,k+bs]+hnom[j+2bs,k+3bs])
                    #                 # else
                    #                     # @inbounds add_to_expression!(rcons[Locb], 2, hnom[j,k+bs])
                    #                 # end
                    #             # end
                    #         # else
                    #         #     Locb = bfind(tsupp, ltsupp, bi[2:-1:1])
                    #         #     if ipart == true
                    #         #         @inbounds add_to_expression!(icons[Locb], -1, hnom[j,k+3bs]-hnom[k+bs,j+2bs])
                    #         #         @inbounds add_to_expression!(rcons[Locb], hnom[j,k+bs]+hnom[j+2bs,k+3bs])
                    #         #     else
                    #         #         @inbounds add_to_expression!(rcons[Locb], hnom[j,k+bs])
                    #         #     end
                    #         # end
                    # end
                    for j = 1:bs, k = 1:bs
                        a=deepcopy(temp[i+n])
                        b=deepcopy(wbasis[j])
                        c=deepcopy(wbasis[k])
                        bi = qtermadd3(b,a,c,n)
                        Locb = bfind(tsupp, ltsupp, bi)
                        @inbounds add_to_expression!(rcons[Locb], hnom[j,k+bs])
                    end
                    for j = 1:bs, k = 1:bs
                        a=deepcopy(temp[i])
                        b=deepcopy(wbasis[j])
                        c=deepcopy(wbasis[k])
                        bi = qtermadd3(b,a,c,n)
                        Locb = bfind(tsupp, ltsupp, bi)
                        @inbounds add_to_expression!(rcons[Locb], hnom[j+bs,k])
                    end
                end
            end
        end
        # pos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m-numeq+1)
        pos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m+1)
        pos[1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, 1)
            @inbounds bs = length(basis[1])
            println(bs)
            if bs == 1
                pos[1][1] = @variable(model, lower_bound=0,base_name="Q")
                a=deepcopy(basis[1][1])
                b=deepcopy(basis[1][1])
                @inbounds bi = qtermadd(a, b,n)
                if nb > 0
                    bi = qreduce_unitnorm(bi;nb=nb)
                end
                # println(bi)
                Locb = bfind(tsupp, ltsupp, bi)
                @inbounds add_to_expression!(rcons[Locb], pos[1][1])
            else
                # pos[1][1] = @variable(model, [1:4bs, 1:4bs], PSD,base_name="Q")
                pos[1][1] = @variable(model, [1:bs, 1:bs], PSD,base_name="Q")
                Temp=[[UInt16[],UInt16[],UInt16[]]]
                for t = 1:bs, r = 1:bs
                    a=deepcopy(basis[1][t])
                    b=deepcopy(basis[1][r])
                    @inbounds bi = qtermadd(a,b,n)
                    if nb > 0
                        bi = qreduce_unitnorm(bi;nb=nb)
                    end
                    # if bi==[UInt16[],UInt16[],UInt16[]]||!(bi in Temp) 
                    # println(bi)
                    #     push!(Temp,bi)
                    #     Locb = bfind(tsupp, ltsupp, bi)
                    #     @inbounds add_to_expression!(rcons[Locb], pos[1][1][t,r])
                    # end
                    # push!(Temp,bi)
                    Locb = bfind(tsupp, ltsupp, bi)
                    # println(Locb)
                    # println(basis[1][t],basis[1][r],bi,Locb)
                    # if t == r
                        @inbounds add_to_expression!(rcons[Locb], pos[1][1][t,r])
                        # @inbounds add_to_expression!(rcons[Locb], pos[1][1][t,r]+pos[1][1][t+bs,r+bs]+pos[1][1][t+2*bs,r+2*bs]+pos[1][1][t+3*bs,r+3*bs])
                        # @inbounds add_to_expression!(icons[Locb], pos[1][1][t+bs,r]-pos[1][1][t,r+bs]+pos[1][1][t+3*bs,r+2*bs]-pos[1][1][t+2*bs,r+3*bs]) 
                        # @inbounds add_to_expression!(jcons[Locb], pos[1][1][t+2*bs,r]-pos[1][1][t,r+2*bs]-pos[1][1][t+3*bs,r+bs]+pos[1][1][t+bs,r+3*bs])
                        # @inbounds add_to_expression!(kcons[Locb], pos[1][1][t+3*bs,r]-pos[1][1][t,r+3*bs]+pos[1][1][t+2*bs,r+bs]-pos[1][1][t+bs,r+2*bs])
                    # else
                        # @inbounds add_to_expression!(rcons[Locb], 2,pos[1][1][t,r])
                        # @inbounds add_to_expression!(rcons[Locb], 2,pos[1][1][t,r]+pos[1][1][t+bs,r+bs]+pos[1][1][t+2*bs,r+2*bs]+pos[1][1][t+3*bs,r+3*bs])
                        # @inbounds add_to_expression!(icons[Locb], 2,pos[1][1][t+bs,r]-pos[1][1][t,r+bs]+pos[1][1][t+3*bs,r+2*bs]-pos[1][1][t+2*bs,r+3*bs]) 
                        # @inbounds add_to_expression!(jcons[Locb], 2,pos[1][1][t+2*bs,r]-pos[1][1][t,r+2*bs]-pos[1][1][t+3*bs,r+bs]+pos[1][1][t+bs,r+3*bs])
                        # @inbounds add_to_expression!(kcons[Locb], 2,pos[1][1][t+3*bs,r]-pos[1][1][t,r+3*bs]+pos[1][1][t+2*bs,r+bs]-pos[1][1][t+bs,r+2*bs])
                    # end               
                end
            end
        for k=1:m
            pos[k+1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, 1)
            bs = length(basis[k+1])
            if bs == 1
                if k <= m-numeq
                    pos[k+1][1] = @variable(model, lower_bound=0)
                else
                    pos[k+1][1] = @variable(model)
                end
                for s = 1:length(supp[k+1])
                    a=deepcopy(basis[k+1][1])
                    b=deepcopy(supp[k+1][s])
                    c=deepcopy(basis[k+1][1])
                    @inbounds bi = qtermadd3(a,b,c,n)
                    if nb > 0
                        bi = qreduce_unitnorm(bi;nb=nb)
                    end
                    # @inbounds bi = qtermadd3(a,c,b,n)
                    # println(bi)
                    # println(tsupp)
                    Locb = bfind(tsupp, ltsupp, bi)
                    @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[k+1][1])
                end 
            else
                if k <= m-numeq
                    pos[k+1][1] = @variable(model, [1:bs, 1:bs], PSD)
                else
                    pos[k+1][1] = @variable(model, [1:bs, 1:bs], Symmetric)
                end
                for t = 1:bs, r = 1:bs
                    Temp=[[UInt16[],UInt16[],UInt16[]]]
                    for s = 1:length(supp[k+1])
                        a=deepcopy(basis[k+1][t])
                        b=deepcopy(supp[k+1][s])
                        c=deepcopy(basis[k+1][r])
                        @inbounds bi = qtermadd3(a,b,c,n)
                        if nb > 0
                            bi = qreduce_unitnorm(bi;nb=nb)
                        end
                        Locb = bfind(tsupp, ltsupp, bi)
                        # if Locb==nothing
                        #     println(basis[k+1][t],supp[k+1][s],basis[k+1][r],bi)
                        # end
                        @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[k+1][1][t,r])
                    end
                end
            end
        end
        # for (j, w) in enumerate(Ic)
        #     pos[j+1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, 1)
        #         bs = length(basis[j+1])
        #         # println(bs)
        #         if bs == 1
        #             pos[j+1][1] = @variable(model, lower_bound=0,base_name="X")
        #             for s = 1:length(supp[w+1])
        #                 a=deepcopy(basis[j+1][1])
        #                 b=deepcopy(supp[w+1][s])
        #                 c=deepcopy(basis[j+1][1])
        #                 @inbounds bi = qtermadd3(a,b,c,n)
        #                 println(bi)
        #                 Locb = bfind(tsupp, ltsupp, bi)
        #                 @inbounds add_to_expression!(rcons[Locb], real(coe[w+1][s]), pos[j+1][1])
        #                 # @inbounds add_to_expression!(icons[Locb], imag_part(coe[w+1][s])[1], pos[j+1][1])
        #                 # @inbounds add_to_expression!(jcons[Locb], imag_part(coe[w+1][s])[2], pos[j+1][1])
        #                 # @inbounds add_to_expression!(kcons[Locb], imag_part(coe[w+1][s])[3], pos[j+1][1])
        #             end 
        #         else
        #             # pos[j+1][1] = @variable(model, [1:4bs, 1:4bs], PSD,base_name="X")
        #             pos[j+1][1] = @variable(model, [1:bs, 1:bs], PSD,base_name="X")
        #             for t = 1:bs, r = 1:bs
        #                 for s = 1:length(supp[w+1])
        #                     a=deepcopy(basis[j+1][t])
        #                     b=deepcopy(supp[w+1][s])
        #                     c=deepcopy(basis[j+1][r])
        #                     @inbounds bi = qtermadd3(a,b,c,n)
        #                     Locb = bfind(tsupp, ltsupp, bi)
        #                     # println(j,',',t,',',w,',',s,',',r,basis[j+1][t],supp[w+1][s],basis[j+1][r],Locb)
        #                     # if t == r
        #                     @inbounds add_to_expression!(rcons[Locb], real(coe[w+1][s]), pos[j+1][1][t,r])
        #                         # @inbounds add_to_expression!(rcons[Locb], real(coe[w+1][s]), pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
        #                         # @inbounds add_to_expression!(rcons[Locb], imag_part(coe[w+1][s])[1], pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
        #                         # @inbounds add_to_expression!(rcons[Locb], imag_part(coe[w+1][s])[2], pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
        #                         # @inbounds add_to_expression!(rcons[Locb], imag_part(coe[w+1][s])[3], pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])

        #                         # @inbounds add_to_expression!(icons[Locb], -imag_part(coe[w+1][s])[1], pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
        #                         # @inbounds add_to_expression!(icons[Locb], real(coe[w+1][s]), pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
        #                         # @inbounds add_to_expression!(icons[Locb], -imag_part(coe[w+1][s])[3], pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
        #                         # @inbounds add_to_expression!(icons[Locb], imag_part(coe[w+1][s])[2], pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
                                        
        #                         # @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[w+1][s])[2], pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
        #                         # @inbounds add_to_expression!(jcons[Locb], imag_part(coe[w+1][s])[3], pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
        #                         # @inbounds add_to_expression!(jcons[Locb], real(coe[w+1][s]), pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
        #                         # @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[w+1][s])[1], pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
                                        
        #                         # @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[w+1][s])[3], pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
        #                         # @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[w+1][s])[2], pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
        #                         # @inbounds add_to_expression!(kcons[Locb], imag_part(coe[w+1][s])[1], pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
        #                         # @inbounds add_to_expression!(kcons[Locb], real(coe[w+1][s]), pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
        #                     # else
        #                     #     @inbounds add_to_expression!(rcons[Locb], 2*real(coe[w+1][s]), pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
        #                     #     @inbounds add_to_expression!(rcons[Locb], 2*imag_part(coe[w+1][s])[1], pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
        #                     #     @inbounds add_to_expression!(rcons[Locb], 2*imag_part(coe[w+1][s])[2], pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
        #                     #     @inbounds add_to_expression!(rcons[Locb], 2*imag_part(coe[w+1][s])[3], pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
                                
        #                     #     @inbounds add_to_expression!(icons[Locb], 2*-imag_part(coe[w+1][s])[1], pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
        #                     #     @inbounds add_to_expression!(icons[Locb], 2*real(coe[w+1][s]), pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
        #                     #     @inbounds add_to_expression!(icons[Locb], 2*-imag_part(coe[w+1][s])[3], pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
        #                     #     @inbounds add_to_expression!(icons[Locb], 2*imag_part(coe[w+1][s])[2], pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
                                        
        #                     #     @inbounds add_to_expression!(jcons[Locb], 2*-imag_part(coe[w+1][s])[2], pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
        #                     #     @inbounds add_to_expression!(jcons[Locb], 2*imag_part(coe[w+1][s])[3], pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
        #                     #     @inbounds add_to_expression!(jcons[Locb], 2*real(coe[w+1][s]), pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
        #                     #     @inbounds add_to_expression!(jcons[Locb], 2*-imag_part(coe[w+1][s])[1], pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
                                        
        #                     #     @inbounds add_to_expression!(kcons[Locb], 2*-imag_part(coe[w+1][s])[3], pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
        #                     #     @inbounds add_to_expression!(kcons[Locb], 2*-imag_part(coe[w+1][s])[2], pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
        #                     #     @inbounds add_to_expression!(kcons[Locb], 2*imag_part(coe[w+1][s])[1], pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
        #                     #     @inbounds add_to_expression!(kcons[Locb], 2*real(coe[w+1][s]), pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
        #                     # end
        #                     # @inbounds bi = [sadd(basis[j+1][t], supp[w+1][s][1]), sadd(basis[j+1][r], supp[w+1][s][2])]
        #                     # if bi[1] <= bi[2]
        #                     #     Locb = bfind(tsupp, ltsupp, bi)
        #                     #     # @inbounds add_to_expression!(rcons[Locb], real(coe[w+1][s]), pos[j+1][1][t,r])
        #                     #         @inbounds add_to_expression!(rcons[Locb], real(coe[w+1][s]), pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
        #                     #         @inbounds add_to_expression!(rcons[Locb], imag_part(coe[w+1][s])[1], pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
        #                     #         @inbounds add_to_expression!(rcons[Locb], imag_part(coe[w+1][s])[2], pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
        #                     #         @inbounds add_to_expression!(rcons[Locb], imag_part(coe[w+1][s])[3], pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
        #                     #         if bi[1] != bi[2]
        #                     #             @inbounds add_to_expression!(icons[Locb], -imag_part(coe[w+1][s])[1], pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
        #                     #             @inbounds add_to_expression!(icons[Locb], real(coe[w+1][s]), pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
        #                     #             @inbounds add_to_expression!(icons[Locb], -imag_part(coe[w+1][s])[3], pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
        #                     #             @inbounds add_to_expression!(icons[Locb], imag_part(coe[w+1][s])[2], pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
                                        
        #                     #             @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[w+1][s])[2], pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
        #                     #             @inbounds add_to_expression!(jcons[Locb], imag_part(coe[w+1][s])[3], pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
        #                     #             @inbounds add_to_expression!(jcons[Locb], real(coe[w+1][s]), pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
        #                     #             @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[w+1][s])[1], pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
                                        
        #                     #             @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[w+1][s])[3], pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
        #                     #             @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[w+1][s])[2], pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
        #                     #             @inbounds add_to_expression!(kcons[Locb], imag_part(coe[w+1][s])[1], pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
        #                     #             @inbounds add_to_expression!(kcons[Locb], real(coe[w+1][s]), pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
        #                     #         end
        #                     # end
        #                 end
        #             end
        #         end
        # end
        # if numeq > 0
        #     free = Vector{Vector{Vector{VariableRef}}}(undef, 1)
        #         if !isempty(Jc)
        #             free[1] = Vector{Vector{VariableRef}}(undef, length(Jc))
        #             for (j, w) in enumerate(Jc)
        #                 mons = hbasis[j][Vector(1:length(hbasis[j]))]
        #                 # temp = mons[[item[1] <= item[2] for item in mons]]
        #                 # lb = length(temp)
        #                 lb=length(mons)
        #                 free[1][j] = @variable(model, [1:lb])
        #                 # free[1][j] = @variable(model, [1:4*lb])
        #                 for k = 1:lb
        #                     for s = 1:length(supp[w+1])
        #                         a=deepcopy(hbasis[j][k])
        #                         b=deepcopy(supp[w+1][s])
        #                         c=deepcopy(hbasis[j][k])
        #                         @inbounds bi = qtermadd3(a,b,c,n)
        #                 # for k in Vector(1:length(hbasis[j])), s = 1:length(supp[w+1])
        #                     # @inbounds bi = [sadd(hbasis[j][k][1], supp[w+1][s][1]), sadd(hbasis[j][k][2], supp[w+1][s][2])]
        #                     # if bi[1] <= bi[2]
        #                         Locb = bfind(tsupp, ltsupp, bi)
        #                         loc=bfind(mons,lb,hbasis[j][k])
        #                         # if hbasis[j][k][1] <= hbasis[j][k][2]
        #                             # loc = bfind(temp, lb, hbasis[j][k])
        #                             # tag = 1
        #                             # if hbasis[j][k][1] == hbasis[j][k][2]
        #                                 # tag = 0
        #                             # end
        #                         # else
        #                             # loc = bfind(temp, lb, hbasis[j][k][2:-1:1])
        #                             # tag = -1
        #                         # end
        #                         @inbounds add_to_expression!(rcons[Locb], real(coe[w+1][s])*free[1][j][loc])
        #                         # @inbounds add_to_expression!(rcons[Locb], real(coe[w+1][s])*free[1][j][loc]-imag_part(coe[w+1][s])[1]*free[1][j][loc+lb]-imag_part(coe[w+1][s])[2]*free[1][j][loc+2*lb]-imag_part(coe[w+1][s])[3]*free[1][j][loc+3*lb])
        #                             # if bi[1] != bi[2]
        #                         # @inbounds add_to_expression!(icons[Locb], real(coe[w+1][s])*free[1][j][loc+lb]+imag_part(coe[w+1][s])[1]*free[1][j][loc]+imag_part(coe[w+1][s])[2]*free[1][j][loc+3*lb]-imag_part(coe[w+1][s])[3]*free[1][j][loc+2*lb])
        #                         # @inbounds add_to_expression!(jcons[Locb], real(coe[w+1][s])*free[1][j][loc+2*lb]-imag_part(coe[w+1][s])[1]*free[1][j][loc+3*lb]+imag_part(coe[w+1][s])[2]*free[1][j][loc]+imag_part(coe[w+1][s])[3]*free[1][j][loc+lb])
        #                         # @inbounds add_to_expression!(kcons[Locb], real(coe[w+1][s])*free[1][j][loc+3*lb]+imag_part(coe[w+1][s])[1]*free[1][j][loc+2*lb]-imag_part(coe[w+1][s])[2]*free[1][j][loc+lb]-imag_part(coe[w+1][s])[3]*free[1][j][loc])
        #                     end    
        #                             # end
        #                     # end
        #                 end
        #             end
        #         end
        # end
        rbc = zeros(ltsupp)
        # ibc = zeros(length(ltsupp))
        # jbc = zeros(length(ltsupp))
        # kbc = zeros(length(ltsupp))
        # ibc = zeros(ltsupp)
        # jbc = zeros(ltsupp)
        # kbc = zeros(ltsupp)
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
                @error "The monomial basis is not enough!"
                return nothing,ksupp,nothing,nothing,nothing
            else
            rbc[Locb] = rbc[Locb]+real(coe[1][i])
                # if supp[1][i][1] != supp[1][i][2]
            # Locb = bfind(itsupp, length(itsupp), supp[1][i])
            # ibc[Locb] = imag_part(coe[1][i])[1]
            # jbc[Locb] = imag_part(coe[1][i])[2]
            # kbc[Locb] = imag_part(coe[1][i])[3]
                # end
            end
        end
        @variable(model, lower)
        rcons[1] += lower
        @constraint(model, rcon[i=1:ltsupp], rcons[i]==rbc[i])
        # @constraint(model, icon[i=1:ltsupp], icons[i]==ibc[i])
        # @constraint(model, jcon[i=1:ltsupp], jcons[i]==jbc[i])
        # @constraint(model, kcon[i=1:ltsupp], kcons[i]==kbc[i])
        # @constraint(model, icon[i=1:length(itsupp)], icons[i]==ibc[i])
        # @constraint(model, jcon[i=1:length(itsupp)], jcons[i]==jbc[i])
        # @constraint(model, kcon[i=1:length(itsupp)], kcons[i]==kbc[i])
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

# function quaternion_to_real(qpop, q; symmetry=true)
#     n = Int(length(q)/2)
#     @polyvar x[1:4n]
#     if symmetry == false
#         pop = Vector{Polynomial{true,Float64}}(undef, 4*length(qpop))
#     else
#         pop = Vector{Polynomial{true,Float64}}(undef, length(qpop))
#     end
#     I=quat(0,1,0,0)
#     J=quat(0,0,1,0)
#     K=quat(0,0,0,1)
#     for (i,qp) in enumerate(qpop)
#         temp = qp(q[1:n]=>x[1:n]+I*x[n+1:2n]+J*x[2n+1:3n]+K*x[3n+1:4n], q[n+1:2n]=>x[1:n]-I*x[n+1:2n]-J*x[2n+1:3n]-K*x[3n+1:4n])
#         pop[4*(i-1)+1] = real.(MultivariatePolynomials.coefficients(temp))'*MultivariatePolynomials.monomials(temp)
#         if symmetry == false
#             ip1 = [imag_part(coe)[1] for coe in MultivariatePolynomials.coefficients(temp)]
#             ip2 = [imag_part(coe)[2] for coe in MultivariatePolynomials.coefficients(temp)]
#             ip3 = [imag_part(coe)[3] for coe in MultivariatePolynomials.coefficients(temp)]
#             pop[4*(i-1)+2] = ip1'*MultivariatePolynomials.monomials(temp)
#             pop[4*(i-1)+3] = ip2'*MultivariatePolynomials.monomials(temp)
#             pop[4*(i-1)+4] = ip3'*MultivariatePolynomials.monomials(temp)
#         end
#     end
#     return pop,x
# end
