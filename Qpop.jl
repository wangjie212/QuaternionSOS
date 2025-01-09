using TSSOS
using DynamicPolynomials
using Quaternions
using MosekTools
using JuMP
using LinearAlgebra
include("function.jl")
function qs_tssos_first(pop::Vector{Polynomial{true, T}}, z, n, d; numeq=0, RemSig=false, nb=0, CS="MF", cliques=[], minimize=false, 
    TS=false, merge=false, md=3, solver="Mosek", reducebasis=false, QUIET=false, solve=true, tune=false, solution=false, ipart=true, 
    dualize=false, balanced=false, MomentOne=false, Gram=false, Mommat=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), 
    writetofile=false, normality=0, NormalSparse=false) where {T<:Number}
    supp,coe = Qpolys_info(pop, z, n)
    opt,cons,I,J,basis,hbasis,icons,jcons= qs_tssos_first(supp, coe, n, d, numeq=numeq, RemSig=RemSig, nb=nb, CS=CS, cliques=cliques, minimize=minimize,
    TS=TS, merge=merge, md=md, solver=solver, reducebasis=reducebasis, QUIET=QUIET, solve=solve, tune=tune, 
    solution=solution, ipart=ipart, dualize=dualize, balanced=balanced, MomentOne=MomentOne, Gram=Gram, Mommat=Mommat, 
    cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, writetofile=writetofile, normality=normality, NormalSparse=NormalSparse, cpop=pop, z=z)
    return opt,cons,I,J,basis,hbasis,icons,jcons
end

function qs_tssos_first(supp::Vector{Vector{Vector{Vector{UInt16}}}}, coe, n, d; numeq=0, RemSig=false, nb=0, CS="MF", cliques=[], 
    minimize=false, TS=false, merge=false, md=3, solver="Mosek", reducebasis=false, QUIET=false, solve=true, tune=false, solution=false, 
    ipart=true, dualize=false, balanced=false, MomentOne=false, Gram=false, Mommat=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), 
    writetofile=false, normality=0, NormalSparse=false, cpop=nothing, z=nothing)
    println("*********************************** TSSOS ***********************************")
    println("TSSOS is launching...")
    m = length(supp) - 1
    ind = [supp[1][i][1]<=supp[1][i][2] for i=1:length(supp[1])]
    supp[1] = supp[1][ind]
    coe[1] = coe[1][ind]
    dc = zeros(Int, m)
    for i = 1:m
        dc[i] = maximum([max(length(supp[i+1][j][1]), length(supp[i+1][j][2])) for j=1:length(supp[i+1])])
    end
    cliques = [UInt16[i for i=1:n]]
    cql=1
    I,J,ncc = assign_constraint(m, numeq, supp, cliques, cql)
    if d == "min"
        rlorder = isempty(I[1]) && isempty(J[1]) ? 1 : maximum(dc[[I[1]; J[1]]])
    else
        rlorder = d
    end
    basis = Vector{Vector{Vector{UInt16}}}(undef, m-numeq+1)
    hbasis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, numeq)
    basis[1] = get_basis(cliques[1], rlorder)
    for s = 1:length(I[1])
            basis[s+1] = get_basis(cliques[1], rlorder-dc[I[1][s]])
    end
    for s = 1:length(J[1])
            temp = get_basis(cliques[1], rlorder-dc[J[1][s]])
            hbasis[s] = vec([[item1, item2] for item1 in temp, item2 in temp])
            sort!(hbasis[s])
    end  
    tsupp = copy(supp[1])
    for i = 2:m+1, j = 1:length(supp[i])
        if supp[i][j][1] <= supp[i][j][2]
            push!(tsupp, supp[i][j])
        end
    end
    sort!(tsupp)
    unique!(tsupp)
    opt,ksupp,SDP_status,cons,icons,jcons= qsolvesdp(n, m, rlorder, supp, coe, basis, hbasis,I, J, ncc,
    numeq=numeq, QUIET=QUIET, TS=TS, solver=solver, solve=solve, tune=tune, solution=solution, ipart=ipart, MomentOne=MomentOne, balanced=balanced,
    Gram=Gram, Mommat=Mommat, nb=nb, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, dualize=dualize, writetofile=writetofile, normality=normality, NormalSparse=NormalSparse)
    return opt,cons,I,J,basis,hbasis,icons,jcons
end

function Qpolys_info(pop, z, n)
    coe = Vector{Vector{QuaternionF64}}(undef, length(pop))
    supp = Vector{Vector{Vector{Vector{UInt16}}}}(undef, length(pop))
    for k in eachindex(pop)
        mon = MultivariatePolynomials.monomials(pop[k])
        coe[k] = MultivariatePolynomials.coefficients(pop[k])
        lm = length(mon)
        supp[k] = [[[], []] for i=1:lm]
        for i = 1:lm
            ind = mon[i].z.> 0
            vars = mon[i].vars[ind]
            exp = mon[i].z[ind]
            for j in eachindex(vars)
                l = ncbfind(z, 2n, vars[j])
                if l <= n
                    append!(supp[k][i][1], l*ones(UInt16, exp[j]))
                else
                    append!(supp[k][i][2], (l-n)*ones(UInt16, exp[j]))
                end
            end
        end
    end
    return supp,coe
end

function qsolvesdp(n, m, rlorder, supp::Vector{Vector{Vector{Vector{UInt16}}}}, coe, basis, hbasis,I, J, ncc; 
    numeq=0, nb=0, QUIET=false, TS=false, solver="Mosek", tune=false, solve=true, dualize=false, solution=false, Gram=false, MomentOne=false, ipart=true, Mommat=false, 
    cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false,balanced=false, normality=0, NormalSparse=false)
    tsupp = Vector{Vector{UInt16}}[]
    println(basis[1])
    for k = 1:length(basis[1]), r = k:length(basis[1])
        @inbounds bi = [basis[1][k], basis[1][r]]
        if bi[1] <= bi[2]
            push!(tsupp, bi)
        else
            push!(tsupp, bi[2:-1:1])
        end
    end
    ksupp = copy(tsupp)
    sort!(tsupp)
    unique!(tsupp)
    # println(tsupp)
    #initial set
    objv = moment = GramMat = SDP_status= nothing
    if solve == true
        ltsupp = length(tsupp)
        println(ltsupp)
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
            if tune == true
                set_optimizer_attributes(model,
                "MSK_DPAR_INTPNT_CO_TOL_MU_RED" => 1e-7,
                "MSK_DPAR_INTPNT_CO_TOL_INFEAS" => 1e-7,
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1e-7,
                "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => 1e-7,
                "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-7,
                "MSK_DPAR_INTPNT_CO_TOL_NEAR_REL" => 1e6,
                "MSK_IPAR_BI_IGNORE_NUM_ERROR" => 1,
                "MSK_DPAR_BASIS_TOL_X" => 1e-3,
                "MSK_DPAR_BASIS_TOL_S" => 1e-3,
                "MSK_DPAR_BASIS_REL_TOL_S" => 1e-5)
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
        pos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m-numeq+1)
        pos[1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, 1)
            @inbounds bs = length(basis[1])
            if bs == 1
                pos[1][1] = @variable(model, lower_bound=0,base_name="Q")
                @inbounds bi = [basis[1][1], basis[1][1]]
                Locb = bfind(tsupp, ltsupp, bi)
                @inbounds add_to_expression!(rcons[Locb], pos[1][1])
            else
                pos[1][1] = @variable(model, [1:4bs, 1:4bs], PSD,base_name="Q")
                # pos[1][1] = @variable(model, [1:bs, 1:bs], PSD,base_name="Q")
                for t = 1:bs, r = t:bs
                    @inbounds bi = [basis[1][t], basis[1][r]]
                    if bi[1] <= bi[2]
                        Locb = bfind(tsupp, ltsupp, bi)
                        if bi[1] != bi[2]
                            @inbounds add_to_expression!(jcons[Locb], pos[1][1][t+2*bs,r]-pos[1][1][t,r+2*bs]-pos[1][1][t+3*bs,r+bs]+pos[1][1][t+bs,r+3*bs])
                            @inbounds add_to_expression!(icons[Locb], pos[1][1][t+bs,r]-pos[1][1][t,r+bs]+pos[1][1][t+3*bs,r+2*bs]-pos[1][1][t+2*bs,r+3*bs]) 
                            @inbounds add_to_expression!(kcons[Locb], pos[1][1][t+3*bs,r]-pos[1][1][t,r+3*bs]+pos[1][1][t+2*bs,r+bs]-pos[1][1][t+bs,r+2*bs])
                        end
                    else
                        Locb = bfind(tsupp, ltsupp, bi[2:-1:1])
                        @inbounds add_to_expression!(icons[Locb], -1, pos[1][1][t+bs,r]-pos[1][1][t,r+bs]+pos[1][1][t+3*bs,r+2*bs]-pos[1][1][t+2*bs,r+3*bs])
                        @inbounds add_to_expression!(jcons[Locb], -1, pos[1][1][t+2*bs,r]-pos[1][1][t,r+2*bs]-pos[1][1][t+3*bs,r+bs]+pos[1][1][t+bs,r+3*bs])
                        @inbounds add_to_expression!(kcons[Locb], -1, pos[1][1][t+3*bs,r]-pos[1][1][t,r+3*bs]+pos[1][1][t+2*bs,r+bs]-pos[1][1][t+bs,r+2*bs])
                    end
                    @inbounds add_to_expression!(rcons[Locb], pos[1][1][t,r]+pos[1][1][t+bs,r+bs]+pos[1][1][t+2*bs,r+2*bs]+pos[1][1][t+3*bs,r+3*bs])
                    # @inbounds add_to_expression!(rcons[Locb], pos[1][1][t,r])
                
                end
            end
        for (j, w) in enumerate(I[1])
            pos[j+1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, 1)
                bs = length(basis[j+1])
                # println(basis[j+1])
                # println(bs)
                if bs == 1
                    pos[j+1][1] = @variable(model, lower_bound=0,base_name="X")
                    for s = 1:length(supp[w+1])
                        @inbounds bi = [sadd(basis[j+1][1], supp[w+1][s][1]), sadd(basis[j+1][1], supp[w+1][s][2])]
                        if bi[1] <= bi[2]
                            Locb = bfind(tsupp, ltsupp, bi)
                            @inbounds add_to_expression!(rcons[Locb], real(coe[w+1][s]), pos[j+1][1])
                            @inbounds add_to_expression!(icons[Locb], imag_part(coe[w+1][s])[1], pos[j+1][1])
                            @inbounds add_to_expression!(jcons[Locb], imag_part(coe[w+1][s])[2], pos[j+1][1])
                            @inbounds add_to_expression!(kcons[Locb], imag_part(coe[w+1][s])[3], pos[j+1][1])
                        end
                    end 
                else
                        pos[j+1][1] = @variable(model, [1:4bs, 1:4bs], PSD,base_name="X")
                        # pos[j+1][1] = @variable(model, [1:bs, 1:bs], PSD,base_name="X")
                    for t = 1:bs, r = 1:bs
                        for s = 1:length(supp[w+1])
                            @inbounds bi = [sadd(basis[j+1][t], supp[w+1][s][1]), sadd(basis[j+1][r], supp[w+1][s][2])]
                            if bi[1] <= bi[2]
                                Locb = bfind(tsupp, ltsupp, bi)
                                # @inbounds add_to_expression!(rcons[Locb], real(coe[w+1][s]), pos[j+1][1][t,r])
                                    @inbounds add_to_expression!(rcons[Locb], real(coe[w+1][s]), pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
                                    @inbounds add_to_expression!(rcons[Locb], imag_part(coe[w+1][s])[1], pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
                                    @inbounds add_to_expression!(rcons[Locb], imag_part(coe[w+1][s])[2], pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
                                    @inbounds add_to_expression!(rcons[Locb], imag_part(coe[w+1][s])[3], pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
                                    if bi[1] != bi[2]
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[w+1][s])[1], pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], real(coe[w+1][s]), pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], -imag_part(coe[w+1][s])[3], pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(icons[Locb], imag_part(coe[w+1][s])[2], pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
                                        
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[w+1][s])[2], pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], imag_part(coe[w+1][s])[3], pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], real(coe[w+1][s]), pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(jcons[Locb], -imag_part(coe[w+1][s])[1], pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
                                        
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[w+1][s])[3], pos[j+1][1][t,r]+pos[j+1][1][t+bs,r+bs]+pos[j+1][1][t+2*bs,r+2*bs]+pos[j+1][1][t+3*bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[w+1][s])[2], pos[j+1][1][t+bs,r]-pos[j+1][1][t,r+bs]+pos[j+1][1][t+3*bs,r+2*bs]-pos[j+1][1][t+2*bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], imag_part(coe[w+1][s])[1], pos[j+1][1][t+2*bs,r]-pos[j+1][1][t,r+2*bs]-pos[j+1][1][t+3*bs,r+bs]+pos[j+1][1][t+bs,r+3*bs])
                                        @inbounds add_to_expression!(kcons[Locb], real(coe[w+1][s]), pos[j+1][1][t+3*bs,r]-pos[j+1][1][t,r+3*bs]+pos[j+1][1][t+2*bs,r+bs]-pos[j+1][1][t+bs,r+2*bs])
                                    end
                            end
                        end
                    end
                end
        end
        if numeq > 0
            free = Vector{Vector{Vector{VariableRef}}}(undef, 1)
                if !isempty(J[1])
                    free[1] = Vector{Vector{VariableRef}}(undef, length(J[1]))
                    for (j, w) in enumerate(J[1])
                        mons = hbasis[j][Vector(1:length(hbasis[j]))]
                        temp = mons[[item[1] <= item[2] for item in mons]]
                        lb = length(temp)
                        free[1][j] = @variable(model, [1:4*lb])
                        for k in Vector(1:length(hbasis[j])), s = 1:length(supp[w+1])
                            @inbounds bi = [sadd(hbasis[j][k][1], supp[w+1][s][1]), sadd(hbasis[j][k][2], supp[w+1][s][2])]
                            if bi[1] <= bi[2]
                                Locb = bfind(tsupp, ltsupp, bi)
                                if hbasis[j][k][1] <= hbasis[j][k][2]
                                    loc = bfind(temp, lb, hbasis[j][k])
                                    tag = 1
                                    if hbasis[j][k][1] == hbasis[j][k][2]
                                        tag = 0
                                    end
                                else
                                    loc = bfind(temp, lb, hbasis[j][k][2:-1:1])
                                    tag = -1
                                end
                                    @inbounds add_to_expression!(rcons[Locb], real(coe[w+1][s])*free[1][j][loc]-tag*imag_part(coe[w+1][s])[1]*free[1][j][loc+lb]-tag*imag_part(coe[w+1][s])[2]*free[1][j][loc+2*lb]-tag*imag_part(coe[w+1][s])[3]*free[1][j][loc+3*lb])
                                    if bi[1] != bi[2]
                                        @inbounds add_to_expression!(icons[Locb], tag*real(coe[w+1][s])*free[1][j][loc+lb]+imag_part(coe[w+1][s])[1]*free[1][j][loc]+imag_part(coe[w+1][s])[2]*free[1][j][loc+3*lb]-imag_part(coe[w+1][s])[3]*free[1][j][loc+2*lb])
                                        @inbounds add_to_expression!(jcons[Locb], tag*real(coe[w+1][s])*free[1][j][loc+2*lb]-imag_part(coe[w+1][s])[1]*free[1][j][loc+3*lb]+imag_part(coe[w+1][s])[2]*free[1][j][loc]+imag_part(coe[w+1][s])[3]*free[1][j][loc+lb])
                                        @inbounds add_to_expression!(kcons[Locb], tag*real(coe[w+1][s])*free[1][j][loc+3*lb]+imag_part(coe[w+1][s])[1]*free[1][j][loc+2*lb]-imag_part(coe[w+1][s])[2]*free[1][j][loc+lb]-imag_part(coe[w+1][s])[3]*free[1][j][loc])
                                    
                                    end
                            end
                        end
                    end
                end
        end
        rbc = zeros(ltsupp)
        ncons = ltsupp
        itsupp = nothing
        ind = [item[1] != item[2] for item in tsupp]
        itsupp = tsupp[ind]
        icons = icons[ind]
        jcons = jcons[ind]
        kcons = kcons[ind]
        ibc = zeros(length(itsupp))
        jbc = zeros(length(itsupp))
        kbc = zeros(length(itsupp))
        ncons += length(itsupp)
        if QUIET == false
            println("There are $ncons affine constraints.")
        end
        for i = 1:length(supp[1])
            Locb = bfind(tsupp, ltsupp, supp[1][i])
            if Locb === nothing
                @error "The monomial basis is not enough!"
                return nothing,ksupp,nothing,nothing,nothing
            else
            rbc[Locb] = real(coe[1][i])
                if supp[1][i][1] != supp[1][i][2]
                Locb = bfind(itsupp, length(itsupp), supp[1][i])
                ibc[Locb] = imag_part(coe[1][i])[1]
                jbc[Locb] = imag_part(coe[1][i])[2]
                kbc[Locb] = imag_part(coe[1][i])[3]
                end
            end
        end
        @variable(model, lower)
        rcons[1] += lower
        @constraint(model, rcon[i=1:ltsupp], rcons[i]==rbc[i])
        @constraint(model, icon[i=1:length(itsupp)], icons[i]==ibc[i])
        @constraint(model, jcon[i=1:length(itsupp)], jcons[i]==jbc[i])
        @constraint(model, kcon[i=1:length(itsupp)], kcons[i]==kbc[i])
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
    return objv,ksupp,SDP_status,cons,icons,jcons
end

function quaternion_to_real(qpop, q)
    n = Int(length(q)/2)
    @polyvar x[1:4n]
    pop = Vector{Polynomial{true,Float64}}(undef, length(qpop))
    I=quat(0,1,0,0)
    J=quat(0,0,1,0)
    K=quat(0,0,0,1)
    for (i,qp) in enumerate(qpop)
        temp = qp(q[1:n]=>x[1:n]+I*x[n+1:2n]+J*x[2n+1:3n]+K*x[3n+1:4n], q[n+1:2n]=>x[1:n]-I*x[n+1:2n]-J*x[2n+1:3n]-K*x[3n+1:4n])
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
