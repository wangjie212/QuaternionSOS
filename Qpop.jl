function qs_tssos_first(pop::Vector{Polynomial{false,T}}, z, n, d; numeq=0, RemSig=false, nb=0,
    TS=false, merge=false, md=3, solver="Mosek", reducebasis=false, QUIET=false, solve=true, solution=false, ipart=true, 
    dualize=false, balanced=false, cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), 
    writetofile=false, normality=0, NormalSparse=false, conjubasis=true) where {T<:Number}
    supp,coe = Qpolys_info(pop, z, n)
    println("*********************************** TSSOS ***********************************")
    println("QSSOS is launching...")
    if nb > 0
        supp[1],coe[1] = qresort(supp[1],coe[1];nb=nb)
    end
    m = length(supp) - 1
    dc = zeros(Int, m)
    for i = 1:m
        dc[i] = maximum([length(supp[i+1][j][1])+length(supp[i+1][j][2])+length(supp[i+1][j][3]) for j=1:length(supp[i+1])])
    end
    basis = Vector{Vector{Vector{Vector{UInt16}}}}(undef, m+1)
    basis[1] = get_qncbasis(n, d; conjubasis=conjubasis)
    if QUIET == false
        lb = length(basis[1])
        println("Maximal PSD block size: $lb.")
    end
    for s = 1:m
        basis[s+1] = get_qncbasis(n, d-Int(ceil(dc[s]/2));conjubasis=conjubasis)
    end
    if nb > 0
        for s = 1:m+1
            basis[s] = qreduce_unitnorm.(basis[s], nb=nb)
            unique!(basis[s])
        end
    end
    opt,ksupp,SDP_status = qsolvesdp(n, m, d, supp, coe, basis, numeq=numeq, QUIET=QUIET, TS=TS, solver=solver, solve=solve, solution=solution, ipart=ipart, balanced=balanced,
    nb=nb, cosmo_setting=cosmo_setting, mosek_setting=mosek_setting, dualize=dualize, writetofile=writetofile, normality=normality, NormalSparse=NormalSparse,conjubasis=conjubasis)
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

function qsolvesdp(n, m, d, supp::Vector{Vector{Vector{Vector{UInt16}}}}, coe, basis; 
    numeq=0, nb=0, QUIET=false, TS=false, solver="Mosek", solve=true, dualize=false, solution=false, ipart=true, 
    cosmo_setting=cosmo_para(), mosek_setting=mosek_para(), writetofile=false, balanced=false, normality=0, NormalSparse=false,conjubasis=false)
    tsupp = Vector{Vector{UInt16}}[]
    for k = 1:length(basis[1]), r = 1:length(basis[1])
        a=deepcopy(basis[1][k])
        b=deepcopy(basis[1][r])
        @inbounds bi = qtermadd(a,b,n)
        if nb > 0
            bi = qreduce_unitnorm(bi, nb=nb)
        end
        push!(tsupp,bi)
    end
    ksupp = deepcopy(tsupp)
    if normality > 0
        wbasis = [[UInt16[],UInt16[],UInt16[]]]
        for t=1:n
            wbasis= bget_qncbasis(n, normality,t)
            if nb > 0
                wbasis= qreduce_unitnorm.(wbasis, nb=nb)
                unique!(wbasis)
            end
            ws = length(wbasis)
            for k = 1:ws, r = 1:ws
                a=deepcopy(wbasis[k])
                b=deepcopy(wbasis[r])
                @inbounds bi = qtermadd(a,b,n)
                if nb > 0
                    bi = qreduce_unitnorm(bi, nb=nb)
                end
                push!(tsupp,bi)
            end
        end
    end
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
        if normality > 0
            for t=1:n
                wbasis= bget_qncbasis(n, normality,t)
                if nb > 0
                    wbasis= qreduce_unitnorm.(wbasis, nb=nb)
                    unique!(wbasis)
                end
                ws = length(wbasis)
                if NormalSparse == false
                    hnom[t] = @variable(model, [1:ws, 1:ws], PSD)
                    for j = 1:ws, k = 1:ws
                        a=deepcopy(wbasis[j])
                        b=deepcopy(wbasis[k])
                        bi =qtermadd(a,b,n)
                        if nb > 0
                            bi = qreduce_unitnorm(bi, nb=nb)
                        end
                        Locb = bfind(tsupp, ltsupp, bi)
                        @inbounds add_to_expression!(rcons[Locb], hnom[t][j,k])
                    end
                end
            end
        end
        # pos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m-numeq+1)
        pos = Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m+1)
        pos[1] = Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, 1)
        @inbounds bs = length(basis[1])
        if bs == 1
            pos[1][1] = @variable(model, lower_bound=0,base_name="Q")
            a=deepcopy(basis[1][1])
            b=deepcopy(basis[1][1])
            @inbounds bi = qtermadd(a, b,n)
            if nb > 0
                bi = qreduce_unitnorm(bi;nb=nb)
            end
            Locb = bfind(tsupp, ltsupp, bi)
            @inbounds add_to_expression!(rcons[Locb], pos[1][1])
        else
            if ipart == true
                pos[1][1] = @variable(model, [1:4bs, 1:4bs], PSD,base_name="Q")
            else
                pos[1][1] = @variable(model, [1:bs, 1:bs], PSD,base_name="Q")
            end
            for t = 1:bs, r = 1:bs
                a=deepcopy(basis[1][t])
                b=deepcopy(basis[1][r])
                @inbounds bi = qtermadd(a,b,n)
                if nb > 0
                    bi = qreduce_unitnorm(bi;nb=nb)
                end
                Locb = bfind(tsupp, ltsupp, bi)
                if ipart == true
                    @inbounds add_to_expression!(rcons[Locb], pos[1][1][t,r]+pos[1][1][t+bs,r+bs]+pos[1][1][t+2*bs,r+2*bs]+pos[1][1][t+3*bs,r+3*bs])
                    @inbounds add_to_expression!(icons[Locb], pos[1][1][t+bs,r]-pos[1][1][t,r+bs]+pos[1][1][t+3*bs,r+2*bs]-pos[1][1][t+2*bs,r+3*bs]) 
                    @inbounds add_to_expression!(jcons[Locb], pos[1][1][t+2*bs,r]-pos[1][1][t,r+2*bs]-pos[1][1][t+3*bs,r+bs]+pos[1][1][t+bs,r+3*bs])
                    @inbounds add_to_expression!(kcons[Locb], pos[1][1][t+3*bs,r]-pos[1][1][t,r+3*bs]+pos[1][1][t+2*bs,r+bs]-pos[1][1][t+bs,r+2*bs])
                else
                    @inbounds add_to_expression!(rcons[Locb], pos[1][1][t,r])
                end            
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
                    Locb = bfind(tsupp, ltsupp, bi)
                    @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[k+1][1])
                    if ipart == true
                        @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[1], pos[k+1][1])
                        @inbounds add_to_expression!(jcons[Locb], imag_part(coe[k+1][s])[2], pos[k+1][1])
                        @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], pos[k+1][1])
                    end
                end 
            else
                if ipart == true
                    if k <= m-numeq
                        pos[k+1][1] = @variable(model, [1:4bs, 1:4bs], PSD)
                    else
                        pos[k+1][1] = @variable(model, [1:4bs, 1:4bs], Symmetric)
                    end
                else
                    if k <= m-numeq
                        pos[k+1][1] = @variable(model, [1:bs, 1:bs], PSD)
                    else
                        pos[k+1][1] = @variable(model, [1:bs, 1:bs], Symmetric)
                    end
                end
                for t = 1:bs, r = 1:bs
                    for s = 1:length(supp[k+1])
                        a=deepcopy(basis[k+1][t])
                        b=deepcopy(supp[k+1][s])
                        c=deepcopy(basis[k+1][r])
                        @inbounds bi = qtermadd3(a,b,c,n)
                        if nb > 0
                            bi = qreduce_unitnorm(bi;nb=nb)
                        end
                        Locb = bfind(tsupp, ltsupp, bi)
                        if ipart == true
                            @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[k+1][1][t,r]+pos[k+1][1][t+bs,r+bs]+pos[k+1][1][t+2*bs,r+2*bs]+pos[k+1][1][t+3*bs,r+3*bs])
                            @inbounds add_to_expression!(rcons[Locb], imag_part(coe[k+1][s])[1], pos[k+1][1][t+bs,r]-pos[k+1][1][t,r+bs]+pos[k+1][1][t+3*bs,r+2*bs]-pos[k+1][1][t+2*bs,r+3*bs])
                            @inbounds add_to_expression!(rcons[Locb], imag_part(coe[k+1][s])[2], pos[k+1][1][t+2*bs,r]-pos[k+1][1][t,r+2*bs]-pos[k+1][1][t+3*bs,r+bs]+pos[k+1][1][t+bs,r+3*bs])
                            @inbounds add_to_expression!(rcons[Locb], imag_part(coe[k+1][s])[3], pos[k+1][1][t+3*bs,r]-pos[k+1][1][t,r+3*bs]+pos[k+1][1][t+2*bs,r+bs]-pos[k+1][1][t+bs,r+2*bs])

                            @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[1], pos[k+1][1][t,r]+pos[k+1][1][t+bs,r+bs]+pos[k+1][1][t+2*bs,r+2*bs]+pos[k+1][1][t+3*bs,r+3*bs])
                            @inbounds add_to_expression!(icons[Locb], real(coe[k+1][s]), pos[k+1][1][t+bs,r]-pos[k+1][1][t,r+bs]+pos[k+1][1][t+3*bs,r+2*bs]-pos[k+1][1][t+2*bs,r+3*bs])
                            @inbounds add_to_expression!(icons[Locb], -imag_part(coe[k+1][s])[3], pos[k+1][1][t+2*bs,r]-pos[k+1][1][t,r+2*bs]-pos[k+1][1][t+3*bs,r+bs]+pos[k+1][1][t+bs,r+3*bs])
                            @inbounds add_to_expression!(icons[Locb], imag_part(coe[k+1][s])[2], pos[k+1][1][t+3*bs,r]-pos[k+1][1][t,r+3*bs]+pos[k+1][1][t+2*bs,r+bs]-pos[k+1][1][t+bs,r+2*bs])
                                        
                            @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[2], pos[k+1][1][t,r]+pos[k+1][1][t+bs,r+bs]+pos[k+1][1][t+2*bs,r+2*bs]+pos[k+1][1][t+3*bs,r+3*bs])
                            @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[3], pos[k+1][1][t+bs,r]-pos[k+1][1][t,r+bs]+pos[k+1][1][t+3*bs,r+2*bs]-pos[k+1][1][t+2*bs,r+3*bs])
                            @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), pos[k+1][1][t+2*bs,r]-pos[k+1][1][t,r+2*bs]-pos[k+1][1][t+3*bs,r+bs]+pos[k+1][1][t+bs,r+3*bs])
                            @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[1], pos[k+1][1][t+3*bs,r]-pos[k+1][1][t,r+3*bs]+pos[k+1][1][t+2*bs,r+bs]-pos[k+1][1][t+bs,r+2*bs])
                                            
                            @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[3], pos[k+1][1][t,r]+pos[k+1][1][t+bs,r+bs]+pos[k+1][1][t+2*bs,r+2*bs]+pos[k+1][1][t+3*bs,r+3*bs])
                            @inbounds add_to_expression!(kcons[Locb], -imag_part(coe[k+1][s])[2], pos[k+1][1][t+bs,r]-pos[k+1][1][t,r+bs]+pos[k+1][1][t+3*bs,r+2*bs]-pos[k+1][1][t+2*bs,r+3*bs])
                            @inbounds add_to_expression!(kcons[Locb], imag_part(coe[k+1][s])[1], pos[k+1][1][t+2*bs,r]-pos[k+1][1][t,r+2*bs]-pos[k+1][1][t+3*bs,r+bs]+pos[k+1][1][t+bs,r+3*bs])
                            @inbounds add_to_expression!(kcons[Locb], real(coe[k+1][s]), pos[k+1][1][t+3*bs,r]-pos[k+1][1][t,r+3*bs]+pos[k+1][1][t+2*bs,r+bs]-pos[k+1][1][t+bs,r+2*bs])
                        else
                            @inbounds add_to_expression!(rcons[Locb], real(coe[k+1][s]), pos[k+1][1][t,r])
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
