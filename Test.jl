using Random
using TSSOS
using DynamicPolynomials
using Graphs
using Quaternions
using MosekTools
using JuMP
using LinearAlgebra
using MultivariatePolynomials
import MultivariatePolynomials as MP
import DynamicPolynomials as DP
import CliqueTrees

# include("D:/Programs/QuaternionSOS/ncutils.jl")
# include("D:/Programs/QuaternionSOS/Qpop.jl")
include("C:/Users/qingchefff/Documents/julia/Quat/QuaternionSOS/ncutils.jl")
include("C:/Users/qingchefff/Documents/julia/Quat/QuaternionSOS/Qpop.jl")



#### Test 1 

###set: f = [q]_1^*Q[q]_1, Q:real, n=20,40,60；

rng = Xoshiro(1)
n = 4
@ncpolyvar q[1:2n]
f,Fsupp,Fcoe = randomsymfunc(q, n, 1, rng, conjugates=false, coelimit=false)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## ball

#QSOS: basis[q]_1
@time opt = qs_tssos_first([f, g], q, n, 1, fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false,ipart=false,conjubasis=false, solution = false, QUIET=true)

#RSOS: d=1
pop,x = quaternion_to_real([f, g], q)
@time tssos_first(pop, x, 1, TS=false,solution=true, QUIET=true)

## unit norm
@time opt = qs_tssos_first([f], q, n, 1, nb=n,fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=false, conjubasis=true, solution = false , QUIET=true)
pop,x = quaternion_to_real([f;gn], q)
@time opt,sol,data = tssos_first(pop, x, 1, numeq=n, TS=false, solution=true, QUIET=true)
@time opt,sol,data = tssos_first(pop, x, 2, numeq=n, TS=false, solution=true, QUIET=true)
# upper bound
ub = local_solution(data.n,data.m,data.supp,data.coe;nb=data.nb,numeq=data.numeq,startpoint=rand(data.n),QUIET=true)[1]
println(ub)
# @time opt,sol,data = tssos_first(pop, x, 2, numeq=n, TS=false, solution=true, QUIET=true)

### set: f = [q]_1^*Q[q]_1, Q:quaternion, n=10,20,30,40；
rng = Xoshiro(3)
n = 10
@ncpolyvar q[1:2n]
f,fr,Fsupp,Fcoe = qrandomsymfunc(q, n, 1, rng; conjugates=false)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## ball
# QSOS: basis[q]_1
@time opt = qs_tssos_first([f, g], q, n, 1, fsupp=Fsupp, fcoe=Fcoe,TS=false, ipart=true,conjubasis=false, solution = false, QUIET=true)

# RSOS:d =1
pop,x = quaternion_to_real([fr, g], q)
@time opt,sol,data = tssos_first(pop, x, 1, TS=false, solution=false, QUIET=true)
@time opt,sol,data = tssos_first(pop, x, 2, TS=false, solution=false, QUIET=true)
ub = local_solution(data.n,data.m,data.supp,data.coe;nb=data.nb,numeq=data.numeq,startpoint=rand(data.n),QUIET=true)[1]
println(ub)

## unit norm
# QSOS: basis[q]_1
@time opt= qs_tssos_first([f], q, n, 1, nb=n, fsupp=Fsupp, fcoe=Fcoe,TS=false, ipart=true,conjubasis=false, solution = false, QUIET=true)
pop,x = quaternion_to_real([fr;gn], q)
@time opt,sol,data =tssos_first(pop, x, 1, numeq=n, TS=false, solution=false, QUIET=true)
ub = local_solution(data.n,data.m,data.supp,data.coe;nb=data.nb,numeq=data.numeq,startpoint=rand(data.n),QUIET=true)[1]
println(ub)




#### Test 2

###set: f = [q,\bar(q)]_1^*Q[q,\bar(q)]_1, Q:real, n=20,40,60

rng = Xoshiro(2)
n = 3
@ncpolyvar q[1:2n]
f,Fsupp,Fcoe = randomsymfunc(q, n, 1, rng, conjugates=true, coelimit=false)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]
comm_constraints = typeof(f)[]
for i in 1:n
    for j in 1:n
        if i != j
            expr = q[i]*q[j] + q[i+n]*q[j] - q[j]*q[i] - q[j]*q[i+n]
            push!(comm_constraints, expr)
        end
    end
end
## ball

# QSOS: basis[q,\bar(q)]_1
@time qs_tssos_first([f, g], q, n, 1, fsupp=Fsupp, fcoe=Fcoe, TS=false, ipart=false, conjubasis=true, QUIET=true)

# RSOS: d=1
pop,x = quaternion_to_real([f, g], q)
@time opt,sol,data = tssos_first(pop, x, 1, TS=false, solution=true, QUIET=true)
@time opt,sol,data = tssos_first(pop, x, 2, TS=false, solution=true, QUIET=true)

## unit norm
@time opt = qs_tssos_first([f], q, n, 1, nb=n,fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=false, conjubasis=true, QUIET=true)
@time opt = qs_tssos_first([f;gn;comm_constraints], q, n, 1, numeq=n+length(comm_constraints),rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=false, conjubasis=true, QUIET=true)
# @time opt = qs_tssos_first([f;gn;comm_constraints], q, n, 1, numeq=n+length(comm_constraints),rncnumeq= length(comm_constraints), normality=1,fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=false, conjubasis=false, QUIET=true)
pop,x = quaternion_to_real([f;gn], q)
@time opt,sol,data = tssos_first(pop, x, 1, numeq=n, TS=false, solution=true, QUIET=true)
@time opt,sol,data = tssos_first(pop, x, 2, numeq=n, TS=false, solution=true, QUIET=true)
ub = local_solution(data.n,data.m,data.supp,data.coe;nb=data.nb,numeq=data.numeq,startpoint=rand(data.n),QUIET=true)[1]
println(ub)


### set Q:quaternion
rng = Xoshiro(3)
n = 4
@ncpolyvar q[1:2n]
f,fr,Fsupp,Fcoe = qrandomsymfunc(q, n, 1, rng; conjugates=true)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]
comm_constraints = typeof(f)[]
for i in 1:n
    for j in 1:n
        if i != j
            expr = q[i]*q[j] + q[i+n]*q[j] - q[j]*q[i] - q[j]*q[i+n]
            push!(comm_constraints, expr)
        end
    end
end
println(f)
println(length(Fsupp))
## ball

# QSOS: basis[q,\bar(q)]_1
@time qs_tssos_first([f, g], q, n, 1, fsupp=Fsupp, fcoe=Fcoe, TS=false, ipart=true, conjubasis=true, QUIET=true)
@time qs_tssos_first([f, g], q, n, 1, TS=false, ipart=true, conjubasis=true, QUIET=true)

@time qs_tssos_first([[f, g];comm_constraints], q, n, 1, numeq = length(comm_constraints), rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=true, conjubasis=true, QUIET=true) 
# RSOS: d=1
pop,x = quaternion_to_real([fr, g], q)
@time opt,sol,data = tssos_first(pop, x, 1, TS=false, solution=true, QUIET=true)
@time opt,sol,data = tssos_first(pop, x, 2, TS=false, solution=true, QUIET=true)

## unit norm
@time opt = qs_tssos_first([f], q, n, 1, nb=n,fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=true, conjubasis=true, QUIET=true)
@time opt = qs_tssos_first([f;gn;comm_constraints], q, n, 1, numeq=n+length(comm_constraints),rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=true, conjubasis=true, QUIET=true)
@time opt = qs_tssos_first([f;gn;comm_constraints], q, n, 1, numeq=n+length(comm_constraints),rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=true, conjubasis=true, QUIET=true)
@time opt = qs_tssos_first([f;comm_constraints], q, n, 1, numeq=length(comm_constraints),rncnumeq= length(comm_constraints),nb=n, fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=true, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([fr;gn], q)
@time opt,sol,data = tssos_first(pop, x, 1, numeq=n, TS=false, solution=true, QUIET=true)
@time opt,sol,data = tssos_first(pop, x, 2, numeq=n, TS=false, solution=true, QUIET=true)
ub = local_solution(data.n,data.m,data.supp,data.coe;nb=data.nb,numeq=data.numeq,startpoint=rand(data.n),QUIET=true)[1]
println(ub)


#### Test 3

###set: f = [q]_2^*Q[q]_2, Q:real, n=1,2,3,4,5,6；

rng = Xoshiro(11)
n = 4
@ncpolyvar q[1:2n]
f,Fsupp,Fcoe = randomsymfunc(q, n, 2, rng, conjugates=false, coelimit=true)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]
comm_constraints = typeof(f)[]
for i in 1:n
    for j in 1:n
        if i != j
            expr = q[i]*q[j] + q[i+n]*q[j] - q[j]*q[i] - q[j]*q[i+n]
            push!(comm_constraints, expr)
        end
    end
end
## ball
# QSOS:basis[q]_2
@time qs_tssos_first([f, g], q, n, 2, fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=false, conjubasis=false, QUIET=true) 
@time qs_tssos_first([[f, g];comm_constraints], q, n, 2, numeq = length(comm_constraints), rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=false, conjubasis=false, QUIET=true) 
# strengthened QSOS :basis [q]_2+\bar(q_i)[q]_2
@time qs_tssos_first([f, g], q, n, 2, fsupp=Fsupp, fcoe=Fcoe, CS=false,TS=false, ipart=false, normality = 1,conjubasis=false, QUIET=true)
@time qs_tssos_first([[f, g];comm_constraints], q, n, 2, numeq = length(comm_constraints), rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=false, normality = 1, conjubasis=false, QUIET=true) 

#QSOS: full basis [q,\bar(q)]_2
@time qs_tssos_first([f, g], q, n, 2, fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=false, conjubasis=true, QUIET=true) 
@time qs_tssos_first([[f, g];comm_constraints], q, n, 2, numeq = length(comm_constraints), rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=false, conjubasis=true, QUIET=true) 

#RSOS
pop,x = quaternion_to_real([f, g], q)
@time tssos_first(pop, x, 2, TS=false, solution=true, QUIET=true)

## unit norm
@time qs_tssos_first([f], q, n, 2, nb=n, fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=false, conjubasis=false, QUIET=true)
@time qs_tssos_first([f;gn;comm_constraints], q, n, 2, numeq=n+length(comm_constraints), rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=false, conjubasis=false, QUIET=true)
@time qs_tssos_first([f], q, n, 2, nb=n, fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=false, normality = 1,conjubasis=false, QUIET=true)
@time qs_tssos_first([f;gn;comm_constraints], q, n, 2, numeq=n+length(comm_constraints), rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=false, normality = 1,conjubasis=false, QUIET=true)
@time qs_tssos_first([f], q, n, 2, nb=n, fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=false, conjubasis=true, QUIET=true)
@time qs_tssos_first([f;gn;comm_constraints], q, n, 2, numeq=n+length(comm_constraints), rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=false, conjubasis=true, QUIET=true)
@time qs_tssos_first([f;comm_constraints], q, n, 2, numeq=length(comm_constraints), rncnumeq= length(comm_constraints),nb=n, fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=false, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([f; gn], q)
@time tssos_first(pop, x, 2, numeq=n, TS=false, solution=true, QUIET=true)


### set Q:quaternion
rng = Xoshiro(3)
n = 3
@ncpolyvar q[1:2n]
f,fr,Fsupp,Fcoe = qrandomsymfunc(q, n, 2, rng; conjugates=false)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]
comm_constraints = typeof(f)[]
for i in 1:n
    for j in 1:n
        if i != j
            expr = q[i]*q[j] + q[i+n]*q[j] - q[j]*q[i] - q[j]*q[i+n]
            push!(comm_constraints, expr)
        end
    end
end

## ball
# QSOS:basis[q]_2
@time qs_tssos_first([f, g], q, n, 2, fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=true, conjubasis=false, QUIET=true) 
@time qs_tssos_first([[f, g];comm_constraints], q, n, 3, numeq = length(comm_constraints), rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=true, conjubasis=false, QUIET=true) 
# strengthened QSOS :basis [q]_2+\bar(q_i)[q]_2
@time qs_tssos_first([f, g], q, n, 3, fsupp=Fsupp, fcoe=Fcoe, CS=false,TS=false, ipart=false, normality = 1,conjubasis=false, QUIET=true)
@time qs_tssos_first([[f, g];comm_constraints], q, n, 3, numeq = length(comm_constraints), rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=true, normality = 1, conjubasis=false, QUIET=true) 

#QSOS: full basis [q,\bar(q)]_2
@time qs_tssos_first([f, g], q, n, 2, fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=false, conjubasis=true, QUIET=true) 
@time qs_tssos_first([[f, g];comm_constraints], q, n, 2, numeq = length(comm_constraints), rncnumeq= length(comm_constraints),fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=true, conjubasis=true, QUIET=true) 

pop,x = quaternion_to_real([fr, g], q)
@time tssos_first(pop, x, 2, TS=false, solution=true, QUIET=true)

## unit norm
# @time qs_tssos_first([f], q, n, 2, nb=n, fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=false, conjubasis=false, QUIET=true)
@time qs_tssos_first([f;comm_constraints], q, n, 2, numeq=length(comm_constraints), rncnumeq= length(comm_constraints),nb=n, fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=true, conjubasis=false, QUIET=true)
@time qs_tssos_first([f;comm_constraints], q, n, 3, numeq=length(comm_constraints), rncnumeq= length(comm_constraints),nb=n, fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=true, normality = 1,conjubasis=false, QUIET=true)
# @time qs_tssos_first([f], q, n, 2, nb=n, fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=false, conjubasis=true, QUIET=true)
@time qs_tssos_first([f;gn;comm_constraints], q, n, 2, numeq=n+length(comm_constraints), rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=true, conjubasis=true, QUIET=true)
@time qs_tssos_first([f;comm_constraints], q, n, 2, numeq=length(comm_constraints), rncnumeq= length(comm_constraints), nb=n, fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=true, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([fr; gn], q)
@time tssos_first(pop, x, 2, numeq=n, TS=false, solution=true, QUIET=true)

####Test4 correlative sparsity

### 按cliques构造二次函数 \sum_{i=1}^k[\q_i]_1^*Q_i[\q_i]_1 

###[n,cn,size] = [100,33,4],[200,67,4],[300,100,4]

rng = Xoshiro(1)
n = 100
cn = 33 #number of cliques
size =4 #cliquesize
@ncpolyvar q[1:2n]
# f,gs = cliques_sparserandomsymfunc(q, n, cn,size,1,rng, 0.1,conjugates=false, coelimit=false)
f,gs,Fsupp,Fcoe = cliques_randomsymfunc(q, n, cn,size,1,rng,conjugates=false, coelimit=false)

## sphere

# QSOS: d=1
@time qs_tssos_first([f;gs], q, n, 1,fsupp=Fsupp, fcoe=Fcoe, numeq=cn, TS=false,CS ="MF",ipart=false,conjubasis=true, QUIET=false)

# RSOS: d=1
pop,x = quaternion_to_real([f; gs], q)
@time opt,sol,data = tssos_first(pop, x, 1, numeq=cn, TS="MD", GroebnerBasis=false, solution=false, QUIET=false)



####Test5 

### 按cliques构造四次函数 \sum_{i=1}^k[\q_i]_2^*Q_i[\q_i]_2

###[n,cn,size] = [60,20,4], [90,30,4], [120,40,4]

rng = Xoshiro(1)
n = 60
cn = 20 #number of cliques
size =4 #cliquesize
@ncpolyvar q[1:2n]
f,gs,Fsupp,Fcoe  = cliques_randomsymfunc(q, n, cn,size,2,rng,conjugates=false, coelimit=false)
# f,g = cliques_sparserandomsymfunc(q, n, cn,size,2,rng, 0.05,conjugates=false, coelimit=false)

## sphere

# QSOS :d=2
@time qs_tssos_first([f;gs], q, n, 2, fsupp=Fsupp, fcoe=Fcoe, numeq=cn, CS="MF",TS=false,ipart=false, conjubasis=false, QUIET=false)

# strengthened QSOS
@time qs_tssos_first([f;gs], q, n, 2, fsupp=Fsupp, fcoe=Fcoe, numeq=cn, CS="MF",TS=false,ipart=false, normality = 1, conjubasis=false, QUIET=false)

# RSOS :d=2
pop,x = quaternion_to_real([f; gs], q)
@time opt,sol,data = cs_tssos_first(pop, x, 2, numeq=cn, TS="MD", solution=true, QUIET=false)


#### Application1 QMMC

### n = 20,30,40,60

### set
using Serialization

## step1 generate data
if !isdir("data")
    mkdir("data")
end

function create_sample_quaternion_vector(seed, n)
    rng = MersenneTwister(seed)   # 新建种子固定的随机数生成器
    return [quat(0, randn(rng), randn(rng), randn(rng)) for _ in 1:n]
end

function generate_data()
    n = 20 # 每个样本的四元数向量长度
    class1 = [create_sample_quaternion_vector(i, n) for i in 1:5]
    class2 = [create_sample_quaternion_vector(i+100, n) for i in 1:5] #20,50,100

    serialize("data/class1.jls", class1)
    serialize("data/class2.jls", class2)
end

generate_data()

## step2
function load_data(path)
    return deserialize(path)
end

function mean_vector(vectors::Vector{Vector{Quaternion{Float64}}})
    m = length(vectors)          # 样本数
    n = length(vectors[1])       # 每个样本向量长度
    result = Vector{Quaternion{Float64}}(undef, n)
    for i in 1:n
        s = quat(0.0, 0.0, 0.0, 0.0)
        for j in 1:m
            s += vectors[j][i]
        end
        result[i] = s / m
    end
    return result
end

function compute_scatter_matrices(class1_path, class2_path)
    X1 = load_data(class1_path)
    X2 = load_data(class2_path)

    n1 = length(X1)
    n2 = length(X2)
    μ1 = mean_vector(X1)
    μ2 = mean_vector(X2)
    μ = mean_vector(vcat(X1,X2))

    n = length(μ)
    S_B = n1*(μ1 - μ)*conj.(transpose(μ1 - μ))+n2*(μ2 - μ)*conj.(transpose(μ2 - μ))
    S_W = [quat(0, 0, 0, 0) for _ in 1:n, _ in 1:n]
    for x in X1
            S_W+= (x- μ1) * conj.(transpose(x - μ1))
    end
    for x in X2
            S_W+= (x- μ2) * conj.(transpose(x - μ2))
    end

    return S_B, S_W
end

S_B , S_W = compute_scatter_matrices("data/class1.jls", "data/class2.jls")

## step3 solve
n = length(S_B[1,:])
@ncpolyvar q[1:2n]  # n变量及共轭

Q = -(S_B - S_W)
Q_sym = (Q + Q') / 2
f = transpose(q[n+1:2n])*Q_sym*q[1:n]
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
P = conj.(Q_sym)
fr = transpose(q[1:n])*P*q[n+1:2n]

# QSOS d=1
@time qs_tssos_first([f,g], q, n, 1, numeq=1, TS=false,CS=false,ipart=true, conjubasis=false,QUIET=false)

# RSOS d=1
pop,x = quaternion_to_real([fr;g], q)
@time opt,sol,data = tssos_first(pop, x, 1, numeq=1, TS="MD",solution=false, QUIET=false)

# 从 moment矩阵里恢复近似最优解
@time opt,sol = qs_tssos_first([f,g], q, n, 1, numeq=1, TS=false,CS=false,ipart=true, conjubasis=false, solution=true,QUIET=true)
# 计算模长
qnorm = norm(sol)
# 如果模长不是1，则归一化
if abs(qnorm - 1) > 1e-6  # 允许一定误差
    sol = sol ./ qnorm
end
upper_bound = real(transpose(sol)*P*conj.(sol))
println("Upper bound:",upper_bound)


#### Application2 旋转同步

# 生成单位四元数
function rand_unit_quaternion(rng)
    v = randn(rng,4)
    v ./= norm(v)
    return Quaternion(v[1], v[2], v[3], v[4])
end

# 设置参数
n = 60                 # 节点个数
rng = Xoshiro(1)        # 使用种子为1的伪随机数生成器（可复现）
prob_density = 0.2     # 图的稀疏程度
noise_level = 0.2       # 噪声大小

# Step 1: 数据生成
qs_true = [rand_unit_quaternion(rng) for i in 1:n]

edges_list = []
for i in 1:n
    for j in i+1:n
        if rand(rng) < prob_density
            push!(edges_list, (i, j))
        end
    end
end

Q = Dict()
for (i,j) in edges_list
    noise = rand_unit_quaternion(rng)
    Q[(i,j)] = qs_true[i] * conj(qs_true[j]) + noise_level * noise
end
# Step 2: 构造目标函数
@ncpolyvar q[1:2n]

f = zero(q[1])
for (i,j) in keys(Q)
    term = q[i+n]*Q[(i,j)]/2*q[j]+ q[j+n]*conj(Q[(i,j)])/2*q[i]
    f -= term
end
gn = [1 - q[i]*q[i+n] for i in 1:n] 

# Step 3: Solve
opt = @time qs_tssos_first([f], q, n, 1, nb=n, TS="MD",CS=false,ipart=true, conjubasis=false, solution=false,QUIET=false)
# opt = @time qs_tssos_first([f], q, n, 1, nb=n, TS=false,CS=false,ipart=true, conjubasis=false, solution=false,QUIET=false)
# Step 4 : qualify

# TSSOS
fr = zero(q[1])
for (i,j) in keys(Q)
    term = q[i]*conj(Q[(i,j)])/2*q[j+n]+ q[j]*Q[(i,j)]/2*q[i+n]
    fr -= term
end
pop,x = quaternion_to_real([fr;gn], q)
@time opt,sol,data = tssos_first(pop, x, 1, numeq=n, TS="MD", GroebnerBasis=false, solution=false, QUIET=false)
# @time opt,sol,data = tssos_first(pop, x, 1, numeq=n, TS=false,solution=false, QUIET=false)
# println(data.blocksize[1])

# true value
function objective(sol, Q)
    f = zero(sol[1])
    for (i,j) in keys(Q)
        term = conj(sol[i])*Q[(i,j)]/2*sol[j] + conj(sol[j])*conj(Q[(i,j)])/2*sol[i]
        f -= term
    end
    return real(f)
end
qs_all = vcat(qs_true, [conj(q) for q in qs_true])
val_true = objective(qs_all, Q)
println("true value:",val_true)


###set: f = [q]_2^*Q[q]_2, Q:real, n=1,2,3,4,5,6；

rng = Xoshiro(3)
n = 3
@ncpolyvar q[1:2n]
f,Fsupp,Fcoe = randomsymfunc(q, n, 2, rng, conjugates=true, coelimit=true)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## ball
# QSOS:basis[q]_2
@time qs_tssos_first([f, g], q, n, 3, fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=false, conjubasis=false, QUIET=true) 

# strengthened QSOS :basis [q]_2+\bar(q_i)[q]_2
@time qs_tssos_first([f, g], q, n, 2, fsupp=Fsupp, fcoe=Fcoe, CS=false,TS=false, ipart=false, normality = 1,conjubasis=false, QUIET=true)

#QSOS: full basis [q,\bar(q)]_2
@time qs_tssos_first([f, g], q, n, 2, fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=false, conjubasis=true, QUIET=true) 

#RSOS
pop,x = quaternion_to_real([f, g], q)
@time tssos_first(pop, x, 2, TS=false, solution=true, QUIET=true)

## unit norm

@time qs_tssos_first([f], q, n, 2, nb=n, fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=false, conjubasis=false, QUIET=true)
@time qs_tssos_first([f], q, n, 2, nb=n, fsupp=Fsupp, fcoe=Fcoe, CS=false,TS=false, ipart=false, normality = 1,conjubasis=false, QUIET=true)
@time qs_tssos_first([f], q, n, 2, nb=n, fsupp=Fsupp, fcoe=Fcoe, TS=false, CS=false, ipart=false, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([f; gn], q)
@time tssos_first(pop, x, 2, numeq=n, TS=false, solution=true, QUIET=true)


n = 3
@ncpolyvar q[1:2n]
# f,Fsupp,Fcoe = randomsymfunc(q, n, 2, rng, conjugates=true, coelimit=true)
# f = (q[1]*q[2]-q[2]*q[1])+(q[4]*q[3]-q[3]*q[4])
# f = q[1]*q[2]*q[3]*q[4]+q[8]*q[7]*q[6]*q[5]
# f = q[3]*q[4]*q[3]*q[2]+q[4]*q[1]*q[2]*q[1]
f = q[1]*q[2]*q[3]+q[6]*q[5]*q[4]-q[1]*q[3]*q[2]-q[5]*q[6]*q[4]
gn = [1 - q[i]*q[i+n] for i = 1:n]
@time qs_tssos_first([f], q, n, 2, nb=n, TS=false, CS=false, ipart=false, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([f; gn], q)
@time tssos_first(pop, x, 2, numeq=n, TS=false, solution=true, QUIET=true)
println( QPolys_info([f],q,n))

n = 2
@ncpolyvar q[1:2n]
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
f_non = q[1]^2+q[2]^2+q[3]^2+q[4]^2+q[1]*q[2]+q[4]*q[3]+q[1]*q[4]+q[2]*q[3]
comm_constraints = typeof(f_non)[]
for i in 1:n
    for j in 1:n
        if i != j
            expr = q[i]*q[j] + q[i+n]*q[j] - q[j]*q[i] - q[j]*q[i+n]
            push!(comm_constraints, expr)
        end
    end
end
@time opt= qs_tssos_first([f_non, g], q, n, 1, TS=false, ipart=false,conjubasis=true, QUIET=true)
@time qs_tssos_first([[f_non,g];comm_constraints], q, n, 1, numeq=length(comm_constraints),rncnumeq=length(comm_constraints), TS=false, ipart=false, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([f_non, g], q)
@time tssos_first(pop, x, 2, TS=false, solution=false, QUIET=false)

n = 3
@ncpolyvar q[1:2n]
gn = [1 - q[i]*q[i+n] for i = 1:n]
gc = [(q[1]+q[4])*(q[2]+q[5])-(q[2]+q[5])*(q[1]+q[4]), (q[2]+q[5])*(q[3]+q[6])-(q[3]+q[6])*(q[2]+q[5]), (q[3]+q[6])*(q[1]+q[4])-(q[1]+q[4])*(q[3]+q[6])]
comm_constraints = typeof(f_non)[]
for i in 1:n
    for j in 1:n
        if i != j
            expr = q[i]*q[j] + q[i+n]*q[j] - q[j]*q[i] - q[j]*q[i+n]
            push!(comm_constraints, expr)
        end
    end
end
println(comm_constraints)
lc = length(comm_constraints)
f_non = 1+1/4*(q[1]+q[4])*(q[2]+q[5])+1/4*(q[2]+q[5])*(q[3]+q[6])+1/4*(q[3]+q[6])*(q[1]+q[4])
f_non = 1+1/8*(q[1]+q[4])*(q[2]+q[5])+1/8*(q[2]+q[5])*(q[1]+q[4])+1/8*(q[2]+q[5])*(q[3]+q[6])+1/8*(q[3]+q[6])*(q[2]+q[5])+1/8*(q[3]+q[6])*(q[1]+q[4])+1/8*(q[1]+q[4])*(q[3]+q[6])
@time qs_tssos_first([f_non], q, n, 2 , nb=n, TS=false, ipart=false,conjubasis=true, QUIET=false)
@time qs_tssos_first([f_non;comm_constraints], q, n, 2, numeq=6, nb=n, TS=false, ipart=false, conjubasis=true, QUIET=true)
@time qs_tssos_first([f_non;gn;comm_constraints], q, n, 2 , numeq=n+6, rncnumeq=6,TS=false, ipart=false,conjubasis=true, QUIET=false)
pop,x = quaternion_to_real([f_non;gn], q)
@time tssos_first(pop, x, 2, numeq=n, TS=false, solution=true, QUIET=false)


n = 3
@ncpolyvar q[1:2n]
gn = [1 - q[i]*q[i+n] for i = 1:n]
f = 5/4 + 7/32*q[1]*q[2] + 9/32*q[1]*q[5] + 7/32*q[4]*q[2] + 7/32*q[4]*q[5] + 1/16*q[2]*q[4] + 7/32*q[1]*q[3] + 9/32*q[1]*q[6] + 7/32*q[4]*q[3] + 7/32*q[4]*q[6] + 1/16*q[3]*q[4] + 7/32*q[2]*q[3] + 9/32*q[2]*q[6] + 7/32*q[5]*q[3] + 7/32*q[5]*q[6] + 1/16*q[3]*q[5]
f = 5/4 + 7/32*(q[1]+q[4])*(q[2]+q[5])+1/16*(q[1]*q[5]+q[2]*q[4])+ 7/32*(q[1]+q[4])*(q[3]+q[6])+1/16*(q[1]*q[6]+q[3]*q[4])+ 7/32*(q[2]+q[5])*(q[3]+q[6])+1/16*(q[2]*q[6]+q[3]*q[5])
comm_constraints = typeof(f)[]
for i in 1:n
    for j in 1:n
        if i != j
            expr = q[i]*q[j] + q[i+n]*q[j] - q[j]*q[i] - q[j]*q[i+n]
            push!(comm_constraints, expr)
        end
    end
end
g3 = 1 - quat(0,1,1,0)*q[1] - quat(0,0,1,0)*q[2]
@time qs_tssos_first([f;gn;comm_constraints], q, n, 2 , numeq=n+6, rncnumeq=6,TS=false, ipart=false,conjubasis=true, QUIET=false)
@time qs_tssos_first([f;gn;comm_constraints;g3], q, n, 2 , numeq=n+7, rncnumeq=6, qncnumeq=1,TS=false, ipart=false,conjubasis=true, QUIET=false)
pop,x = quaternion_to_real([f;gn], q)
@time tssos_first(pop, x, 2, numeq=n, TS=false, solution=true, QUIET=false)
pop,x = quaternion_to_real([f;gn;g3], q)
@time tssos_first(pop, x, 2, numeq=n+1, TS=false, solution=true, QUIET=false)


rng = Xoshiro(5)
n = 3
@ncpolyvar q[1:2n]
# f,fr,Fsupp,Fcoe = qrandomsymfunc(q, n, 1, rng; conjugates=true)
# f = quat(1,0,0,0) + q[4]*quat(2,0,0,0)*q[1] + q[5]*quat(3,0,0,0)*q[2] + q[6]*quat(4,0,0,0)*q[3] + quat(0,1,0,0)*q[1] + q[4]*quat(0,-1,0,0) + quat(0,0,1,0)*q[2] + q[5]*quat(0,0,-1,0) + quat(0,0,0,1)*q[3] + q[6]*quat(0,0,0,-1) + q[4]*quat(1,0,0,0)*q[2] + q[5]*quat(1,0,0,0)*q[1] + q[6]*quat(0,1,0,0)*q[4] + q[1]*quat(0,-1,0,0)*q[3]
# 定义矩阵 Q (7×7，元素为四元数)
Q = zeros(Quaternion{Float64}, 7, 7)

# 对角线 (实数)
Q[1,1] = quat(1,0,0,0)   # 对应 1
Q[2,2] = quat(2,0,0,0)   # q1* q1
Q[3,3] = quat(3,0,0,0)   # q2* q2
Q[4,4] = quat(4,0,0,0)   # q3* q3
Q[5,5] = quat(5,0,0,0)   # q1 q1*
Q[6,6] = quat(6,0,0,0)   # q2 q2*
Q[7,7] = quat(7,0,0,0)   # q3 q3*

# 非对角元 (包含虚部)
Q[1,2] = quat(1,1,0,0);   Q[2,1] = conj(Q[1,2])  # i
Q[1,3] = quat(0,2,1,1);   Q[3,1] = conj(Q[1,3])  # j
Q[1,4] = quat(0,5,0,1);   Q[4,1] = conj(Q[1,4])  # k
Q[2,3] = quat(1,0,0,0);   Q[3,2] = Q[2,3]        # 1 (实数)
Q[4,5] = quat(0,6,2,1);   Q[5,4] = conj(Q[4,5])

# f,fr,Fsupp,Fcoe = qrandomsymfunc2(q,n,1,rng,Q;conjugates=true,given=true)
f,fr,Fsupp,Fcoe = qrandomsymfunc2(q,n,1,rng,Q;conjugates=true,given=false)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]
comm_constraints = typeof(f)[]
for i in 1:n
    for j in 1:n
        if i != j
            expr = q[i]*q[j] + q[i+n]*q[j] - q[j]*q[i] - q[j]*q[i+n]
            push!(comm_constraints, expr)
        end
    end
end
println(f)
## ball

# QSOS: basis[q,\bar(q)]_1
@time qs_tssos_first([f, g], q, n, 1, fsupp=Fsupp, fcoe=Fcoe, TS=false, ipart=true, conjubasis=true, QUIET=true)
@time qs_tssos_first([f, g], q, n, 1, TS=false, ipart=true, conjubasis=true, QUIET=true)
@time qs_tssos_first([[f, g];comm_constraints], q, n, 1, numeq = length(comm_constraints), rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=true, conjubasis=true, QUIET=true) 
# RSOS: d=1
pop,x = quaternion_to_real([fr, g], q)
@time opt,sol,data = tssos_first(pop, x, 1, TS=false, solution=true, QUIET=true)
@time opt,sol,data = tssos_first(pop, x, 2, TS=false, solution=true, QUIET=true)



### set Q:quaternion
rng = Xoshiro(1)
n = 2
@ncpolyvar q[1:2n]
# f,fr,Fsupp,Fcoe = qrandomsymfunc(q, n, 2, rng; conjugates=false)
f,fr,Fsupp,Fcoe = qrandomsymfunc2(q,n,2,rng,Q;conjugates=false,given=false)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]
gc1 = 2*q[1]*q[3]*q[2]*q[4]+q[1]*q[4]*q[1]*q[4]+q[2]*q[3]*q[2]*q[3]-q[3]*q[3]*q[2]*q[2]-q[4]*q[4]*q[1]*q[1]-q[4]*q[3]*q[2]*q[1]-q[3]*q[4]*q[1]*q[2]
gc = 2*q[1]*q[3]*q[2]*q[4]+q[1]*q[4]*q[1]*q[4]+q[2]*q[3]*q[2]*q[3]-q[3]*q[3]*q[2]*q[2]-q[4]*q[4]*q[1]*q[1]-q[1]*q[4]*q[3]*q[2]-q[2]*q[3]*q[4]*q[1]
gc2 = q[4]*q[3]*q[2]-q[3]*q[4]*q[2]
gc3 = q[4]*q[1]*q[2]-q[1]*q[4]*q[2]
comm_constraints = typeof(f)[]
for i in 1:n
    for j in 1:n
        if i != j
            expr = q[i]*q[j] + q[i+n]*q[j] - q[j]*q[i] - q[j]*q[i+n]
            push!(comm_constraints, expr)
        end
    end
end
Fsupp = Vector{Vector{UInt16}}[[[0x0002], [0x0004], [0x0001]], [[0x0002], [0x0004], [0x0001]], [[0x0002], [0x0004], [0x0003]], [[], [], [0x0001, 0x0001, 0x0004, 0x0004]], [[0x0002], [0x0004], [0x0003]], [[], [], [0x0002, 0x0002, 0x0003, 0x0003]]]
println(Fsupp)
println(Fcoe)
println(f)
println(length(Fsupp))

## ball
# QSOS:basis[q]_2
@time qs_tssos_first([f, g], q, n, 2, fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=true, conjubasis=false, QUIET=true) 
@time qs_tssos_first([[f, g];comm_constraints], q, n, 2, numeq = length(comm_constraints), rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=true, conjubasis=false, QUIET=true) 
# strengthened QSOS :basis [q]_2+\bar(q_i)[q]_2
@time qs_tssos_first([f, g], q, n, 2, fsupp=Fsupp, fcoe=Fcoe, CS=false,TS=false, ipart=false, normality = 1,conjubasis=false, QUIET=true)
@time qs_tssos_first([[f, g];comm_constraints], q, n, 2, numeq = length(comm_constraints), rncnumeq= length(comm_constraints), fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=true, normality = 1, conjubasis=false, QUIET=true) 

#QSOS: full basis [q,\bar(q)]_2
@time qs_tssos_first([f, g], q, n, 2, fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=false, conjubasis=true, QUIET=true) 
@time qs_tssos_first([[f, g];comm_constraints], q, n, 2, numeq = length(comm_constraints), rncnumeq= length(comm_constraints),fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=true, conjubasis=true, QUIET=true) 
@time qs_tssos_first([[f, g];gc1;gc3;comm_constraints], q, n, 2, numeq = length(comm_constraints)+2, rncnumeq= length(comm_constraints)+2,fsupp=Fsupp, fcoe=Fcoe, CS=false, TS=false, ipart=true, conjubasis=true, QUIET=true) 

pop,x = quaternion_to_real([fr, g], q)
@time opt,sol,data = tssos_first(pop, x, 2, TS=false, solution=true, QUIET=true)
ub = local_solution(data.n,data.m,data.supp,data.coe;nb=data.nb,numeq=data.numeq,startpoint=rand(data.n),QUIET=true)[1]
println(ub)


