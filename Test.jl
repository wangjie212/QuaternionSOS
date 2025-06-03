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

# Test 1
## n = 20 seed = 4,2,3; n = 40 seed = 1,2,3; n = 60 seed = 11,12,13
rng = Xoshiro(11)
n = 20
@ncpolyvar q[1:2n]
f = randomsymfunc(q, n, 1, rng, conjugates=false, coelimit=false)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## ball
@time qs_tssos_first([f, g], q, n, 1, TS=false, ipart=false,conjubasis=false, QUIET=true)
pop,x = quaternion_to_real([f, g], q)
@time tssos_first(pop, x, 1, TS=false, solution=true, QUIET=true)


## unit norm
@time qs_tssos_first([f], q, n, 1, nb=n, TS=false, ipart=false, conjubasis=false, QUIET=true)
pop,x = quaternion_to_real([f;gn], q)
@time tssos_first(pop, x, 1, numeq=n, TS=false, solution=true, QUIET=true)


# Test 2
# seed = 1,2,3
rng = Xoshiro(1)
n = 20
@ncpolyvar q[1:2n]
f = randomsymfunc(q, n, 1, rng, conjugates=true, coelimit=false)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## ball
@time qs_tssos_first([f, g], q, n, 1, TS=false, ipart=false, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([f, g], q)
@time tssos_first(pop, x, 1, TS=false, solution=true, QUIET=true)

## unit norm
@time qs_tssos_first([f], q, n, 1, nb=n, TS=false, ipart=false, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([f; gn], q)
@time tssos_first(pop, x, 1, numeq=n, TS=false, solution=true, QUIET=true)

# Test 3
# seed = 11,12,13
rng = Xoshiro(13)
n = 5
@ncpolyvar q[1:2n]
f = randomsymfunc(q, n, 2, rng, conjugates=true, coelimit=true)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## ball
@time qs_tssos_first([f, g], q, n, 2, QUIET=true, TS=false, ipart=false, conjubasis=true)
pop,x = quaternion_to_real([f, g], q)
@time tssos_first(pop, x, 2, QUIET=true, TS=false, solution=true)

## unit norm
@time qs_tssos_first([f], q, n, 2, nb=n, TS=false, ipart=false, QUIET=true, conjubasis=true)
pop,x = quaternion_to_real([f; gn], q)
@time tssos_first(pop, x, 2, numeq=n, TS=false, QUIET=true)


# Test 4
# seed = 30,60,90
rng = Xoshiro(90)
n = 5
@ncpolyvar q[1:2n]
f = randomsymfunc(q, n, 2, rng, conjugates=false, coelimit=true)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## ball
@time qs_tssos_first([f, g], q, n, 2, CS=false,TS=false, ipart=false, conjubasis=true, QUIET=true)
@time qs_tssos_first([f, g], q, n, 2, CS=false,TS=false, ipart=false, normality = 1,conjubasis=false, QUIET=true)
pop,x = quaternion_to_real([f, g], q)
@time tssos_first(pop, x, 2, TS=false, solution=true, QUIET=true)

## unit norm
@time qs_tssos_first([f], q, n, 2, nb=n, TS=false, ipart=false, conjubasis=false, QUIET=true)
pop,x = quaternion_to_real([f; gn], q)
@time tssos_first(pop, x, 2, numeq=n, TS=false, solution=true, QUIET=true)

# Test 5
## n = 10,20,30 seed = 1,2,3
rng = Xoshiro(1)
n = 50
@ncpolyvar q[1:2n]
f,fr = qrandomsymfunc(q, n, 1, rng, conjugates=false)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## ball
@time qs_tssos_first([f, g], q, n, 1, TS=false, ipart=true,conjubasis=false, QUIET=true)
pop,x = quaternion_to_real([fr, g], q)
@time tssos_first(pop, x, 1, TS=false, solution=true, QUIET=true)

## unit norm
@time qs_tssos_first([f], q, n, 1, nb=n, TS=false, ipart=true, conjubasis=false, QUIET=true)
pop,x = quaternion_to_real([fr;gn], q)
@time tssos_first(pop, x, 1, numeq=n, TS=false, solution=true, QUIET=true)

# sparsity

# Test 1
## seed = 1,2,3
rng = Xoshiro(1)
n = 10
@ncpolyvar q[1:2n]
f = sparserandomsymfunc(q, n, 1, rng, 0.1, conjugates=false, coelimit=false)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## ball
@time qs_tssos_first([f, g], q, n, 1, TS="MD", ipart=false,conjubasis=false, QUIET=true)
pop,x = quaternion_to_real([f, g], q)
@time tssos_first(pop, x, 1, TS=false, solution=true, QUIET=true)
## unit norm
# @time qs_tssos_first([f;gn], q, n, 1, numeq=n, TS="MD", ipart=false, conjubasis=false, QUIET=true)
@time qs_tssos_first([f], q, n, 1, nb=n, TS="MD",CS=false,ipart=false, conjubasis=false, QUIET=false)
pop,x = quaternion_to_real([f;gn], q)
@time opt,sol,data = tssos_first(pop, x, 1, numeq=n, TS="MD", GroebnerBasis=false, solution=false, QUIET=false)


# Test 2
# seed = 1,2,3
rng = Xoshiro(1)
n = 40
@ncpolyvar q[1:2n]
f = sparserandomsymfunc(q, n, 1,rng, 0.1,conjugates=true, coelimit=false)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## ball
@time qs_tssos_first([f, g], q, n, 1, TS="MD", ipart=false, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([f, g], q)
@time tssos_first(pop, x, 1, TS=false, solution=true, QUIET=true)

## unit norm
@time qs_tssos_first([f], q, n, 1, nb=n, TS=false, ipart=false, conjubasis=true, QUIET=true)
@time qs_tssos_first([f], q, n, 1, nb=n, TS="MD",CS=false,ipart=false, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([f; gn], q)
# @time tssos_first(pop, x, 1, numeq=n, TS=false, solution=true, QUIET=true)
@time opt,sol,data = tssos_first(pop, x, 1, numeq=n, TS="MD", GroebnerBasis=false, solution=false, QUIET=false)

# Test 3
# seed = 11,12,13
rng = Xoshiro(1)
n = 5
@ncpolyvar q[1:2n]
f = sparserandomsymfunc(q, n, 2, rng, 0.2, conjugates=true, coelimit=true)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## ball
@time qs_tssos_first([f, g], q, n, 2, QUIET=true, TS="MD", ipart=false, conjubasis=true)
pop,x = quaternion_to_real([f, g], q)
@time tssos_first(pop, x, 2, QUIET=true, TS="MD", solution=true)

## unit norm
@time qs_tssos_first([f], q, n, 2, nb=n, TS="MD",ipart=false, conjubasis=true, QUIET=false)
pop,x = quaternion_to_real([f; gn], q)
@time opt,sol,data = tssos_first(pop, x, 2, numeq=n, TS="MD", GroebnerBasis=false, solution=false, QUIET=false)
@time opt,sol,data = cs_tssos_first(pop, x, 2, numeq=n, TS="MD", solution=false, QUIET=false)


# Test 4
# seed = 1,2,3 n= 
rng = Xoshiro(3)
n = 5
@ncpolyvar q[1:2n]
f = sparserandomsymfunc(q, n, 2, rng, 0.2, conjugates=false, coelimit=true)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]
## ball
@time qs_tssos_first([f, g], q, n, 2, TS="MD", ipart=false, conjubasis=true, QUIET=true) 
@time qs_tssos_first([f, g], q, n, 2, normality = 1, TS="MD", ipart=false, conjubasis=false, QUIET=true) 
pop,x = quaternion_to_real([f, g], q)
@time tssos_first(pop, x, 2, TS="MD", solution=true, QUIET=true)

## unit norm
@time qs_tssos_first([f], q, n, 2, nb=n, TS="MD",ipart=false, conjubasis=true, QUIET=false)
# @time qs_tssos_first([f], q, n, 2, nb=n, TS="MD",normality = 1,ipart=false, conjubasis=false, QUIET=false)
pop,x = quaternion_to_real([f; gn], q)
@time opt,sol,data = cs_tssos_first(pop, x, 2, numeq=n, TS="MD", solution=false, QUIET=false)


# 按cliques构造函数

#Test5 [n,cn,size] = [50,7,8] [100,33,4],[150,16,10],[200,11,20]
#Test5 [n,cn,size] = [50,17,4] [100,33,4],[200,67,4],[300,100,4]

rng = Xoshiro(1)
n = 100
cn = 33 #number of cliques
size =4 #cliquesize
@ncpolyvar q[1:2n]
# f,gs = cliques_sparserandomsymfunc(q, n, cn,size,1,rng, 0.1,conjugates=false, coelimit=false)
f,gs = cliques_randomsymfunc(q, n, cn,size,1,rng,conjugates=false, coelimit=false)

## sphere
@time qs_tssos_first([f;gs], q, n, 1, numeq=cn, TS="MD",CS =false,ipart=false, conjubasis=false, QUIET=false)
pop,x = quaternion_to_real([f; gs], q)
@time opt,sol,data = tssos_first(pop, x, 1, numeq=cn, TS="MD", GroebnerBasis=false, solution=false, QUIET=false)



#Test6 [n,cn,size] = [10,3,4],[20,5,5],[30,6,6]
#Test6 [n,cn,size] = [10,3,4],[20,7,4],[30,10,4]
rng = Xoshiro(1)
n = 30
cn = 6 #number of cliques
size =6 #cliquesize
@ncpolyvar q[1:2n]
f,gs = cliques_randomsymfunc(q, n, cn,size,2,rng,conjugates=false, coelimit=false)
# f,g = cliques_sparserandomsymfunc(q, n, cn,size,2,rng, 0.05,conjugates=false, coelimit=false)
## sphere
@time qs_tssos_first([f;gs], q, n, 2, numeq=cn, TS="MD",ipart=false, conjubasis=true, QUIET=false)
@time qs_tssos_first([f;gs], q, n, 2, numeq=cn, TS="MD",ipart=false, normality = 1, conjubasis=false, QUIET=false)
pop,x = quaternion_to_real([f; gs], q)
@time opt,sol,data = cs_tssos_first(pop, x, 2, numeq=cn, TS="MD", solution=true, QUIET=false)

#Test7

rng = Xoshiro(3)
n = 100
cn = 10
size =11
@ncpolyvar q[1:2n]
f,g = cliques_sparserandomsymfunc(q, n, cn,size,1,rng, 0.1,conjugates=true, coelimit=false)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## ball
@time qs_tssos_first([f, g], q, n, 1, TS="MD", ipart=false, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([f, g], q)
@time tssos_first(pop, x, 1, TS=false, solution=true, QUIET=true)

## sphere
@time qs_tssos_first([f;g], q, n, 1, numeq=cn, TS="MD",CS=false,ipart=false, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([f; g], q)
@time opt,sol,data = tssos_first(pop, x, 1, numeq=cn, TS="MD", GroebnerBasis=false, solution=false, QUIET=false)



#Test8
rng = Xoshiro(2)
n = 7
cn = 3
size =3
@ncpolyvar q[1:2n]
f,g = cliques_sparserandomsymfunc(q, n, cn,size,2,rng, 0.05,conjugates=true, coelimit=false)
println(g)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## sphere
@time qs_tssos_first([f;g], q, n, 2, numeq=cn, TS="MD",ipart=false, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([f; g], q)
@time opt,sol,data = cs_tssos_first(pop, x, 2, numeq=cn, TS="MD", solution=false, QUIET=false)

@time opt,sol,data = tssos_first(pop, x, 2, numeq=cn, TS="MD", GroebnerBasis=false, solution=false, QUIET=false)

# actual

# Test1 
using Serialization
#step1 generate data
if !isdir("data")
    mkdir("data")
end
function create_sample_quaternion_vector(seed, n)
    rng = MersenneTwister(seed)   # 新建种子固定的随机数生成器
    return [quat(0, randn(rng), randn(rng), randn(rng)) for _ in 1:n]
end

function generate_data()
    n = 15 # 每个样本的四元数向量长度
    class1 = [create_sample_quaternion_vector(i, n) for i in 1:5]
    class2 = [create_sample_quaternion_vector(i+100, n) for i in 1:5]

    serialize("data/class1.jls", class1)
    serialize("data/class2.jls", class2)
end

generate_data()

# step2
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

#step3 solve
n = length(S_B[1,:])
@ncpolyvar q[1:2n]  # n变量及共轭
Q = -(S_B - S_W)
Q_sym = (Q + Q') / 2
f = transpose(q[n+1:2n])*Q_sym*q[1:n]
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
P = conj.(Q_sym)
fr = transpose(q[1:n])*P*q[n+1:2n]
@time qs_tssos_first([f,g], q, n, 1, numeq=1, TS=false,CS=false,ipart=true, conjubasis=false, QUIET=false)
pop,x = quaternion_to_real([fr;g], q)
@time opt,sol,data = tssos_first(pop, x, 1, numeq=1, TS=false,solution=true, QUIET=false)


#Test2

# 生成单位四元数
function rand_unit_quaternion()
    v = randn(4)
    v ./= norm(v)
    return Quaternion(v[1], v[2], v[3], v[4])
end

# 设置参数
n = 10 
prob_density = 0.5  # 边的稀疏程度
noise_level = 0.01  # 噪声大小

# Step 1: 数据生成
qs_true = [rand_unit_quaternion() for i in 1:n]

edges_list = []
for i in 1:n
    for j in i+1:n
        if rand() < prob_density
            push!(edges_list, (i, j))
        end
    end
end

Qij = Dict()
for (i,j) in edges_list
    noise = rand_unit_quaternion()
    Qij[(i,j)] = qs_true[i] * conj(qs_true[j]) + noise_level * noise
end

# Step 2: 构造目标函数
@ncpolyvar q[1:2n]

f = zero(q[1])
for (i,j) in keys(Qij)
    term = q[i+n]*Qij[(i,j)]/2*q[j]+ q[j+n]*conj(Qij[(i,j)])/2*q[i]
    f -= term
end
gn = [1 - q[i]*q[i+n] for i in 1:n] 

# Step 4: 
@time qs_tssos_first([f], q, n, 1, nb=n, TS=false,CS=false,ipart=true, conjubasis=false, QUIET=false)
fr = zero(q[1])
for (i,j) in keys(Qij)
    term = q[i]*conj(Qij[(i,j)])/2*q[j+n]+ q[j]*Qij[(i,j)]/2*q[i+n]
    fr -= term
end
pop,x = quaternion_to_real([fr;gn], q)
@time opt,sol,data = tssos_first(pop, x, 1, numeq=n, TS=false,solution=true, QUIET=false)


