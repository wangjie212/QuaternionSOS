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
n = 4
@ncpolyvar q[1:2n]
f = randomsymfunc(q, n, 2, rng, conjugates=false, coelimit=true)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## ball
@time qs_tssos_first([f, g], q, n, 2, TS=false, ipart=false, conjubasis=true, QUIET=true)
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
n = 140
@ncpolyvar q[1:2n]
f = sparserandomsymfunc(q, n, 1, rng, 0.2, conjugates=false, coelimit=false)
g = 1 - sum(q[i]*q[i+n] for i = 1:n)
gn = [1 - q[i]*q[i+n] for i = 1:n]

## ball
@time qs_tssos_first([f, g], q, n, 1, TS="MD", ipart=false,conjubasis=false, QUIET=true)
pop,x = quaternion_to_real([f, g], q)
@time tssos_first(pop, x, 1, TS=false, solution=true, QUIET=true)
## unit norm
@time qs_tssos_first([f;gn], q, n, 1, numeq=n, TS="MD", ipart=false, conjubasis=false, QUIET=true)
@time qs_tssos_first([f], q, n, 1, nb=n, TS="MD",CS=false,ipart=false, conjubasis=false, QUIET=true)
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
# @time qs_tssos_first([f], q, n, 2, nb=n, TS=false, ipart=false, QUIET=true, conjubasis=true)
# pop,x = quaternion_to_real([f; gn], q)
# @time tssos_first(pop, x, 2, numeq=n, TS=false, QUIET=true)
@time qs_tssos_first([f], q, n, 2, nb=n, TS="MD",ipart=false, conjubasis=true, QUIET=false)
pop,x = quaternion_to_real([f; gn], q)
@time opt,sol,data = tssos_first(pop, x, 2, numeq=n, TS="MD", GroebnerBasis=false, solution=false, QUIET=false)
@time opt,sol,data = cs_tssos_first(pop, x, 2, numeq=n, TS="MD", solution=false, QUIET=false)


# Test 4
# seed = 30,60,90
rng = Xoshiro(1)
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
# @time qs_tssos_first([f], q, n, 2, nb=n, TS="MD", CS=false,ipart=false, conjubasis=false, QUIET=true)
# pop,x = quaternion_to_real([f; gn], q)
# @time cs_tssos_first(pop, x, 2, numeq=n, TS="MD", CS=false,solution=true, QUIET=true)
@time qs_tssos_first([f], q, n, 2, nb=n, TS="MD",ipart=false, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([f; gn], q)
# @time tssos_first(pop, x, 1, numeq=n, TS=false, solution=true, QUIET=true)
@time opt,sol,data = cs_tssos_first(pop, x, 2, numeq=n, TS="MD", solution=false, QUIET=false)


# 按cliques构造函数

#Test5 [n,cn,size] = [50,7,8] [100,33,4],[150,16,10],[200,11,20]

rng = Xoshiro(1)
n = 50
cn = 7 #number of cliques
size =8 #cliquesize
@ncpolyvar q[1:2n]
# f,gs = cliques_sparserandomsymfunc(q, n, cn,size,1,rng, 0.1,conjugates=false, coelimit=false)
f,gs = cliques_randomsymfunc(q, n, cn,size,1,rng,conjugates=false, coelimit=false)

## ball
@time qs_tssos_first([f, g], q, n, 1, TS="MD", ipart=false, conjubasis=false, QUIET=true)
pop,x = quaternion_to_real([f, g], q)
@time tssos_first(pop, x, 1, TS=false, solution=true, QUIET=true)

## sphere
@time qs_tssos_first([f;gs], q, n, 1, numeq=cn, TS="MD",CS=false,ipart=false, conjubasis=false, QUIET=true)
pop,x = quaternion_to_real([f; gs], q)
@time opt,sol,data = tssos_first(pop, x, 1, numeq=cn, TS="MD", GroebnerBasis=false, solution=false, QUIET=false)



#Test6 [n,cn,size] = [10,3,4],[20,5,5],[30,6,6]
rng = Xoshiro(3)
n = 30
cn = 6 #number of cliques
size =6 #cliquesize
@ncpolyvar q[1:2n]
f,gs = cliques_randomsymfunc(q, n, cn,size,2,rng,conjugates=false, coelimit=false)
# f,g = cliques_sparserandomsymfunc(q, n, cn,size,2,rng, 0.05,conjugates=false, coelimit=false)
## sphere
@time qs_tssos_first([f;gs], q, n, 2, numeq=cn, TS="MD",ipart=false, conjubasis=true, QUIET=true)
pop,x = quaternion_to_real([f; gs], q)
@time opt,sol,data = cs_tssos_first(pop, x, 2, numeq=cn, TS="MD", solution=true, QUIET=false)

@time opt,sol,data = tssos_first(pop, x, 2, numeq=cn, TS="MD", GroebnerBasis=false, solution=false, QUIET=false)

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