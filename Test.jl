using Random
using TSSOS
using DynamicPolynomials
using Quaternions
using MosekTools
using JuMP
using LinearAlgebra

include("D:/Programs/QuaternionSOS/ncutils.jl")
include("D:/Programs/QuaternionSOS/Qpop.jl")

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
n = 2
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
n = 2
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
## n = 20,30 seed = 1,2,3
rng = Xoshiro(1)
n = 20
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
