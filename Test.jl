using Random
include("Qpop.jl")
# Example1 
# Transform to symmetric equation
n = 2 # set the number of quaternion variables
@polyvar q[1:2n]
f = quat(1,0,0,0)*q[1]*q[3]+quat(1,0,0,0)*q[2]*q[4]
g1 = quat(3,0,0,0)-quat(1,0,0,0)*q[1]*q[3]-quat(1,0,0,0)*q[2]*q[4]
g2 = quat(2,0,0,0)-quat(1,0,0,0)*q[1]-quat(1,0,0,0)*q[2]-quat(1,0,0,0)*q[3]-quat(1,0,0,0)*q[4]
g3 = quat(0,1,1,1)*q[1]-quat(0,1,1,1)*q[3]+quat(0,1,1,1)*q[2]-quat(0,1,1,1)*q[4]
qpop = [f,g1,g2,g3]
order = 2 # set the relaxation order
opt,cons,I,J= qs_tssos_first(qpop, q, n, order, numeq=2, TS=false)
sqrt(opt)
pop,x = quaternion_to_real(qpop,q)
opt2 = tssos_first(pop, x, order, numeq=2)
sqrt(opt2[1])

# Example1 delete one equation
n = 1
@polyvar q[1:2n]
f = quat(2,0,0,0)*q[1]*q[2]+quat(1,0,0,0)-quat(1,0,0,0)*q[1]-quat(1,0,0,0)*q[2]
g1 = quat(3,0,0,0)-quat(2,0,0,0)*q[1]*q[2]-quat(1,0,0,0)+quat(1,0,0,0)*q[1]+quat(1,0,0,0)*q[2]
qpop = [f, g1]
order = 2
opt,cons,I,basis= qs_tssos_first(qpop, q, n, order, numeq=0, TS=false)
sqrt(opt)
#qualitify
pop,x = quaternion_to_real(qpop,q)
opt2 = tssos_first(pop, x, order, numeq=0)

# Example2 infeasible for nonsymmetric
n = 2
@polyvar q[1:2n]
f = quat(1,0,0,0)*q[1]*q[3]+quat(1,0,0,0)*q[2]*q[4]+quat(1/2,-3/2,-3/2,-3/2)*q[1]+quat(1/2,3/2,3/2,3/2)*q[2]+quat(1/2,3/2,3/2,3/2)*q[3]+quat(1/2,-3/2,-3/2,-3/2)*q[4]
g1 = quat(-7,0,0,0)-quat(1,0,0,0)*q[1]*q[3]+quat(1,-1,-1,-1)*q[1]+quat(1,1,1,1)*q[3]-quat(1,0,0,0)*q[2]*q[4]+quat(1,-1,-1,-1)*q[2]+quat(1,1,1,1)*q[4]
g2 = quat(1,1,1,1)*q[1]+quat(1,-1,-1,-1)*q[3]+quat(2,-1,-1,-1)*q[2]+quat(2,1,1,1)*q[4]
g3 = quat(-2,-2,0,2)-quat(-3,1,1,1)*q[1]-quat(3,2,2,2)*q[2]-quat(-3,-1,-1,-1)*q[3]-quat(3,-2,-2,-2)*q[4]
qpop = [f,g1,g2,g3]
order = 2
opt,cons,I,J,basis,hbasis= qs_tssos_first(qpop, q, n, order, numeq=2, TS=false)
pop,x = quaternion_to_real(qpop,q)
opt2 = tssos_first(pop, x, 2, numeq=2,TS=false)

# Example 2 transform to one variable
#infeasible?
n = 1 
@polyvar q[1:2n]
f = quat(11/7,0,0,0)*q[1]*q[2]+quat(-6/7,0,0,0)+quat(33/14,-145/98,-157/98,-136/98)*q[1]+quat(33/14,145/98,157/98,136/98)*q[2]
g = quat(-48/7,0,0,0)-quat(11/7,0,0,0)*q[1]*q[2]+quat(0,78/49,72/49,88/49)*q[2]+quat(0,-78/49,-72/49,-88/49)*q[1]
qpop = [f,g]
order = 3
opt,cons,I,J,basis,hbasis= qs_tssos_first(qpop, q, n, order, numeq=0, TS=false)

pop,x = quaternion_to_real(qpop,q)
opt2 = tssos_first(pop, x, 4, TS=false, numeq=0)


# # construct some 
# #1
# n = 2
# @polyvar q[1:2n]
# f = quat(1,0,0,0)*q[1]*q[3]+quat(1,0,0,0)*q[2]*q[4]
# g1 = quat(10,0,0,0)-quat(1,0,0,0)*q[1]*q[3]-quat(1,0,0,0)*q[2]*q[4]
# g2 = quat(1,0,0,0)-quat(1,1,0,0)*q[1]*q[2]-quat(1,-1,0,0)*q[3]*q[4]
# qpop = [f,g2,g1]
# order = 3
# opt,cons,I,J,basis,hbasis= qs_tssos_first(qpop, q, n, order, numeq=1, TS=false)
# #qualitify
# pop,x = quaternion_to_real(qpop,q)
# opt2 = tssos_first(pop, x, 2, numeq=1,solution=true)
# #2 ??
# n = 3
# @polyvar q[1:2n]
# f = quat(1,1,0,0)*q[1]*q[3]+quat(1,-1,0,0)*q[4]*q[6]+quat(1,1,1,1)*q[2]*q[3]+quat(1,-1,-1,-1)*q[5]*q[6]
# g1 = quat(1,0,0,0)-quat(1,0,0,0)*q[1]*q[4]
# g2 = quat(1,0,0,0)-quat(1,0,0,0)*q[2]*q[5]
# g3 = quat(1,0,0,0)-quat(1,0,0,0)*q[3]*q[6]
# g = quat(1,0,0,0)-quat(1,0,0,0)*q[1]*q[4]-quat(1,0,0,0)*q[2]*q[5]-quat(1,0,0,0)*q[3]*q[6]
# qpop = [f,g1,g2,g3]
# qpop = [f,g]
# order = 7
# opt,cons,I,J,basis,hbasis= qs_tssos_first(qpop, q, n, order, numeq=0, TS=false)
# #qualitify
# pop,x = quaternion_to_real(qpop,q)
# opt2 = tssos_first(pop, x, 3, numeq=0,solution=true)
# #3
# n = 3
# @polyvar q[1:2n]
# f = quat(1,0,0,0)*q[1]+quat(1,0,0,0)*q[4]+quat(1,0,0,0)*q[2]+quat(1,0,0,0)*q[5]+quat(1,0,0,0)*q[3]+quat(1,0,0,0)*q[6]
# g1 = quat(1,0,0,0)-quat(1,0,0,0)*q[1]*q[4]
# g2 = quat(1,0,0,0)-quat(1,0,0,0)*q[2]*q[5]
# g3 = quat(1,0,0,0)-quat(1,0,0,0)*q[3]*q[6]
# g = quat(1,0,0,0)-quat(1,0,0,0)*q[1]*q[4]-quat(1,0,0,0)*q[2]*q[5]-quat(1,0,0,0)*q[3]*q[6]
# qpop = [f,g1,g2,g3]
# order = 3
# opt,cons,I,J,basis,hbasis= qs_tssos_first(qpop, q, n, order, numeq=0, TS=false)
# #qualitify
# pop,x = quaternion_to_real(qpop,q)
# opt2 = tssos_first(pop, x, 3, numeq=0,solution=true)
# #4
# rng = Xoshiro(12)
# fcoef=rand(rng,Float64,2)
# n = 1
# @polyvar q[1:2n]
# f = quat(1,0,0,0)*q[1]*q[2]
# g = quat(1,0,0,0)-quat(1,fcoef[2],fcoef[1],0)*q[1]-quat(1,-fcoef[2],-fcoef[1],0)*q[2]
# qpop = [f,g]
# order = 2
# opt,cons,I,J,basis,hbasis= qs_tssos_first(qpop, q, n, order, numeq=1, TS=false)
# #qualitify
# pop,x = quaternion_to_real(qpop,q)
# opt2 = tssos_first(pop, x, 2, numeq=1)

#Random
#n=1
rng = Xoshiro(12)
qcoef=rand(rng,QuaternionF64,1)
fcoef=rand(rng,Float64,3)
n = 1
@polyvar q[1:2n]
f = qcoef[1]*q[1]+conj(qcoef[1])*q[2]+quat(fcoef[1],0,0,0)*q[1]*q[2]
g = quat(1,0,0,0)-quat(1,0,0,0)*q[1]*q[2]
qpop = [f,g]
order = 1
opt,cons,I,J,basis,hbasis= qs_tssos_first(qpop, q, n, order, numeq=1, TS=false,QUIET=true)
#qualify
pop,x = quaternion_to_real(qpop,q)
opt2 = tssos_first(pop, x, 1, numeq=1,solution=true,QUIET=true)


#n=2
rng = Xoshiro(88)
qcoef=rand(rng,QuaternionF64,4)
fcoef=rand(rng,Float64,8)
n = 2
@polyvar q[1:2n]
# General form//wrong
f = qcoef[1]*q[1]+conj(qcoef[1])*q[3]+qcoef[2]*q[2]+conj(qcoef[2])*q[4]+qcoef[3]*q[1]*q[2]+conj(qcoef[3])*q[3]*q[4]
+qcoef[4]*q[1]*q[4]+conj(qcoef[4])*q[3]*q[2]+quat(fcoef[1],0,0,0)*q[1]*q[3]+quat(fcoef[2],0,0,0)*q[2]*q[4]
g = quat(1,0,0,0)-quat(1,0,0,0)*q[1]*q[3]-quat(1,0,0,0)*q[2]*q[4]
g1 = quat(1,0,0,0)-quat(1,0,0,0)*q[1]*q[3]
g2 = quat(1,0,0,0)-quat(1,0,0,0)*q[2]*q[4]
# qpop = [f,g]
# order = 2
# opt,cons,I,J,basis,hbasis= qs_tssos_first(qpop, q, n, order, numeq=0, TS=false)
# #qualify
qpop =[f,g1,g2]
order = 7
opt,cons,I,J,basis,hbasis= qs_tssos_first(qpop, q, n, order, numeq=0, TS=false, QUIET=true)
#qualify
pop,x = quaternion_to_real(qpop,q)
opt2 = tssos_first(pop, x, 2, TS=false, numeq=0,solution=true, QUIET=true)

#n=2
rng = Xoshiro(21)
qcoef=rand(rng,QuaternionF64,4)
fcoef=rand(rng,Float64,8)
n = 2
@polyvar q[1:2n]
# only r-part and i-part//right
f = quat(fcoef[1],fcoef[2],0,0)*q[1]+quat(fcoef[1],-fcoef[2],0,0)*q[3]
+quat(fcoef[3],fcoef[4],0,0)*q[2]+quat(fcoef[3],-fcoef[4],0,0)*q[4]
+quat(fcoef[5],fcoef[6],0,0)*q[1]*q[2]+quat(fcoef[5],-fcoef[6],0,0)*q[3]*q[4]
+quat(fcoef[7],fcoef[8],0,0)*q[1]*q[4]+quat(fcoef[7],-fcoef[8],0,0)*q[3]*q[2]
g = quat(1,0,0,0)-quat(1,0,0,0)*q[1]*q[3]-quat(1,0,0,0)*q[2]*q[4]
qpop = [f,g]
order = 7
opt,cons,I,J,basis,hbasis= qs_tssos_first(qpop, q, n, order, numeq=0, TS=false, QUIET=true)
#qualify
pop,x = quaternion_to_real(qpop,q)
opt2 = tssos_first(pop, x, 1, numeq=0,solution=true, QUIET=true)

#n=2
rng = Xoshiro(1234)
qcoef=rand(rng,QuaternionF64,4)
fcoef=rand(rng,Float64,8)
n = 2
@polyvar q[1:2n]
# only i-part//right
f = quat(0,1,0,0)*q[1]+quat(0,-1,0,0)*q[3]+quat(0,1,0,0)*q[2]+quat(0,-1,0,0)*q[4]
g = quat(1,0,0,0)-quat(1,0,0,0)*q[1]*q[3]-quat(1,0,0,0)*q[2]*q[4]
qpop = [f,g]
order = 2
opt,cons,I,J,basis,hbasis,icons,jcons= qs_tssos_first(qpop, q, n, order, numeq=0, TS=false, QUIET=true)
#qualify
pop,x = quaternion_to_real(qpop,q)
opt2 = tssos_first(pop, x, 1, numeq=0,solution=true, QUIET=true)

# only j-part//right
f = quat(0,0,1,0)*q[1]+quat(0,0,-1,0)*q[3]+quat(0,0,1,0)*q[2]+quat(0,0,-1,0)*q[4]
g = quat(1,0,0,0)-quat(1,0,0,0)*q[1]*q[3]-quat(1,0,0,0)*q[2]*q[4]
qpop = [f,g]
order = 2
opt,cons,I,J,basis,hbasis,icons,jcons= qs_tssos_first(qpop, q, n, order, numeq=0, TS=false, QUIET=true)
#qualify
pop,x = quaternion_to_real(qpop,q)
opt2 = tssos_first(pop, x, 1, numeq=0,solution=true, QUIET=true)

# i-part and j-part//wrong
f = quat(0,1,0,0)*q[1]+quat(0,-1,0,0)*q[3]+quat(0,0,1,0)*q[2]+quat(0,0,-1,0)*q[4]
g = quat(1,0,0,0)-quat(1,0,0,0)*q[1]*q[3]-quat(1,0,0,0)*q[2]*q[4]
qpop = [f,g]
order = 2
opt,cons,I,J,basis,hbasis,icons,jcons= qs_tssos_first(qpop, q, n, order, numeq=0, TS=false, QUIET=true)
#qualify
pop,x = quaternion_to_real(qpop,q)
opt2 = tssos_first(pop, x, 1, numeq=0,solution=true, QUIET=true)
