using Random
include("Qpop.jl")

#Test1
##n=20 seed=4,2,3;n=40 seed=1,2,3;n=60 seed=11,12,13

##generate objective function and constraints
rng = Xoshiro(11)
n = 60
@ncpolyvar q[1:2n]
f=randomsymfunc(q,n,1,rng;conjugates=false,coelimit=false)
g=1
for i=1:n
    g=g+(-1)*q[i]*q[i+n]
end
gn=[]
for j=1:n
    temp=1+(-1)*q[j]*q[j+n]
    push!(gn,temp)
end
order = 1

###ball
qpop=[f,g]
opt= qs_tssos_first(qpop, q, n, order, numeq=0,TS=false,conjubasis=false) #conjubasis=false means using basis [q]_d
# opt= qs_tssos_first(qpop, q, n, order, numeq=0,TS=false,conjubasis=true) #conjubasis=true means using basis [q,\bar(q)]_d

#### tssos qualify
pop,x= quaternion_to_real(qpop,q)
opt2= tssos_first(pop, x, order, numeq=0,TS=false,solution=true)
# opt2,sol,data = tssos_first(pop, x, order, numeq=0,TS=false,solve=false)
# println(data.blocksize[1][1])

# ###sphere
# qpop = [f,g]
# opt= qs_tssos_first(qpop, q, n,order, numeq=1,TS=false,conjubasis=false)
# ####qualify
# pop,x = quaternion_to_real(qpop,q)
# opt2 = tssos_first(pop, x, order, numeq=1,TS=false,solution=true)

###unit norm
qpop = append!([f],gn)
opt= qs_tssos_first(qpop, q, n, order,numeq=n,TS=false,conjubasis=false)
#turn on nb
opt= qs_tssos_first(qpop, q, n, order,numeq=n,nb=n,TS=false,conjubasis=false)

####qualify
pop,x = quaternion_to_real(qpop,q)
opt2 = tssos_first(pop, x, order,numeq=n ,TS=false,solution=true)



#Test2
# n=20,40,60 seed=1,2,3;
rng = Xoshiro(1)
n = 60
@ncpolyvar q[1:2n]
f=randomsymfunc(q,n,1,rng;conjugates=true,coelimit=false)
g=1
for i=1:n
    g=g+(-1)*q[i]*q[i+n]
end
gn=[]
for j=1:n
    temp=1+(-1)*q[j]*q[j+n]
    push!(gn,temp)
end
order = 1

###ball
qpop=[f,g]
opt= qs_tssos_first(qpop, q, n, order,numeq=0,TS=false,conjubasis=true)

####qualify
pop,x= quaternion_to_real(qpop,q)
opt2= tssos_first(pop, x, order, numeq=0,TS=false,solution=true)

###unit norm
qpop = append!([f],gn)
opt= qs_tssos_first(qpop, q, n, order,numeq=n,TS=false,conjubasis=true)
#turn on nb
opt= qs_tssos_first(qpop, q, n, order,numeq=n,nb=n,TS=false,conjubasis=true)

####qualify
pop,x = quaternion_to_real(qpop,q)
opt2= tssos_first(pop, x, order, numeq=n,TS=false,solution=true)



#Test3
# n=1,2,3,4,5,6 seed=11,12,13;
rng = Xoshiro(11)
n = 1
@ncpolyvar q[1:2n]
f=randomsymfunc(q,n,2,rng;conjugates=true,coelimit=true)
g=1
for i=1:n
    g=g+(-1)*q[i]*q[i+n]
end
gn=[]
for j=1:n
    temp=1+(-1)*q[j]*q[j+n]
    push!(gn,temp)
end
order = 2

###ball
qpop=[f,g]
opt= qs_tssos_first(qpop, q, n, order,numeq=0,TS=false,conjubasis=true)

####qualify
pop,x= quaternion_to_real(qpop,q)
opt2 = tssos_first(pop, x, order, numeq=0,TS=false,solution=true)
# opt2,sol,data = tssos_first(pop, x, order, numeq=0,TS=false,solve=false)
# println(data.blocksize[1][1])

# ###sphere
# qpop = [f,g]
# opt= qs_tssos_first(qpop, q, n, order, numeq=1, TS=false,conjubasis=true)
# ####qualify
# pop,x = quaternion_to_real(qpop,q)
# opt2 = tssos_first(pop, x, order, numeq=1,TS=false,solution=true)

###unit norm
qpop = append!([f],gn)
opt= qs_tssos_first(qpop, q, n, order, numeq=n, TS=false,conjubasis=true)
#turn on nb
opt= qs_tssos_first(qpop, q, n, order, numeq=n,nb=n, TS=false,conjubasis=true)

####qualify
pop,x = quaternion_to_real(qpop,q)
opt2 = tssos_first(pop, x, order,numeq=n ,TS=false,solution=true)



#Test4
# n=1,2,3,4,5,6 seed=30,60,90;
rng = Xoshiro(90)
n = 5
@ncpolyvar q[1:2n]
f=randomsymfunc(q,n,2,rng;conjugates=false,coelimit=true)
g=1
for i=1:n
    g=g+(-1)*q[i]*q[i+n]
end
gn=[]
for j=1:n
    temp=1+(-1)*q[j]*q[j+n]
    push!(gn,temp)
end
order = 2

###ball
qpop=[f,g]
opt= qs_tssos_first(qpop, q, n, order,numeq=0,TS=false,conjubasis=true)
opt= qs_tssos_first(qpop, q, n, order,numeq=0,TS=false,conjubasis=false)

####qualify
pop,x= quaternion_to_real(qpop,q)
opt2 = tssos_first(pop, x, order, numeq=0,TS=false,solution=true)

# ###sphere
# qpop = [f,g]
# opt= qs_tssos_first(qpop, q, n, order, numeq=1, TS=false)
# ####qualify
# pop,x = quaternion_to_real(qpop,q)
# opt2 = tssos_first(pop, x, 2, numeq=1,TS=false,solution=true)

###unit norm
qpop = append!([f],gn)
opt= qs_tssos_first(qpop, q, n, order,numeq=n,TS=false,conjubasis=true)
#turn on nb
opt= qs_tssos_first(qpop, q, n, order,numeq=n,nb=n,TS=false,conjubasis=false)

####qualify
pop,x = quaternion_to_real(qpop,q)
opt2,sol,data = tssos_first(pop, x, order,numeq=n ,TS=false,solution=true)
println(data.blocksize[1][1])

