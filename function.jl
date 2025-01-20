function bfind(A, l, a)
    if l == 0
        return 0
    end
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if ndims(A) == 2
            temp = A[:, mid]
        else
            temp = A[mid]
        end
        if temp == a
            return mid
        elseif temp < a
            low = mid + 1
        else
            high = mid - 1
        end
    end
    return nothing
end
function ncbfind(A, l, a)
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if A[mid] == a
           return mid
        elseif A[mid] < a
            high = mid - 1
        else
            low = mid + 1
        end
    end
    return nothing
end
# function get_basis(n, d; nb=0, lead=[])
#     lb = binomial(n+d, d)
#     basis = zeros(UInt8, n, lb)
#     i = 0
#     t = 1
#     while i < d+1
#         t += 1
#         if basis[n,t-1] == i
#            if i < d
#               basis[1,t] = i+1
#            end
#            i += 1
#         else
#             j = findfirst(x->basis[x,t-1]!=0, 1:n)
#             basis[:,t] = basis[:,t-1]
#             if j == 1
#                basis[1,t] -= 1
#                basis[2,t] += 1
#             else
#                basis[1,t] = basis[j,t] - 1
#                basis[j,t] = 0
#                basis[j+1,t] += 1
#             end
#         end
#     end
#     if nb > 0
#         basis_bin = basis[1:nb,:]
#         basis_valid = all.(x->x<=1, eachcol(basis_bin))
#         basis = basis[:, basis_valid]
#     end
#     if !isempty(lead)
#         basis_valid = map(a->!divide(a, lead, n, size(lead,2)), eachcol(basis))
#         basis = basis[:, basis_valid]
#     end
#     return basis
# end

# function get_basis(var::Vector{UInt16}, d)
#     n = length(var)
#     lb = binomial(n+d, d)
#     basis = Vector{Vector{UInt16}}(undef, lb)
#     basis[1] = UInt16[]
#     i = 0
#     t = 1
#     while i < d+1
#         t += 1
#         if length(basis[t-1])>=i && basis[t-1][end-i+1:end] == var[n]*ones(UInt16, i)
#             if i < d
#                 basis[t] = var[1]*ones(UInt16, i+1)
#             end
#             i += 1
#         else
#             j = bfind(var, n, basis[t-1][1])
#             basis[t] = copy(basis[t-1])
#             ind = findfirst(x->basis[t][x]!=var[j], 1:length(basis[t]))
#             if ind === nothing
#                 ind = length(basis[t])+1
#             end
#             if j != 1
#                 basis[t][1:ind-2] = var[1]*ones(UInt16, ind-2)
#             end
#             basis[t][ind-1] = var[j+1]
#         end
#     end
#     return basis
# end
# function get_basis(n, d; nb=0, lead=[])
#     lb = binomial(n+d, d)
#     basis = zeros(UInt8, n, lb)
#     i = 0
#     t = 1
#     while i < d+1
#         t += 1
#         if basis[n,t-1] == i
#            if i < d
#               basis[1,t] = i+1
#            end
#            i += 1
#         else
#             j = findfirst(x->basis[x,t-1]!=0, 1:n)
#             basis[:,t] = basis[:,t-1]
#             if j == 1
#                basis[1,t] -= 1
#                basis[2,t] += 1
#             else
#                basis[1,t] = basis[j,t] - 1
#                basis[j,t] = 0
#                basis[j+1,t] += 1
#             end
#         end
#     end
#     if nb > 0
#         basis_bin = basis[1:nb,:]
#         basis_valid = all.(x->x<=1, eachcol(basis_bin))
#         basis = basis[:, basis_valid]
#     end
#     if !isempty(lead)
#         basis_valid = map(a->!divide(a, lead, n, size(lead,2)), eachcol(basis))
#         basis = basis[:, basis_valid]
#     end
#     return basis
# end

function get_ncbasis(n, d; ind=Vector{UInt16}(1:n), binary=false)
    basis = [UInt16[]]
    for i = 1:d
        append!(basis, _get_ncbasis_deg(n, i, ind=ind, binary=binary))
    end
    return basis
end

function _get_ncbasis_deg(n, d; ind=Vector{UInt16}(1:n), binary=false)
    if d > 0
        basis = Vector{UInt16}[]
        for i = 1:n
            temp = _get_ncbasis_deg(n, d-1, ind=ind, binary=binary)
            if binary == false || d == 1
                push!.(temp, ind[i])
                append!(basis, temp)
            else
                for item in temp
                    if item[end] != ind[i]
                        push!(basis, [item;ind[i]])
                    end
                end
            end
        end
        return basis
    else
        return [UInt16[]]
    end
end

function assign_constraint(m, numeq, supp::Vector{Vector{Vector{Vector{UInt16}}}}, cliques, cql)
    I = [UInt32[] for i=1:cql]
    J = [UInt32[] for i=1:cql]
    ncc = UInt32[]
    for i = 1:m
        ind = findall(k->issubset(unique(reduce(vcat, [item[1] for item in supp[i+1]])), cliques[k]), 1:cql)
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
function sadd(a, b; nb=0)
    c = [a; b]
    # sort!(c)
    # if nb > 0
    #     i = 1
    #     while i < length(c)
    #         if c[i] <= nb
    #             if c[i] == c[i+1]
    #                 deleteat!(c, i:i+1)
    #             else
    #                 i += 1
    #             end
    #         else
    #             break
    #         end
    #     end
    # end
    return c
end