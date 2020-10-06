#
# implementation of walter's method
#

module WalterMethod

export WalterTree
export move, update!

mutable struct WalterTree
    k::Int
    p::Float64
    sum::Float64
    left::Union{Nothing, WalterTree}
    right::Union{Nothing, WalterTree}

    function WalterTree(ps::Vector, depth=0, maxdepth=nothing)
        if depth == 0
            N = length(ps)
            maxdepth = ceil(Int, log2(N + 1)) - 1
            ps = collect(enumerate(vcat(ps, zeros(2^(maxdepth + 1) - 1 - N))))
        end
        k = div(2^(maxdepth - depth + 1) - 1, 2) + 1
        total = sum([p for (_, p) in ps])
        if depth < maxdepth
            left = WalterTree(ps[1:k-1], depth+1, maxdepth)
            right = WalterTree(ps[k+1:end], depth+1, maxdepth)
        else
            left = nothing
            right = nothing
        end
        return new(ps[k][1], ps[k][2], total, left, right)
    end
end

function move(tree::WalterTree, x=nothing)
    x = tree.sum * rand()
    if tree.left == nothing
        return tree.k
    else
        if x < tree.left.sum
            return move(tree.left, x)
        elseif x > tree.left.sum + tree.p
            return move(tree.right, x - tree.left.sum - tree.p)
        else
            return tree.k
        end
    end
end

function update!(tree::WalterTree, Δps::Vector{Tuple{Int, Float64}})
    for (k, Δp) in Δps
        update!(tree, k, Δp)
    end
end

function update!(tree::WalterTree, k::Int, Δp::Float64)
    tree.sum += Δp
    if k < tree.k
        update!(tree.left, k, Δp)
    elseif k > tree.k
        update!(tree.right, k, Δp)
    else
        tree.p += Δp
    end
end


end
