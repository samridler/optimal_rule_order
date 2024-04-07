include("utils.jl")

function swapheuristic(ruleevaltimes::Matrix{Float64}, ruleevalpass::BitMatrix)::Vector{Int64}
    # todo:
    # - make recalculating of rule eval time more efficient

    (m, n) = size(ruleevaltimes)
    ruleevalorder = [1:n;]
    firstrulefailindex = calcfirstrulefailindex(ruleevalpass, ruleevalorder)

    # create copies of vectors to be re-used
    ruleevalordertemp = copy(ruleevalorder)
    firstrulefailindextemp = copy(firstrulefailindex)

    # perform swapping of rules that decreases rule evaluation time
    ruleevaltimedecreased = true
    ruleevalorderhashes = Set(hash(ruleevalorder))
    while ruleevaltimedecreased
        ruleevaltimedecreased = false
        for k = 1:n
            for l = 1:n
                if k == l
                    continue
                end

                copy!(ruleevalordertemp, ruleevalorder)
                copy!(firstrulefailindextemp, firstrulefailindex)

                # swap rule at position k with rule at position l
                ruleevalordertemp[k], ruleevalordertemp[l] = ruleevalordertemp[l], ruleevalordertemp[k]

                # calculate rule evaluation time unless rule eval order has been tried before
                hashvalue = hash(ruleevalordertemp)
                if hashvalue in ruleevalorderhashes
                    dt = Inf
                else
                    dt = calcdeltaruleevaltime_swaprules!(
                        ruleevaltimes, ruleevalpass, ruleevalorder, firstrulefailindextemp, k, l)
                    push!(ruleevalorderhashes, hashvalue)
                end

                if dt < 0
                    copy!(ruleevalorder, ruleevalordertemp)
                    copy!(firstrulefailindex, firstrulefailindextemp)
                    ruleevaltimedecreased = true
                end
            end
        end
    end

    return ruleevalorder
end

function insertheuristic(ruleevaltimes::Matrix{Float64}, ruleevalpass::BitMatrix)::Vector{Int64}
    # todo:
    # - make recalculating of rule eval time more efficient

    (m, n) = size(ruleevaltimes)
    ruleevalorder = [1:n;]
    ruleevaltime = calcruleevaltime(ruleevaltimes, ruleevalpass, ruleevalorder)

    # perform re-insertion of a rule to a different index to decrease rule evaluation time
    ruleevaltimedecreased = true
    ruleevalorderhashes = Set(hash(ruleevalorder))
    while ruleevaltimedecreased
        ruleevaltimedecreased = false
        for k = 1:n
            for l = 1:n
                if k == l
                    continue
                end
                reinsert!(ruleevalorder, k, l)

                # calculate rule evaluation time unless rule eval order has been tried before
                hashvalue = hash(ruleevalorder)
                if hashvalue in ruleevalorderhashes
                    ruleevaltimenew = Inf
                else
                    ruleevaltimenew = calcruleevaltime(ruleevaltimes, ruleevalpass, ruleevalorder)
                    push!(ruleevalorderhashes, hashvalue)
                end

                if ruleevaltimenew < ruleevaltime
                    ruleevaltime = ruleevaltimenew
                    ruleevaltimedecreased = true
                else
                    # change rule eval order back
                    reinsert!(ruleevalorder, l, k)
                end
            end
        end
    end

    return ruleevalorder
end

function reinsert!(vector::Vector{T}, k::Int64, l::Int64) where {T<:Any}
    # move items in vector in place so that item in index k is moved to index l, shifting other items as needed
    if k == l
        return
    end

    v = vector[k]
    if k < l
        # shift items left by 1
        vector[k:l-1] = vector[k+1:l]
    else # k > l
        # shift items right by 1
        vector[l+1:k] = vector[l:k-1]
    end
    vector[l] = v
end

function insertheuristic2(ruleevaltimes::Matrix{Float64}, ruleevalpass::BitMatrix)::Vector{Int64}
    # note:
    # This code runs quite slow, and I'm not sure if the algorithm is correct,
    # may be more effective to improve the original insertheuristic method but
    # keep this anyway for reference.

    (m, n) = size(ruleevaltimes)
    ruleevalorder = [1:n;]
    ruleevaltime = calcruleevaltime(ruleevaltimes, ruleevalpass, ruleevalorder)

    firstrulefailindex = calcfirstrulefailindex(ruleevalpass, ruleevalorder)

    # perform re-insertion of a rule to a different index to decrease rule evaluation time
    ruleevaltimedecreased = true
    while ruleevaltimedecreased
        ruleevaltimedecreased = false
        for i = 1:n
            dtremove, ruleevalordertemp, firstrulefailindextemp = calcdeltaruleevaltime_removerule(
                ruleevaltimes, ruleevalpass, ruleevalorder, firstrulefailindex, i)
            for j = 1:n
                if i == j
                    continue
                end
                dtinsert, ruleevalordernew, firstrulefailindexnew = calcdeltaruleevaltime_insertrule(
                    ruleevaltimes, ruleevalpass, ruleevalordertemp, firstrulefailindextemp, ruleevalorder[j], j)
                ruleevaltimenew = ruleevaltime + dtremove + dtinsert
                if ruleevaltimenew < ruleevaltime
                    ruleevalorder = ruleevalordernew
                    ruleevaltime = ruleevaltimenew
                    firstrulefailindex = firstrulefailindexnew
                    ruleevaltimedecreased = true

                    dtremove, ruleevalordertemp, firstrulefailindextemp = calcdeltaruleevaltime_removerule(
                        ruleevaltimes, ruleevalpass, ruleevalorder, firstrulefailindex, i)
                end
            end
        end
    end

    return ruleevalorder
end
