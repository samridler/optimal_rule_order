include("utils.jl")

function ruleofthumb(ruleevaltimes::Matrix{Float64}, ruleevalpass::BitMatrix)::Vector{Int64}
    # find rule with highest ratio of fail rate to evaluation time, this rule should be evaluated first
    # then remove candidates that fail that rule and repeat

    (m, n) = size(ruleevaltimes)
    ruleevalorder = []

    ruleevalfails = map(!, ruleevalpass)

    rows = ones(Bool, m)
    cols = ones(Bool, n)

    # find rule with highest ratio of fail rate to evaluation time
    for iter = 1:n
        ruleevalfail = sum(view(ruleevalfails, rows, cols), dims=1)
        ruleevaltime = sum(view(ruleevaltimes, rows, cols), dims=1)
        _, k = findmax(i -> ruleevalfail[i] / ruleevaltime[i], 1:n-iter+1)
        # k is the index in the view of the remaining rules, need to get back the original rule index
        j = view(1:n, cols)[k]
        push!(ruleevalorder, j)
        # remove the rule from consideration, and remove remaining candidates that fail that rule
        cols[j] = false
        for (i, row) in enumerate(rows)
            if row && ruleevalfails[i, j]
                rows[i] = false
            end
        end
    end

    @assert(sort(ruleevalorder) == 1:n)

    return ruleevalorder
end

function insertheuristic(
    ruleevaltimes::Matrix{Float64},
    ruleevalpass::BitMatrix;
    ruleinsertorder::Vector{Int64}=Int64[]
)::Vector{Int64}
    (m, n) = size(ruleevaltimes)
    if isempty(ruleinsertorder)
        ruleinsertorder = [1:n;]
    end

    # insert rules one at a time based on ruleinsertorder, greedily minimising rule evaluation time

    ruleevalorder = []
    if n == 0
        return ruleevalorder
    end

    ruleevalorder = [ruleinsertorder[1]]
    firstrulefailindex = calcfirstrulefailindex(ruleevalpass, ruleevalorder)

    # create copy of vector to be re-used
    firstrulefailindextemp = copy(firstrulefailindex)

    for k = 2:n
        # put new rule last then try re-insert it in each other position
        push!(ruleevalorder, ruleinsertorder[k])
        firstrulefailindex = calcfirstrulefailindex(ruleevalpass, ruleevalorder)
        bestdt = 0.0
        bestindex = k
        for l = 1:k-1
            copy!(firstrulefailindextemp, firstrulefailindex)
            dt = calcdeltaruleevaltime_reinsertrule!(
                ruleevaltimes, ruleevalpass, ruleevalorder, firstrulefailindextemp, k, l)
            if dt < bestdt
                bestdt = dt
                bestindex = l
            end
        end

        reinsert!(ruleevalorder, k, bestindex)
    end

    return ruleevalorder
end

function repeatedinsertheuristic(
    ruleevaltimes::Matrix{Float64},
    ruleevalpass::BitMatrix;
    ruleinsertorder::Vector{Int64}=Int64[]
)::Vector{Int64}
    if isempty(ruleinsertorder)
        n = size(ruleevaltimes, 2)
        ruleinsertorder = [1:n;]
    end

    # repeat the insertheuristic until no change is made to rule order
    ruleevalorder = insertheuristic(ruleevaltimes, ruleevalpass, ruleinsertorder=ruleinsertorder)
    while true
        ruleevalordernew = insertheuristic(ruleevaltimes, ruleevalpass, ruleinsertorder=ruleevalorder)
        if ruleevalorder == ruleevalordernew
            return ruleevalorder
        end
        ruleevalorder = ruleevalordernew
    end
end

function swapheuristic(
    ruleevaltimes::Matrix{Float64},
    ruleevalpass::BitMatrix;
    ruleevalorder::Vector{Int64}=Int64[]
)::Vector{Int64}
    (m, n) = size(ruleevaltimes)
    ruleevalorder = isempty(ruleevalorder) ? [1:n;] : copy(ruleevalorder)
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
                    continue
                end

                dt = calcdeltaruleevaltime_swaprules!(
                    ruleevaltimes, ruleevalpass, ruleevalorder, firstrulefailindextemp, k, l)
                push!(ruleevalorderhashes, hashvalue)

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

function reinsertheuristic(
    ruleevaltimes::Matrix{Float64},
    ruleevalpass::BitMatrix;
    ruleevalorder::Vector{Int64}=Int64[]
)::Vector{Int64}
    (m, n) = size(ruleevaltimes)
    ruleevalorder = isempty(ruleevalorder) ? [1:n;] : copy(ruleevalorder)
    firstrulefailindex = calcfirstrulefailindex(ruleevalpass, ruleevalorder)

    # create copies of vectors to be re-used
    ruleevalordertemp = copy(ruleevalorder)
    firstrulefailindextemp = copy(firstrulefailindex)

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

                copy!(ruleevalordertemp, ruleevalorder)
                copy!(firstrulefailindextemp, firstrulefailindex)

                reinsert!(ruleevalordertemp, k, l)

                # calculate rule evaluation time unless rule eval order has been tried before
                hashvalue = hash(ruleevalordertemp)
                if hashvalue in ruleevalorderhashes
                    continue
                end

                dt = calcdeltaruleevaltime_reinsertrule!(
                    ruleevaltimes, ruleevalpass, ruleevalorder, firstrulefailindextemp, k, l)
                push!(ruleevalorderhashes, hashvalue)

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

function reinsertheuristic2(ruleevaltimes::Matrix{Float64}, ruleevalpass::BitMatrix)::Vector{Int64}
    # note:
    # This code runs quite slow, and I'm not sure if the algorithm is correct,
    # may be more effective to improve the original reinsertheuristic method but
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
