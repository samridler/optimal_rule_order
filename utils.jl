# calculate total rule evaluation time for a given rule evaluation order
function calcruleevaltime(
    ruleevaltimes::Matrix{Float64},
    ruleevalpass::BitMatrix,
    ruleevalorder::Vector{Int64}
)::Float64
    (m, n) = size(ruleevaltimes)
    time = 0.0
    for i = 1:m
        # at times until the first rule that fails (if any) for candidate i
        for j = 1:n
            k = ruleevalorder[j]
            time += ruleevaltimes[i, k]
            if !ruleevalpass[i, k]
                break
            end
        end
    end

    return time
end

function calcfirstrulefailindex(
    ruleevalpass::BitMatrix,
    ruleevalorder::Vector{Int64}
)::Vector{Int64}
    (m, n) = size(ruleevaltimes)
    # firstrulefailindex[i] gives index of first rule that fails in ruleevalorder,
    # index of n + 1 means that no rule fails
    firstrulefailindex = fill(n + 1, m)
    for i = 1:m
        for k = 1:n
            j = ruleevalorder[k]
            if !ruleevalpass[i, j]
                firstrulefailindex[i] = k
                break
            end
        end
    end
    return firstrulefailindex
end

function calcdeltaruleevaltime_removerule!(
    ruleevaltimes::Matrix{Float64},
    ruleevalpass::BitMatrix,
    ruleevalorder::Vector{Int64},
    firstrulefailindex::Vector{Int64},
    removeruleindex::Int64,
)::Tuple{Float64,Vector{Int64},Vector{Int64}}
    # Calculate change in rule evaluation time if rule at index firstrulefailindex is removed from ruleevalorder.

    ruleevalorder = copy(ruleevalorder)
    firstrulefailindex = copy(firstrulefailindex)

    # calculate change in rule eval time
    dt = 0.0
    (m, n) = size(ruleevaltimes)
    for i = 1:m
        if firstrulefailindex[i] == removeruleindex
            # find next rule that fails
            firstrulefailindex[i] = n + 1
            for k = removeruleindex+1:n
                j = ruleevalorder[k]
                dt += ruleevaltimes[i, j]
                if !ruleevalpass[i, j]
                    firstrulefailindex[i] = k
                    break
                end
            end
        end

        if firstrulefailindex[i] >= removeruleindex
            firstrulefailindex[i] -= 1
            dt -= ruleevaltimes[i, ruleevalorder[removeruleindex]]
        end
    end

    popat!(ruleevalorder, removeruleindex)

    return dt, ruleevalorder, firstrulefailindex
end

function calcdeltaruleevaltime_insertrule!(
    ruleevaltimes::Matrix{Float64},
    ruleevalpass::BitMatrix,
    ruleevalorder::Vector{Int64},
    firstrulefailindex::Vector{Int64},
    ruleindex::Int64, # index of rule to be inserted
    insertindex::Int64, # index where rule will be inserted in ruleevalorder
)::Tuple{Float64,Vector{Int64},Vector{Int64}}
    # Calculate change in rule evaluation time if rule (with index ruleindex)
    # is inserted into ruleevalorder at index insertindex.

    ruleevalorder = copy(ruleevalorder)
    firstrulefailindex = copy(firstrulefailindex)

    # calculate change in rule eval time
    dt = 0.0
    (m, n) = size(ruleevaltimes)
    for i = 1:m
        if firstrulefailindex[i] >= insertindex
            if ruleevalpass[i, ruleindex]
                firstrulefailindex[i] += 1
            else
                # subtract time for rules that no longer need evaluation
                k = firstrulefailindex[i] <= n ? firstrulefailindex[i] : n
                dt -= sum(ruleevaltimes[i, ruleevalorder[insertindex:k]])

                firstrulefailindex[i] = insertindex
            end
            dt += ruleevaltimes[i, ruleindex]
        end
    end

    insert!(ruleevalorder, insertindex, ruleindex)

    return dt, ruleevalorder, firstrulefailindex
end