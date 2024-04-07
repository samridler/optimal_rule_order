using MathOptInterface

"""
solveip(ruleevaltimes, ruleevalpass, optimizer)

Formulate and solve an instance of the integer programming problem to find
a rule evaluation order minimising total rule evaluation time, given
the matrix of rule evaluation times for each pair of rules and candidates,
and a binary matrix indicating which candidates pass each rule.

Inputs:

* `ruleevaltimes`: 2D matrix of rule evaluation times (seconds)
* `ruleevalpass`: 2D binary matrix indicating for each rule and candidate
  whether the candidate passes the rule (1) or not (0)
* `optimizer`: a MathOptInterface optimizer instance

Returns:

* An array of integers giving the rule evaluation order
"""
function solveip(
    ruleevaltimes::Matrix{Float64},
    ruleevalpass::BitMatrix,
    optimizer::MathOptInterface.OptimizerWithAttributes;
    objmax::Float64=Inf
)::Vector{Int64}

    # todo:
    # - Simplify problem when 2+ rules have the same rule results for all candidates,
    # or when 2+ candidates have the same rule results.
    # - Simplify problem when one rule fails on a subset of candidates compared to another rule
    # which takes less time per candidate; can remove the first rule.

    (m, n) = size(ruleevaltimes)

    model = Model(optimizer)

    # x[j, k] = 1 if rule j is evaluated kth, 0 otherwise
    @variable(model, x[1:n, 1:n], Bin)

    # sum of x in each row is 1 and in each column is 1
    @constraint(model, [i = 1:n], sum(x[i, k] for k = 1:n) == 1)
    @constraint(model, [k = 1:n], sum(x[i, k] for i = 1:n) == 1)

    # t[i, k] = time it takes to evaluate kth rule on candidate i
    @variable(model, t[1:m, 1:n] >= 0)

    # ruleevaltimes[candidate i, rule j] * x[rule j, position k] = rule eval time for candidate i on kth rule
    @constraint(model, time[i=1:m, k=1:n], t[i, k] == sum(ruleevaltimes[i, j] * x[j, k] for j = 1:n))

    # pass[i, k] = 1 if kth rule passes for candidate i, 0 otherwise
    @variable(model, p[1:m, 1:n], Bin)

    # ruleevalpass[candidate i, rule j] * p[rule j, position k] = 1 if rule j is passed for candidate i on kth rule
    @constraint(model, pass[i=1:m, k=1:n], p[i, k] == sum(ruleevalpass[i, j] * x[j, k] for j = 1:n))

    # e[i, k] = 1 if kth rule is evaluated for candidate i, 0 otherwise
    @variable(model, e[1:m, 1:n], Bin)

    # first rule must always be evaluated for all candidates
    @constraint(model, [i = 1:m], e[i, 1] == 1)

    # subsequent rules only need to be evaulated if the previous rules pass
    @constraint(model, [i = 1:m, k = 2:n], e[i, k] >= e[i, k-1] + p[i, k-1] - 1)

    # et[i,k] = time it takes to evaluate kth rule on candidate i if it is evaluated
    @variable(model, et[1:m, 1:n] >= 0)
    T = maximum(ruleevaltimes) # need a big-m value
    @constraint(model, etc[i=1:m, k=1:n], et[i, k] >= t[i, k] - (1 - e[i, k]) * T)

    @variable(model, obj >= 0)
    @constraint(model, obj == sum(et))
    if objmax == Inf
        objmax = sum(ruleevaltimes)
    end
    @constraint(model, obj <= objmax)

    @objective(model, Min, obj)

    optimize!(model)

    # return ordering of rules, first value gives index of first rule, etc.
    xs = value.(x)
    ruleevalorder = [findfirst(xs[:, j] .== 1) for j = 1:n]
    return ruleevalorder

    # # check:
    # xs = Int64.(value.(x) .== 1)
    # f(x) = maximum(abs.(x))
    # f(ruleevaltimes * xs - value.(t)) # should be 0
    # f(ruleevalpass * xs - value.(p)) # should be 0
    # # todo: also check that e is correct, given p
    # f(value.(et) - value.(t) .* value.(e)) # should be 0
end
