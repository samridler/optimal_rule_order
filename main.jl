using JuMP
using GLPK
# using Cbc
include("utils.jl")
include("heuristic.jl")
include("ip.jl")

# problem: find ordering of rules that minimises total rule evaluation time,
# if a candidate fails a rule then the following rules are not evaluated.

# set up inputs
n = 10 # number of rules
m = 10 # number of candidates
ruleevaltimes = rand(m, n)
ruleevalpass = rand(m, n) .<= 0.9

# todo:
# - If a rule passes for all candidates then it can be removed,
#   as rule should be evaluated last.
# - If a candidate passes all rules then it can be removed,
#   as rule order will not affect rule eval time for the candidate.

@enum SolutionMethod none thumb insert repeatinsert swap reinsert ip remmeantime

ruleevalorders = Dict{SolutionMethod,Vector{Int64}}()

ruleevalorders[none] = [1:n;]
@time ruleevalorders[thumb] = ruleofthumb(ruleevaltimes, ruleevalpass)
@time ruleevalorders[insert] = insertheuristic(ruleevaltimes, ruleevalpass)
@time ruleevalorders[repeatinsert] = repeatedinsertheuristic(ruleevaltimes, ruleevalpass)
@time ruleevalorders[swap] = swapheuristic(ruleevaltimes, ruleevalpass)
@time ruleevalorders[reinsert] = reinsertheuristic(ruleevaltimes, ruleevalpass)
@time ruleevalorders[remmeantime] = remainingmeantimeheuristic(ruleevaltimes, ruleevalpass)

# GLPK seems faster than Cbc for this problem
# optimizer = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
optimizer = optimizer_with_attributes(GLPK.Optimizer, "mip_gap" => 0.0)
# @time ruleevalorders[ip] = solveip(ruleevaltimes, ruleevalpass, optimizer)
# objmax = calcruleevaltime(ruleevaltimes, ruleevalpass, ruleevalorders[thumb])
# @time ruleevalorders[ip] = solveip(ruleevaltimes, ruleevalpass, optimizer, objmax=objmax) # solving can be slower or faster with objmax

print("Solution method, total rule evaluation time:\n")
for solutionmethod in instances(SolutionMethod)
    if !haskey(ruleevalorders, solutionmethod)
        continue
    end
    t = calcruleevaltime(ruleevaltimes, ruleevalpass, ruleevalorders[solutionmethod])
    print(solutionmethod, ": ", round(t, digits=2), "\n")
end
