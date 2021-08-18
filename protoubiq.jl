#!/usr/bin/env julia
using ArgParse
using Glob
using Distributed, ClusterManagers
import JSON
using Base.Filesystem

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--input", "-i"
            help = ""
            nargs = '+'
            arg_type = String
            action = :store_arg
            required = true
        "--exec", "-e"
            help = ""
            required = true
            arg_type = String
        "--suffix", "-s"
            arg_type = String
            required = true
    end
    parse_args(s)
end
args = parse_commandline()
files = vcat([glob(i) for i in args["input"]]...)

@everywhere function subst(str, dir)
    c = str
    for (k, v) in dir
        c = replace(c, k => v)
    end
    c
end

@everywhere SUFFIX = $args["suffix"]
@everywhere EXEC = $args["exec"]

@sync @distributed for f = files
    @info f
    outputf = "$f.$SUFFIX"
    metaf = "$outputf.m.json"
    subst_map = Dict(
        "{{}}" => f,
        "<output>" => outputf
    )
    basecmd = Cmd([subst(e, subst_map) for e in split(EXEC, " ")])
    memfile = tempname()
    pgm = Sys.isapple() ? "gtime" : "/usr/bin/time"
    finalcmd = `$pgm -f "%M" -o $memfile $basecmd`
    @info finalcmd
    stats = @timed run(
        pipeline(finalcmd, stdout="$outputf.stdout", stderr="$outputf.stderr")
    )
    mem = parse(Int, read(memfile, String))
    open(metaf, "w+") do f
        JSON.print(f, Dict(
            "mem" => mem,
            "time" => stats.time
            ))
    end
end