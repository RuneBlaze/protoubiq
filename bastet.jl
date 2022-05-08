#!/usr/bin/env julia
const shell_special = "#{}()[]<>|&*?~;"
using ArgParse, Glob
using Base.Filesystem
using DataFrames
using CSV

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--input", "-i"
            help = ""
            nargs = '+'
            arg_type = String
            action = :store_arg
            required = true
        "--range", "-r"
            help=""
            nargs = '+'
            arg_type = Int64
            nargs = 3
            default = [1, 1, 1]
        "--partition", "-p"
            help="partition"
            arg_type = String
        default = "secondary-Eth"
        "--cores", "-c"
        help="cores"
        default = 1
        "--exec", "-e"
            help = ""
            # required = true
            arg_type = String
        "--suffix", "-s"
            arg_type = String
        "--bins", "-b"
            arg_type = Int
            required = false
            default = -1
        "--output", "-o"
            arg_type = String
            # required = true
        "--across", "-a"
            arg_type = String
            nargs = '+'
        "--against"
            arg_type = String
    end
    parse_args(s)
end

args = parse_commandline()
files = vcat([glob(i) for i in args["input"]]...)

function subst(str, dir)
    c = str
    for (k, v) in dir
        c = replace(c, k => v)
    end
    c
end

function distributebins(balls, nbins = length(balls))
    bins = [[] for _=1:nbins]
    while !isempty(balls)
        for i=1:nbins
            if isempty(balls) break end
            push!(bins[i], pop!(balls))
        end
    end
    bins
end


function print_shell_word(io::IO, word::AbstractString, special::AbstractString = "")
    has_single = false
    has_special = false
    for c in word
        if isspace(c) || c=='\\' || c=='\'' || c=='"' || c=='$' || c in special
            has_special = true
            if c == '\''
                has_single = true
            end
        end
    end
    if isempty(word)
        print(io, "''")
    elseif !has_special
        print(io, word)
    elseif !has_single
        print(io, '\'', word, '\'')
    else
        print(io, '"')
        for c in word
            if c == '"' || c == '$'
                print(io, '\\')
            end
            print(io, c)
        end
        print(io, '"')
    end
    nothing
end

function showclean(io::IO, cmd::Cmd)
    print_env = cmd.env !== nothing
    print_dir = !isempty(cmd.dir)
    (print_env || print_dir) && print(io, "setenv(")
    # print(io, '`')
    join(io, map(cmd.exec) do arg
        replace(sprint(context=io) do io
            print_shell_word(io, arg, shell_special)
        end, '`' => "\\`")
    end, ' ')
    # print(io, '`')
    print_env && (print(io, ","); show(io, cmd.env))
    print_dir && (print(io, "; dir="); show(io, cmd.dir))
    (print_dir || print_env) && print(io, ")")
    nothing
end

function mk_slurm(io, fname, jobname, cmds, part, cores)
    template = read(joinpath(@__DIR__, "slurmtpl.sh"), String)
    sub = Dict(
        "<jobname>" => jobname,
        "<fullpath>" => fname,
        "<partition>" => part,
        "<cores>" => cores,
    )
    println(io, subst(template, sub))
    for c in cmds
        println(io, c)
    end
end


function levelrange(x,s,y)
    res = []
    v = collect(StepRange(x, s, y))
    for i in 1:(length(v))
        push!(res, i == length(v) ? (v[i], y) : (v[i], v[i+1]-1))
    end
    res
end

using Base.Iterators

function execmode()
    SUFFIX = args["suffix"]
    EXEC = args["exec"]
    PART = args["partition"]
    BINS = haskey(args, "bins") ? args["bins"] : -1
    RANGE = levelrange(args["range"]...)
    CORES = args["cores"]
    cmds = String[]
    @info "running on partition $PART"
    for (f, (st, ed)) = product(files, RANGE)
        @info "processing $f"
        outputf = "$f.$SUFFIX"
        noext = replace(splitext(f)[1], "." => "_")
        metaf = "$outputf.meta"
        if isfile(outputf) && filesize(outputf) > 0
            @warn "$f has already finished processing"
            continue
        end
        condition_level = 1
        cond = splitpath(f)[condition_level]
        splitted = splitpath(dirname(f))
        repname = basename(dirname(f))
        bname = basename(f)
        subst_map = Dict(
            "{{}}" => f,
            "<output>" => outputf,
            "<s>" => st,
            "<e>" => ed,
            "<repname>" => repname,
            "<repname-1>" => splitted[max(end - 1, 1)],
            "<cond>" => cond,
            "<noext>" => noext,
            "<bname>" => bname,
            "<bname->" => replace(bname, "." => "_"),
            "{//}" => dirname(f),
            "{.}" => splitext(f)[1],
        )
        basecmd = subst(EXEC, subst_map) #Cmd([subst(e, subst_map) for e in split(EXEC, " ")])
        # memfile = tempname()
        pgm = Sys.isapple() ? "gtime" : "/usr/bin/time"
        finalcmd = "$pgm -f '%e;%M' -o $metaf $basecmd"
        @info "final command: $finalcmd"

        push!(cmds, finalcmd)
    end

    bins = distributebins(cmds, BINS == -1 ? length(cmds) : BINS)
    mkpath(args["output"])
    for (i,b)=enumerate(bins)
        of = joinpath(args["output"], "$SUFFIX.$i.sh")
        @info "written to $of"
        open(of, "w+") do f
            mk_slurm(f, of, "$SUFFIX.$i", b, PART, CORES)
        end
    end

    function commandhelper(output)
        "find $output -name \"*.sh\" | xargs -I '{}' sbatch {}"
    end

    println("Generation finished!")
    c = commandhelper(args["output"])
    println("Run `$c` to queue your job")
end

function againstmode()
    # we now compare phylogenetic trees, again
    # we first make a temporary dataframe of all the parameters
    # and leave the heavy-lifting to Python
    condition_level = 1
    #tmp = splitpath(files[1])[condition_level]
    #for f = files
    #    if splitpath(f)[condition_level] != tmp
    #        condition_level += 1
    #        break
    #    end
    #end

    df = DataFrame(
        inputpath = String[],
        condition = String[],
        streepath = String[],
        method = String[],
        k = Int[],
        gtreepath = String[],
    )
    for f = files
        cond = splitpath(f)[condition_level]
        splitted = splitpath(dirname(f))
        repname = basename(dirname(f))
        for a = args["across"]
            ip = "$f.$a"
            if !isfile(ip) || filesize(ip) <= 0
                @warn "$ip does not exist for method $a. Embrace yourself."
                continue
            end
            sp = subst(args["against"], Dict("<repname>" => repname, "<repname-1>" => splitted[max(end - 1, 1)], "<cond>" => cond))
            
            #k = parse(Int, split(f, ".")[end-1])
            k = parse(Int, split(f, ".")[end][2:end])
            push!(df, (ip, cond, sp, a, k, f))
        end
    end
    csvp = tempname()
    CSV.write(csvp, df)
    reporter = joinpath(@__DIR__,"reporter.py")
    if !isfile(reporter)
        reporter = joinpath(@__DIR__,"protoubiq/reporter.py")
    end
    run(`python3 $reporter -i $csvp`)
end

if !isnothing(args["exec"])
    execmode()
elseif !isnothing(args["against"])
    againstmode()
end
