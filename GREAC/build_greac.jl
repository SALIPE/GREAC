using Pkg
using PackageCompiler

Pkg.activate(".")

output_dir = joinpath(@__DIR__, "greac")

println("Iniciando a compilação do executável GREAC com create_app...")
create_app(
    @__DIR__,
    output_dir;
    force=true,
    cpu_target="x86-64",
    incremental=true,
    # sysimage_build_args=`-O2`,
    precompile_execution_file=joinpath(@__DIR__, "src", "GREAC_warmup.jl"))
println("Compilação concluída! O executável está em: ", output_dir)

